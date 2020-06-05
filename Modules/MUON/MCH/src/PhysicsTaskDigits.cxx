///
/// \file   PhysicsTaskDigits.cxx
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <algorithm>

#include "Headers/RAWDataHeader.h"
#include "DPLUtils/DPLRawParser.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
//#include "MCHPreClustering/PreClusterFinder.h"
#include "QualityControl/QcInfoLogger.h"
#include "MCH/PhysicsTaskDigits.h"

using namespace std;

#define QC_MCH_SAVE_TEMP_ROOTFILE 1

static FILE* flog = NULL;

struct CRUheader {
  uint8_t header_version;
  uint8_t header_size;
  uint16_t block_length;
  uint16_t fee_id;
  uint8_t priority_bit;
  uint8_t reserved_1;
  uint16_t next_packet_offset;
  uint16_t memory_size;
  uint8_t link_id;
  uint8_t packet_counter;
  uint16_t source_id;
  uint32_t hb_orbit;
  //uint16_t cru_id;
  //uint8_t dummy1;
  //uint64_t dummy2;
};

enum decode_state_t {
  DECODE_STATE_UNKNOWN,
  DECODE_STATE_SYNC_FOUND,
  DECODE_STATE_HEADER_FOUND,
  DECODE_STATE_CSIZE_FOUND,
  DECODE_STATE_CTIME_FOUND,
  DECODE_STATE_SAMPLE_FOUND
};

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{
PhysicsTaskDigits::PhysicsTaskDigits() : TaskInterface(), count(1) {}

PhysicsTaskDigits::~PhysicsTaskDigits() { fclose(flog); }

void PhysicsTaskDigits::initialize(o2::framework::InitContext& /*ctx*/)
{
  QcInfoLogger::GetInstance() << "initialize PhysicsTaskDigits" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(stdout, "initialize PhysicsTaskDigits\n");

  mDecoder.initialize();

  mPrintLevel = 1;

  flog = stdout; //fopen("/root/qc.log", "w");
  fprintf(stdout, "PhysicsTaskDigits initialization finished\n");

  uint32_t dsid;
  for (int cruid = 0; cruid < 3; cruid++) {
    QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE CRUID " << cruid << AliceO2::InfoLogger::InfoLogger::endm;
    for (int linkid = 0; linkid < 24; linkid++) {
      QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE LINKID " << linkid << AliceO2::InfoLogger::InfoLogger::endm;

      {
        int index = 24 * cruid + linkid;
        mHistogramNhits[index] = new TH2F(TString::Format("QcMuonChambers_NHits_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - Number of hits (CRU link %02d)", index), 40, 0, 40, 64, 0, 64);
        //mHistogramPedestals->SetDrawOption("col");
        //getObjectsManager()->startPublishing(mHistogramNhits[index]);

        mHistogramADCamplitude[index] = new TH1F(TString::Format("QcMuonChambers_ADC_Amplitude_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - ADC amplitude (CRU link %02d)", index), 5000, 0, 5000);
        //mHistogramPedestals->SetDrawOption("col");
        //getObjectsManager()->startPublishing(mHistogramADCamplitude[index]);
      }

      int32_t link_id = mDecoder.getMapCRU(cruid, linkid);
      QcInfoLogger::GetInstance() << "  LINK_ID " << link_id << AliceO2::InfoLogger::InfoLogger::endm;
      if (link_id == -1)
        continue;
      for (int ds_addr = 0; ds_addr < 40; ds_addr++) {
        QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE DS_ADDR " << ds_addr << AliceO2::InfoLogger::InfoLogger::endm;
        uint32_t de;
        int32_t result = mDecoder.getMapFEC(link_id, ds_addr, de, dsid);
        QcInfoLogger::GetInstance() << "C'EST LA LIGNE APRÈS LE GETMAPFEC, DE " << de << "  RESULT " << result << AliceO2::InfoLogger::InfoLogger::endm;
        if(result < 0) continue;

        if (std::find(DEs.begin(), DEs.end(), de) == DEs.end()) {
          DEs.push_back(de);

          TH1F* h = new TH1F(TString::Format("QcMuonChambers_ADCamplitude_DE%03d", de),
              TString::Format("QcMuonChambers - ADC amplitude (DE%03d)", de), 5000, 0, 5000);
          mHistogramADCamplitudeDE.insert(make_pair(de, h));
          //getObjectsManager()->startPublishing(h);

          float Xsize = 40 * 5;
          float Xsize2 = Xsize / 2;
          float Ysize = 50;
          float Ysize2 = Ysize / 2;

          TH2F* h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d", de),
              TString::Format("QcMuonChambers - Number of hits (DE%03d)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsDE.insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);
          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d", de),
              TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsHighAmplDE.insert(make_pair(de, h2));
          //getObjectsManager()->startPublishing(h2);
        }
      }
    }
  }

  for(int de = 1; de <= 1030; de++) {
    const o2::mch::mapping::Segmentation* segment = &(o2::mch::mapping::segmentation(de));
    if (segment == nullptr) continue;

    float Xsize = 40 * 5;
    float Xsize2 = Xsize / 2;
    float Ysize = 50;
    float Ysize2 = Ysize / 2;
    float scale = 0.5;
      {
    TH2F* hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_B_XY_%03d", de),
          TString::Format("QcMuonChambers - Occupancy XY (DE%03d B)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramOccupancyXY[0].insert(make_pair(de, hXY));

//      hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_NB_XY_%03d", de),
//           TString::Format("QcMuonChambers - Occupancy XY (DE%03d NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
//      mHistogramOccupancyXY[1].insert(make_pair(de, hXY));
//
//      hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_BNB_XY_%03d", de),
//           TString::Format("QcMuonChambers - Occupancy XY (DE%03d B+NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
//      mHistogramOccupancyXY[2].insert(make_pair(de, hXY));
    }
  }
    
  mHistogramOccupancy[0] = new GlobalHistogram("QcMuonChambers_Occupancy_den", "Occupancy");
  mHistogramOccupancy[0]->init();
//  mHistogramOccupancy[1] = new GlobalHistogram("QcMuonChambers_Occupancy", "Occupancy");
//  mHistogramOccupancy[1]->init();
//  mHistogramOccupancy[2] = new GlobalHistogram("QcMuonChambers_Occupancy_BNB", "Occupancy - B+NB");
//  mHistogramOccupancy[2]->init();
}

void PhysicsTaskDigits::startOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "startOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::startOfCycle()
{
  QcInfoLogger::GetInstance() << "startOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::monitorDataReadout(o2::framework::ProcessingContext& ctx)
{
  // Reset the containers
  mDecoder.reset();

  // For some reason the input selection doesn't work, to be investigated...
  o2::framework::DPLRawParser parser(ctx.inputs() /*, o2::framework::select("readout:MCH/RAWDATA")*/);

  for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
    // retrieving RDH v4
    auto const* rdh = it.get_if<o2::header::RAWDataHeaderV4>();
    if (!rdh)
      continue;
    // retrieving the raw pointer of the page
    auto const* raw = it.raw();
    // size of payload
    size_t payloadSize = it.size();
    if (payloadSize == 0)
      continue;

    if (mPrintLevel >= 1) {
      std::cout<<"payloadSize: "<<payloadSize<<std::endl;
      //std::cout<<"raw:     "<<(void*)raw<<std::endl;
      //std::cout<<"payload: "<<(void*)payload<<std::endl;
    }

    // Run the decoder on the CRU buffer
    mDecoder.processData((const char*)raw, (size_t)(payloadSize + sizeof(o2::header::RAWDataHeaderV4)));
  }

  std::vector<SampaHit>& hits = mDecoder.getHits();
  if (mPrintLevel >= 1)
    fprintf(flog, "hits.size()=%d\n", (int)hits.size());
  for (uint32_t i = 0; i < hits.size(); i++) {
    //continue;
    SampaHit& hit = hits[i];
    if (mPrintLevel >= 1)
      fprintf(stdout, "hit[%d]: link_id=%d, ds_addr=%d, chan_addr=%d\n",
          i, hit.link_id, hit.ds_addr, hit.chan_addr);
    if (hit.link_id >= 24 || hit.ds_addr >= 40 || hit.chan_addr >= 64) {
      fprintf(stdout, "hit[%d]: link_id=%d, ds_addr=%d, chan_addr=%d\n",
          i, hit.link_id, hit.ds_addr, hit.chan_addr);
      continue;
    }
    //if( hit.csum > 500 ) {
    mHistogramNhits[hit.link_id]->Fill(hit.ds_addr, hit.chan_addr);
    mHistogramADCamplitude[hit.link_id]->Fill(hit.csum);
    //}
  }

  std::vector<o2::mch::Digit>& digits = mDecoder.getDigits();
  if (mPrintLevel >= 1)
    fprintf(flog, "digits.size()=%d\n", (int)digits.size());
  for (uint32_t i = 0; i < digits.size(); i++) {
    o2::mch::Digit& digit = digits[i];
    plotDigit(digit);
  }
}


void PhysicsTaskDigits::monitorDataDigits(const o2::framework::DataRef& input)
{
  if (input.spec->binding != "digits")
    return;

  const auto* header = o2::header::get<header::DataHeader*>(input.header);
  if (mPrintLevel >= 1)
    fprintf(flog, "Header: %p\n", (void*)header);
  if (!header)
    return;
  //QcInfoLogger::GetInstance() << "payloadSize: " << header->payloadSize << AliceO2::InfoLogger::InfoLogger::endm;
  if (mPrintLevel >= 1)
    fprintf(flog, "payloadSize: %d\n", (int)header->payloadSize);
  if (mPrintLevel >= 1)
    fprintf(flog, "payload: %p\n", input.payload);

  std::vector<o2::mch::Digit> digits{ 0 };
  o2::mch::Digit* digitsBuffer = NULL;
  digitsBuffer = (o2::mch::Digit*)input.payload;
  size_t ndigits = ((size_t)header->payloadSize / sizeof(o2::mch::Digit));

  if (mPrintLevel >= 1)
    std::cout << "There are " << ndigits << " digits in the payload" << std::endl;

  o2::mch::Digit* ptr = (o2::mch::Digit*)digitsBuffer;
  for (size_t di = 0; di < ndigits; di++) {
    digits.push_back(*ptr);
    ptr += 1;
  }

  for (auto& d : digits) {
    if (mPrintLevel >= 1) {
      std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
          d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTime().sampaTime);
      std::cout << std::endl;
    }
    plotDigit(d);
  }
}


void PhysicsTaskDigits::monitorData(o2::framework::ProcessingContext& ctx)
{
#ifdef QC_MCH_SAVE_TEMP_ROOTFILE_
  if ((count % 100) == 0) {

    TFile f("/tmp/qc.root", "RECREATE");
    for (int i = 0; i < 3 * 24; i++) {
      //mHistogramNhits[i]->Write();
      //mHistogramADCamplitude[i]->Write();
    }
    //std::cout<<"mHistogramADCamplitudeDE.size() = "<<mHistogramADCamplitudeDE.size()<<"  DEs.size()="<<DEs.size()<<std::endl;
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      //std::cout<<"  de="<<de<<std::endl;
      {
        auto h = mHistogramADCamplitudeDE.find(de);
        if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsDE.find(de);
        if ((h != mHistogramNhitsDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsHighAmplDE.find(de);
        if ((h != mHistogramNhitsHighAmplDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
    }

    f.ls();
    f.Close();
    if (mPrintLevel == 0) {
      printf("count: %d\n", count);
    }
  }
  if (mPrintLevel >= 1) {
    printf("count: %d\n", count);
  }
  count += 1;
#endif

  QcInfoLogger::GetInstance() << "monitorData" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(flog, "\n================\nmonitorData\n================\n");
  monitorDataReadout(ctx);
  for (auto&& input : ctx.inputs()) {
    if (mPrintLevel >= 1) {
      QcInfoLogger::GetInstance() << "run PhysicsTaskDigits: input " << input.spec->binding << AliceO2::InfoLogger::InfoLogger::endm;
    }
    if (input.spec->binding == "digits") {
      monitorDataDigits(input);
    }
  }
//  monitorDataReadout(ctx);
}


void PhysicsTaskDigits::plotDigit(const o2::mch::Digit& digit)
{
  int ADC = digit.getADC();
  int de = digit.getDetID();
  int padid = digit.getPadID();

  //fprintf(stdout, "digit[%d]: ADC=%d, DetId=%d, PadId=%d\n",
  //        i, ADC, de, padid);
  if (ADC < 0 || de < 0 || padid < 0) {
    return;
  }

  try {
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(de);

    double padX = segment.padPositionX(padid);
    double padY = segment.padPositionY(padid);
    float padSizeX = segment.padSizeX(padid);
    float padSizeY = segment.padSizeY(padid);
    int cathode = segment.isBendingPad(padid) ? 0 : 1;


    if (mPrintLevel >= 1)
      fprintf(flog, "de=%d pad=%d x=%f y=%f\n", de, padid, padX, padY);
    //if(pad.fX>=32 && pad.fX<=34 && pad.fY>=1.1 && pad.fY<=1.4)
    //  fprintf(flog, "mapping: link_id=%d ds_addr=%d chan_addr=%d  ==>  de=%d x=%f y=%f A=%d\n",
    //    hit.link_id, hit.ds_addr, hit.chan_addr, pad.fDE, pad.fX, pad.fY, hit.csum);

    auto h = mHistogramADCamplitudeDE.find(de);
    if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
      h->second->Fill(ADC);
    }

    if (cathode == 0 && ADC > 0) {
      auto h2 = mHistogramNhitsDE.find(de);
      if ((h2 != mHistogramNhitsDE.end()) && (h2->second != NULL)) {
        int binx_min = h2->second->GetXaxis()->FindBin(padX - padSizeX / 2 + 0.1);
        int binx_max = h2->second->GetXaxis()->FindBin(padX + padSizeX / 2 - 0.1);
        int biny_min = h2->second->GetYaxis()->FindBin(padY - padSizeY / 2 + 0.1);
        int biny_max = h2->second->GetYaxis()->FindBin(padY + padSizeY / 2 - 0.1);
        for (int by = biny_min; by <= biny_max; by++) {
          float y = h2->second->GetYaxis()->GetBinCenter(by);
          for (int bx = binx_min; bx <= binx_max; bx++) {
            float x = h2->second->GetXaxis()->GetBinCenter(bx);
            h2->second->Fill(x, y);
          }
        }
      }
    }
    if (cathode == 0 && ADC > 500) {
      auto h2 = mHistogramNhitsHighAmplDE.find(de);
      if ((h2 != mHistogramNhitsHighAmplDE.end()) && (h2->second != NULL)) {
        int binx_min = h2->second->GetXaxis()->FindBin(padX - padSizeX / 2 + 0.1);
        int binx_max = h2->second->GetXaxis()->FindBin(padX + padSizeX / 2 - 0.1);
        int biny_min = h2->second->GetYaxis()->FindBin(padY - padSizeY / 2 + 0.1);
        int biny_max = h2->second->GetYaxis()->FindBin(padY + padSizeY / 2 - 0.1);
        for (int by = biny_min; by <= biny_max; by++) {
          float y = h2->second->GetYaxis()->GetBinCenter(by);
          for (int bx = binx_min; bx <= binx_max; bx++) {
            float x = h2->second->GetXaxis()->GetBinCenter(bx);
            h2->second->Fill(x, y);
          }
        }
      }
    }
  } catch (const std::exception& e) {
    QcInfoLogger::GetInstance() << "[MCH] Detection Element " << de << " not found in mapping." << AliceO2::InfoLogger::InfoLogger::endm;
    return;
  }
}


void PhysicsTaskDigits::endOfCycle()
{
  QcInfoLogger::GetInstance() << "endOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;

  for(int de = 100; de <= 1030; de++) {
    for(int i = 0; i < 1; i++) {
        auto ih = mHistogramOccupancyXY[i+1].find(de);
        ih = mHistogramOccupancyXY[i].find(de);
        if (ih == mHistogramOccupancyXY[i].end()) {
          continue;
        }
    }
  }

    mHistogramOccupancy[0]->add(mHistogramNhitsDE, mHistogramNhitsDE);
//    mHistogramOccupancy[1]->add(mHistogramPreclustersXY[1], mHistogramPreclustersXY[2]);
//    mHistogramOccupancy[2]->add(mHistogramPreclustersXY[3], mHistogramPreclustersXY[3]);

#ifdef QC_MCH_SAVE_TEMP_ROOTFILE
    TFile f("/tmp/qc.root", "RECREATE");
    for (int i = 0; i < 3 * 24; i++) {
      //mHistogramNhits[i]->Write();
      //mHistogramADCamplitude[i]->Write();
    }
    //std::cout<<"mHistogramADCamplitudeDE.size() = "<<mHistogramADCamplitudeDE.size()<<"  DEs.size()="<<DEs.size()<<std::endl;
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      //std::cout<<"  de="<<de<<std::endl;
      {
        auto h = mHistogramADCamplitudeDE.find(de);
        if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsDE.find(de);
        if ((h != mHistogramNhitsDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsHighAmplDE.find(de);
        if ((h != mHistogramNhitsHighAmplDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
    }
    
    mHistogramOccupancy[0]->Write();
//    mHistogramOccupancy[1]->Write();
//    mHistogramOccupancy[2]->Write();

    //f.ls();
    f.Close();
#endif
}

void PhysicsTaskDigits::endOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "endOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::reset()
{
  // clean all the monitor objects here

  QcInfoLogger::GetInstance() << "Reseting the histogram" << AliceO2::InfoLogger::InfoLogger::endm;
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

