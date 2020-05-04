///
/// \file   PhysicsTask.cxx
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
#include "MCHBase/PreCluster.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingFactory/CreateSegmentation.h"
//#include "MCHPreClustering/PreClusterFinder.h"
#include "QualityControl/QcInfoLogger.h"
#include "MCH/PhysicsTask.h"

using namespace std;

static int gPrintLevel;

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
PhysicsTask::PhysicsTask() : TaskInterface(), count(1) {}

PhysicsTask::~PhysicsTask() { fclose(flog); }

void PhysicsTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  QcInfoLogger::GetInstance() << "initialize PhysicsTask" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(stdout, "initialize PhysicsTask\n");

  mDecoder.initialize();

  uint32_t dsid;
  std::vector<int> DEs;
  for (int cruid = 0; cruid < 3; cruid++) {
    QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE CRUID " << cruid << AliceO2::InfoLogger::InfoLogger::endm;
    for (int linkid = 0; linkid < 24; linkid++) {
      QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE LINKID " << linkid << AliceO2::InfoLogger::InfoLogger::endm;

      {
        int index = 24 * cruid + linkid;
        mHistogramNhits[index] = new TH2F(TString::Format("QcMuonChambers_NHits_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - Number of hits (CRU link %02d)", index), 40, 0, 40, 64, 0, 64);
        //mHistogramPedestals->SetDrawOption("col");
        getObjectsManager()->startPublishing(mHistogramNhits[index]);

        mHistogramADCamplitude[index] = new TH1F(TString::Format("QcMuonChambers_ADC_Amplitude_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - ADC amplitude (CRU link %02d)", index), 5000, 0, 5000);
        //mHistogramPedestals->SetDrawOption("col");
        getObjectsManager()->startPublishing(mHistogramADCamplitude[index]);
      }

      int32_t link_id = mDecoder.getMapCRU(cruid, linkid);
      if (link_id == -1)
        continue;
      for (int ds_addr = 0; ds_addr < 40; ds_addr++) {
        QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE DS_ADDR " << ds_addr << AliceO2::InfoLogger::InfoLogger::endm;
        uint32_t de = mDecoder.getMapFEC(link_id, ds_addr, de, dsid);
        QcInfoLogger::GetInstance() << "C'EST LA LIGNE APRÈS LE GETMAPFEC, DE " << de << AliceO2::InfoLogger::InfoLogger::endm;

        if (!(std::find(DEs.begin(), DEs.end(), de) != DEs.end())) {
          DEs.push_back(de);

          TH1F* h = new TH1F(TString::Format("QcMuonChambers_ADCamplitude_DE%03d", de),
              TString::Format("QcMuonChambers - ADC amplitude (DE%03d)", de), 5000, 0, 5000);
          mHistogramADCamplitudeDE.insert(make_pair(de, h));
          getObjectsManager()->startPublishing(h);

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
          getObjectsManager()->startPublishing(h2);
        }
      }
    }
  }

  mPrintLevel = 0;

  flog = stdout; //fopen("/root/qc.log", "w");
  fprintf(stdout, "PhysicsTask initialization finished\n");
}

void PhysicsTask::startOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "startOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::startOfCycle()
{
  QcInfoLogger::GetInstance() << "startOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::monitorDataReadout(o2::framework::ProcessingContext& ctx)
{
#ifdef QC_MCH_SAVE_TEMP_ROOTFILE
  if ((count % 1) == 0) {

    TFile f("/tmp/qc.root", "RECREATE");
    for (int i = 0; i < 3 * 24; i++) {
      mHistogramNhits[i]->Write();
      mHistogramADCamplitude[i]->Write();
    }
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      auto h = mHistogramADCamplitudeDE.find(de);
      if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
        h->second->Write();
      }
      auto h2 = mHistogramNhitsDE.find(de);
      if ((h2 != mHistogramNhitsDE.end()) && (h2->second != NULL)) {
        h2->second->Write();
      }

      f.ls();
      f.Close();
    }
  }
  printf("count: %d\n", count);
  count += 1;
#endif

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

    //std::cout<<"\n\npayloadSize: "<<payloadSize<<std::endl;
    //std::cout<<"raw:     "<<(void*)raw<<std::endl;
    //std::cout<<"payload: "<<(void*)payload<<std::endl;


    // Run the decoder on the CRU buffer
    mDecoder.processData((const char*)raw, (size_t)(payloadSize + sizeof(o2::header::RAWDataHeaderV4)));
  }

  std::vector<SampaHit>& hits = mDecoder.getHits();
  if (gPrintLevel >= 1)
    fprintf(flog, "hits.size()=%d\n", (int)hits.size());
  for (uint32_t i = 0; i < hits.size(); i++) {
    //continue;
    SampaHit& hit = hits[i];
    if (gPrintLevel >= 1)
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
  if (gPrintLevel >= 1)
    fprintf(flog, "digits.size()=%d\n", (int)digits.size());
  for (uint32_t i = 0; i < digits.size(); i++) {
    o2::mch::Digit& digit = digits[i];
    plotDigit(digit);
  }
}

void PhysicsTask::monitorDataDigits(const o2::framework::DataRef& input)
{
#ifdef QC_MCH_SAVE_TEMP_ROOTFILE
  if ((count % 1) == 0) {

    TFile f("/tmp/qc.root", "RECREATE");
    for (int i = 0; i < 3 * 24; i++) {
      mHistogramNhits[i]->Write();
      mHistogramADCamplitude[i]->Write();
    }
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      auto h = mHistogramADCamplitudeDE.find(de);
      if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
        h->second->Write();
      }
      auto h2 = mHistogramNhitsDE.find(de);
      if ((h2 != mHistogramNhitsDE.end()) && (h2->second != NULL)) {
        h2->second->Write();
      }

      f.ls();
      f.Close();
    }
  }
  printf("count: %d\n", count);
  count += 1;
#endif

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

  for (uint32_t i = 0; i < digits.size(); i++) {
    o2::mch::Digit& digit = digits[i];
    plotDigit(digit);
  }
}


void PhysicsTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  monitorDataReadout(ctx);
  for (auto&& input : ctx.inputs()) {
    QcInfoLogger::GetInstance() << "run PedestalsTask: input " << input.spec->binding << AliceO2::InfoLogger::InfoLogger::endm;
    if (input.spec->binding == "digits")
      monitorDataDigits(input);
  }
}


void PhysicsTask::plotDigit(const o2::mch::Digit& digit)
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


    if (gPrintLevel >= 1)
      fprintf(flog, "de=%d pad=%d x=%f y=%f\n", de, padid, padX, padY);
    //if(pad.fX>=32 && pad.fX<=34 && pad.fY>=1.1 && pad.fY<=1.4)
    //  fprintf(flog, "mapping: link_id=%d ds_addr=%d chan_addr=%d  ==>  de=%d x=%f y=%f A=%d\n",
    //    hit.link_id, hit.ds_addr, hit.chan_addr, pad.fDE, pad.fX, pad.fY, hit.csum);

    auto h = mHistogramADCamplitudeDE.find(de);
    if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
      h->second->Fill(ADC);
    }

    if (ADC > 0) {
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
    if (ADC > 500) {
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


void PhysicsTask::endOfCycle()
{
  QcInfoLogger::GetInstance() << "endOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::endOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "endOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::reset()
{
  // clean all the monitor objects here

  QcInfoLogger::GetInstance() << "Reseting the histogram" << AliceO2::InfoLogger::InfoLogger::endm;
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
