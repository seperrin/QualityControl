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

  mPrintLevel = 2;

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
          mHistogramNhitsDE[0].insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);

          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d_B", de),
              TString::Format("QcMuonChambers - Number of hits (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsDE[1].insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);
            
          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d_NB", de),
              TString::Format("QcMuonChambers - Number of hits (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsDE[2].insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);
            
          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d_BNB", de),
              TString::Format("QcMuonChambers - Number of hits (DE%03d BNB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsDE[3].insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);
            
          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d_B", de),
              TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsHighAmplDE[0].insert(make_pair(de, h2));
          //getObjectsManager()->startPublishing(h2);
            
            h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d_NB", de),
                TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
            mHistogramNhitsHighAmplDE[1].insert(make_pair(de, h2));
            //getObjectsManager()->startPublishing(h2);
            h2 = new TH2F(TString::Format("QcMuonChambers_Norbits_DE%03d_B", de),
                TString::Format("QcMuonChambers - Number of orbits (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
            mHistogramNorbitsDE[0].insert(make_pair(de, h2));
            h2 = new TH2F(TString::Format("QcMuonChambers_Norbits_DE%03d_NB", de),
                TString::Format("QcMuonChambers - Number of orbits (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
            mHistogramNorbitsDE[1].insert(make_pair(de, h2));
        }
      }
    }
  }

  for(int de = 1; de <= 1030; de++) {
      norbits[de] = 0;
      firstorbitseen[de] = 0;
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

      hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_NB_XY_%03d", de),
           TString::Format("QcMuonChambers - Occupancy XY (DE%03d NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramOccupancyXY[1].insert(make_pair(de, hXY));
//
//      hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_BNB_XY_%03d", de),
//           TString::Format("QcMuonChambers - Occupancy XY (DE%03d B+NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
//      mHistogramOccupancyXY[2].insert(make_pair(de, hXY));
    }
  }
    
  mHistogramOccupancy[0] = new GlobalHistogram("QcMuonChambers_Occupancy_denB", "Occupancy B");
  mHistogramOccupancy[0]->init();
  mHistogramOccupancy[1] = new GlobalHistogram("QcMuonChambers_Occupancy_NB", "Occupancy NB");
  mHistogramOccupancy[1]->init();
  mHistogramOccupancy[2] = new GlobalHistogram("QcMuonChambers_Occupancy_BNB", "Occupancy - B+NB");
  mHistogramOccupancy[2]->init();
    
    mHistogramOrbits[0] = new GlobalHistogram("QcMuonChambers_Orbits_denB", "Orbits");
    mHistogramOrbits[0]->init();
    mHistogramOrbits[1] = new GlobalHistogram("QcMuonChambers_Orbits_NB", "Orbits");
    mHistogramOrbits[1]->init();
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


void PhysicsTaskDigits::monitorDataDigits(o2::framework::ProcessingContext& ctx)
{
    fprintf(flog, "\n================\nmonitorDataDigits\n================\n");
  // get the input preclusters and associated digits
  auto digits = ctx.inputs().get<gsl::span<o2::mch::Digit>>("digits");
  auto orbits = ctx.inputs().get<gsl::span<uint64_t>>("orbits");

  if (mPrintLevel >= 1) {
  std::cout<<"digits.size()="<<digits.size()<<std::endl;
  }
    
    for (auto& orb : orbits){ //Normalement une seule fois
        
        if (mPrintLevel >= 0) {
            std::cout << fmt::format(" ORBIT {}", orb);
        }
        
        for (auto& d : digits) {
            
            if(firstorbitseen[d.getDetID()] == 0){
                firstorbitseen[d.getDetID()] = orb;
                std::cout<< "First orbit of DE " << d.getDetID() << " is set to " << orb <<std::endl;
            }
            norbits[d.getDetID()] = (orb-firstorbitseen[d.getDetID()]+1);
            std::cout<< "Number of orbits of DE " << d.getDetID() << " is set to " << norbits[d.getDetID()] <<std::endl;
            
            if (mPrintLevel >= 0) {
              std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
                  d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTime().sampaTime);
              std::cout << std::endl;
            }
            
            plotDigit(d);
        }
   }
    
}

void PhysicsTaskDigits::monitorData(o2::framework::ProcessingContext& ctx)
{
  QcInfoLogger::GetInstance() << "monitorData" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(flog, "\n================\nmonitorData\n================\n");
  //monitorDataReadout(ctx);
  for (auto&& input : ctx.inputs()) {
    if (mPrintLevel >= 1) {
      QcInfoLogger::GetInstance() << "run PhysicsTaskDigits: input " << input.spec->binding << AliceO2::InfoLogger::InfoLogger::endm;
    }
    if (input.spec->binding == "digits") {
      monitorDataDigits(ctx);
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
      
    if (ADC > 0) {
      auto h2 = mHistogramNhitsDE[0].find(de);
      if ((h2 != mHistogramNhitsDE[0].end()) && (h2->second != NULL)) {
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
      
    if (cathode == 0 && ADC > 0) {
      auto h2 = mHistogramNhitsDE[1].find(de);
      if ((h2 != mHistogramNhitsDE[1].end()) && (h2->second != NULL)) {
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
    
    if (cathode == 1 && ADC > 0) {
      auto h2 = mHistogramNhitsDE[2].find(de);
      if ((h2 != mHistogramNhitsDE[2].end()) && (h2->second != NULL)) {
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
      
    if ((cathode == 0 || cathode == 1) && ADC > 0) {
      auto h2 = mHistogramNhitsDE[3].find(de);
      if ((h2 != mHistogramNhitsDE[3].end()) && (h2->second != NULL)) {
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
      
      if (cathode == 0 && ADC > 0) {
       auto h2 = mHistogramNorbitsDE[0].find(de);
        if ((h2 != mHistogramNorbitsDE[0].end()) && (h2->second != NULL)) {
          int NYbins = h2->second->GetYaxis()->GetNbins();
          int NXbins = h2->second->GetXaxis()->GetNbins();
          for (int by = 0; by <= NYbins; by++) {
            //float y = h2->second->GetYaxis()->GetBinCenter(by);
            for (int bx = 0; bx <= NXbins; bx++) {
              //float x = h2->second->GetXaxis()->GetBinCenter(bx);
              h2->second->SetBinContent(bx, by, norbits[de]);
            }
          }
        }
      }
      
      if (cathode == 1 && ADC > 0) {
       auto h2 = mHistogramNorbitsDE[1].find(de);
        if ((h2 != mHistogramNorbitsDE[1].end()) && (h2->second != NULL)) {
          int NYbins = h2->second->GetYaxis()->GetNbins();
          int NXbins = h2->second->GetXaxis()->GetNbins();
          for (int by = 0; by <= NYbins; by++) {
            //float y = h2->second->GetYaxis()->GetBinCenter(by);
            for (int bx = 0; bx <= NXbins; bx++) {
              //float x = h2->second->GetXaxis()->GetBinCenter(bx);
              h2->second->SetBinContent(bx, by, norbits[de]);
            }
          }
        }
      }

      
    if (cathode == 0 && ADC > 500) {
      auto h2 = mHistogramNhitsHighAmplDE[0].find(de);
      if ((h2 != mHistogramNhitsHighAmplDE[0].end()) && (h2->second != NULL)) {
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
      
    if (cathode == 1 && ADC > 500) {
         auto h2 = mHistogramNhitsHighAmplDE[1].find(de);
         if ((h2 != mHistogramNhitsHighAmplDE[1].end()) && (h2->second != NULL)) {
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
    for(int i = 0; i < 2; i++) {
        auto ih = mHistogramOccupancyXY[i+1].find(de);
        ih = mHistogramOccupancyXY[i].find(de);
        if (ih == mHistogramOccupancyXY[i].end()) {
          continue;
        }
    }
  }

     mHistogramOrbits[0]->add(mHistogramNorbitsDE[0], mHistogramNorbitsDE[0]);
     mHistogramOrbits[1]->add(mHistogramNorbitsDE[1], mHistogramNorbitsDE[1]);
    
    mHistogramOccupancy[0]->add(mHistogramNhitsDE[1], mHistogramNhitsDE[1]);
    mHistogramOccupancy[0]->Divide(mHistogramOrbits[0]);
    mHistogramOccupancy[1]->add(mHistogramNhitsDE[2], mHistogramNhitsDE[2]);
    mHistogramOccupancy[1]->Divide(mHistogramOrbits[1]);
    
    mHistogramOccupancy[2]->add(mHistogramNhitsDE[3], mHistogramNhitsDE[3]);
    mHistogramOccupancy[2]->Divide(mHistogramOrbits[0]);

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
      
        for(int i=0; i<4;i++){
          {
            auto h = mHistogramNhitsDE[i].find(de);
            if ((h != mHistogramNhitsDE[i].end()) && (h->second != NULL)) {
              h->second->Write();
            }
          }
        }
        
        for(int i=0; i<2;i++){
                 {
                   auto h = mHistogramNorbitsDE[i].find(de);
                   if ((h != mHistogramNorbitsDE[i].end()) && (h->second != NULL)) {
                     h->second->Write();
                   }
                 }
               }
        
        for(int i=0; i<2;i++){
          {
            auto h = mHistogramNhitsHighAmplDE[i].find(de);
            if ((h != mHistogramNhitsHighAmplDE[i].end()) && (h->second != NULL)) {
              h->second->Write();
            }
          }
        }
    }
    
    mHistogramOrbits[0]->Write();
    mHistogramOrbits[1]->Write();
    mHistogramOccupancy[0]->Write();
    mHistogramOccupancy[1]->Write();
    mHistogramOccupancy[2]->Write();

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
