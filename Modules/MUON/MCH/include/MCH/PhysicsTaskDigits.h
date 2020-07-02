///
/// \file   PhysicsTaskDigits.h
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PHYSICSTASKDIGITS_H
#define QC_MODULE_MUONCHAMBERS_PHYSICSTASKDIGITS_H

#include <TRandom3.h>

#include "QualityControl/TaskInterface.h"
#include "MCH/Mapping.h"
#include "MCH/Decoding.h"
#include "MCH/GlobalHistogram.h"
#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"

class TH1F;
class TH2F;

using namespace o2::quality_control::core;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

/// \brief Quality Control Task for the analysis of MCH physics data
/// \author Andrea Ferrero
/// \author Sebastien Perrin
class PhysicsTaskDigits /*final*/ : public TaskInterface // todo add back the "final" when doxygen is fixed
{
 public:
  /// \brief Constructor
  PhysicsTaskDigits();
  /// Destructor
  ~PhysicsTaskDigits() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorDataReadout(o2::framework::ProcessingContext& ctx);
  void monitorDataDigits(o2::framework::ProcessingContext& ctx);
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

  ssize_t getNumberOfDigits();
  void storeDigits(void* bufferPtr);

  void plotDigit(const o2::mch::Digit& digit);

 private:
  int count;
  Decoder mDecoder;
  uint64_t nhits[24][40][64];
  uint32_t norbits[1030];
  uint32_t firstorbitseen[1030];
    //Assumed that if there is an orbit with the B side of a DE, there is also the NB side. So Norbits B and NB are the same for the moment. 

  std::vector<std::unique_ptr<mch::Digit>> digits;
  mch::Digit* digitsBuffer;
  int nDigits;

  TH2F* mHistogramNhits[72];
  TH1F* mHistogramADCamplitude[72];
  std::vector<int> DEs;
  std::map<int, TH1F*> mHistogramADCamplitudeDE;
  std::map<int, TH2F*> mHistogramNhitsDE[4]; // 1 B, 1 NB, 2 B+NB
  std::map<int, TH2F*> mHistogramNorbitsDE[2]; //1 B dt 1 NB
  std::map<int, TH2F*> mHistogramNhitsHighAmplDE[2]; // 1 B, 1 NB
    
  std::map<int, TH2F*> mHistogramOccupancyXY[3]; //Inutile
  TRandom3 rnd;

  GlobalHistogram* mHistogramOccupancy[3]; // 1 B, 1 NB, 1 B+NB
  GlobalHistogram* mHistogramOrbits[2];// 1 B, 1 NB, assume B pour B+NB
    //mHistogramOrbits[1 and 2] are the same for the moment (assumed B and NB are always taken together in the same orbit. 

  int mPrintLevel;
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_MUONCHAMBERS_PhysicsTaskDIGITS_H
