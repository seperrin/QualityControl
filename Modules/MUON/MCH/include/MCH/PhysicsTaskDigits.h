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
  void monitorDataDigits(const o2::framework::DataRef& input);
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

  std::vector<std::unique_ptr<mch::Digit>> digits;
  mch::Digit* digitsBuffer;
  int nDigits;

  TH2F* mHistogramNhits[72];
  TH1F* mHistogramADCamplitude[72];
  std::vector<int> DEs;
  std::map<int, TH1F*> mHistogramADCamplitudeDE;
  std::map<int, TH2F*> mHistogramNhitsDE;
  std::map<int, TH2F*> mHistogramNhitsHighAmplDE;
    
  std::map<int, TH2F*> mHistogramOccupancyXY[3];
  TRandom3 rnd;

  GlobalHistogram* mHistogramOccupancy[1];

  int mPrintLevel;
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_MUONCHAMBERS_PhysicsTaskDIGITS_H

