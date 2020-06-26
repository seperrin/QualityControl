///
/// \file   PhysicsTask.h
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PHYSICSTASK_H
#define QC_MODULE_MUONCHAMBERS_PHYSICSTASK_H

#include <TRandom3.h>
#include <vector>

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
class PhysicsTask /*final*/ : public TaskInterface // todo add back the "final" when doxygen is fixed
{
 public:
  /// \brief Constructor
  PhysicsTask();
  /// Destructor
  ~PhysicsTask() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorDataReadout(o2::framework::ProcessingContext& ctx);
  void monitorDataDigits(o2::framework::ProcessingContext& ctx);
  void monitorDataPreclusters(o2::framework::ProcessingContext& ctx);
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

  ssize_t getNumberOfDigits();
  void storeDigits(void* bufferPtr);

  void plotDigit(const o2::mch::Digit& digit);
  bool plotPrecluster(const o2::mch::PreCluster& preCluster, gsl::span<const o2::mch::Digit> digits);
  void checkPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits);
  void printPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits);

 private:
  int count;
  Decoder mDecoder;
  uint64_t nhits[24][40][64];
  uint32_t norbits[24];
  uint32_t firstorbitseen[24];

  std::vector<std::unique_ptr<mch::Digit>> digits;
  mch::Digit* digitsBuffer;
  int nDigits;

  TH2F* mHistogramNhits[72];
  TH1F* mHistogramADCamplitude[72];
  std::vector<int> DEs;
  std::map<int, TH1F*> mHistogramADCamplitudeDE;
  std::map<int, TH2F*> mHistogramNhitsDE;
  std::map<int, TH2F*> mHistogramNorbitsDE;
  std::map<int, TH2F*> mHistogramNhitsHighAmplDE;

  std::map<int, TH1F*> mHistogramClchgDE;
  std::map<int, TH1F*> mHistogramClsizeDE;

  std::map<int, TH2F*> mHistogramPreclustersXY[4];
  std::map<int, TH2F*> mHistogramPseudoeffXY[3];
  std::map<int, TH2F*> mHistogramOccupancyXY[3];
  TRandom3 rnd;

  GlobalHistogram* mHistogramPseudoeff[3];
  GlobalHistogram* mHistogramOccupancy[1];
  GlobalHistogram* mHistogramOrbits[1];

  int mPrintLevel;
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_MUONCHAMBERS_PHYSICSDATAPROCESSOR_H
