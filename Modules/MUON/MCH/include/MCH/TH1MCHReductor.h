// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   TH1MCHReductor.h
/// \author Piotr Konopka, Sebastien Perrin
///
#ifndef QUALITYCONTROL_TH1MCHREDUCTOR_H
#define QUALITYCONTROL_TH1MCHREDUCTOR_H

#include "QualityControl/Reductor.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief A Reductor which obtains the most popular characteristics of TH1.
///
/// A Reductor which obtains the most popular characteristics of TH1.

class TH1MCHReductor : public quality_control::postprocessing::Reductor
{
 public:
  TH1MCHReductor() = default;
  ~TH1MCHReductor() = default;

  void* getBranchAddress() override;
  const char* getBranchLeafList() override;
  void update(TObject* obj) override;

 private:
 struct {
   union {
     struct {
         Double_t occ819;
         Double_t occ8O;
     } named;
     Double_t array[2];
   } occs;
   Double_t entries;
 } mStats;
};

} // namespace o2::quality_control_modules::muonchambers

#endif //QUALITYCONTROL_TH1mchREDUCTOR_H


