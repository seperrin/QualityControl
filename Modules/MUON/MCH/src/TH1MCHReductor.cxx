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
/// \file   TH1MCHReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH1.h>
#include "MCH/TH1MCHReductor.h"
#include <iostream>
#include <string>
#include <regex>

namespace o2::quality_control_modules::muonchambers
{

void* TH1MCHReductor::getBranchAddress()
{
  return &mStats;
}

const char* TH1MCHReductor::getBranchLeafList()
{
  return "occ819/D:occ8O:entries";
}

void TH1MCHReductor::update(TObject* obj)
{
 auto histo = dynamic_cast<TH1*>(obj);
  if (histo) {
      double mean = 0;
    mStats.entries = histo->GetEntries();
    mStats.occs.array[0] = histo->GetBinContent(819+1);
      
      for(int i=807; i<=819; i++){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.occs.array[1] = mean;
    
      std::cout << "MeanOccupancy DE819 obtained from reductor " << mStats.occs.array[0] << std::endl;
      std::cout << "MeanOccupancy Ch8O obtained from reductor " << mStats.occs.array[1] << std::endl;
      
  }
}

} // namespace o2::quality_control_modules::muonchambers


