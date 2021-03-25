//
//    rfnoc-hls-neuralnet: Vivado HLS code for neural-net building blocks
//
//    Copyright (C) 2017 EJ Kreinar
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef TAUNNIDHW_H_
#define TAUNNIDHW_H_

#include <complex>
#include "ap_int.h"
#include "ap_fixed.h"
//#include "parameters.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/utils/parameters.h"

typedef ap_ufixed<16, 14> pt_t;
typedef ap_fixed<10, 4> etaphi_t;

//typedef ap_fixed<14, 6> pt_t;
//typedef ap_fixed<14, 6> etaphi_t;

//typedef float pt_t;
//typedef float input_t;
//typedef float result_t;
//typedef float etaphi_t;

class TauNNIdHW {
public:
  TauNNIdHW();
  ~TauNNIdHW();

  void initialize(const std::string &iName, const std::string &iWeightFile, int iNParticles);
  void SetNNVectorVar();
  result_t EvaluateNN();
  result_t compute(l1t::PFCandidate &iSeed, std::vector<edm::Ptr<l1t::PFCandidate>> &iParts);

  std::string fInput_;
  unsigned fNParticles_;
  unique_ptr<pt_t[]> fPt_;
  unique_ptr<etaphi_t[]> fEta_;
  unique_ptr<etaphi_t[]> fPhi_;
  unique_ptr<id_t[]> fId_;

private:
  std::vector<input_t> NNvectorVar_;
};


#endif
