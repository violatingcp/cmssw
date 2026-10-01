#include <algorithm>
#include <array>
#include <utility>
#include "L1Trigger/Phase2L1ParticleFlow/interface/taus/TauNNIdHW.h"

static const ap_fixed<18, 1> ETAPHI_LSB_FX = 571.0 / 131072.0;

// The firmware multiplies that constant by SortPart::deta, which is the seeded
// cone's detaphi_t = ap_int<9> (seededcone/firmware/data.h), NOT the wider
// ap_fixed<13,13> that shares the name in L1TauEmu. Mirror the narrow type so
// the out-of-range behaviour matches too, not just the in-cone arithmetic.
typedef ap_int<9> fw_detaphi_t;

TauNNIdHW::TauNNIdHW(const std::shared_ptr<hls4mlEmulator::Model> model) : modelRef_(model) { NNvectorVar_.clear(); }

void TauNNIdHW::initialize(const std::string &iInput, int iNParticles) {
  fNParticles_ = iNParticles;
  fPt_ = std::make_unique<pt_t[]>(fNParticles_);
  fEta_ = std::make_unique<etaphi_t[]>(fNParticles_);
  fPhi_ = std::make_unique<etaphi_t[]>(fNParticles_);
  fId_ = std::make_unique<id_t[]>(fNParticles_);
  fInput_ = iInput;
}

//Prepare the inputs for the Tau NN
void TauNNIdHW::SetNNVectorVar() {
  NNvectorVar_.clear();
  for (unsigned i0 = 0; i0 < fNParticles_; i0++) {
    L1TauEmu::tauinput_t pPt = L1TauEmu::tauinput_t(fPt_.get()[i0]);
    L1TauEmu::tauinput_t pEta = L1TauEmu::tauinput_t(fEta_.get()[i0]);
    L1TauEmu::tauinput_t pPhi = L1TauEmu::tauinput_t(fPhi_.get()[i0]);

    NNvectorVar_.push_back(pPt);
    NNvectorVar_.push_back(pEta);
    NNvectorVar_.push_back(pPhi);
    if (fPt_.get()[i0] == 0) {
      for (unsigned i1 = 0; i1 < 5; i1++)
        NNvectorVar_.push_back(0);
      continue;
    }
    NNvectorVar_.push_back(fId_.get()[i0] == l1t::PFCandidate::Photon);         // Photon
    NNvectorVar_.push_back(fId_.get()[i0] == l1t::PFCandidate::Electron);       // Electron
    NNvectorVar_.push_back(fId_.get()[i0] == l1t::PFCandidate::Muon);           // Muon
    NNvectorVar_.push_back(fId_.get()[i0] == l1t::PFCandidate::NeutralHadron);  // Neutral Had
    NNvectorVar_.push_back(fId_.get()[i0] == l1t::PFCandidate::ChargedHadron);  // Charged Had
  }
}

// Main architecture of the NN here: delegates to the externally-built
// NNPuppiTauModel package (loaded via hls4mlEmulator::Model) instead of
// in-tree weight arrays, avoiding ELF symbol clashes between hls4ml models
// loaded into the same process (cms-sw/cmssw#49632).
Tau_NN_Result TauNNIdHW::EvaluateNN() {
  constexpr unsigned kNInputs = 80;
  L1TauEmu::tauinput_t model_input[kNInputs];
  for (unsigned int i = 0; i < NNvectorVar_.size(); i++) {
    model_input[i] = L1TauEmu::tauinput_t(NNvectorVar_[i]);
  }

  typedef std::pair<std::array<L1TauEmu::tauresult_t, 1>, std::array<L1TauEmu::tauresult_t, 1>> pairtype;
  pairtype modelResult;
  modelRef_->prepare_input(model_input);
  modelRef_->predict();
  modelRef_->read_result(&modelResult);

  // Return both pT correction and the NN ID
  Tau_NN_Result nn_out;
  nn_out.nn_pt_correction = modelResult.first[0];
  nn_out.nn_id = modelResult.second[0];

  return nn_out;
}

/*
// Uncomment for debugging purposes
void TauNNIdHW::print() { 
  for (unsigned i0 = 0; i0 < fNParticles_; i0++) {
    tauinput_t pPt  = tauinput_t(fPt_.get()[i0]);
    tauinput_t pEta = tauinput_t(fEta_.get()[i0]);
    tauinput_t pPhi = tauinput_t(fPhi_.get()[i0]);
    tauinput_t pId  = tauinput_t(fId_.get()[i0]);
    fprintf(file_, " %08x", pPt.to_uint());
    fprintf(file_, " %08x", pEta.to_uint());
    fprintf(file_, " %08x", pPhi.to_uint());
    fprintf(file_, " %08x", pId.to_uint());
  }
  fprintf(file_, "\n");
}
*/

Tau_NN_Result TauNNIdHW::compute(const l1t::PFCandidate &iSeed, std::vector<l1t::PFCandidate> &iParts) {
  // Initialize the input vector
  for (unsigned i0 = 0; i0 < fNParticles_; i0++) {
    fPt_.get()[i0] = 0.;
    fEta_.get()[i0] = 0.;
    fPhi_.get()[i0] = 0.;
    fId_.get()[i0] = 0.;
  }

  // Sort the candidates by pT.
  // stable_sort, not sort: the comparator ranks on pt_t, quantised to 0.25 GeV,
  // so equal keys are common and std::sort leaves their order unspecified.
  // NOTE: this makes the tie-break DETERMINISTIC but not yet equal to the
  // firmware's, which is the order jet_loop emits constituents in (deregionizer
  // link/slot order). Aligning that ordering is the outstanding item; a
  // deta tie-break was measured in hardware and cost +16% LUT and +41 cycles
  // on tau_select, so the merge-the-bits route was not taken.
  std::stable_sort(iParts.begin(), iParts.end(), [](const l1t::PFCandidate &i, const l1t::PFCandidate &j) {
    return (pt_t(i.pt()) > pt_t(j.pt()));
  });


  // Compute the values w.r.t to the seeds
  for (unsigned int i0 = 0; i0 < iParts.size(); i0++) {
    if (i0 >= fNParticles_)
      break;

    fPt_.get()[i0] = pt_t(iParts[i0].pt());
    // Take the firmware's path for deta: it differences the pi/720 integer eta
    // values the Puppi objects already carry, then converts once. Differencing
    // the raw floats instead means the two implementations quantise at
    // different points and disagree by up to an LSB near the boundaries.
    // makeGlbEta is round(eta / (pi/720)) - the same quantisation that produced
    // hwEta upstream - so this reproduces the hardware value rather than
    // approximating it.
    l1ct::glbeta_t lSeedEta = l1ct::Scales::makeGlbEta(iSeed.eta());
    l1ct::glbeta_t lPartEta = l1ct::Scales::makeGlbEta(iParts[i0].eta());
    // Scale exactly as tau_nn.cpp does - l1ct::detaphi_t integer difference
    // times the ap_fixed constant, truncated once into ap_fixed<10,4> - rather
    // than going out to float and back.
    fEta_.get()[i0] = etaphi_t(fw_detaphi_t(lSeedEta - lPartEta) * ETAPHI_LSB_FX);
    // Same treatment as deta, and here it matters. This used to quantise EACH
    // absolute phi to etaphi_t (LSB 0.015625 rad, coarser than the pi/720 grid
    // the candidates sit on) and only then subtract, so up to two LSB of
    // rounding entered a quantity the firmware gets exactly. Difference the
    // integers first, wrap there, and convert once - which also makes the wrap
    // a whole turn (1440 = 2*pi) rather than the half turn it used to be.
    l1ct::glbphi_t lSeedPhi = l1ct::Scales::makeGlbPhi(iSeed.phi());
    l1ct::glbphi_t lPartPhi = l1ct::Scales::makeGlbPhi(iParts[i0].phi());
    int lDPhiI = lSeedPhi.to_int() - lPartPhi.to_int();
    if (lDPhiI > l1ct::Scales::INTPHI_PI)
      lDPhiI -= l1ct::Scales::INTPHI_TWOPI;
    if (lDPhiI < -l1ct::Scales::INTPHI_PI)
      lDPhiI += l1ct::Scales::INTPHI_TWOPI;

    fPhi_.get()[i0] = etaphi_t(fw_detaphi_t(lDPhiI) * ETAPHI_LSB_FX);
    fId_.get()[i0] = id_t(iParts[i0].id());

  }

  // Set the inputs
  SetNNVectorVar();

  // Return the N outputs with the inputs
  return EvaluateNN();
}
