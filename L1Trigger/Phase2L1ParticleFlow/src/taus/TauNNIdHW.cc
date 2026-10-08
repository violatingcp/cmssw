#include <cstdint>
#include <algorithm>
#include <array>
#include <utility>
#include "L1Trigger/Phase2L1ParticleFlow/interface/taus/TauNNIdHW.h"

// pi/720 rad per integer unit of glbeta_t/glbphi_t, as the FIRMWARE holds it.
// tau_nn.cpp scales deta/dphi with an ap_fixed<18,1> copy of this constant, not
// with a float, and ap_fixed<18,1> only carries multiples of 2^-17. Scaling by
// the float l1ct::Scales::ETAPHI_LSB here and by the fixed-point value there
// put the two on scales 0.17% apart, which is 0.18 LSB of the ap_fixed<10,4>
// the network is fed - enough to push ~9% of the integers across a truncation
// boundary and show up as the 1-LSB deta/dphi disagreements in 07.
// Written as the explicit ratio because `= M_PI / 720` would be TRUNCATED into
// those 17 bits (571/2^17, -0.159%) rather than rounded (572/2^17, +0.016%);
// this must track tau_nn.cpp::ETAPHI_LSB_FX exactly, whichever value it holds.
// The 10 constituents actually handed to the network, after compute() sorts by
// pt and truncates. Off unless TAU_DUMP_NNINPUTS names a file. Species is
// printed as a name rather than a raw id because the two sides encode the
// enum differently (l1t::PFCandidate here, l1ct 3-bit ids in the firmware).
namespace {
  FILE *nnDumpFile() {
    static FILE *f = [] {
      const char *path = getenv("TAU_DUMP_NNINPUTS");
      return path ? fopen(path, "w") : nullptr;
    }();
    return f;
  }
  const char *speciesName(int pfid) {
    switch (pfid) {
      case l1t::PFCandidate::Photon: return "gamma";
      case l1t::PFCandidate::Electron: return "e";
      case l1t::PFCandidate::Muon: return "mu";
      case l1t::PFCandidate::NeutralHadron: return "h0";
      case l1t::PFCandidate::ChargedHadron: return "h+-";
      default: return "?";
    }
  }
}  // namespace

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

  // Sort the candidates on the MERGED KEY: pt in the high bits, then the low 3
  // bits of deta and dphi. This mirrors makeSortKey() in the firmware's
  // tau_types.h bit for bit - keep the two in step.
  //
  // pt alone is quantised to 0.25 GeV, so ties are common, and the two
  // implementations broke them differently: the firmware kept whichever particle
  // the deregionizer happened to present first (link/slot order), which nothing
  // here reproduces. Appending position bits makes the key almost always unique,
  // so both sides pick the same 10 particles without the emulator having to model
  // the hardware's traversal order.
  //
  // Three details make the two sides agree exactly:
  //   - deta/dphi are seed-minus-particle in BOTH implementations, so there is no
  //     sign convention to track;
  //   - the firmware truncates them to ap_int<9>, which cannot change the low 3
  //     bits of the same integer;
  //   - the firmware wraps dphi by +-INTPHI_TWOPI = 1440, a multiple of 8, so the
  //     low 3 bits are the same wrapped or not and this does not repeat the wrap.
  auto sortKey = [&iSeed](const l1t::PFCandidate &p) -> uint32_t {
    uint32_t ptbits = pt_t(p.pt()).range().to_uint();
    int deta = l1ct::Scales::makeGlbEta(iSeed.eta()).to_int() - l1ct::Scales::makeGlbEta(p.eta()).to_int();
    int dphi = l1ct::Scales::makeGlbPhi(iSeed.phi()).to_int() - l1ct::Scales::makeGlbPhi(p.phi()).to_int();
    // TIEBREAK_BITS in the firmware's tau_types.h - keep the two equal.
    constexpr int kTieBits = 3;
    constexpr uint32_t kMask = (1u << kTieBits) - 1u;
    return (ptbits << (2 * kTieBits)) | ((uint32_t(deta) & kMask) << kTieBits) | (uint32_t(dphi) & kMask);
  };
  // stable_sort keeps the residual collisions (equal pt and deta == dphi == 0
  // mod 8) deterministic rather than unspecified.
  std::stable_sort(
      iParts.begin(), iParts.end(), [&sortKey](const l1t::PFCandidate &i, const l1t::PFCandidate &j) {
        return sortKey(i) > sortKey(j);
      });

  if (FILE *df = nnDumpFile())
    fprintf(df, "NNSEED eta=%+6d phi=%+6d seedpt=%8.3f nparts=%u\n",
            l1ct::Scales::makeGlbEta(iSeed.eta()).to_int(),
            l1ct::Scales::makeGlbPhi(iSeed.phi()).to_int(),
            iSeed.pt(), unsigned(iParts.size()));

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

    if (FILE *df = nnDumpFile())
      fprintf(df, "  NNC k=%u pt=%8.3f deta=%+5d dphi=%+5d id=%s\n", i0,
              double(fPt_.get()[i0]), (lSeedEta - lPartEta).to_int(), lDPhiI,
              speciesName(int(iParts[i0].id())));
  }
  if (FILE *df = nnDumpFile()) { fprintf(df, "\n"); fflush(df); }

  // Set the inputs
  SetNNVectorVar();

  // Return the N outputs with the inputs
  return EvaluateNN();
}
