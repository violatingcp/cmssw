#ifndef L1Trigger_Phase2L1ParticleFlow_TAUNNIDHW_H_
#define L1Trigger_Phase2L1ParticleFlow_TAUNNIDHW_H_

#include <cstdio>
#include <complex>
#include <memory>

#include "ap_fixed.h"
#include "ap_int.h"
#include "hls4ml/emulator.h"

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

typedef ap_ufixed<16, 14> pt_t;
typedef ap_fixed<10, 4> etaphi_t;

namespace L1TauEmu {
  // Input/output precision of the NNPuppiTauModel hls4ml network, matching the
  // model's defines.h (kept here since the weights/layer code now live in the
  // externally-built NNPuppiTauModel package, not in-tree). Named distinctively
  // (not input_t/result_t) and namespaced to avoid clashing with the generic
  // typedefs other hls4ml model headers define at global scope.
  typedef ap_fixed<16, 10> tauinput_t;
  typedef ap_fixed<16, 6> tauresult_t;

  // Data types and constants used in the FPGA and FPGA-optimized functions
  // etaphi_base maps physical eta/phi onto integer units.
  //
  // 720/pi, so one unit is exactly pi/720 rad - the same step the firmware and
  // l1ct::Scales use. It was 100/64, giving a 0.01 rad step, which made the
  // emulator's cone coarser than the firmware's and left particles within
  // ~0.005 rad of the cone edge on different sides in the two implementations.
  //
  // It also makes the file self-consistent: deltaPhi() already compares against
  // l1ct::Scales::INTPHI_PI (= 720) and INTPHI_TWOPI (= 1440), which are in
  // pi/720 units. Under the old base pi was 4.909, so those comparisons were
  // against a value beyond the reachable range and the phi wrap never fired.
  // Under this base they are exactly pi and 2pi.
  //
  // The tau cone follows: 0.1^2 * etaphi_base^2 = 525.25, which is the 525 the
  // firmware uses.
  static constexpr float etaphi_base = 720. / 3.14159265358979323846;
  typedef ap_ufixed<14, 12, AP_TRN, AP_SAT> pt_t;  // 1 unit = 0.25 GeV;
  // Integer-valued: one unit IS the quantum, so no fractional bits. eta spans
  // +-5.11 -> +-1171 units, phi +-pi -> +-720.
  typedef ap_fixed<12, 12> etaphi_t;               // 1 unit = pi/720 rad
  typedef ap_fixed<13, 13> detaphi_t;              // difference between etas or phis

  // A cone radius squared, and the squared distance it is compared against.
  // Worst case is (2*1171)^2 + 1440^2 ~ 7.6e6, so 25 bits.
  typedef ap_fixed<26, 26> detaphi2_t;             // type for detaphi_t squared
  typedef ap_fixed<22, 16> pt_etaphi_t;            // type for product of pt with deta or phi
  typedef ap_int<8> dxy_t;
  typedef ap_int<10> z0_t;
  typedef ap_uint<5> count_t;  // type for multiplicity
  typedef ap_uint<5> id_t;     // type for multiplicity

  // constants for the axis update
  typedef ap_ufixed<18, -2> inv_pt_t;
  static constexpr int N_table_inv_pt = 1024;
  static const detaphi_t TWOPI = 3.14159 * 2. * etaphi_base;
  static const detaphi_t PI = 3.14159 * etaphi_base;
  static const detaphi_t HALFPI = 3.14159 / 2 * etaphi_base;
  static const detaphi_t RCONE = 0.4 * 100 / 128;
  static const detaphi_t R2CONE = RCONE * RCONE;
  //
  static const etaphi_t FIDUCIAL_ETA_PHI = 5.11 * etaphi_base;

  constexpr int ceillog2(int x) { return (x <= 2) ? 1 : 1 + ceillog2((x + 1) / 2); }
  constexpr int floorlog2(int x) { return (x < 2) ? 0 : 1 + floorlog2(x / 2); }
  constexpr int pow2(int x) { return x == 0 ? 1 : 2 * pow2(x - 1); }

  template <class data_T, int N>
  inline float real_val_from_idx(unsigned i) {
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    data_T x(0);
    // The MSB of 1 is implicit in the table
    x[x.width - 1] = 1;
    // So we can use the next NB bits for real data
    x(x.width - 2, x.width - NB - 1) = i;
    return (float)x;
  }

  template <class data_T, int N>
  inline unsigned idx_from_real_val(data_T x) {
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    // Slice the top-1 NB bits of the value
    // the MSB of '1' is implicit, so only slice below that
    ap_uint<NB> y = x(x.width - 2, x.width - NB - 1);
    return (unsigned)y(NB - 1, 0);
  }

  template <class data_T, class table_T, int N>
  void init_invert_table(table_T table_out[N]) {
    // The template data_T is the data type used to address the table
    for (unsigned i = 0; i < N; i++) {
      float x = real_val_from_idx<data_T, N>(i);
      table_T inv_x = 1 / x;
      table_out[i] = inv_x;
    }
  }

  template <class in_t, class table_t, int N>
  table_t invert_with_shift(in_t in, bool debug = false) {
    table_t inv_table[N];
    init_invert_table<in_t, table_t, N>(inv_table);

    // find the first '1' in the denominator
    int msb = 0;
    for (int b = 0; b < in.width; b++) {
      if (in[b])
        msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width - msb - 1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t, N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width - msb - 1);

    return out;
  }

  inline detaphi_t deltaPhi(l1t::PFCandidate a, l1t::PFCandidate b) {
    // Reconstruct the integers the HARDWARE actually holds. The firmware never
    // sees a float: layer-1 quantised phi with l1ct::Scales::makeGlbPhi, which
    // ROUNDS, and those integers are what the deregionizer hands the seeded
    // cone. etaphi_t(phi * etaphi_base) instead TRUNCATES (ap_fixed<12,12> is
    // AP_TRN), landing on a different unit often enough to flip the cone test
    // for ~2.6% of pairs near the boundary - which is how particles the
    // firmware never had ended up inside CMSSW's cone.
    detaphi_t dphi = detaphi_t(l1ct::Scales::makeGlbPhi(a.phi()).to_int() -
                               l1ct::Scales::makeGlbPhi(b.phi()).to_int());
    // phi wrap
    detaphi_t dphi0 =
        dphi > detaphi_t(l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI - dphi) : detaphi_t(dphi);
    detaphi_t dphi1 =
        dphi < detaphi_t(-l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI + dphi) : detaphi_t(dphi);
    //dphi > PI ? detaphi_t(TWOPI - dphi) : detaphi_t(dphi);
    //dphi < -PI ? detaphi_t(TWOPI + dphi) : detaphi_t(dphi);
    detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
    return dphiw;
  }

  inline bool inCone(l1t::PFCandidate seed, l1t::PFCandidate part, detaphi2_t cone2) {
    // makeGlbEta, not a truncating cast - see deltaPhi above.
    detaphi_t deta = detaphi_t(l1ct::Scales::makeGlbEta(seed.eta()).to_int() -
                               l1ct::Scales::makeGlbEta(part.eta()).to_int());
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = (deta * deta + dphi * dphi) < cone2;
    return ret;
  }

};  // namespace L1TauEmu

// Tau NN returns two values
struct Tau_NN_Result {
  L1TauEmu::tauresult_t nn_pt_correction;
  L1TauEmu::tauresult_t nn_id;
};

class TauNNIdHW {
public:
  TauNNIdHW(const std::shared_ptr<hls4mlEmulator::Model> model);
  ~TauNNIdHW() = default;

  void initialize(const std::string &iName, int iNParticles);
  void SetNNVectorVar();
  L1TauEmu::tauinput_t *NNVectorVar() { return NNvectorVar_.data(); }
  Tau_NN_Result EvaluateNN();
  Tau_NN_Result compute(const l1t::PFCandidate &iSeed, std::vector<l1t::PFCandidate> &iParts);
  //void print();

  std::string fInput_;
  unsigned fNParticles_;
  unique_ptr<pt_t[]> fPt_;
  unique_ptr<etaphi_t[]> fEta_;
  unique_ptr<etaphi_t[]> fPhi_;
  unique_ptr<id_t[]> fId_;
  //FILE *file_;

private:
  std::vector<L1TauEmu::tauinput_t> NNvectorVar_;
  std::shared_ptr<hls4mlEmulator::Model> modelRef_;
};

#endif
