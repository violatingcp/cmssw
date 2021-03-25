#include "L1Trigger/Phase2L1ParticleFlow/interface/L1NNTauProducer.hh"
#include <TLorentzVector.h>
#include <cmath>

#include "ap_int.h"
#include "ap_fixed.h"

L1NNTauProducer::L1NNTauProducer(const edm::ParameterSet& cfg)
    : fSeedPt_(cfg.getParameter<double>("seedpt")),
      fConeSize_(cfg.getParameter<double>("conesize")),
      fTauSize_(cfg.getParameter<double>("tausize")),
      fMaxTaus_(cfg.getParameter<unsigned>("maxtaus")),
      fNParticles_(cfg.getParameter<int>("nparticles")),
      fHW(cfg.getParameter<bool>("HW")),
      fDebug(cfg.getParameter<bool>("debug")),
      fL1PFToken_(consumes<vector<l1t::PFCandidate> >(cfg.getParameter<edm::InputTag>("L1PFObjects"))) {
  
  std::string lNNFile = cfg.getParameter<std::string>("NNFileName");  //,"L1Trigger/Phase2L1Taus/data/tau_3layer.pb");
  
  if(fHW) { 
    fTauNNIdHW_ = std::make_unique<TauNNIdHW>();
    fTauNNIdHW_->initialize("input_1:0", lNNFile, fNParticles_);
  } else { 
    fTauNNId_   = std::make_unique<TauNNId>();
    if (lNNFile.find("v0") == std::string::npos)
      fTauNNId_->initialize("input_1:0", lNNFile, fNParticles_);
    else if (lNNFile.find("v0") != std::string::npos)
      fTauNNId_->initialize("dense_1_input:0", lNNFile, fNParticles_);
  }
  produces<l1t::PFTauCollection>("L1PFTausNN");
}

void L1NNTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(fL1PFToken_, l1PFCandidates);

  std::vector<edm::Ptr<l1t::PFCandidate>> particles;
  for (unsigned i = 0; i < (*l1PFCandidates).size(); i++) {
    particles.push_back(edm::Ptr<l1t::PFCandidate>(l1PFCandidates, i));
  }
  auto lTaus = std::make_unique<l1t::PFTauCollection>();

  if (fHW) {
    process_HW(particles,*lTaus);
  } else {
    process_SW(particles,*lTaus);
  }

  if (lTaus->empty()) {
    PFTau dummy;
    lTaus->push_back(dummy);
  }
  std::sort(lTaus->begin(), lTaus->end(), [](l1t::PFTau i, l1t::PFTau j) { return (i.pt() > j.pt()); });
  iEvent.put(std::move(lTaus), "L1PFTausNN");
}

void L1NNTauProducer::process_SW(std::vector<edm::Ptr<l1t::PFCandidate>> parts, l1t::PFTauCollection &iTaus) { 

  std::vector<edm::Ptr<l1t::PFCandidate>> pfChargedHadrons;
  std::vector<edm::Ptr<l1t::PFCandidate>> pfChargedHadrons_sort;
  std::vector<edm::Ptr<l1t::PFCandidate>> pfChargedHadrons_seeds;
  for (const auto& l1PFCand : parts)
    if ((l1PFCand->id() == l1t::PFCandidate::ChargedHadron || l1PFCand->id() == l1t::PFCandidate::Electron) &&
        std::abs(l1PFCand->eta()) < 2.5)
      pfChargedHadrons_sort.push_back(l1PFCand);
  std::sort(pfChargedHadrons_sort.begin(), pfChargedHadrons_sort.end(), [](edm::Ptr<l1t::PFCandidate> i, edm::Ptr<l1t::PFCandidate> j) {
    return (i->pt() > j->pt());
    });

  pfChargedHadrons_seeds.push_back(pfChargedHadrons_sort[0]);
  for (unsigned int i0 = 1; i0 < pfChargedHadrons_sort.size(); i0++) {
    bool pMatch = false;
    l1t::PFCandidate part0 = *(pfChargedHadrons_sort[i0]);
    for (unsigned int i1 = 0; i1 < pfChargedHadrons_seeds.size(); i1++) {
      edm::Ptr<l1t::PFCandidate> part1 = pfChargedHadrons_seeds[i1];
      if (deltaR(part0, *part1) < fConeSize_)
        pMatch = true;
    }
    if (pMatch)continue;
    pfChargedHadrons_seeds.push_back(pfChargedHadrons_sort[i0]);
    if (pfChargedHadrons_seeds.size() > fMaxTaus_ - 1)break;
  }
  for (unsigned int i0 = 0; i0 < pfChargedHadrons_seeds.size(); i0++) {
    addTau(*(pfChargedHadrons_seeds[i0]), parts, iTaus);
  }
}

// create taus based on grid structure
void L1NNTauProducer::addTau(const l1t::PFCandidate& iCand,
                             std::vector<edm::Ptr<l1t::PFCandidate>> iParts,
                             l1t::PFTauCollection &outputTaus) {
  l1t::PFCandidateCollection pfTauCands;
  TLorentzVector lTot;
  lTot.SetPtEtaPhiM(0, 0, 0, 0);
  TLorentzVector lCand;
  lCand.SetPtEtaPhiM(0, 0, 0, 0);
  int lId = 0;
  for (auto l1PFCand : iParts) {
    if (deltaR(iCand, *l1PFCand) > fConeSize_) continue;
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(l1PFCand->pt(), l1PFCand->eta(), l1PFCand->phi(), 0);
    lTot += pVec;
    if (deltaR(iCand, *l1PFCand) < fTauSize_ &&
        (l1PFCand->id() == l1t::PFCandidate::Electron || l1PFCand->id() == l1t::PFCandidate::ChargedHadron ||
         l1PFCand->id() == l1t::PFCandidate::Photon)) {
      lId++;
      lCand += pVec;
    }
    pfTauCands.push_back(*l1PFCand);
  }
  if (lTot.Pt() < fSeedPt_)
    return;
  std::sort(
      pfTauCands.begin(), pfTauCands.end(), [](l1t::PFCandidate i, l1t::PFCandidate j) { return (i.pt() > j.pt()); });
  float NN = fTauNNId_->compute(iCand, pfTauCands);
  math::PtEtaPhiMLorentzVector tempP4(lCand.Pt(), lCand.Eta(), lCand.Phi(), lCand.M());
  l1t::PFTau l1PFTau(tempP4, NN, 0, lId);
  outputTaus.push_back(l1PFTau);
}
float L1NNTauProducer::deltaR(const l1t::PFCandidate& iPart1, const l1t::PFCandidate& iPart2) {
  float delta_r = 20;
  float pDPhi = fabs(iPart1.phi() - iPart2.phi());
  if (pDPhi > 2. * M_PI - pDPhi)
    pDPhi = 2. * M_PI - pDPhi;
  delta_r = sqrt((iPart1.eta() - iPart2.eta()) * (iPart1.eta() - iPart2.eta()) + pDPhi * pDPhi);
  return delta_r;
}


namespace L1TauEmu {
  // Data types and constants used in the FPGA and FPGA-optimized functions

  //etaphi_base maps physical eta phi units onto bits
  //This way, the least significant bit of etaphi_t is exactly 0.01
  //Even though 0.01 is not a power of 2
  static float etaphi_base = 100. / 64;

  typedef ap_ufixed<16, 14> pt_t;        // 1 unit = 0.25 GeV;
  typedef ap_fixed<10, 4> etaphi_t;      // 1 unit = 0.01;
  typedef ap_fixed<12, 6> detaphi_t;     // type for the difference between etas or phis
  typedef ap_fixed<18, 9> detaphi2_t;    // type for detaphi_t squared
  typedef ap_fixed<22, 16> pt_etaphi_t;  // type for product of pt with deta or phi
  typedef ap_uint<5> count_t;            // type for multiplicity
  typedef ap_uint<5> id_t;               // type for multiplicity

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
  //static const etaphi_t FIDUCIAL_ETA_PHI = 5.11;

  //typedef float  pt_t;
  //typedef float  etaphi_t;
  //typedef float  detaphi_t;
  //typedef float  detaphi2_t;    // type for detaphi_t squared                                                                                                                                                                      
  //typedef float  pt_etaphi_t;  // type for product of pt with deta or phi                                                                                                                                                         
  //typedef float  count_t;            // type for multiplicity                                                                                                                                                                           
  //typedef float  id_t;               // type for multiplicity   

  //typedef float inv_pt_t;
  //static constexpr int N_table_inv_pt = 1024;
  //static const float TWOPI = 3.14159 * 2. * etaphi_base;
  //static const float PI = 3.14159 * etaphi_base;
  //static const float HALFPI = 3.14159 / 2 * etaphi_base;

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
    /*
    if (debug) {
      std::cout << "           x " << in << ", msb = " << msb << ", shift = " << (in.width - msb) << ", idx = " << idx
                << std::endl;
      std::cout << "     pre 1 / " << in_shifted << " = " << inv_in << "(" << 1 / (float)in_shifted << ")" << std::endl;
      std::cout << "    post 1 / " << in << " = " << out << "(" << 1 / (float)in << ")" << std::endl;
    }
    */
    return out;
  }
 
  detaphi_t deltaPhi(edm::Ptr<l1t::PFCandidate> a, edm::Ptr<l1t::PFCandidate> b) {
    // scale the particle eta, phi to hardware units
    etaphi_t aphi = etaphi_t(a->phi() * etaphi_base);
    etaphi_t bphi = etaphi_t(b->phi() * etaphi_base);
    detaphi_t dphi = detaphi_t(aphi) - detaphi_t(bphi);
    // phi wrap
    detaphi_t dphi0 = dphi > PI ? detaphi_t(TWOPI - dphi) : detaphi_t(dphi);
    detaphi_t dphi1 = dphi < -PI ? detaphi_t(TWOPI + dphi) : detaphi_t(dphi);
    detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
    return dphiw;
  }

  bool inCone(edm::Ptr<l1t::PFCandidate> seed, edm::Ptr<l1t::PFCandidate> part, detaphi_t cone2, bool debug) {
    // scale the particle eta, phi to hardware units
    etaphi_t seta = etaphi_t(seed->eta() * etaphi_base);
    etaphi_t peta = etaphi_t(part->eta() * etaphi_base);
    detaphi_t deta = detaphi_t(seta) - detaphi_t(peta);
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = deta * deta + dphi * dphi < cone2;
    //bool ret = r2 < cone2;
    if (debug) {
      /*
      detaphi2_t r2 = detaphi2_t(deta) * detaphi2_t(deta) + detaphi2_t(dphi) * detaphi2_t(dphi);
      std::cout << "  part eta, seed eta: " << etaphi_t(part->eta() * etaphi_base) << ", "
                << etaphi_t(seed->eta() * etaphi_base) << std::endl;
      std::cout << "  part phi, seed phi: " << etaphi_t(part->phi() * etaphi_base) << ", "
                << etaphi_t(seed->phi() * etaphi_base) << std::endl;
      std::cout << "  pt, deta, dphi, r2, lt: " << part->pt() << ", " << deta << ", " << dphi << ", " << r2 << ", "
                << ret << std::endl;
      */
    }
    return ret;
  }

  void OutinCone(edm::Ptr<l1t::PFCandidate> seed, edm::Ptr<l1t::PFCandidate> part, detaphi_t cone2, bool debug) {
    // scale the particle eta, phi to hardware units
    etaphi_t seta = etaphi_t(seed->eta() * etaphi_base);
    etaphi_t peta = etaphi_t(part->eta() * etaphi_base);
    detaphi_t deta = detaphi_t(seta) - detaphi_t(peta);
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = deta * deta + dphi * dphi < cone2;
    std::cout << deta*deta + dphi*dphi << " -- " << cone2 << std::endl; 
    return;
  }
};  // namespace L1SCJetEmu

void L1NNTauProducer::makeTau_HW(edm::Ptr<l1t::PFCandidate> seed,std::vector<edm::Ptr<l1t::PFCandidate>>& parts,l1t::PFTauCollection &iTaus) const {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  //edm::Ptr<l1t::PFCandidate> seed = parts.at(0);
  
  l1t::PFCandidate tmpSeed = *seed;
  L1TauEmu::detaphi_t rCone2 = L1TauEmu::detaphi_t(fTauSize_ * fTauSize_ * L1TauEmu::etaphi_base * L1TauEmu::etaphi_base);
  unsigned lId = 0; 
  L1TauEmu::pt_t pt_tot = 0; 
  L1TauEmu::pt_t pt = 0; 
  for(unsigned i0 = 0; i0 < parts.size(); i0++) { 
    //std::cout << " ----> Adding Candidate: " << pt_t(parts[i0]->pt()) << "  -- " << parts[i0]->id() << " -- " << fTauSize_ << " -- ";
    //L1TauEmu::OutinCone(seed, (parts[i0]), rCone2, false);
    pt_tot = pt_tot +  L1TauEmu::pt_t(parts[i0]->pt());
    if (L1TauEmu::inCone(seed, (parts[i0]), rCone2, false)) { 
      if((parts[i0])->id() == l1t::PFCandidate::Electron      || 
	 (parts[i0])->id() == l1t::PFCandidate::ChargedHadron ||
	 (parts[i0])->id() == l1t::PFCandidate::Photon) {
	lId++;
	pt = pt + L1TauEmu::pt_t(parts[i0]->pt());
      }
    }
  }
  //std::cout << "---> pt " << pt << " -- " << fSeedPt_ << std::endl;
  if(pt < fSeedPt_) return;
   //L1TauEmu::pt_etaphi_t pt  = pt_t(seed.pt());
  result_t NN = fTauNNIdHW_->compute(tmpSeed,parts);
  L1TauEmu::etaphi_t    eta = etaphi_t(seed->eta() * L1TauEmu::etaphi_base);
  L1TauEmu::etaphi_t    phi = etaphi_t(seed->phi() * L1TauEmu::etaphi_base);

  math::PtEtaPhiMLorentzVector tempP4(float(pt), float(eta)/L1TauEmu::etaphi_base, float(phi)/L1TauEmu::etaphi_base, 0.);
  l1t::PFTau l1PFTau(tempP4, NN, 0, lId);
  iTaus.push_back(l1PFTau);
  /*
  if (fDebug) {
    std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x) { std::cout << "pt_dphi: " << x << std::endl; });
    std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x) { std::cout << "pt_deta: " << x << std::endl; });
    std::cout << " sum_pt_eta: " << sum_pt_eta << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    std::cout << " sum_pt_phi: " << sum_pt_phi << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
    std::cout << " uncorr eta: " << etaphi_t(seed->eta() * etaphi_base)
              << ", phi: " << etaphi_t(seed->phi() * etaphi_base) << std::endl;
    std::cout << "   corr eta: " << eta << ", phi: " << phi << std::endl;
    std::cout << "         pt: " << pt << std::endl;
  }
  */

  //return l1PFTau;
}

static constexpr float track_trigger_eta_max = 2.5;

void L1NNTauProducer::process_HW(std::vector<edm::Ptr<l1t::PFCandidate>> parts, l1t::PFTauCollection &iTaus) const { 
  // The fixed point algorithm emulation
  using namespace L1TauEmu;
  std::vector<edm::Ptr<l1t::PFCandidate>> work;
  work.resize(parts.size());
  

  std::transform(parts.begin(), parts.end(), work.begin(), [](const edm::Ptr<l1t::PFCandidate>& part) { return part; });
  std::sort(work.begin(), work.end(), [](edm::Ptr<l1t::PFCandidate> i, edm::Ptr<l1t::PFCandidate> j) { return (pt_t(i->pt()) > pt_t(j->pt())); });

  std::vector<edm::Ptr<l1t::PFCandidate>> seeds;
  std::copy_if(
	       work.begin(), work.end(), std::back_inserter(seeds), [&](const edm::Ptr<l1t::PFCandidate>& part) {
		 return ((part->id() ==  l1t::PFCandidate::Electron ||  part->id() ==  l1t::PFCandidate::ChargedHadron)  && std::abs(part->eta()) < track_trigger_eta_max);
	       });
  // It would be nice to transform the inputs to the etaphi_base of the FW here, as in the line below
  // However the phi may wrap around if the etaphi_base > 1, so don't do it...
  //std::for_each(work.begin(), work.end(), [](l1t::PFCandidate& x){x.setP4(math::PtEtaPhiMLorentzVector(pt_t(x.pt()), etaphi_t(x.eta()*etaphi_base), etaphi_t(x.phi()*etaphi_base), x.mass()));});
  detaphi_t rCone2 = detaphi_t(fConeSize_ * fConeSize_ * etaphi_base * etaphi_base);

  iTaus.reserve(fMaxTaus_);
  while (!seeds.empty() && iTaus.size() < fMaxTaus_) {
    // Take the first (highest pt) candidate as a seed
    edm::Ptr<l1t::PFCandidate> seed = seeds.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<edm::Ptr<l1t::PFCandidate>> particlesInCone;
    std::copy_if(
        work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const edm::Ptr<l1t::PFCandidate>& part) {
          return inCone(seed, part, rCone2, false);
        });
    /*
    if (fDebug) {
      std::cout << "Seed: " << pt_t(seed->pt()) << ", " << etaphi_t(seed->eta() * etaphi_base) << ", "
                << etaphi_t(seed->phi() * etaphi_base) << std::endl;
      std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](edm::Ptr<l1t::PFCandidate>& part) {
        std::cout << "  Part: " << pt_t(part->pt()) << ", " << etaphi_t(part->eta() * etaphi_base) << ", "
                  << etaphi_t(part->phi() * etaphi_base) << std::endl;
        inCone(seed, part, rCone2, true);
      });
    }
    */
    makeTau_HW(seed,particlesInCone,iTaus);
    // remove the clustered particles
    work.erase(
        std::remove_if(work.begin(),
                       work.end(),
                       [&](const edm::Ptr<l1t::PFCandidate>& part) { return inCone(seed, part, rCone2, false); }),
        work.end());

    seeds.erase(
        std::remove_if(seeds.begin(),
                       seeds.end(),
                       [&](const edm::Ptr<l1t::PFCandidate>& part) { return inCone(seed, part, rCone2, false); }),
        seeds.end());
  }
}


L1NNTauProducer::~L1NNTauProducer() {}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1NNTauProducer);
