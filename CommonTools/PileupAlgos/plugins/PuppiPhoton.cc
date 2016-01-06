// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Common/interface/Association.h"
//Main File
#include "CommonTools/PileupAlgos/plugins/PuppiPhoton.h"

// ------------------------------------------------------------------------------------------
PuppiPhoton::PuppiPhoton(const edm::ParameterSet& iConfig) {
  tokenPFCandidates_     = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("candName"));
  tokenPhotonCandidates_ = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("photonName"));
  tokenWeights_          = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("weightsName"));
  tokenPhotonId_         = consumes<edm::ValueMap<bool>  >(iConfig.getParameter<edm::InputTag>("photonId")); 
  pt_                    = iConfig.getUntrackedParameter<double>("pt");
  dRMatch_               = iConfig.getUntrackedParameter<std::vector<double> > ("dRMatch");
  pdgIds_                = iConfig.getUntrackedParameter<std::vector<int32_t> >("pdgids");
  usePFRef_              = iConfig.getUntrackedParameter<bool>("useRefs");

  produces<edm::ValueMap<float> > ();
  produces<edm::ValueMap<LorentzVector> > ();
  produces< edm::ValueMap<reco::CandidatePtr> >(); 
  produces<PFOutputCollection>();
}
// ------------------------------------------------------------------------------------------
PuppiPhoton::~PuppiPhoton(){}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<CandidateView> hPhoProduct;
  iEvent.getByToken(tokenPhotonCandidates_,hPhoProduct);
  const CandidateView *phoCol = hPhoProduct.product();

  edm::Handle<edm::ValueMap<bool> > photonId;
  iEvent.getByToken(tokenPhotonId_,photonId);
  int iC = -1;
  std::vector<const reco::Candidate*> phoCands;
  std::vector<uint16_t> phoIndx;

  // Get PFCandidate Collection
  edm::Handle<CandidateView> hPFProduct;
  iEvent.getByToken(tokenPFCandidates_,hPFProduct);
  const CandidateView *pfCol = hPFProduct.product();

  for(CandidateView::const_iterator itPho = phoCol->begin(); itPho!=phoCol->end(); itPho++) {
    iC++;
    bool passObject = false;
    if(itPho->isPhoton())   passObject =  (*photonId)  [phoCol->ptrAt(iC)];
    if(itPho->pt() < pt_) continue;
    if(!passObject) continue;
    phoCands.push_back(&(*itPho)); 
    if(!usePFRef_) continue;
    try{
      const pat::Photon *pPho = dynamic_cast<const pat::Photon*>(&(*itPho));
      if(pPho != 0) for( const edm::Ref<pat::PackedCandidateCollection> & ref : pPho->associatedPackedPFCandidates() ) 
		      if(matchPFCandidate(&(*(pfCol->ptrAt(ref.key()))),&(*itPho))) phoIndx.push_back(ref.key());
      if(pPho != 0) continue;
      const pat::Electron *pElectron = dynamic_cast<const pat::Electron*>(&(*itPho));
      if(pElectron != 0) for( const edm::Ref<pat::PackedCandidateCollection> & ref : pElectron->associatedPackedPFCandidates() ) 
			   if(matchPFCandidate(&(*(pfCol->ptrAt(ref.key()))),&(*itPho)))  phoIndx.push_back(ref.key());
    } catch(...) { 
      throw edm::Exception(edm::errors::LogicError,"PuppiPhoton: photon references not found"); 
    }
  }

  //Get Weights
  edm::Handle<edm::ValueMap<float> > pupWeights; 
  iEvent.getByToken(tokenWeights_,pupWeights);

  std::auto_ptr<edm::ValueMap<LorentzVector> > p4PupOut(new edm::ValueMap<LorentzVector>());
  LorentzVectorCollection puppiP4s;
  std::vector<reco::CandidatePtr> values(hPFProduct->size());
  int iPF = 0; 
  std::vector<float> lWeights;
  static const reco::PFCandidate dummySinceTranslateIsNotStatic;
  corrCandidates_.reset( new PFOutputCollection );
  for(CandidateView::const_iterator itPF = pfCol->begin(); itPF!=pfCol->end(); itPF++) {  
    auto id = dummySinceTranslateIsNotStatic.translatePdgIdToType(itPF->pdgId());
    const reco::PFCandidate *pPF = dynamic_cast<const reco::PFCandidate*>(&(*itPF));
    reco::PFCandidate pCand( pPF ? *pPF : reco::PFCandidate(itPF->charge(), itPF->p4(), id) );
    LorentzVector pVec = itPF->p4();
    float pWeight = (*pupWeights)[pfCol->ptrAt(iPF)];     
    if(!usePFRef_) { 
      for(std::vector<const reco::Candidate*>::iterator itPho = phoCands.begin(); itPho!=phoCands.end(); itPho++) {
	if(matchPFCandidate(&(*itPF),*itPho)) pWeight = 1.;
      }
    } else { 
      for(std::vector<uint16_t>::const_iterator itPho = phoIndx.begin(); itPho!=phoIndx.end(); itPho++) {
	if(pfCol->refAt(iPF).key() != *itPho) continue;
	pWeight = 1.;
      }
    }
    pVec.SetPxPyPzE(itPF->px()*pWeight,itPF->py()*pWeight,itPF->pz()*pWeight,itPF->energy()*pWeight);
    lWeights.push_back(pWeight);
    pCand.setP4(pVec);
    puppiP4s.push_back( pVec );
    pCand.setSourceCandidatePtr( itPF->sourceCandidatePtr(0) );
    corrCandidates_->push_back(pCand);
    iPF++;
  }
  //Fill it into the event
  std::auto_ptr<edm::ValueMap<float> > lPupOut(new edm::ValueMap<float>());

  edm::ValueMap<float>::Filler  lPupFiller(*lPupOut);
  lPupFiller.insert(hPFProduct,lWeights.begin(),lWeights.end());
  lPupFiller.fill();
  
  //Compute the modified p4s
  edm::ValueMap<LorentzVector>::Filler  p4PupFiller(*p4PupOut);
  p4PupFiller.insert(hPFProduct,puppiP4s.begin(), puppiP4s.end() );
  p4PupFiller.fill();
  
  iEvent.put(lPupOut);
  iEvent.put(p4PupOut);
  edm::OrphanHandle<reco::PFCandidateCollection> oh = iEvent.put( corrCandidates_ );
  for(unsigned int ic=0, nc = oh->size(); ic < nc; ++ic) {
      reco::CandidatePtr pkref( oh, ic );
      values[ic] = pkref;
  }  
  std::auto_ptr<edm::ValueMap<reco::CandidatePtr> > pfMap_p(new edm::ValueMap<reco::CandidatePtr>());
  edm::ValueMap<reco::CandidatePtr>::Filler filler(*pfMap_p);
  filler.insert(hPFProduct, values.begin(), values.end());
  filler.fill();
  iEvent.put(pfMap_p);
}
// ------------------------------------------------------------------------------------------
bool PuppiPhoton::matchPFCandidate(const reco::Candidate *iPF,const reco::Candidate *iPho) { 
  double pDEta = iPF->eta()-iPho->eta();
  double pDPhi = std::abs(iPF->phi()-iPho->phi());
  if(pDPhi > 2.*M_PI-pDPhi) pDPhi =  2.*M_PI-pDPhi;
  double pDR2 = pDEta*pDEta+pDPhi*pDPhi;
  for(unsigned int i0 = 0; i0 < pdgIds_.size(); i0++) {
    if(std::abs(iPF->pdgId()) == pdgIds_[i0] && pDR2 < dRMatch_[i0])  return true;
  }
  return false;
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::beginJob() {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::endJob() {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::beginRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::endRun(edm::Run&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}
// ------------------------------------------------------------------------------------------
void PuppiPhoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PuppiPhoton);
