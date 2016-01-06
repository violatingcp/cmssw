#ifndef CommonTools_Puppi_PuppiPhoton_h_
#define CommonTools_Puppi_PuppiPhoton_h_
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// ------------------------------------------------------------------------------------------
class PuppiPhoton : public edm::stream::EDProducer<> {

public:
	explicit PuppiPhoton(const edm::ParameterSet&);
	~PuppiPhoton();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	typedef math::XYZTLorentzVector                        LorentzVector;
	typedef std::vector<LorentzVector>                     LorentzVectorCollection;
	typedef edm::View<reco::Candidate>                     CandidateView;
	typedef std::vector< reco::PFCandidate >               PFOutputCollection;
	typedef edm::View<reco::PFCandidate>                   PFView;

private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
      
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	bool matchPFCandidate(const reco::Candidate *iPF,const reco::Candidate *iPho);
	edm::EDGetTokenT< CandidateView >        tokenPFCandidates_;
	edm::EDGetTokenT< CandidateView >        tokenPhotonCandidates_;
	edm::EDGetTokenT< edm::ValueMap<float> > tokenWeights_; 
	edm::EDGetTokenT< edm::ValueMap<bool>  > tokenPhotonId_; 
	double pt_;
	bool   usePFRef_;
	std::vector<double>  dRMatch_;
	std::vector<int32_t> pdgIds_; 
        std::auto_ptr< PFOutputCollection > corrCandidates_;
};
#endif
