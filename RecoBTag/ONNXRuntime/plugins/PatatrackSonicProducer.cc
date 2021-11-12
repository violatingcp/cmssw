
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/SiPixelDigi/interface/SiPixelDigisSoA.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/BeamSpot/interface/BeamSpotPOD.h"
#include "DataFormats/SiPixelRawData/interface/SiPixelErrorsSoA.h"
#include "CUDADataFormats/Track/interface/PixelTrackHeterogeneous.h"
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEFast.h"

#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"
#include "EventFilter/SiPixelRawToDigi/interface/ErrorChecker.h"
#include "EventFilter/SiPixelRawToDigi/interface/PixelDataFormatter.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitsSoA.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <nlohmann/json.hpp>

class PatatrackSonicProducer : public TritonEDProducer<> {
public:
  explicit PatatrackSonicProducer(const edm::ParameterSet &);
  void acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) override;
  void produce(edm::Event &iEvent, edm::EventSetup const &iSetup, Output const &iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:

  edm::ESWatcher<SiPixelFedCablingMapRcd> recordWatcher_;
  edm::ESGetToken<SiPixelFedCablingMap, SiPixelFedCablingMapRcd> cablingMapToken_;

  edm::EDGetTokenT<FEDRawDataCollection> rawGetToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsGetToken_;

  edm::EDPutTokenT<SiPixelRecHitsSoA> hitsSOA_;
  edm::EDPutTokenT<SiPixelDigisSoA> digiPutToken_;
  edm::EDPutTokenT<SiPixelErrorsSoA> digiErrorPutToken_;

  edm::EDPutTokenT<ZVertexHeterogeneous> vertexSOA_;
  edm::EDPutTokenT<PixelTrackHeterogeneous> trackSOA_;

  PixelDataFormatter::Errors errors_;
  std::vector<unsigned int> fedIds_;
  const SiPixelFormatterErrors* formatterErrors_ = nullptr;
  bool debug_ = false;
};

PatatrackSonicProducer::PatatrackSonicProducer(const edm::ParameterSet &iConfig)
    : TritonEDProducer<>(iConfig, "PatatrackProducer"),
      cablingMapToken_(esConsumes<SiPixelFedCablingMap, SiPixelFedCablingMapRcd>(
										 edm::ESInputTag("", iConfig.getParameter<std::string>("CablingMapLabel")))),
      rawGetToken_(consumes<FEDRawDataCollection>(iConfig.getParameter<edm::InputTag>("InputLabel"))),
      bsGetToken_{consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))},
      hitsSOA_(produces<SiPixelRecHitsSoA>()),
      digiPutToken_(produces<SiPixelDigisSoA>()),
      digiErrorPutToken_(produces<SiPixelErrorsSoA>()),
      vertexSOA_(produces<ZVertexHeterogeneous>()),
      trackSOA_(produces<PixelTrackHeterogeneous>()),
      debug_(iConfig.getUntrackedParameter<bool>("debugMode", false)) {
	formatterErrors_ = new SiPixelFormatterErrors();
      }

void PatatrackSonicProducer::acquire(edm::Event const &iEvent, edm::EventSetup const &iSetup, Input &iInput) {
  const reco::BeamSpot& bs = iEvent.get(bsGetToken_);
  const auto& buffers = iEvent.get(rawGetToken_);
  // initialize cabling map or update if necessary
  if (recordWatcher_.check(iSetup)) {
    //cabling map, which maps online address (fed->link->ROC->local pixel) to offline (DetId->global pixel)
    auto cablingMap = iSetup.getTransientHandle(cablingMapToken_);
    fedIds_ = cablingMap->fedIds();
  }
  //Note Fed data quality  checks ae curently done on the CPU server, and could be moved here
  auto& input = iInput.at("input");
  auto  feds  = input.allocate<uint32_t>();
  auto& vin   = (*feds)[0];
  unsigned int pSize = 0; 

  BeamSpotPOD bsHost;
  bsHost.x = bs.x0();
  bsHost.y = bs.y0();
  bsHost.z = bs.z0();
  bsHost.sigmaZ = bs.sigmaZ();
  bsHost.beamWidthX = bs.BeamWidthX();
  bsHost.beamWidthY = bs.BeamWidthY();
  bsHost.dxdz = bs.dxdz();
  bsHost.dydz = bs.dydz();
  bsHost.emittanceX = bs.emittanceX();
  bsHost.emittanceY = bs.emittanceY();
  bsHost.betaStar = bs.betaStar();
  unsigned bsSize = 11;
  vin.resize(vin.size()+bsSize); pSize+=bsSize;
  std::memcpy(&vin[vin.size() - bsSize],&bsHost,bsSize*sizeof(float));
  //Now deal with FEDs
  vin.push_back(fedIds_.size()); pSize++;
  //ErrorChecker errorcheck;
  //bool errorsInEvent = false;
  //errors_.clear();
  for (unsigned int fedId : fedIds_) {
    vin.push_back(fedId); pSize++;
    const FEDRawData& rawData = buffers.FEDData(fedId);
    /*
    int nWords = rawData.size() / sizeof(uint64_t);
    if (nWords == 0) {
      std::cout << " !!!! Continuing " << std::endl;
      continue;
    }
    const cms_uint64_t* trailer = reinterpret_cast<const cms_uint64_t*>(rawData.data()) + (nWords - 1);
    if (not errorcheck.checkCRC(errorsInEvent, fedId, trailer, errors_)) {
      continue;
    }
    */
    unsigned int rawsize=rawData.size()/4;
    vin.push_back(rawsize); pSize++;
    vin.resize(vin.size()+rawsize);
    std::memcpy(&vin[vin.size() - rawsize],rawData.data(), rawData.size());
    pSize += rawsize;
  }
  input.toServer(feds);
} 

void PatatrackSonicProducer::produce(edm::Event &iEvent,
				     const edm::EventSetup &iSetup,
				     Output const &iOutput) {
  
  uint32_t pdigi_[150000];
  uint32_t rawIdArr_[150000];
  uint16_t adc_ [150000];
  int32_t  clus_[150000];
  uint32_t hits_[2001];
  float    pos_ [4*35000];
  SiPixelErrorCompact  pixerrors_[20];
  
  auto hits   = std::make_unique<SiPixelRecHitsSoA>();

  //PixelTrackHeterogeneous tracks;
  const auto &output1 = iOutput.begin()->second;
  const auto &outputs_from_server = output1.fromServer<int8_t>();
  auto output = (outputs_from_server[0]);  
  unsigned int pCount = 0;
  uint32_t nHits      = 0; //output[pCount]; pCount++;
  std::memcpy(&nHits,&(output.front())+pCount,sizeof(uint32_t)); pCount += 4;
  static const unsigned nMax = 2000; 
  //if(nHits_ < 2000) nMax = nHits_;
  std::memcpy(hits_,&(output.front())+pCount,(nMax+1)*sizeof(uint32_t));    pCount += 4*(nMax+1);
  std::memcpy(pos_, &(output.front())+pCount,4*nHits*sizeof(float));         pCount += 4*4*nHits;
  iEvent.emplace(hitsSOA_,      nHits, hits_, pos_); 

  uint32_t nDigis    = 0; //output[pCount]; pCount++;
  std::memcpy(&nDigis,&(output.front())+pCount,sizeof(uint32_t)); pCount += 4;
  std::memcpy(pdigi_,   &(output.front())+pCount,nDigis*sizeof(uint32_t)); pCount += 4*nDigis;
  std::memcpy(rawIdArr_,&(output.front())+pCount,nDigis*sizeof(uint32_t)); pCount += 4*nDigis;
  std::memcpy(adc_,     &(output.front())+pCount,nDigis*sizeof(uint16_t)); pCount += 2*nDigis;
  std::memcpy(clus_,    &(output.front())+pCount,nDigis*sizeof(int32_t));  pCount += 4*nDigis;
  iEvent.emplace(digiPutToken_, nDigis, pdigi_, rawIdArr_, adc_, clus_);

  uint32_t nErrors = 0; 
  std::memcpy(&nErrors,&(output.front())+pCount,sizeof(uint32_t)); pCount += 4;
  std::memcpy(pixerrors_, &(output.front())+pCount,10*nErrors);     pCount += 10*nErrors;
  iEvent.emplace(digiErrorPutToken_, nErrors, pixerrors_, formatterErrors_);

  static constexpr uint32_t MAXTRACKS = 32 * 1024;
  unsigned int nTracks = 0;
  auto tracks = std::make_unique<pixelTrack::TrackSoA>();
  std::memcpy(&nTracks,&(output.front())+pCount,sizeof(uint32_t)); pCount += 4;
  tracks->ntFinal = nTracks;
  std::memcpy((*tracks).chi2.data(),      &(output.front())+pCount,nTracks*sizeof(float));                 pCount+=4*nTracks;
  std::memcpy((*tracks).qualityData(),    &(output.front())+pCount,nTracks*sizeof(uint8_t));               pCount+=1*nTracks;
  std::memcpy((*tracks).eta.data(),       &(output.front())+pCount,nTracks*sizeof(float));                 pCount+=4*nTracks;
  std::memcpy((*tracks).pt.data(),        &(output.front())+pCount,nTracks*sizeof(float));                 pCount+=4*nTracks;
  for(unsigned i1 = 0; i1 < 5; i1++) { 
    std::memcpy(tracks->stateAtBS.state(0).data()+MAXTRACKS*i1,     &(output.front())+pCount,nTracks*sizeof(float));  pCount+=4*(nTracks);
  }
  for(unsigned i1 = 0; i1 < 15; i1++) { 
    std::memcpy(tracks->stateAtBS.covariance(0).data()+MAXTRACKS*i1,&(output.front())+pCount,nTracks*sizeof(float)); pCount+=4*(nTracks);
  }
  std::memcpy((void*)(*tracks).hitIndices.content.data(),&(output.front())+pCount,nTracks*sizeof(uint32_t)*5);     pCount+=4*(nTracks*5);
  std::memcpy((*tracks).hitIndices.off.data(),           &(output.front())+pCount,(nTracks+1)*sizeof(int32_t));    pCount+=4*(nTracks+1);
  std::memcpy((void*)(*tracks).detIndices.content.data(),&(output.front())+pCount,nTracks*sizeof(uint32_t)*5);     pCount+=4*(nTracks*5);
  std::memcpy((*tracks).detIndices.off.data(),           &(output.front())+pCount,(nTracks+1)*sizeof(int32_t));    pCount+=4*(nTracks+1);
  iEvent.emplace(trackSOA_,  PixelTrackHeterogeneous(std::move(tracks)));

  auto vertices = std::make_unique<ZVertexSoA>();
  std::memcpy(&(vertices->nvFinal),&(output.front())+pCount,sizeof(uint32_t)); pCount += 4;
  unsigned int nVtx = vertices->nvFinal;
  std::memcpy((vertices)->idv    , &(output.front())+pCount,nTracks*sizeof(int16_t));   pCount+=2*nTracks;
  std::memcpy((vertices)->zv     , &(output.front())+pCount,nVtx*sizeof(float));        pCount+=4*nVtx;
  std::memcpy((vertices)->wv     , &(output.front())+pCount,nVtx*sizeof(float));        pCount+=4*nVtx;
  std::memcpy((vertices)->chi2   , &(output.front())+pCount,nVtx*sizeof(float));        pCount+=4*nVtx;
  std::memcpy((vertices)->ptv2   , &(output.front())+pCount,nVtx*sizeof(float));        pCount+=4*nVtx;
  std::memcpy((vertices)->ndof   , &(output.front())+pCount,nVtx*sizeof(int32_t));      pCount+=4*nVtx;
  std::memcpy((vertices)->sortInd, &(output.front())+pCount,nVtx*sizeof(uint16_t));     pCount+=2*nVtx;
  iEvent.emplace(vertexSOA_, ZVertexHeterogeneous(std::move(vertices)));
}

void PatatrackSonicProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("InputLabel");
  desc.add<edm::InputTag>("beamSpot");
  desc.add<std::string>("CablingMapLabel");
  desc.addOptionalUntracked<bool>("debugMode", false);
  descriptions.add("PatatrackSonicProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PatatrackSonicProducer);
