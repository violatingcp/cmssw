#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include <vector>

class FacileHcalReconstructor : public TritonEDProducer<> {
public:
  explicit FacileHcalReconstructor(const edm::ParameterSet&);
  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  edm::InputTag fChannelInfoName_;
  edm::EDGetTokenT<std::vector<HBHEChannelInfo>> fTokChannelInfo_;
  std::vector<HcalDetId> hcalIds_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> htopoToken_;

};

FacileHcalReconstructor::FacileHcalReconstructor(edm::ParameterSet const& cfg)
    : TritonEDProducer<>(cfg,"FacileHcalReconstructor"),
      fChannelInfoName_(cfg.getParameter<edm::InputTag>("ChannelInfoName")),
      fTokChannelInfo_(consumes<std::vector<HBHEChannelInfo>>(fChannelInfoName_)),
      htopoToken_(esConsumes<HcalTopology, HcalRecNumberingRecord>()) {
    produces<HBHERecHitCollection>();
    //setDebugName("FacileHcalReconstructor");
}


void FacileHcalReconstructor::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) {
    const auto& hChannelInfo = iEvent.get(fTokChannelInfo_);

    const HcalTopology* htopo = &iSetup.getData(htopoToken_);

    auto& input = iInput.at("continuousinputs");
    auto& input_depth = iInput.at("depth");
    auto& input_ieta = iInput.at("ieta");
    //auto& input1 = iInput.begin()->second;
    client_->setBatchSize(hChannelInfo.size());

    auto tdata = input.allocate<float>(true);
    auto tdata_depth = input_depth.allocate<int>(true);
    auto tdata_ieta = input_ieta.allocate<int>(true);

    hcalIds_.clear();
    hcalIds_.reserve(hChannelInfo.size());
    unsigned int i = 0;
    for (const auto& pChannel : hChannelInfo) {
      //std::vector<float> input;
      const HcalDetId pDetId = pChannel.id();
      hcalIds_.push_back(pDetId);
      auto& idata = (*tdata)[i];
      auto& idata_depth = (*tdata_depth)[i];
      auto& idata_ieta = (*tdata_ieta)[i];

      //inputs for Facile: iphi, gain, raw[8], depth (categorical), ieta (categorical)
      //idata["iphi"].push_back(pDetId.iphi());
      //idata["gain"].push_back(pChannel.tsGain(0.));
      for (unsigned int itTS = 0; itTS < pChannel.nSamples(); ++itTS) {
        idata.push_back(pChannel.tsRawCharge(itTS) - pChannel.tsPedestal(itTS));
      }
      idata.push_back(pChannel.tsGain(0.));
      idata.push_back(pDetId.iphi());

      idata_ieta.push_back(pDetId.ietaAbs());
      idata_depth.push_back(pDetId.depth());
      i = i+1;
      //for (int itDepth = 1; itDepth <= htopo->maxDepth(); itDepth++) {
      //  input.push_back(pDetId.depth() == itDepth);
      //}

      //for (int itIeta = 1; itIeta <= htopo->lastHERing(); itIeta++) {
      //  input.push_back(pDetId.ietaAbs() == itIeta);
      //}

      //data1->push_back(input);
    }
    //for (int i = hChannelInfo.size(); i < 10000; i++){
    //  auto& input = data1[i];
      //for (int ii = 0; ii < 12; ii++){
      //  input.push_back(0.f);
      //}
    //  data1->push_back(input);
    //}
    input.toServer(tdata);
    input_ieta.toServer(tdata_ieta);
    input_depth.toServer(tdata_depth);
}

void FacileHcalReconstructor::produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {
    std::unique_ptr<HBHERecHitCollection> out;
    out = std::make_unique<HBHERecHitCollection>();
    out->reserve(hcalIds_.size());
    //std::cout <<"doing facile" << std::endl;
    const auto& output1 = iOutput.at("output"); //begin()->second;
    const auto& outputs = output1.fromServer<float>();
    for (std::size_t iB = 0; iB < hcalIds_.size(); iB++) {
      float rhE = outputs[iB][0];
      if (rhE < 0.f or std::isnan(rhE) or std::isinf(rhE))
        rhE = 0;
      HBHERecHit rh(hcalIds_[iB], rhE, 0.f, 0.f);
      //std::cout << "Hcal ieta/iphi" << rh.id().ietaAbs() << "/" << rh.id().iphi() << "\n" <<
      //    "\tenergy:" << rhE << std::endl;
      std::cout << "Running FACILE." << std::endl;
      out->push_back(rh);
    }
    iEvent.put(std::move(out));
}

void FacileHcalReconstructor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    TritonClient::fillPSetDescription(desc);
    desc.add<edm::InputTag>("ChannelInfoName");
    descriptions.add("FacileHcalReconstructor", desc);
}


DEFINE_FWK_MODULE(FacileHcalReconstructor);
