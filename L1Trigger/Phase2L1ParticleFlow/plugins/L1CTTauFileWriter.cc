#include <memory>
#include <numeric>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"

//
// Pattern file writer for the NN taus, following L1CTJetFileWriter.
//
// l1gt::Tau is 128 bits and packs into two 64-bit words, exactly as l1gt::Jet
// does, so the link encoding is the same shape as the jet writer's: two words
// per object, zero padded up to nTaus.
//
// Taus carry no accompanying sums, so unlike the jet writer there is no second
// token per collection - a collection is just a tau InputTag and a multiplicity.
//

class L1CTTauFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1CTTauFileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  std::vector<edm::ParameterSet> collections_;

  size_t nFramesPerBX_;
  size_t ctl2BoardTMUX_;
  size_t gapLengthOutput_;
  size_t maxLinesPerFile_;
  std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>> channelSpecsOutputToGT_;

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::vector<ap_uint<64>> encodeTaus(const std::vector<l1t::PFTau> taus, unsigned nTaus);

  l1t::demo::BoardDataWriter fileWriterOutputToGT_;
  std::vector<edm::EDGetTokenT<edm::View<l1t::PFTau>>> tokens_;
  std::vector<unsigned> nTaus_;
};

L1CTTauFileWriter::L1CTTauFileWriter(const edm::ParameterSet& iConfig)
    : collections_(iConfig.getParameter<std::vector<edm::ParameterSet>>("collections")),
      nFramesPerBX_(iConfig.getParameter<unsigned>("nFramesPerBX")),
      ctl2BoardTMUX_(iConfig.getParameter<unsigned>("TMUX")),
      gapLengthOutput_(iConfig.getParameter<unsigned>("gapLengthOutput")),
      maxLinesPerFile_(iConfig.getParameter<unsigned>("maxLinesPerFile")),
      channelSpecsOutputToGT_{{{"taus", 0}, {{ctl2BoardTMUX_, gapLengthOutput_}, {0}}}},
      fileWriterOutputToGT_(l1t::demo::parseFileFormat(iConfig.getParameter<std::string>("format")),
                            iConfig.getParameter<std::string>("outputFilename"),
                            iConfig.getParameter<std::string>("outputFileExtension"),
                            nFramesPerBX_,
                            ctl2BoardTMUX_,
                            maxLinesPerFile_,
                            channelSpecsOutputToGT_) {
  for (const auto& pset : collections_) {
    unsigned nTaus = pset.getParameter<unsigned>("nTaus");
    nTaus_.push_back(nTaus);
    // An empty token is never read: analyze() skips any collection asking for
    // zero taus, which keeps a zero-multiplicity entry from consuming a branch
    // that may not exist in the event.
    edm::EDGetTokenT<edm::View<l1t::PFTau>> tauToken;
    if (nTaus > 0) {
      tauToken = consumes<edm::View<l1t::PFTau>>(pset.getParameter<edm::InputTag>("taus"));
    }
    tokens_.push_back(tauToken);
  }
}

void L1CTTauFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // 1) Pack collections in the order they are specified
  std::vector<ap_uint<64>> link_words;
  for (unsigned iCollection = 0; iCollection < collections_.size(); iCollection++) {
    if (nTaus_.at(iCollection) == 0)
      continue;

    // 2) Encode tau information onto vectors containing link data
    const edm::View<l1t::PFTau>& taus = iEvent.get(tokens_.at(iCollection));
    std::vector<l1t::PFTau> sortedTaus;
    sortedTaus.reserve(taus.size());
    std::copy(taus.begin(), taus.end(), std::back_inserter(sortedTaus));

    // stable_sort so that taus of equal pt keep the producer's ordering, which
    // is what the firmware comparison relies on
    std::stable_sort(
        sortedTaus.begin(), sortedTaus.end(), [](l1t::PFTau i, l1t::PFTau j) { return (i.hwPt() > j.hwPt()); });
    const auto outputTaus(encodeTaus(sortedTaus, nTaus_.at(iCollection)));
    link_words.insert(link_words.end(), outputTaus.begin(), outputTaus.end());
  }

  // 3) Pack tau information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataTaus;
  eventDataTaus.add({"taus", 0}, link_words);
  fileWriterOutputToGT_.addEvent(eventDataTaus);
}

// ------------ method called once each job just after ending the event loop  ------------
void L1CTTauFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterOutputToGT_.flush();
}

std::vector<ap_uint<64>> L1CTTauFileWriter::encodeTaus(const std::vector<l1t::PFTau> taus, const unsigned nTaus) {
  // Encode up to nTaus taus, padded with 0s
  std::vector<ap_uint<64>> tau_words(2 * nTaus, 0);  // allocate 2 words per tau
  for (unsigned i = 0; i < std::min(nTaus, (uint)taus.size()); i++) {
    const l1t::PFTau& t = taus.at(i);
    // encodedTau() is the l1gt::Tau already packed by the producer, so the
    // words written here are bit-for-bit what the GT link carries
    tau_words[2 * i] = t.encodedTau()[0];
    tau_words[2 * i + 1] = t.encodedTau()[1];
  }
  return tau_words;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1CTTauFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.addOptional<edm::InputTag>("taus");
    vpsd1.add<uint>("nTaus", 0);
    desc.addVPSet("collections", vpsd1);
  }
  desc.add<std::string>("outputFilename");
  desc.add<std::string>("outputFileExtension", "txt");
  desc.add<uint32_t>("nFramesPerBX", 9);
  desc.add<uint32_t>("gapLengthOutput", 4);
  desc.add<uint32_t>("TMUX", 6);
  desc.add<uint32_t>("maxLinesPerFile", 1024);
  desc.add<std::string>("format", "EMPv2");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CTTauFileWriter);
