import FWCore.ParameterSet.Config as cms

# Pattern file writer for the NN taus, matching l1tJetFileWriter_cfi.
#
# L1NNTauProducer publishes its collection under the "L1PFTausNN" instance
# label, so the InputTag needs both the module and the instance.
#
# nTaus defaults to 5, which is maxtaus in L1NNTauProducer_cff - writing more
# than the producer can emit would only pad with zeros.

l1tNNTauFileWriter = cms.EDAnalyzer('L1CTTauFileWriter',
  collections = cms.VPSet(cms.PSet(taus = cms.InputTag("l1tNNTauProducerPuppi", "L1PFTausNN"),
                                   nTaus = cms.uint32(5))),
  nFramesPerBX = cms.uint32(9), # 360 MHz clock or 25 Gb/s link
  gapLengthOutput = cms.uint32(4),
  TMUX = cms.uint32(6),
  maxLinesPerFile = cms.uint32(1024),
  outputFilename = cms.string("L1CTTausPatterns"),
  format = cms.string("EMPv2"),
  outputFileExtension = cms.string("txt.gz")
)
