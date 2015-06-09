import FWCore.ParameterSet.Config as cms

process = cms.Process('TestPUMods')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.GlobalTag.globaltag = 'START53_V7G::All'

process.load('CommonTools/PileupAlgos/Puppi_cff')
process.load('CommonTools/PileupAlgos/PUPuppi_cff')
process.load('CommonTools/PileupAlgos/softKiller_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.source = cms.Source("PoolSource",
                            fileNames  = cms.untracked.vstring('/store/relval/CMSSW_7_4_1/RelValTTbar_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9_gensim71X-v1/00000/24367CFB-9BEC-E411-90A3-0025905A60C6.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)


process.puSequence = cms.Sequence(process.puppi*process.softKiller*process.pupuppi)
process.p = cms.Path(process.puSequence)
process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_particleFlow_*_*',
                                                                         'keep *_*_*_TestPUMods'),
                                  fileName       = cms.untracked.string ("Output.root")
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.output) 
