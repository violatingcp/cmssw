import FWCore.ParameterSet.Config as cms
import os,sys

sys.path = sys.path + [os.path.expandvars("$CMSSW_BASE/src/HLTrigger/Configuration/test/"), os.path.expandvars("$CMSSW_RELEASE_BASE/src/HLTrigger/Configuration/test/")]

from OnLine_HLT_GRun import process

process.hltHbherecopre = process.hltHbhereco.clone(
    makeRecHits = cms.bool(False),
    saveInfos = cms.bool(True),
)
process.load("HeterogeneousCore.SonicTriton.TritonService_cff")

process.TritonService.verbose = True
process.TritonService.fallback.verbose = True
process.TritonService.servers.append(
    cms.PSet(
            name = cms.untracked.string("default"),
            address = cms.untracked.string("ailab01.fnal.gov"),
            port = cms.untracked.uint32(8001),
            useSsl = cms.untracked.bool(False),
            rootCertificates = cms.untracked.string(""),
            privateKey = cms.untracked.string(""),
            certificateChain = cms.untracked.string(""),
    )
)

from RecoLocalCalo.HcalRecProducers.facileHcalReconstructor_cfi import sonic_hbheprereco
process.hltHbhereco = sonic_hbheprereco.clone(
    ChannelInfoName = cms.InputTag("hltHbherecopre")
)

process.maxEvents.input = cms.untracked.int32(100)

process.HLTDoLocalHcalSequence = cms.Sequence( process.hltHcalDigis + process.hltHbherecopre + process.hltHbhereco + process.hltHfprereco + process.hltHfreco + process.hltHoreco )
process.HLTStoppedHSCPLocalHcalReco = cms.Sequence( process.hltHcalDigis + process.hltHbherecopre + process.hltHbhereco)

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '120X_mcRun3_2021_realistic_v2')

process.source.fileNames = cms.untracked.vstring("/store/relval/CMSSW_11_2_0_pre6_ROOT622/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/112X_mcRun3_2021_realistic_v7-v1/20000/FED4709C-569E-0A42-8FF7-20E565ABE999.root")

keepMsgs = ['TritonClient','TritonService','FastTimerService']
for producer in process._Process__producers.values():
    if hasattr(producer,'Client'):
        if hasattr(producer.Client,'verbose'):
            producer.Client.verbose = True 
            keepMsgs.extend([producer._TypedParameterizable__type,producer._TypedParameterizable__type+":TritonClient"])
        if hasattr(producer.Client,'compression'):
            producer.Client.compression = True 
        if hasattr(producer.Client,'useSharedMemory'):
            producer.Client.useSharedMemory = True

if 1:
    process.load('FWCore/MessageService/MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 500
    for msg in keepMsgs:
        setattr(process.MessageLogger.cerr,msg,
            cms.untracked.PSet(
                limit = cms.untracked.int32(10000000),
            )
        )

# remove any instance of the FastTimerService
if 'FastTimerService' in process.__dict__:
    del process.FastTimerService

# instrument the menu with the FastTimerService
process.load( "HLTrigger.Timer.FastTimerService_cfi" )

# print a text summary at the end of the job
process.FastTimerService.printEventSummary        = True
process.FastTimerService.printRunSummary          = True
process.FastTimerService.printJobSummary          = True

