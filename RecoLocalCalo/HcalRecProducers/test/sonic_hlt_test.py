import FWCore.ParameterSet.Config as cms
import os,sys

sys.path = sys.path + [os.path.expandvars("$CMSSW_BASE/src/HLTrigger/Configuration/test/"), os.path.expandvars("$CMSSW_RELEASE_BASE/src/HLTrigger/Configuration/test/")]

from OnLine_HLT_GRun import process
process.maxEvents.input = cms.untracked.int32(1000)

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '120X_mcRun3_2021_realistic_v2')
process.hltHbherecopre = process.hltHbhereco.clone(
    makeRecHits = cms.bool(True),
    saveInfos = cms.bool(True),
    #saveDroppedInfos = cms.bool(True)
    saveInfosVector = cms.bool(True),
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

process.source.fileNames = cms.untracked.vstring(
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/8BD20F29-96F9-7C44-9078-E641186F0B19.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/53FC7D98-3F75-8148-81B2-9E866BA325BC.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/44D71998-2DC0-6744-AF35-C8781E72499E.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/08FB950E-6A55-374B-AA96-C43C992B55AD.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/D60BECD5-C9B7-054D-8569-1E71ADE201DA.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/7D870A3B-64E2-CB40-94DA-434B1F4F2E41.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/95F73111-5194-2F45-9116-C91EDBFD1D82.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/FF9081F1-5F25-9F41-9BA9-44F447E24BA3.root",
"/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/4D5F2C2B-331A-554E-BDC6-A1826A3AF6C9.root",
)
'''
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'keep *_*_*_*',
        'keep *_hltHbhereco_*_*', 
        'keep *_hltAK8PFJetsCorrected_*_*',
        'keep *_hltAK4PFJetsCorrected_*_*',
        'keep *_ak4GenJets_*_*',
        'keep *_ak8GenJets_*_*',
        'keep *_hltMuons_*_*',
        'keep *_genMetTrue_*_*',
        'keep *_hltFastPrimaryVertex_*_*',
 
    ),
    fileName = cms.untracked.string("facile_rechits.root")
)
'''
#process.finalize = cms.EndPath(process.out)

process.HLTDoLocalHcalSequence = cms.Sequence( process.hltHcalDigis + process.hltHbherecopre + process.hltHbhereco + process.hltHfprereco + process.hltHfreco + process.hltHoreco )
process.HLTStoppedHSCPLocalHcalReco = cms.Sequence( process.hltHcalDigis + process.hltHbherecopre + process.hltHbhereco)
print("test",process.HLTDoLocalHcalSequence)
#process.HLTSchedule.append(process.finalize)



#keepMsgs = [] # ['TritonClient','TritonService','FastTimerService']
#for producer in process._Process__producers.values():
#    if hasattr(producer,'Client'):
#        if hasattr(producer.Client,'verbose'):
#            producer.Client.verbose = True 
#            keepMsgs.extend([producer._TypedParameterizable__type,producer._TypedParameterizable__type+":TritonClient"])
#        if hasattr(producer.Client,'compression'):
#            producer.Client.compression = True 
#        if hasattr(producer.Client,'useSharedMemory'):
#            producer.Client.useSharedMemory = True

keepMsgs = ['FastTimerService','FastReport',]
if 1:
    process.load('FWCore/MessageService/MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 10
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
#process.FastTimerService.printEventSummary        = True
#process.FastTimerService.printRunSummary          = True
#process.FastTimerService.printJobSummary          = True

