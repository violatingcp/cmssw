import FWCore.ParameterSet.Config as cms
import os,sys

sys.path = sys.path + [os.path.expandvars("$CMSSW_BASE/src/HLTrigger/Configuration/test/"), os.path.expandvars("$CMSSW_RELEASE_BASE/src/HLTrigger/Configuration/test/")]

from OnLine_HLT_GRun import process
process.maxEvents.input = cms.untracked.int32(10000)

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '120X_mcRun3_2021_realistic_v2')

process.source.fileNames = cms.untracked.vstring("/store/relval/CMSSW_11_2_0_pre7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_112X_mcRun3_2021_realistic_v8-v1/20000/08FB950E-6A55-374B-AA96-C43C992B55AD.root")#/store/relval/CMSSW_11_2_0_pre6_ROOT622/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/112X_mcRun3_2021_realistic_v7-v1/20000/FED4709C-569E-0A42-8FF7-20E565ABE999.root")


keepMsgs = ['FastTimerService']
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

