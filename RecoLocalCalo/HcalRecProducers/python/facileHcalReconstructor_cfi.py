import FWCore.ParameterSet.Config as cms

sonic_hbheprereco = cms.EDProducer("FacileHcalReconstructor",
    Client = cms.PSet(
        #batchSize = cms.untracked.uint32(16000),
        #port = cms.untracked.uint32(8001),
        mode = cms.string("Async"),
        preferredServer = cms.untracked.string(""),
        timeout = cms.untracked.uint32(300),
        modelName = cms.string("facile_tf"),
        modelVersion = cms.string(""),
        modelConfigPath = cms.FileInPath("HeterogeneousCore/SonicTriton/data/models/facile_plan_10k/config.pbtxt"),        
        verbose = cms.untracked.bool(True),
        allowedTries = cms.untracked.uint32(0),
        #outputs = cms.untracked.vstring("Identity:0"),
    ),
    ChannelInfoName = cms.InputTag("hbhechannelinfo")
)
