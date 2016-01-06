import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

puppiPhoton = cms.EDProducer("PuppiPhoton",
                             candName       = cms.InputTag('packedPFCandidates'),
                             photonName     = cms.InputTag('slimmedPhotons'),
                             weightsName    = cms.InputTag('puppi'),
                             photonId       = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                             pt             = cms.untracked.double(10),
                             useRefs        = cms.untracked.bool(True),
                             dRMatch        = cms.untracked.vdouble(0.03,10, 10),
                             pdgids         = cms.untracked.vint32 ( 22,11,211)
                             )


def setupPuppiPhoton(process):
    my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


#puppiPhotonSeq = cms.Sequence(egmPhotonIDSequence*puppiPhoton)
