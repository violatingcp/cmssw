import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia8GeneratorFilter",
    crossSection = cms.untracked.double(71.39e+09),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring(
'Main:timesAllowErrors = 10000',
            'ParticleDecays:limitTau0 = on',
'ParticleDecays:tauMax = 10',
            'SoftQCD:nonDiffractive = on',
            'SoftQCD:singleDiffractive = on',
            'SoftQCD:doubleDiffractive = on',
            'Tune:pp 2',
            'Tune:ee 3'),
        parameterSets = cms.vstring('processParameters')
    )
)

ProductionFilterSequence = cms.Sequence(generator)

