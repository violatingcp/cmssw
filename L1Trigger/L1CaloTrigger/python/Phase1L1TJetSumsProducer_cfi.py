import FWCore.ParameterSet.Config as cms
from math import pi

sinPhi = cms.vdouble(
  -0.0353352962792, -0.122533930843, -0.208795013406, -0.293458528818, -0.375876685504, -0.455418871948, -0.531476481737, -0.603467570232, -0.670841307236, -0.733082191603, -0.789713995522, -0.840303408309, -0.884463351833, -0.921855942186, -0.952195074957, -0.975248614326, -0.990840169216, -0.998850442928, -0.999218145922, -0.99194046477, -0.977073083675, -0.954729758418, -0.925081445966, -0.888354996422, -0.844831417308, -0.794843723474, -0.738774389082, -0.677052421152, -0.610150077076, -0.538579251202, -0.462887558141, -0.383654142772, -0.301485248985, -0.217009581095, -0.130873493387, -0.0437360446299, 0.0437360446299, 0.130873493387, 0.217009581095, 0.301485248985, 0.383654142772, 0.462887558141, 0.538579251202, 0.610150077076, 0.677052421152, 0.738774389082, 0.794843723474, 0.844831417308, 0.888354996422, 0.925081445966, 0.954729758418, 0.977073083675, 0.99194046477, 0.999218145922, 0.998850442928, 0.990840169216, 0.975248614326, 0.952195074957, 0.921855942186, 0.884463351833, 0.840303408309, 0.789713995522, 0.733082191603, 0.670841307236, 0.603467570232, 0.531476481737, 0.455418871948, 0.375876685504, 0.293458528818, 0.208795013406, 0.122533930843, 0.0353352962792 )

cosPhi = cms.vdouble( -0.999375513427, -0.992464324695, -0.977959427777, -0.955971804952, -0.926669691581, -0.890277288868, -0.847073048421, -0.797387541713, -0.741600930761, -0.680140059366, -0.613475187173, -0.542116391547, -0.466609664777, -0.387532736497, -0.305490653258, -0.2211111491, -0.135039842524, -0.0479352966351, 0.039536019772, 0.126704831606, 0.212904178348, 0.297474517214, 0.379768769555, 0.459157271892, 0.535032593708, 0.606814185113, 0.673952818851, 0.735934792636, 0.792285859677, 0.842574857312, 0.886417005995, 0.923476853383, 0.953470841004, 0.976169473869, 0.991399076421, 0.999043121392, 0.999043121392, 0.991399076421, 0.976169473869, 0.953470841004, 0.923476853383, 0.886417005995, 0.842574857312, 0.792285859677, 0.735934792636, 0.673952818851, 0.606814185113, 0.535032593708, 0.459157271892, 0.379768769555, 0.297474517214, 0.212904178348, 0.126704831606, 0.039536019772, -0.0479352966351, -0.135039842524, -0.2211111491, -0.305490653258, -0.387532736497, -0.466609664777, -0.542116391547, -0.613475187173, -0.680140059366, -0.741600930761, -0.797387541713, -0.847073048421, -0.890277288868, -0.926669691581, -0.955971804952, -0.977959427777, -0.992464324695, -0.999375513427 )


Phase1L1TJetSumsProducer = cms.EDProducer('Phase1L1TJetSumsProducer',
  inputJetCollectionTag = cms.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates"),
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-pi),
  phiUp = cms.double(pi),
  etaLow = cms.double(-5),
  etaUp = cms.double(5),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  htPtThreshold = cms.double(30),
  htAbsEtaCut = cms.double(2.4),
  mhtPtThreshold = cms.double(30),
  mhtAbsEtaCut = cms.double(2.4),
  ptlsb = cms.double(0.25),
  outputCollectionName = cms.string("Sums"),
)
