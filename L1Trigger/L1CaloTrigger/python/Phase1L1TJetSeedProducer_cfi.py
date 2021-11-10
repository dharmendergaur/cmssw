import FWCore.ParameterSet.Config as cms
from math import pi

from Phase1L1TJets_sincosLUT_cff import sinPhi, cosPhi

caloEtaSegmentation = cms.vdouble(
-3 , -2.91491519955 , -2.83201206065 , -2.74910892175 , -2.66620578285 , -2.58330264395 ,
-2.5 , -2.41491519955 , -2.33201206065 , -2.24910892175 , -2.16620578285 , -2.08330264395 , -2.00039950505 , -1.91749636615 , -1.83459322725 , -1.75169008835 , -1.66878694945 , -1.58588381055 ,
-1.5 , -1.41491519955 , -1.33201206065 , -1.24910892175 , -1.16620578285 , -1.08330264395 ,
-1.0 , -0.91491519955 , -0.83201206065 , -0.74910892175 , -0.66620578285 , -0.58330264395 ,
-0.5 , -0.41491519955 , -0.33201206065 , -0.24910892175 , -0.16620578285 , -0.08330264395 ,
0 , 0.08508480045 , 0.16798793935 , 0.25089107825 , 0.33379421715 , 0.41669735605 ,
0.5 , 0.58508480045 , 0.66798793935 , 0.75089107825 , 0.83379421715 , 0.91669735605 ,
1 , 1.08508480045 , 1.16798793935 , 1.25089107825 , 1.33379421715 , 1.41669735605 ,
1.5 , 1.58508480045 , 1.66798793935 , 1.75089107825 , 1.83379421715 , 1.91669735605 , 1.99960049495 , 2.08250363385 , 2.16540677275 , 2.24830991165 , 2.33121305055 , 2.41411618945 ,
2.5 , 2.58508480045 , 2.66798793935 , 2.75089107825 , 2.83379421715 , 2.91669735605 , 3
)

Phase1L1TJetSeedProducer9x9trimmed = cms.EDProducer('Phase1L1TJetSeedProducer',
  inputCollectionTag = cms.InputTag("l1ctLayer1", "Puppi"),
  etaBinning = caloEtaSegmentation,
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.15),
  phiUp = cms.double(3.15),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  trimmedGrid = cms.bool(True),
  seedPtThreshold = cms.double(5), # GeV
  ptlsb = cms.double(0.25),
  philsb = cms.double(0.0043633231),
  etalsb = cms.double(0.0043633231),
  outputCollectionName = cms.string("histoJetSeeds9x9trimmed"),
  etaRegions = cms.vdouble( -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3 ),
  phiRegions = cms.vdouble( -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 ),#, 4.2, 4.9, 5.6, 6.3 ),
  maxInputsPerRegion = cms.uint32( 18 )
)
