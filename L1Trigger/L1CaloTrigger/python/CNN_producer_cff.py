import FWCore.ParameterSet.Config as cms
from math import pi


CNNProducer = cms.EDProducer('CNNProducer',
  inputCollectionTag = cms.InputTag("l1ctLayer1", "Puppi"),
  phiLow = cms.double(-pi),
  phiUp = cms.double(pi),
  etaLow = cms.double(-3),
  etaUp = cms.double(3),
  ptlsb = cms.double(0.25),
  philsb = cms.double(0.0043633231),
  etalsb = cms.double(0.0043633231),
  outputCollectionName = cms.string("Keras"),
  etaRegions = cms.vdouble( -3., -2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5, 3. ),
  phiRegions = cms.vdouble( -3.5, -2.8, -2.1, -1.4, -0.7, 0, 0.7, 1.4, 2.1, 2.8, 3.5 ),#, 4.2, 4.9, 5.6, 6.3 ),
  maxInputsPerRegion = cms.uint32( 18 ),
  runInference = cms.bool(True),
  cnnGraph = cms.FileInPath('L1Trigger/L1CaloTrigger/data/Qcnn_graph_bigNetwork.pb'),
  cnnImageSizePhi = cms.uint32(12),
  cnnImageSizeEta = cms.uint32(12),
  cnnNPadPhi = cms.uint32(3),
  cnnNPadEta = cms.uint32(3),
  maxPixelValue = cms.double(63),

  writeImages = cms.bool(True),
  etaNBins_writeToFile = cms.uint32(120),
  etaMin_writeToFile = cms.double(-5),
  etaMax_writeToFile = cms.double(5),
  phiNBins_writeToFile = cms.uint32(72),
  phiMin_writeToFile = cms.double(-pi),
  phiMax_writeToFile = cms.double(pi),
  imagesFileName = cms.string("images.root")
)


CNNProducerSequence = cms.Sequence(
  CNNProducer 
)