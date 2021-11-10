import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer9x9trimmed
#9x9 sample with 1x1 seed size
l1tHistoSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone( useExternalSeeds = True,
                                                                                   allowDoubleCounting = False,
                                                                                  JetSeeds = ('l1tPhase1JetSeedProducer9x9trimmed', 'histoJetSeeds9x9trimmed'),
                                                                                  nJets = 12,
                                                                                  L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
                                                                                  doCorrections = cms.bool(False),
                                                                                  correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                                  correctorDir = cms.string('L1PuppiSC4EmuJets')
                                                                                  )
#task for sample (9x9)
L1TPFHistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer9x9trimmed, l1tLayer2Deregionizer, l1tHistoSeedsSCPFL1PuppiCorrectedEmulator
)
#Double counting allowed
l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC = l1tSeedConePFJetEmulatorProducer.clone( useExternalSeeds = True,
                                                                                   allowDoubleCounting = True,
                                                                                  JetSeeds = ('l1tPhase1JetSeedProducer9x9trimmed', 'histoJetSeeds9x9trimmed'),
                                                                                  nJets = 12,
                                                                                  L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
                                                                                  doCorrections = cms.bool(False),
                                                                                  correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                                  correctorDir = cms.string('L1PuppiSC4EmuJets')
                                                                                  )
#task for sample (9x9) with Double counting
L1TPFHistoSeedJetsTaskDC = cms.Task(
    l1tPhase1JetSeedProducer9x9trimmed, l1tLayer2Deregionizer, l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC
)
#----------------------------------------------------------------------------------------------
#for 9x9  (trimmed)
l1tPhase1JetSeedProducer9x9 = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 9,
                                                                        jetIPhiSize = 9,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds9x9"
                                                                        )

l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tHistoSeedsSCPFL1PuppiCorrectedEmulator.clone( JetSeeds = ('l1tPhase1JetSeedProducer9x9', 'histoJetSeeds9x9') )

#Task for 9x9 (trimmed)
L1TPF9x9HistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer9x9, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator
)

#for 9x9 with double counting (trimmed)
l1tPhase1JetSeedProducer9x9DC = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 9,
                                                                        jetIPhiSize = 9,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds9x9DC"
                                                                        )

l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDC = l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC.clone( JetSeeds = ('l1tPhase1JetSeedProducer9x9DC', 'histoJetSeeds9x9DC') )

#Task for 9x9 with double counting (trimmed)
L1TPF9x9HistoSeedJetsTaskDC = cms.Task(
    l1tPhase1JetSeedProducer9x9DC, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDC
)

#----------------------------------------------------------------------------------------------
#for 7x7(trimmed)
l1tPhase1JetSeedProducer7x7T = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 7,
                                                                        jetIPhiSize = 7,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds7x7T"
                                                                        )

l1t7x7THistoSeedsSCPFL1PuppiCorrectedEmulator = l1tHistoSeedsSCPFL1PuppiCorrectedEmulator.clone( JetSeeds = ('l1tPhase1JetSeedProducer7x7T', 'histoJetSeeds7x7T') )

#Task for 7x7(trimmed)
L1TPF7x7THistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer7x7T, l1tLayer2Deregionizer, l1t7x7THistoSeedsSCPFL1PuppiCorrectedEmulator
)

#for 7x7T with double counting(trimmed)
l1tPhase1JetSeedProducer7x7TDC = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 7,
                                                                        jetIPhiSize = 7,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds7x7TDC"
                                                                        )

l1t7x7THistoSeedsSCPFL1PuppiCorrectedEmulatorDC = l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC.clone( JetSeeds = ('l1tPhase1JetSeedProducer7x7TDC', 'histoJetSeeds7x7TDC') )

#Task for 7x7T with double counting(trimmed)
L1TPF7x7THistoSeedJetsTaskDC = cms.Task(
    l1tPhase1JetSeedProducer7x7TDC, l1tLayer2Deregionizer, l1t7x7THistoSeedsSCPFL1PuppiCorrectedEmulatorDC
)

#----------------------------------------------------------------------------------------------
#for 7x7
l1tPhase1JetSeedProducer7x7 = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 7,
                                                                        jetIPhiSize = 7,
                                                                        trimmedGrid = False,
                                                                        outputCollectionName = "histoJetSeeds7x7"
                                                                        )

l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tHistoSeedsSCPFL1PuppiCorrectedEmulator.clone( JetSeeds = ('l1tPhase1JetSeedProducer7x7', 'histoJetSeeds7x7') )

#Task for 7x7
L1TPF7x7HistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer7x7, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator
)

#for 7x7 with double counting
l1tPhase1JetSeedProducer7x7DC = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 7,
                                                                        jetIPhiSize = 7,
                                                                        trimmedGrid = False,
                                                                        outputCollectionName = "histoJetSeeds7x7DC"
                                                                        )

l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulatorDC = l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC.clone( JetSeeds = ('l1tPhase1JetSeedProducer7x7DC', 'histoJetSeeds7x7DC') )

#Task for 7x7 with double counting
L1TPF7x7HistoSeedJetsTaskDC = cms.Task(
    l1tPhase1JetSeedProducer7x7DC, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulatorDC
)
#----------------------------------------------------------------------------------------------
#for 5x5 
l1tPhase1JetSeedProducer5x5 = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 5,
                                                                        jetIPhiSize = 5,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds5x5"
                                                                        )

l1t5x5HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tHistoSeedsSCPFL1PuppiCorrectedEmulator.clone( JetSeeds = ('l1tPhase1JetSeedProducer5x5', 'histoJetSeeds5x5') )
#Task for 5x5
L1TPF5x5HistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer5x5, l1tLayer2Deregionizer, l1t5x5HistoSeedsSCPFL1PuppiCorrectedEmulator
)
#for 5x5 with double counting
l1tPhase1JetSeedProducer5x5DC = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 5,
                                                                        jetIPhiSize = 5,
                                                                        trimmedGrid = True,
                                                                        outputCollectionName = "histoJetSeeds5x5DC"
                                                                        )

l1t5x5HistoSeedsSCPFL1PuppiCorrectedEmulatorDC = l1tHistoSeedsSCPFL1PuppiCorrectedEmulatorDC.clone( JetSeeds = ('l1tPhase1JetSeedProducer5x5DC', 'histoJetSeeds5x5DC') )
#Task for 5x5 with DC
L1TPF5x5HistoSeedJetsTaskDC = cms.Task(
    l1tPhase1JetSeedProducer5x5DC, l1tLayer2Deregionizer, l1t5x5HistoSeedsSCPFL1PuppiCorrectedEmulatorDC
)
