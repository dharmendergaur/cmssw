import FWCore.ParameterSet.Config as cms

l1tSeedConePFJetProducer = cms.EDProducer("L1SeedConePFJetProducer",
                           L1PFObjects = cms.InputTag("l1tLayer1","Puppi"),
                           nJets       = cms.uint32(16),
                           coneSize    = cms.double(0.4),
                           HW          = cms.bool(False),
                           debug       = cms.bool(False),
                           doCorrections = cms.bool(False),
                           useExternalSeeds = cms.bool(False),
                           allowDoubleCounting = cms.bool(False),
                           JetSeeds = cms.InputTag("l1tPhase1JetSeedProducer9x9trimmed","histoJetSeeds9x9trimmed"),
                           correctorFile = cms.string(""),
                           correctorDir  = cms.string("")
                         )

l1tSeedConePFJetEmulatorProducer = l1tSeedConePFJetProducer.clone(HW = True)

