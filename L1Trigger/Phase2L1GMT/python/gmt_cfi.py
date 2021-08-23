import FWCore.ParameterSet.Config as cms

gmtStubs = cms.EDProducer("Phase2L1TGMTStubProducer",
    verbose = cms.int32(0),
    srcCSC = cms.InputTag("simCscTriggerPrimitiveDigis"),
    srcDT = cms.InputTag("dtTriggerPhase2PrimitiveDigis"),
    srcDTTheta = cms.InputTag("simDtTriggerPrimitiveDigis"),
    srcRPC = cms.InputTag("simMuonRPCDigis"),
    Endcap =cms.PSet(                            
        verbose              = cms.uint32(0),
        minBX                = cms.int32(0),                           
        maxBX                = cms.int32(0),         
        coord1LSB            = cms.double(0.00076660156*32), 
        eta1LSB              = cms.double(7.68334e-04*32), 
        coord2LSB            = cms.double(0.00076660156*32), 
        eta2LSB              = cms.double(7.68334e-04*32),
        phiMatch             = cms.double(0.05),
        etaMatch             = cms.double(0.1)
    ),
    Barrel = cms.PSet(                         
        verbose            = cms.int32(0),
        minPhiQuality      = cms.int32(0),#was 5
        minThetaQuality    = cms.int32(0),
        minBX              = cms.int32(0),                           
        maxBX              = cms.int32(0),                           
        phiLSB             = cms.double(0.00076660156*32),
        phiBDivider        = cms.int32(16),
        etaLSB             = cms.double(7.68334e-04*32), 
        eta_1              = cms.vint32(-1503/32,-1446/32,-1387/32,-1327/32,-1266/32,-1194/32,-1125/32,-985/32,-916/32,-839/32,-752/32,-670/32,-582/32,-489/32,-315/32,-213/32,-115/32,-49/32,49/32, 115/32, 213/32, 315/32, 489/32, 582/32, 670/32, 752/32, 839/32, 916/32, 985/32, 1125/32, 1194/32, 1266/32, 1327/32, 1387/32, 1446/32, 1503),
        eta_2              = cms.vint32(-1334/32,-1279/32,-1227/32,-1168/32,-1109/32,-1044/32,-982/32,-861/32,-793/32,-720/32,-648/32,-577/32,-496/32,-425/32,-268/32,-185/32,-97/32,-51/32,51/32, 97/32, 185/32, 268/32, 425/32, 496/32, 577/32, 648/32, 720/32, 793/32, 861/32, 982/32, 1044/32, 1109/32, 1168/32, 1227/32, 1279/32, 1334),
        eta_3              = cms.vint32(-1148/32,-1110/32,-1051/32,-1004/32,-947/32,-895/32,-839/32,-728/32,-668/32,-608/32,-546/32,-485/32,-425/32,-366/32,-222/32,-155/32,-87/32,-40/32,40/32, 87/32, 155/32, 222/32, 366/32, 425/32, 485/32, 546/32, 608/32, 668/32, 728/32, 839/32, 895/32, 947/32, 1004/32, 1051/32, 1110/32, 1148),
        coarseEta_1        = cms.vint32(0/32,758/32,1336/32),
        coarseEta_2        = cms.vint32(0/32,653/32,1168/32),
        coarseEta_3        = cms.vint32(0/32,552/32,1001/32),
        coarseEta_4        = cms.vint32(0/32,478/32,878/32),
        phiOffset          = cms.vint32(33/32,-8/32,+14/32,0)    
   )

)





gmtMuons = cms.EDProducer('Phase2L1TGMTProducer',
                     srcTracks = cms.InputTag("TTTracksFromTrackletEmulation:Level1TTTracks"),
                     srcStubs  = cms.InputTag('gmtStubs'),
                     srcBMTF   = cms.InputTag('simBmtfDigis','BMTF'),
                     srcEMTF   = cms.InputTag('simEmtfDigis','EMTF'),
                     srcOMTF   = cms.InputTag('simOmtfDigis','OMTF'),
                     minTrackStubs = cms.int32(4),     
                     muonBXMin = cms.int32(0),
                     muonBXMax = cms.int32(0),
                          verbose   = cms.int32(0),     
                     trackConverter  = cms.PSet(
                         verbose = cms.int32(0)
                     ),
                     roiTrackAssociator  = cms.PSet(
                         verbose=cms.int32(0)
                     ),
                     trackMatching  = cms.PSet(
                         verbose=cms.int32(0)
                     ),
                     isolation  = cms.PSet(
                       AbsIsoThresholdL = cms.int32(160),
                       AbsIsoThresholdM = cms.int32(120),
                       AbsIsoThresholdT = cms.int32(80),
                       RelIsoThresholdL = cms.double(0.1),
                       RelIsoThresholdM = cms.double(0.05),
                       RelIsoThresholdT = cms.double(0.01),
                       verbose       = cms.int32(0),
                       IsodumpForHLS = cms.int32(1),
                     ),
                    tauto3mu = cms.PSet()

)

standaloneMuons = cms.EDProducer('Phase2L1TGMTSAMuonProducer',
                                 muonToken  = cms.InputTag('simGmtStage2Digis'),
                                 Nprompt    = cms.uint32(12),
                                 Ndisplaced = cms.uint32(12)
                                )
