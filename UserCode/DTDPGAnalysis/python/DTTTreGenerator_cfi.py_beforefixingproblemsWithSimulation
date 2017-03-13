import FWCore.ParameterSet.Config as cms

myDTNtuple = cms.EDAnalyzer('TTreeGenerator',
                            outputFile =cms.string("DTTree.root"),
                            dtDigiLabel = cms.InputTag("muonDTDigis"),
                            dtSegmentLabel = cms.InputTag("dt4DSegments"),
                            cscSegmentLabel = cms.InputTag("cscSegments"),
                            dtTrigTwinMuxInLabel= cms.InputTag("twinMuxStage2Digis","PhIn"),
                            dtTrigTwinMuxOutLabel = cms.InputTag("twinMuxStage2Digis","PhOut"),
                            dtTrigTwinMuxThLabel= cms.InputTag("twinMuxStage2Digis","ThIn"),
                            staMuLabel     = cms.InputTag("muons"),
                            gmtLabel     = cms.InputTag("gtDigis"),  # legacy
                            gtLabel      = cms.InputTag("gtDigis"),  # legacy
                            rpcRecHitLabel = cms.InputTag("rpcRecHits"),
                            dtDigiSize = cms.int32(300),
                            dtSegmentSize = cms.int32(50),
                            cscSegmentSize = cms.int32(50),
                            dtTrigTwinMuxInSize = cms.int32(50),
                            dtTrigTwinMuxOutSize = cms.int32(50),
                            dtTrigTwinMuxThSize = cms.int32(50),
                            gmtSize       = cms.int32(50), # legacy
                            STAMuSize     = cms.int32(20),
                            rpcRecHitSize = cms.int32(300),
                            PrimaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                            TriggerTag = cms.InputTag("TriggerResults::HLT"),
                            beamSpotTag = cms.InputTag("offlineBeamSpot"),
                            scalersResults = cms.InputTag("scalersRawToDigi"),
                            lumiInputTag   = cms.InputTag("lumiProducer"),
                            runOnRaw = cms.bool(True),
                            runOnSimulation = cms.bool(False),
                            localDTmuons    = cms.bool(False),

           ######## Parameters for retrieving the ttrig to correct the recHit times
                            tTrigModeConfig = cms.untracked.PSet(
                                  vPropWire = cms.double(24.4),
                                  doTOFCorrection = cms.bool(False),
                                  ##tofCorrType = cms.int32(2),  ## old
                                  tofCorrType = cms.int32(2),
                                  ##kFactor = cms.double(-1.3),
                                  wirePropCorrType = cms.int32(0),
                                  doWirePropCorrection = cms.bool(False),
                                  doT0Correction = cms.bool(True),
                                  tTrigLabel = cms.string(''),
                                  debug = cms.untracked.bool(False)
                              ),
                              tTrigMode = cms.untracked.string('DTTTrigSyncFromDB')
           ######## END Parameters for retrieving the ttrig to correct the recHit times

                              
)
