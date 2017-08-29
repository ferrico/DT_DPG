import FWCore.ParameterSet.Config as cms

process = cms.Process("DTNT")
# process = cms.Process("DTNTandRPC")
#process = cms.Process("RECLUSTERIZATION")

 
##process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Services_cff')

process.load('Configuration/StandardSequences/GeometryRecoDB_cff')  ##  solve STA problem
process.load('Configuration/EventContent/EventContent_cff')
process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

process.load("EventFilter.DTRawToDigi.dtunpackerDDUGlobal_cfi")
process.dtunpacker.readOutParameters.debug = False
process.dtunpacker.readOutParameters.rosParameters.debug = False

# for DTTF (Not used from 2016)
process.load("EventFilter.DTTFRawToDigi.dttfunpacker_cfi")
process.load("EventFilter.DTTFRawToDigi.dttfpacker_cfi")

# for TWINMUX (Start to use in 2016)
process.load("EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi")

#for RAW
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")  #DB v2, at least since GR_E_V42

process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v9'

# for the emulator
process.load("L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff")
process.load("L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi")
process.dtTriggerPrimitiveDigis.debug = False
process.L1DTConfigFromDB.debug = False

process.load('EventFilter.ScalersRawToDigi.ScalersRawToDigi_cfi')
process.load('RecoLuminosity.LumiProducer.lumiProducer_cfi')

# process.load('EventFilter.L1TRawToDigi.l1tRawtoDigiBMTF_cfi')
process.load('EventFilter.L1TRawToDigi.bmtfDigis_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

  fileNames = cms.untracked.vstring
  (
'/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/283/830/00000/001637C2-489C-E611-A108-02163E0141B5.root',
'/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/283/830/00000/04439A9C-489C-E611-B3F7-02163E01254D.root',

#    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/000BE88F-ED97-E611-B962-02163E011D7E.root',
#    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/0092B080-0798-E611-AD7E-02163E01184D.root',
#    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00958F3E-FE97-E611-A000-02163E014255.root',
#    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00B04838-E597-E611-882E-02163E0120A4.root',
#    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00C5C0B0-EE97-E611-A10B-FA163EE76B26.root',

  ),
  secondaryFileNames = cms.untracked.vstring(
  )
)

#this is to select collisions
process.load("DPGAnalysis.Skims.ZMuSkim_cff")
process.looseMuonsForZMuSkim.cut = cms.string('pt > 20 && abs(eta)<2.4 && isGlobalMuon = 1 && isTrackerMuon = 1 && abs(innerTrack().dxy) < 2.0')
process.dimuonsZMuSkim.cut = cms.string('mass > 60')
process.dimuonsZMuSkim.decay = cms.string("tightMuonsForZMuSkim@+ looseMuonsForZMuSkim@-")

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4"), # && abs(z) <= 15 && position.Rho <= 2" # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
)

process.DTMuonSelection = cms.EDFilter("DTMuonSelection",
                                 src = cms.InputTag('muons'),
                                 Muons = cms.InputTag('muons'),
                                 SAMuons = cms.InputTag('standAloneMuons'),
                                 dtSegmentLabel = cms.InputTag('dt4DSegments'),
                                 etaMin = cms.double(-1.25),
                                 etaMax = cms.double(1.25),
                                 ptMin = cms.double(0.),#3.),
                                 tightness = cms.int32(2) # 0 = loose (e.g. unstable collisions, minimum bias, requires a DT segment)
                                                          # 1 = medium (e.g. cosmics, requires a stand alone muon)
                                                          # 2 = tight (good collisions, requires a global muon)
)


process.load("UserCode/DTDPGAnalysis/DTTTreGenerator_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.myDTNtuple.dtExtrapolation = cms.bool(True) #Set False if DT extrapolation on RPC layer is not needed
process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
process.myDTNtuple.outputFile = "DTNtuple.root"
process.myDTNtuple.dtTrigSimDCCLabel = cms.InputTag("dtTriggerPrimitiveDigis")
process.myDTNtuple.dtDigiLabel = cms.InputTag("dtunpacker")

process.myDTNtuple.bmtfOutputDigis = cms.InputTag("bmtfDigis:BMTF")

## RPC unpacking
process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

## RPC recHit
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
process.rpcRecHits.rpcDigiLabel = cms.InputTag('rpcUnpackingModule')


process.p = cms.Path(process.diMuonSelSeq * 
                     process.dtunpacker * 
                     process.twinMuxStage2Digis  * 
                     process.scalersRawToDigi * 
                     process.lumiProducer * 
                     process.dtTriggerPrimitiveDigis + 
                     process.bmtfDigis + 
                     process.rpcUnpackingModule + 
                     process.rpcRecHits + 
process.myDTNtuple)

