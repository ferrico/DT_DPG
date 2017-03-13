import FWCore.ParameterSet.Config as cms

process = cms.Process("DTNT")
# process = cms.Process("DTNTandRPC")
#process = cms.Process("RECLUSTERIZATION")

 
##process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Services_cff')

###process.load('Configuration/StandardSequences/Geometry_cff')  ##  Deprecated in new versions > 53X
##process.load('Configuration.Geometry.GeometryIdeal_cff')  ## In versions 7X problems with few STA events crashing the runing # but it works in the first tests of 750pre5 
##process.load('Configuration/StandardSequences/GeometryDB_cff')  
process.load('Configuration/StandardSequences/GeometryRecoDB_cff')  ##  solve STA problem
process.load('Configuration/EventContent/EventContent_cff')
process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

###process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
##process.load("Geometry.DTGeometry.dtGeometry_cfi")
##process.load("Geometry.CSCGeometry.cscGeometry_cfi")
##process.load("Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff")
##process.load("Geometry.DTGeometryBuilder.idealForDigiDtGeometry_cff")


process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("EventFilter.DTRawToDigi.dtunpackerDDUGlobal_cfi")
process.dtunpacker.readOutParameters.debug = False
process.dtunpacker.readOutParameters.rosParameters.debug = False

# for DTTF (Not used from 2016)
process.load("EventFilter.DTTFRawToDigi.dttfunpacker_cfi")
process.load("EventFilter.DTTFRawToDigi.dttfpacker_cfi")
process.dttfunpacker.DTTF_FED_Source = "rawDataCollector"

# for TWINMUX (Start to use in 2016)
process.load("EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi")

#for RAW
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
##process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
##process.load('Configuration/StandardSequences/GeometryExtended_cff')  ## For Run1
###process.load('Configuration.Geometry.GeometryExtended2016_cff')  ## 
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

##process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
##process.GlobalTag.globaltag = "GR_E_V48::All"
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")  #DB v2, at least since GR_E_V42
##process.GlobalTag.globaltag = '80X_dataRun2_Express_v0' ##for CMSSW_8_0_X X>0
##process.GlobalTag.globaltag = '80X_dataRun2_Express_v1' ##for CMSSW_8_0_X X>0

#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14' ##for CMSSW_8_0_X X>0
process.GlobalTag.globaltag = '80X_dataRun2_Express_v15' ## for express

# for the emulator
process.load("L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff")
process.load("L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi")
process.dtTriggerPrimitiveDigis.debug = False
process.L1DTConfigFromDB.debug = False

process.load('EventFilter.ScalersRawToDigi.ScalersRawToDigi_cfi')
process.load('RecoLuminosity.LumiProducer.lumiProducer_cfi')
##process.load('RecoLuminosity.LumiProducer.lumiProducer_cff')


### Copied from DQM/DTMonitorModule/python/dt_standalonedatabases_cfi.py trying to fix the problem of crashing (perhaps due to a wrong geom)
###from CondCore.DBCommon.CondDBSetup_cfi import *
###geometryESSource = cms.ESSource("PoolDBESSource",
###      CondDBSetup,
###      toGet = cms.VPSet(cms.PSet(record = cms.string('GlobalPositionRcd'),
###                       tag = cms.string('IdealGeometry')
###                       ),
###              cms.PSet(record = cms.string('DTAlignmentRcd'),
###                       tag = cms.string('DTIdealGeometry200_mc')
###                       ),
###              cms.PSet(record = cms.string('DTAlignmentErrorRcd'),
###                       tag = cms.string('DTIdealGeometryErrors200_mc')
###                       ),
###              cms.PSet(record = cms.string('CSCAlignmentRcd'),
###                       tag = cms.string('CSCIdealGeometry200_mc')
###                       ),
###              cms.PSet(record = cms.string('CSCAlignmentErrorRcd'),
###                       tag = cms.string('CSCIdealGeometryErrors200_mc')
###                       ),
###              cms.PSet(record = cms.string('TrackerAlignmentRcd'),
###                       tag = cms.string('TrackerIdealGeometry210_mc')
###                       ),
###              cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
###                       tag = cms.string('TrackerIdealGeometryErrors210_mc')
###                      )
###                      ),
###     connect = cms.string('frontier://cms_conditions_data/CMS_COND_21X_ALIGNMENT')
###)

process.load('EventFilter.L1TRawToDigi.l1tRawtoDigiBMTF_cfi')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

  fileNames = cms.untracked.vstring
  (
    ##'/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v4/000/258/159/00000/0EFE474F-D26B-E511-9618-02163E011F4B.root'
    ##'/store/express/Commissioning2016/ExpressPhysics/FEVT/Express-v1/000/270/389/00000/A0D5994E-9506-E611-A4A7-02163E0141CE.root'
    # '/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/282/917/00000/12D4E6E2-F891-E611-8EC1-02163E012A22.root'
	#'root://xrootd-cms.infn.it///store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/282/917/00000/06BA6E22-EE91-E611-87B6-02163E012AE5.root'
# 'file:/afs/cern.ch/work/f/ferrico/private/06BA6E22-EE91-E611-87B6-02163E012AE5.root',
#        '/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/282/917/00000/12D4E6E2-F891-E611-8EC1-02163E012A22.root',
#        '/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/282/917/00000/1A2F648B-E891-E611-8262-02163E011995.root',
#        '/store/data/Run2016H/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/282/917/00000/3AA76B2A-F091-E611-B97F-FA163EFD691C.root',

    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/000BE88F-ED97-E611-B962-02163E011D7E.root',
    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/0092B080-0798-E611-AD7E-02163E01184D.root',
    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00958F3E-FE97-E611-A000-02163E014255.root',
    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00B04838-E597-E611-882E-02163E0120A4.root',
    '/store/express/Run2016H/ExpressPhysics/FEVT/Express-v2/000/283/820/00000/00C5C0B0-EE97-E611-A10B-FA163EE76B26.root',

  ),
  secondaryFileNames = cms.untracked.vstring(
  )
)

#this is to select collisions
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   #cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
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
                                 tightness = cms.int32(1) # 0 = loose (e.g. unstable collisions, minimum bias, requires a DT segment)
                                                          # 1 = medium (e.g. cosmics, requires a stand alone muon)
                                                          # 2 = tight (good collisions, requires a global muon)
)


process.load("UserCode/DTDPGAnalysis/DTTTreGenerator_cfi")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
process.myDTNtuple.outputFile = "DTNtuple.root"
process.myDTNtuple.dtTrigSimDCCLabel = cms.InputTag("dtTriggerPrimitiveDigis")
process.myDTNtuple.dtDigiLabel = cms.InputTag("dtunpacker")

process.myDTNtuple.bmtfOutputDigis = cms.InputTag("BMTFStage2Digis:BMTF")
# process.myDTNtuple.bmtfInputPhDigis = cms.InputTag("BMTFStage2Digis:PhiDigis")
# process.myDTNtuple.bmtfInputThDigis = cms.InputTag("BMTFStage2Digis:TheDigis")
##process.myDTNtuple.staMuLabel = cms.InputTag("standAloneMuons")

## RPC unpacking
process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

## RPC recHit
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
process.rpcRecHits.rpcDigiLabel = cms.InputTag('rpcUnpackingModule')


## pp collisions before 2016 (DTTF - DCC)
###process.p = cms.Path(process.DTMuonSelection *  process.dtunpacker * process.dttfunpacker * process.scalersRawToDigi * process.muonDTDigis * process.dtTriggerPrimitiveDigis + process.myDTNtuple)

## pp collisions from 2016 (TM)
process.p = cms.Path(process.DTMuonSelection * process.dtunpacker * process.twinMuxStage2Digis  * process.scalersRawToDigi * process.lumiProducer * process.muonDTDigis * process.dtTriggerPrimitiveDigis + process.BMTFStage2Digis + process.rpcUnpackingModule + process.rpcRecHits + process.myDTNtuple) 

# Output
process.out = cms.OutputModule("PoolOutputModule"
                               , outputCommands = cms.untracked.vstring(
                               											"keep *",
                                                                         "keep *_*_*_testRPCTwinMuxRawToDigi"
                                                                       , "keep *_*_*_DTNTandRPC"
																		)
#                                , fileName = cms.untracked.string("file:cia.root")
                               , SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)

### Cosmics before 2016 (DTTF - DCC) 
###process.p = cms.Path(process.dtunpacker * process.dttfunpacker + process.myDTNtuple)

### Cosmics since 2016 (TM) 
### For cosmics MWGR2016 first is needed to put false the part of Global muons because the global cosmics muon tracks are not available
############ process.myDTNtuple.AnaTrackGlobalMu = cms.untracked.bool(False)    ## Comment this when the global cosmics muon tracks are available again
###process.p = cms.Path(process.dtunpacker * process.twinMuxStage2Digis + process.myDTNtuple)

## RE-RECO with CMSSW712: RECO only in the dataset!
#process.myDTNtuple.runOnRaw = cms.bool(False)
#process.p = cms.Path(process.myDTNtuple)

