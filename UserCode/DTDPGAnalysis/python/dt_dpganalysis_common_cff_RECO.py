from EventFilter.DTTFRawToDigi.dttfunpacker_cfi import *
dttfunpacker.DTTF_FED_Source = 'rawDataCollector'  ## MWGR Feb12 

##dtunpacker = cms.EDFilter("DTUnpackingModule",
dtunpacker = cms.EDProducer("DTUnpackingModule",
                          dataType = cms.string('DDU'),
                          ###useStandardFEDid = cms.untracked.bool(True),
                          useStandardFEDid = cms.bool(True),
                          inputLabel = cms.InputTag('rawDataCollector'), ## MWGR Feb12
                          ##fedbyType = cms.untracked.bool(True),
                          fedbyType = cms.bool(True),
                          dqmOnly = cms.bool(False),   ## needed for new versions, at least >356  
                          readOutParameters = cms.PSet(debug = cms.untracked.bool(False),
                                                       rosParameters = cms.PSet(writeSC = cms.untracked.bool(True),
                                                                                readingDDU = cms.untracked.bool(True),
                                                                                performDataIntegrityMonitor = cms.untracked.bool(True),
                                                                                readDDUIDfromDDU = cms.untracked.bool(True),
                                                                                debug = cms.untracked.bool(False),
                                                                                localDAQ = cms.untracked.bool(False)
                                                                                ),
                                                       localDAQ = cms.untracked.bool(False),
                                                       performDataIntegrityMonitor = cms.untracked.bool(True)
                                                       )
                          )


###from Configuration.StandardSequences.Geometry_cff import *  ##  Deprecated in new versions > 53X
from Configuration.Geometry.GeometryIdeal_cff import *
from RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff import *
###dt1DRecHits.dtDigiLabel = 'dtunpacker'
dt1DRecHits.dtDigiLabel = 'muonDTDigis' #RECO



from RecoTracker.Configuration.RecoTracker_cff import *  ## Needed at least in  710pre8 to avoid an error in RecoMuonCosmics file (GroupedCkfTrajectoryBuilder)
from RecoMuon.Configuration.RecoMuonCosmics_cff import *
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi import *

cosmicMuonsBarrelOnly.TrajectoryBuilderParameters.EnableRPCMeasurement = False

from CondCore.DBCommon.CondDBSetup_cfi import *

###from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import * # DB v1
from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *  #DB v2, at least since GR_E_V42

####GlobalTag.globaltag = 'GR09_31X_V4P::All'  ## During CRAFT
####GlobalTag.globaltag = 'CRAFT09_R_V3::All'  ## For reprocessing with 327 (Sep09)
####GlobalTag.globaltag = 'GR09_P_V7::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR10_P_V2::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR_E_V23::All'  ## For Express 50X 2012 data
##GlobalTag.globaltag = 'GR_E_V25::All'  ## For Express 52X 2012 data
##GlobalTag.globaltag = 'GR_E_V26::All'  ## For Express 52X 2012 data
##GlobalTag.globaltag = 'GR_E_V31::All'  ## For Express 53X 2012 data -need extra conf in the python file for Elec not used for the moment 
##GlobalTag.globaltag = 'GR_E_V33A::All'  ## For Express 53X>=538HI 
##GlobalTag.globaltag = 'GR_E_V47::All' ## For Express >=741
GlobalTag.globaltag = 'GR_E_V47' ## With DB V2






unpackers  = cms.Sequence(dtunpacker + dttfunpacker)
reco       = cms.Sequence(dt1DRecHits * dt2DSegments * dt4DSegments)
###globalreco = cms.Sequence(CosmicMuonSeedBarrelOnly * offlineBeamSpot * cosmicMuonsBarrelOnly)
globalreco = cms.Sequence(CosmicMuonSeed*offlineBeamSpot*cosmicMuons)

#######################################################################################
# DT DPG DQM modules follow

###from UserCode.DTDPGAnalysis.DTOfflineAnalyzer_cfi import *
###from UserCode.DTDPGAnalysis.STAOfflineAnalyzer_cfi import *
###from UserCode.DTDPGAnalysis.DTEffOfflineAnalyzer_cfi import *

from UserCode.DTDPGAnalysis.DTOfflineAnalyzer_RECO_cfi import *
from UserCode.DTDPGAnalysis.STAOfflineAnalyzer_RECO_cfi import *
from UserCode.DTDPGAnalysis.DTEffOfflineAnalyzer_RECO_cfi import *


from DQMServices.Components.MEtoEDMConverter_cfi import *
from DQMServices.Core.DQM_cfg import *

from DQM.DTMonitorModule.dtDataIntegrityTask_cfi import *
DTDataIntegrityTask.dtDDULabel = 'dtunpacker'   ## Needed from, at least, 710  
DTDataIntegrityTask.dtROS25Label = 'dtunpacker' ## Needed from, at least, 710  

from DQM.DTMonitorModule.dtDigiTask_cfi import *
dtDigiMonitor.readDB = True
dtDigiMonitor.doNoiseOccupancies = True
dtDigiMonitor.doInTimeOccupancies = True
dtDigiMonitor.dtDigiLabel = 'muonDTDigis' #RECO

from DQM.DTMonitorModule.dtTriggerTask_cfi import *
dtTriggerMonitor.process_dcc = True
###dtTriggerMonitor.dcc_label   = 'dttfunpacker'
dtTriggerMonitor.dcc_label   = 'dttfDigis'  ## RECO

dtTriggerMonitor.process_seg = True

from DQM.DTMonitorModule.dtEfficiencyTask_cfi import *
from DQM.DTMonitorModule.dtChamberEfficiencyTask_cfi import *
from DQM.DTMonitorModule.dtResolutionTask_cfi import *

from DQM.DTMonitorModule.dtSegmentTask_cfi import *
dtSegmentAnalysisMonitor.detailedAnalysis = True

dummyProducer = cms.EDProducer("ThingWithMergeProducer")


from DQM.L1TMonitor.L1TGMT_cfi import *

##sources = cms.Sequence( dummyProducer + dtDigiMonitor + dtTriggerMonitor + dtEfficiencyMonitor + dtChamberEfficiencyMonitor + dtSegmentAnalysisMonitor + dtResolutionAnalysisMonitor + l1tGmt )
## From, at least, 710 DTDataIntegrityTask MUST be included also in the sources if not the folders 00-DataIntegrity and FEDIntegrity are missing
sources = cms.Sequence( dummyProducer + DTDataIntegrityTask + dtDigiMonitor + dtTriggerMonitor + dtEfficiencyMonitor + dtChamberEfficiencyMonitor + dtSegmentAnalysisMonitor + dtResolutionAnalysisMonitor + l1tGmt)

analysis = cms.Sequence(DTOfflineAnalyzer + DTEffOfflineAnalyzer + STAOfflineAnalyzer)

