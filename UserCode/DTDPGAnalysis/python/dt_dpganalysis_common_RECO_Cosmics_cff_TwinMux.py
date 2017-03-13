import FWCore.ParameterSet.Config as cms

###from EventFilter.DTTFRawToDigi.dttfunpacker_cfi import *
###dttfunpacker.DTTF_FED_Source =  'rawDataCollector'  ## MWGR Feb12 

#### TwinMuxUnpacker  ############################################################
from EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi import *

#### TwinMuxUnpacker  ############################################################
## Substituting DTTF by TwinMux but mantaining the dttfunpacker name to avoid changing 
## all the other modules
###dttfunpacker  = cms.EDProducer("L1TTwinMuxRawToDigi",
###twinMuxStage2Digis = cms.EDProducer("L1TTwinMuxRawToDigi",
###                               DTTM7_FED_Source = cms.InputTag("rawDataCollector"),
###                               feds     = cms.untracked.vint32( 1395,           1391,           1390,           1394,           1393           ),
###                               wheels   = cms.untracked.vint32( -2,             -1,             0,              +1,             +2             ),
###                               # Sector = '1' to '12' in HEX ('1' -> 'C'); 'F' if the AmcId is not associated to any Sector
###                               # AmcId                         (  123456789...)
###                               # Mapping                       (  FFFFF3FFF...)
###                               # Sector                        (  -----3---...)
###                               amcsecmap= cms.untracked.vint64( 0x123456789ABC, 0x123456789ABC, 0x123456789ABC, 0x123456789ABC, 0x123456789ABC ),
###                               debug    = cms.untracked.bool(False),
###                               passbc0  = cms.untracked.bool(False)
###)


########################################################################################
dtunpacker = cms.EDProducer("DTUnpackingModule",
                            dataType = cms.string('DDU'),
                            ##inputLabel = cms.InputTag('source'), ## needed for new versions, at least >356  
                            inputLabel = cms.InputTag('rawDataCollector'), ## MWGR Feb12
                            ###useStandardFEDid = cms.untracked.bool(True),
                            useStandardFEDid = cms.bool(True),
                            ###fedbyType = cms.untracked.bool(True),
                            ##fedbyType = cms.bool(False), ## tracked for new versions, at least >356, and also false ???? 
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



###############################################################################
###from Configuration.StandardSequences.Geometry_cff import *  ##  Deprecated in new versions > 53X
##from Configuration.Geometry.GeometryIdeal_cff import * ## In versions 7X problems with few STA events crashing the runing # but it works in the first tests of 750pre5 
##from Configuration.StandardSequences.GeometryRecoDB_cff import *  # solve STA problem, Doesn't work with 750pre5 (perhaps the globaltag?, still not defining any default when doing the test
from Configuration.StandardSequences.GeometryRecoDB_cff import *  # solve STA problem, works again in >805 
###from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *  # DB v1 
from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *  #DB v2, at least since GR_E_V42


from RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff import *
dt1DRecHits.dtDigiLabel = 'dtunpacker'

from RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi import *
CosmicMuonSeed.EnableCSCMeasurement=False

from RecoTracker.Configuration.RecoTracker_cff import *  ## Needed at least in  710pre8 to avoid an error in RecoMuonCosmics file (GroupedCkfTrajectoryBuilder)
from RecoMuon.Configuration.RecoMuonCosmics_cff import *
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagator_cfi import *

#cosmicMuonsBarrelOnly.TrajectoryBuilderParameters.EnableRPCMeasurement = False
## cosmicMuonsBarrelOnly doesn't exist on 44X, the only barrel used is 1leg
## see RecoLocalMuon.Configuration.RecoMuonCosmics_EventContent_cff and RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff
## We are not using our own reconstruction with this file, we comment it.
## to be fixed if we need to do the reconstruction for ourselves

ancientMuonSeed.EnableCSCMeasurement = False
###ancientMuonSeed.EnableRPCMeasurement = False

standAloneMuons.STATrajBuilderParameters.FilterParameters.EnableCSCMeasurement = False
standAloneMuons.STATrajBuilderParameters.FilterParameters.EnableRPCMeasurement = False
standAloneMuons.STATrajBuilderParameters.BWFilterParameters.EnableCSCMeasurement = False
standAloneMuons.STATrajBuilderParameters.BWFilterParameters.EnableRPCMeasurement = False

cosmicMuons.TrajectoryBuilderParameters.EnableRPCMeasurement = False
cosmicMuons.TrajectoryBuilderParameters.EnableCSCMeasurement = False

###from CondCore.DBCommon.CondDBSetup_cfi import *  ## Obsolete after 8X 
from CondCore.CondDB.CondDB_cfi import *                       # 


##from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *  #DB v2, at least since GR_E_V42
####GlobalTag.globaltag = 'GR09_31X_V4P::All'  ## During CRAFT
####GlobalTag.globaltag = 'CRAFT09_R_V3::All'  ## For reprocessing with 327 (Sep09)
####GlobalTag.globaltag = 'GR09_P_V7::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR10_P_V2::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR10_P_V5::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR10_P_V6::All'  ## For prompt reco on all PD's except Cosmics
##GlobalTag.globaltag = 'GR_H_V13::All'  ## For online (with 3110pre5 date 27/1/11 doesn't work)
##GlobalTag.globaltag = 'GR_R_311_V1::All'  ## For online
##GlobalTag.globaltag = 'GR_P_V20::All'  ## For prompt data 42X 
##GlobalTag.globaltag = 'GR_P_V26::All'  ## For prompt data 44X 
##GlobalTag.globaltag = 'GR_E_V23::All'  ## For Express 50X 2012 data
##GlobalTag.globaltag = 'GR_E_V25::All'  ## For Express 52X 2012 data
##GlobalTag.globaltag = 'GR_E_V26::All'  ## For Express 53X 2012 data
##GlobalTag.globaltag = 'GR_E_V31::All'  ## For Express 53X 2012 data -need extra conf in the python file for Elec not used for the moment 
##GlobalTag.globaltag = 'GR_E_V33A::All'  ## For Express 53X>=538HI 
##GlobalTag.globaltag = 'GR_E_V42::All'
##GlobalTag.globaltag = 'GR_E_V47::All' ## For Express >=741
##GlobalTag.globaltag = 'GR_E_V47' ## With DB V2
##GlobalTag.globaltag = 'GR_E_V49' ## With DB V2

##GlobalTag.globaltag = 'GR_H_V58' ## With DB V2 No tiene el tTrig de cosmicos y cuando lo pongo me da problemas
##GlobalTag.globaltag = '80X_dataRun2_Express_v0' ##for CMSSW_8_0_X X>0
##GlobalTag.globaltag = '80X_dataRun2_Express_v4' ##for CMSSW_8_0_X X>0
##GlobalTag.globaltag = '80X_dataRun2_Express_v7' ##for CMSSW_8_0_X X>0
##GlobalTag.globaltag = '80X_dataRun2_Express_v10' ##for CMSSW_8_0_10 After insertion of TOTEM
GlobalTag.globaltag = '90X_dataRun2_Express_v1' ## Since CMSSW_8_3_0 (commisioning 2017)

##unpackers  = cms.Sequence(dtunpacker + dttfunpacker)
unpackers  = cms.Sequence(dtunpacker + twinMuxStage2Digis)
reco       = cms.Sequence(dt1DRecHits * dt2DSegments * dt4DSegments)
#globalreco = cms.Sequence(CosmicMuonSeedBarrelOnly * offlineBeamSpot * cosmicMuonsBarrelOnly)
## cosmicMuonsBarrelOnly doesn't exist on 44X, the only barrel used is 1leg
## see RecoLocalMuon.Configuration.RecoMuonCosmics_EventContent_cff and RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff
## We are not using our own reconstruction with this file, we comment it.
## to be fixed if we need to do the reconstruction for ourselves
##globalreco = cms.Sequence(standAloneMuonSeeds * offlineBeamSpot * standAloneMuons)  ## old
globalreco = cms.Sequence(CosmicMuonSeed*offlineBeamSpot*cosmicMuons)

#######################################################################################
# DT DPG DQM modules follow

from UserCode.DTDPGAnalysis.DTOfflineAnalyzer_Cosmics_cfi import *
##DTOfflineAnalyzer.SALabel = 'standAloneMuons'
DTOfflineAnalyzer.SALabel = 'cosmicMuons'
from UserCode.DTDPGAnalysis.STAOfflineAnalyzer_Cosmics_cfi import *
##STAOfflineAnalyzer.SALabel = 'standAloneMuons'
STAOfflineAnalyzer.SALabel = 'cosmicMuons'

from UserCode.DTDPGAnalysis.DTEffOfflineAnalyzer_cfi import *


from DQMServices.Components.MEtoEDMConverter_cfi import *
from DQMServices.Core.DQM_cfg import *

from DQM.DTMonitorModule.dtDataIntegrityTask_cfi import *
DTDataIntegrityTask.dtDDULabel = 'dtunpacker'   ## Needed from, at least, 710  
DTDataIntegrityTask.dtROS25Label = 'dtunpacker' ## Needed from, at least, 710  

from DQM.DTMonitorModule.dtDigiTask_cfi import *
dtDigiMonitor.readDB = True
dtDigiMonitor.doNoiseOccupancies = True
dtDigiMonitor.doInTimeOccupancies = True

from DQM.DTMonitorModule.dtTriggerTask_cfi import *
##dtTriggerMonitor.process_dcc = True
##dtTriggerMonitor.dcc_label   = 'dttfunpacker'
dtTriggerMonitor.process_ros = False ## New data has not DDU information
###process.dtTriggerMonitor.process_tm = True (already set in the "official" configuration, and called process_dcc in versions 80X,i tshould be fixed in 81X)
dtTriggerMonitor.process_seg = True

from DQM.DTMonitorModule.dtEfficiencyTask_cfi import *
from DQM.DTMonitorModule.dtChamberEfficiencyTask_cfi import *
from DQM.DTMonitorModule.dtResolutionTask_cfi import *

from DQM.DTMonitorModule.dtSegmentTask_cfi import *
dtSegmentAnalysisMonitor.detailedAnalysis = True

dummyProducer = cms.EDProducer("ThingWithMergeProducer")

from DQM.L1TMonitor.L1TGMT_cfi import *

##sources = cms.Sequence( dummyProducer + dtDigiMonitor + dtTriggerMonitor + dtEfficiencyMonitor + dtChamberEfficiencyMonitor + dtSegmentAnalysisMonitor + dtResolutionAnalysisMonitor + l1tGmt)
## From, at least, 710 DTDataIntegrityTask MUST be included also in the sources if not the folders 00-DataIntegrity and FEDIntegrity are missing
sources = cms.Sequence( dummyProducer + DTDataIntegrityTask + dtDigiMonitor + dtTriggerMonitor + dtEfficiencyMonitor + dtChamberEfficiencyMonitor + dtSegmentAnalysisMonitor + dtResolutionAnalysisMonitor + l1tGmt)

analysis = cms.Sequence(DTOfflineAnalyzer + DTEffOfflineAnalyzer + STAOfflineAnalyzer)

