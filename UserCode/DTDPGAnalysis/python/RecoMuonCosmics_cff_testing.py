import FWCore.ParameterSet.Config as cms

######################## Cosmic Reco #############################

## Full detector ##

# Seed generator
from RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi import *

# Stand alone muon track producer
from RecoMuon.CosmicMuonProducer.cosmicMuons_cff import *

# Global muon track producer
from RecoMuon.CosmicMuonProducer.globalCosmicMuons_cff import *
globalCosmicMuons.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'ctfWithMaterialTracksP5LHCNavigation'


### CAMBIO              
##from RecoLuminosity.LumiProducer.lumiProducer_cff import *
##from RecoLocalCalo.Configuration.RecoLocalCalo_cff import *
from RecoTracker.Configuration.RecoTracker_cff import *
##from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
##from TrackingTools.Configuration.TrackingTools_cff import *
##from RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi import *
##from RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi import *
### CAMBIO              

# Muon Id producer
from RecoMuon.MuonIdentification.muonIdProducerSequence_cff import *

### CAMBIO              muons = muons1stStep.clone()

### CAMBIO              muons.inputCollectionLabels = ['ctfWithMaterialTracksP5LHCNavigation', 'globalCosmicMuons', 'cosmicMuons', "tevMuons:firstHit","tevMuons:picky","tevMuons:dyt"]
### CAMBIO              muons.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks', 'tev firstHit', 'tev picky', 'tev dyt']
### CAMBIO              muons.fillIsolation = True
### CAMBIO              muons.fillGlobalTrackQuality = True
# need to modify track selection as well (not clear to what)
### CAMBIO              muons.TrackExtractorPSet.inputTrackCollection = 'ctfWithMaterialTracksP5LHCNavigation'
### CAMBIO              muons.CaloExtractorPSet.CenterConeOnCalIntersection = True

### CAMBIO              from RecoMuon.MuonIdentification.calomuons_cfi import *
### CAMBIO              calomuons.inputTracks = 'ctfWithMaterialTracksP5LHCNavigation'
### CAMBIO              calomuons.inputCollection = 'muons'
### CAMBIO              calomuons.inputMuons = 'muons'

## Sequences

# Stand Alone Tracking
### CAMBIO              STAmuontrackingforcosmics = cms.Sequence(CosmicMuonSeed*cosmicMuons)
# Stand Alone Tracking plus muon ID
### CAMBIO              STAmuonrecoforcosmics = cms.Sequence(STAmuontrackingforcosmics)

# Stand Alone Tracking plus global tracking
### CAMBIO              muontrackingforcosmics = cms.Sequence(STAmuontrackingforcosmics*globalCosmicMuons)

# Muon Isolation sequence
### CAMBIO              from RecoMuon.MuonIsolationProducers.muIsolation_cff import *
# muisodeposits based on "muons"
# we are using copy extractors now
### CAMBIO              muIsoDepositTk.inputTags = cms.VInputTag(cms.InputTag("muons:tracker"))
### CAMBIO              muIsoDepositJets. inputTags = cms.VInputTag(cms.InputTag("muons:jets"))
### CAMBIO              muIsoDepositCalByAssociatorTowers.inputTags = cms.VInputTag(cms.InputTag("muons:ecal"), cms.InputTag("muons:hcal"), cms.InputTag("muons:ho"))



# TeV refinement
### CAMBIO              from RecoMuon.GlobalMuonProducer.tevMuons_cfi import *
### CAMBIO              tevMuons.MuonCollectionLabel = "globalCosmicMuons"
### CAMBIO              tevMuons.RefitterParameters.PropDirForCosmics = cms.bool(True)

# Glb Track Quality
### CAMBIO              from RecoMuon.GlobalTrackingTools.GlobalTrackQuality_cfi import *
### CAMBIO              glbTrackQual.InputCollection = "globalCosmicMuons"

# all muons id
### CAMBIO              allmuons = cms.Sequence(glbTrackQual*tevMuons*muons*muIsolation*calomuons)

# Final sequence
### CAMBIO              muonrecoforcosmics = cms.Sequence(muontrackingforcosmics*allmuons)
### CAMBIO              muonRecoAllGR = cms.Sequence(muonrecoforcosmics)

# 1 leg mode

# Stand alone muon track producer
### CAMBIO              cosmicMuons1Leg = cosmicMuons.clone()
### CAMBIO              cosmicMuons1Leg.TrajectoryBuilderParameters.BuildTraversingMuon = True
### CAMBIO              cosmicMuons1Leg.MuonSeedCollectionLabel = 'CosmicMuonSeed'

# Global muon track producer
### CAMBIO              globalCosmicMuons1Leg = globalCosmicMuons.clone()
### CAMBIO              globalCosmicMuons1Leg.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'ctfWithMaterialTracksP5'
### CAMBIO              globalCosmicMuons1Leg.MuonCollectionLabel = 'cosmicMuons1Leg'

# Muon Id producer
### CAMBIO              muons1Leg = muons1stStep.clone()
### CAMBIO              muons1Leg.inputCollectionLabels = ['ctfWithMaterialTracksP5', 'globalCosmicMuons1Leg', 'cosmicMuons1Leg']
### CAMBIO              muons1Leg.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
### CAMBIO              muons1Leg.fillIsolation = False
### CAMBIO              muons1Leg.fillGlobalTrackQuality = False
### CAMBIO              muons1Leg.fillGlobalTrackRefits = False
# Sequences

# Stand Alone Tracking
### CAMBIO              STAmuontrackingforcosmics1Leg = cms.Sequence(CosmicMuonSeed*cosmicMuons1Leg)

# Stand Alone Tracking plus global tracking
### CAMBIO              muontrackingforcosmics1Leg = cms.Sequence(STAmuontrackingforcosmics1Leg*globalCosmicMuons1Leg)

# all muons id
### CAMBIO              allmuons1Leg = cms.Sequence(muons1Leg)

# Stand Alone Tracking plus muon ID
### CAMBIO              STAmuonrecoforcosmics1Leg = cms.Sequence(STAmuontrackingforcosmics1Leg)

# Final sequence
### CAMBIO              muonrecoforcosmics1Leg = cms.Sequence(muontrackingforcosmics1Leg*allmuons1Leg)

#####################################################

# t0 Corrections

# Seed generator
### CAMBIO              CosmicMuonSeedWitht0Correction = CosmicMuonSeed.clone()
### CAMBIO              CosmicMuonSeedWitht0Correction.DTRecSegmentLabel = 'dt4DSegmentsT0Seg'

# Stand alone muon track producer
### CAMBIO              cosmicMuonsWitht0Correction = cosmicMuons.clone()
### CAMBIO              cosmicMuonsWitht0Correction.TrajectoryBuilderParameters.BuildTraversingMuon = False
### CAMBIO              cosmicMuonsWitht0Correction.MuonSeedCollectionLabel = 'CosmicMuonSeedWitht0Correction'
### CAMBIO              cosmicMuonsWitht0Correction.TrajectoryBuilderParameters.DTRecSegmentLabel = 'dt4DSegmentsT0Seg'

# Global muon track producer
### CAMBIO              globalCosmicMuonsWitht0Correction = globalCosmicMuons.clone()
### CAMBIO              globalCosmicMuonsWitht0Correction.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'ctfWithMaterialTracksP5'
### CAMBIO              globalCosmicMuonsWitht0Correction.MuonCollectionLabel = 'cosmicMuonsWitht0Correction'

# Muon Id producer
### CAMBIO              muonsWitht0Correction = muons1stStep.clone()
### CAMBIO              muonsWitht0Correction.inputCollectionLabels = ['ctfWithMaterialTracksP5', 'globalCosmicMuonsWitht0Correction', 'cosmicMuonsWitht0Correction']
### CAMBIO              muonsWitht0Correction.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
### CAMBIO              muonsWitht0Correction.fillIsolation = True
### CAMBIO              muonsWitht0Correction.fillGlobalTrackQuality = False
### CAMBIO              muonsWitht0Correction.TimingFillerParameters.DTTimingParameters.UseSegmentT0 = True
### CAMBIO              muonsWitht0Correction.TimingFillerParameters.DTTimingParameters.DTsegments = 'dt4DSegmentsT0Seg'
### CAMBIO              muonsWitht0Correction.TimingFillerParameters.DTTimingParameters.MatchParameters.DTsegments = 'dt4DSegmentsT0Seg'
### CAMBIO              muonsWitht0Correction.TrackExtractorPSet.inputTrackCollection = 'ctfWithMaterialTracksP5'
### CAMBIO              muonsWitht0Correction.CaloExtractorPSet.CenterConeOnCalIntersection = True
### CAMBIO              muonsWitht0Correction.fillGlobalTrackRefits = False
#Sequences


# Stand Alone Tracking
### CAMBIO              STAmuontrackingforcosmicsWitht0Correction = cms.Sequence(CosmicMuonSeedWitht0Correction*cosmicMuonsWitht0Correction)

# Stand Alone Tracking plus global tracking
### CAMBIO              muontrackingforcosmicsWitht0Correction = cms.Sequence(STAmuontrackingforcosmicsWitht0Correction*globalCosmicMuonsWitht0Correction)

# Stand Alone Tracking plus muon ID
### CAMBIO              STAmuonrecoforcosmicsWitht0Correction = cms.Sequence(STAmuontrackingforcosmicsWitht0Correction)

# all muons id
### CAMBIO              allmuonsWitht0Correction = cms.Sequence(muonsWitht0Correction)

# Final sequence
### CAMBIO              muonrecoforcosmicsWitht0Correction = cms.Sequence(muontrackingforcosmicsWitht0Correction*allmuonsWitht0Correction)

### Final sequence ###
### CAMBIO              muonRecoGR = cms.Sequence(muonrecoforcosmics1Leg+muonrecoforcosmicsWitht0Correction)

#####################################################

# Beam halo in Encaps only. GLB reco only is needed

# Seed generator 
### CAMBIO              CosmicMuonSeedEndCapsOnly = CosmicMuonSeed.clone()
### CAMBIO              CosmicMuonSeedEndCapsOnly.EnableDTMeasurement = False

# Stand alone muon track producer
### CAMBIO              cosmicMuonsEndCapsOnly = cosmicMuons.clone()
### CAMBIO              cosmicMuonsEndCapsOnly.TrajectoryBuilderParameters.EnableDTMeasurement = False
### CAMBIO              cosmicMuonsEndCapsOnly.TrajectoryBuilderParameters.MuonNavigationParameters.Barrel = False
### CAMBIO              cosmicMuonsEndCapsOnly.MuonSeedCollectionLabel = 'CosmicMuonSeedEndCapsOnly'

# Global muon track producer
### CAMBIO              globalBeamHaloMuonEndCapslOnly = globalCosmicMuons.clone()
### CAMBIO              globalBeamHaloMuonEndCapslOnly.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'beamhaloTracks'
### CAMBIO              globalBeamHaloMuonEndCapslOnly.MuonCollectionLabel = 'cosmicMuonsEndCapsOnly'


# Muon Id producer
### CAMBIO              muonsBeamHaloEndCapsOnly = muons1stStep.clone()           
### CAMBIO              muonsBeamHaloEndCapsOnly.inputCollectionLabels = ['beamhaloTracks', 'globalBeamHaloMuonEndCapslOnly', 'cosmicMuonsEndCapsOnly']
### CAMBIO              muonsBeamHaloEndCapsOnly.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
### CAMBIO              muonsBeamHaloEndCapsOnly.fillIsolation = True
### CAMBIO              muonsBeamHaloEndCapsOnly.fillGlobalTrackQuality = False
### CAMBIO              muonsBeamHaloEndCapsOnly.TrackExtractorPSet.inputTrackCollection = 'ctfWithMaterialTracksP5'
### CAMBIO              muonsBeamHaloEndCapsOnly.CaloExtractorPSet.CenterConeOnCalIntersection = True
### CAMBIO              muonsBeamHaloEndCapsOnly.fillGlobalTrackRefits = False

# Sequences
### CAMBIO              muonrecoBeamHaloEndCapsOnly = cms.Sequence(CosmicMuonSeedEndCapsOnly*cosmicMuonsEndCapsOnly*globalBeamHaloMuonEndCapslOnly*muonsBeamHaloEndCapsOnly)

########

## Full detector but NO RPC ##

# Stand alone muon track producer
### CAMBIO              cosmicMuonsNoRPC = cosmicMuons.clone()
### CAMBIO              cosmicMuonsNoRPC.TrajectoryBuilderParameters.EnableRPCMeasurement = False

# Global muon track producer
### CAMBIO              globalCosmicMuonsNoRPC = globalCosmicMuons.clone()
### CAMBIO              globalCosmicMuonsNoRPC.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'ctfWithMaterialTracksP5'
### CAMBIO              globalCosmicMuonsNoRPC.MuonCollectionLabel = 'cosmicMuonsNoRPC'

# Muon Id producer
### CAMBIO              muonsNoRPC = muons1stStep.clone()
### CAMBIO              muonsNoRPC.inputCollectionLabels = ['ctfWithMaterialTracksP5', 'globalCosmicMuonsNoRPC', 'cosmicMuonsNoRPC']
### CAMBIO              muonsNoRPC.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
### CAMBIO              muonsNoRPC.fillIsolation = True
### CAMBIO              muonsNoRPC.fillGlobalTrackQuality = False
### CAMBIO              muonsNoRPC.TrackExtractorPSet.inputTrackCollection = 'ctfWithMaterialTracksP5'
### CAMBIO              muonsNoRPC.CaloExtractorPSet.CenterConeOnCalIntersection = True
### CAMBIO              muonsNoRPC.fillGlobalTrackRefits = False

#Sequences

# Stand Alone Tracking
### CAMBIO              STAmuontrackingforcosmicsNoRPC = cms.Sequence(cosmicMuonsNoRPC)

# Stand Alone Tracking plus global tracking
### CAMBIO              muontrackingforcosmicsNoRPC = cms.Sequence(STAmuontrackingforcosmicsNoRPC*globalCosmicMuonsNoRPC)

# all muons id
### CAMBIO              allmuonsNoRPC = cms.Sequence(muonsNoRPC)

# Final sequence
### CAMBIO              muonrecoforcosmicsNoRPC = cms.Sequence(muontrackingforcosmicsNoRPC*allmuonsNoRPC)

##############################################

## Split Tracks  ##

# Global muon track producer
### CAMBIO              globalCosmicSplitMuons = globalCosmicMuons.clone()
### CAMBIO              globalCosmicSplitMuons.TrajectoryBuilderParameters.TkTrackCollectionLabel = 'splittedTracksP5'
### CAMBIO              globalCosmicSplitMuons.MuonCollectionLabel = 'cosmicMuons'

# Muon Id producer

### CAMBIO              splitMuons = muons1stStep.clone()
### CAMBIO              splitMuons.inputCollectionLabels = ['splittedTracksP5', 'globalCosmicSplitMuons', 'cosmicMuons']
### CAMBIO              splitMuons.inputCollectionTypes = ['inner tracks', 'links', 'outer tracks']
### CAMBIO              splitMuons.fillIsolation = True
### CAMBIO              splitMuons.fillGlobalTrackQuality = False
### CAMBIO              splitMuons.TrackExtractorPSet.inputTrackCollection = 'splittedTracksP5'
### CAMBIO              splitMuons.CaloExtractorPSet.CenterConeOnCalIntersection = True
### CAMBIO              splitMuons.fillGlobalTrackRefits = False

#Sequences

# Final sequence
### CAMBIO              muonrecoforsplitcosmics = cms.Sequence(globalCosmicSplitMuons*splitMuons)

##############################################

######################## LHC like Reco #############################

# Standard reco
### CAMBIO              from RecoMuon.Configuration.RecoMuonPPonly_cff import *

# Muon Id producer
### CAMBIO              lhcSTAMuons = muons1stStep.clone()
### CAMBIO              lhcSTAMuons.inputCollectionLabels = ['standAloneMuons']
### CAMBIO              lhcSTAMuons.inputCollectionTypes = ['outer tracks']
### CAMBIO              lhcSTAMuons.fillIsolation = True
### CAMBIO              lhcSTAMuons.fillGlobalTrackQuality = False
### CAMBIO              lhcSTAMuons.TrackExtractorPSet.inputTrackCollection = 'ctfWithMaterialTracksP5LHCNavigation'
### CAMBIO              lhcSTAMuons.CaloExtractorPSet.CenterConeOnCalIntersection = True
### CAMBIO              lhcSTAMuons.fillGlobalTrackRefits = False

# Final sequence
### CAMBIO              muonRecoLHC = cms.Sequence(ancientMuonSeed*standAloneMuons*lhcSTAMuons)



########################### SEQUENCE TO BE ADDED in ReconstructionGR_cff ##############################################

### CAMBIO              muonRecoGR = cms.Sequence(muonRecoAllGR*muonRecoGR*muonrecoBeamHaloEndCapsOnly*muonrecoforcosmicsNoRPC*muonrecoforsplitcosmics*muonRecoLHC)

#######################################################################################################################





