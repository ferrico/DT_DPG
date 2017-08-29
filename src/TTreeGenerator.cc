///
// Package:    TTreeGenerator
// Class:      TTreeGenerator
// 
/**\class TTreeGenerator TTreeGenerator.cc MyTools/TTreeGenerator/src/TTreeGenerator.cc
 Description: <one line class summary>
 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Carlo BATTILANA, Mario PELLICCIONI
//         Created:  Mon Jan 11 14:59:51 CET 2010
// $Id: TTreeGenerator.cc,v 1.33 2012/07/02 16:43:36 guiducci Exp $
//
// Modificated M.C Fouz March/2016 for TwinMux and to include tracks extrapolation and times variable
// Modifications L. Guiducci July 2016 to include TwinMux output data and clean up legacy trigger information

// user include files
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "DataFormats/DTDigi/interface/DTDigi.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/DTDigi/interface/DTLocalTriggerCollection.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"  // New trying to avoid crashes in the topology functions
#include "Geometry/DTGeometry/interface/DTLayer.h" // New trying to avoid crashes in the topology functions
#include "Geometry/DTGeometry/interface/DTTopology.h" // New trying to avoid crashes in the topology functions
#include <DataFormats/MuonDetId/interface/DTLayerId.h> // New trying to avoid crashes in the topology functions

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include <DataFormats/RPCDigi/interface/RPCDigiCollection.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHitCollection.h>
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCRoll.h" 
#include "UserCode/DTDPGAnalysis/interface/DTSegmentToRPC.h"

#include "CalibMuon/DTDigiSync/interface/DTTTrigBaseSync.h"
#include "CalibMuon/DTDigiSync/interface/DTTTrigSyncFactory.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
#include "UserCode/DTDPGAnalysis/interface/TTreeGenerator.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include <DataFormats/MuonData/interface/MuonDigiCollection.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>

#include "RecoLocalMuon/RPCRecHit/interface/RPCPointProducer.h"

#include <iostream>
#include <vector>
#include "TMath.h"


using namespace std;
using namespace reco;

ObjectMap* ObjectMap::mapInstance = NULL;
 
 ObjectMap* ObjectMap::GetInstance(const edm::EventSetup& iSetup){
   if (mapInstance == NULL){
     mapInstance = new ObjectMap(iSetup);
   }
   return mapInstance;
 }
 
 ObjectMap::ObjectMap(const edm::EventSetup& iSetup)
 {
   edm::ESHandle<RPCGeometry> rpcGeo;
   edm::ESHandle<DTGeometry> dtGeo;

  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
   iSetup.get<MuonGeometryRecord>().get(dtGeo);
   
   for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
     if(dynamic_cast< const RPCChamber* >( *it ) != 0 ){
       auto ch = dynamic_cast< const RPCChamber* >( *it ); 
       std::vector< const RPCRoll*> roles = (ch->rolls());
       for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
       		RPCDetId rpcId = (*r)->id();
       		int region=rpcId.region();
       		if(region==0){
       				int wheel=rpcId.ring();
       				int sector=rpcId.sector();
       				int station=rpcId.station();
       				DTStationIndex ind(region,wheel,sector,station);
       				std::set<RPCDetId> myrolls;
       				if (rollstoreDT.find(ind)!=rollstoreDT.end()) myrolls=rollstoreDT[ind];
       				myrolls.insert(rpcId);
       				rollstoreDT[ind]=myrolls;
     		}
       }
 	}
   }
}
 

TTreeGenerator::TTreeGenerator(const edm::ParameterSet& pset):
      rpcToken_(consumes<MuonDigiCollection<RPCDetId,RPCDigi> > (pset.getParameter<edm::InputTag>("rpcLabel"))),
      UnpackingRpcRecHitToken_(consumes<RPCRecHitCollection> (pset.getParameter<edm::InputTag>("UnpackingRpcRecHitLabel")))
{
  // get the tTrigDBInfo
  theSync =
      DTTTrigSyncFactory::get()->create(pset.getUntrackedParameter<std::string>("tTrigMode"),
                                pset.getUntrackedParameter<edm::ParameterSet>("tTrigModeConfig"));


  // get run configuration options
  runOnRaw_        = pset.getParameter<bool>("runOnRaw");
  runOnSimulation_ = pset.getParameter<bool>("runOnSimulation");
  runOnSimulationWithDigis_ = pset.getParameter<bool>("runOnSimulationWithDigis");

  //get parameters from the configuration file
  //names of the different event collections
  dtDigiLabel_     = pset.getParameter<edm::InputTag>("dtDigiLabel");

  dtDigiToken_     = consumes<DTDigiCollection>(edm::InputTag(dtDigiLabel_));
  dtDigiTokenSim_  = consumes<DTDigiSimLinkCollection>(edm::InputTag(dtDigiLabel_));

  dtSegmentLabel_  = pset.getParameter<edm::InputTag>("dtSegmentLabel");
  dtSegmentToken_  = consumes<DTRecSegment4DCollection>(edm::InputTag(dtSegmentLabel_));

  cscSegmentLabel_ = pset.getParameter<edm::InputTag>("cscSegmentLabel");
  cscSegmentToken_ = consumes<CSCSegmentCollection>(edm::InputTag(cscSegmentLabel_));

  dtTrigTwinMuxInLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxInLabel");
  dtTrigTwinMuxThLabel_    = pset.getParameter<edm::InputTag>("dtTrigTwinMuxThLabel");
  dtTrigTwinMuxInToken_    = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxInLabel_)); 
  dtTrigTwinMux_ThToken_ = consumes<L1MuDTChambThContainer>(edm::InputTag(dtTrigTwinMuxThLabel_));

  dtTrigTwinMuxOutLabel_  = pset.getParameter<edm::InputTag>("dtTrigTwinMuxOutLabel");
  dtTrigTwinMuxOutToken_  = consumes<L1MuDTChambPhContainer>(edm::InputTag(dtTrigTwinMuxOutLabel_));

  staMuLabel_      = pset.getParameter<edm::InputTag>("staMuLabel");
  staMuToken_      = consumes<reco::MuonCollection>(edm::InputTag(staMuLabel_));

  gmtLabel_        = pset.getParameter<edm::InputTag>("gmtLabel");
  gmtToken_        = consumes<l1t::MuonBxCollection>(edm::InputTag(gmtLabel_));

  triggerResultTag_   = pset.getParameter<edm::InputTag>("TriggerResultsTag");
  triggerResultToken_ = consumes<edm::TriggerResults>(edm::InputTag(triggerResultTag_));

  triggerEventTag_   = pset.getParameter<edm::InputTag>("TriggerEventTag");
  triggerEventToken_ = consumes<trigger::TriggerEvent>(edm::InputTag(triggerEventTag_));

  gtLabel_         = pset.getParameter<edm::InputTag>("gtLabel"); // legacy
  gtToken_         = consumes<L1GlobalTriggerReadoutRecord>(edm::InputTag(gtLabel_)); //legacy

  rpcRecHitLabel_  = pset.getParameter<edm::InputTag>("rpcRecHitLabel");
  rpcRecHitToken_  = consumes<RPCRecHitCollection>(edm::InputTag(rpcRecHitLabel_));

  //max size of the different saved objects (per event)
  digisSize_       = pset.getParameter<int>("dtDigiSize");
  dtsegmentsSize_  = pset.getParameter<int>("dtSegmentSize");
  cscsegmentsSize_ = pset.getParameter<int>("cscSegmentSize");
  dtltTwinMuxOutSize_     = pset.getParameter<int>("dtTrigTwinMuxOutSize");
  dtltTwinMuxInSize_ = pset.getParameter<int>("dtTrigTwinMuxInSize");
  dtltTwinMuxThSize_ = pset.getParameter<int>("dtTrigTwinMuxThSize");
  gmtSize_         = pset.getParameter<int>("gmtSize");
  hltSize_         = pset.getParameter<int>("hltSize");
  STAMuSize_       = pset.getParameter<int>("STAMuSize");
  rpcRecHitSize_   = pset.getParameter<int>("rpcRecHitSize"); 

  PrimaryVertexTag_ = pset.getParameter<edm::InputTag>("PrimaryVertexTag");
  PrimaryVertexToken_ =consumes<reco::VertexCollection>(edm::InputTag(PrimaryVertexTag_));

  beamSpotTag_      = pset.getParameter<edm::InputTag>("beamSpotTag");
  beamSpotToken_    = consumes<reco::BeamSpot>(edm::InputTag(beamSpotTag_));

  scalersSource_    = pset.getParameter<edm::InputTag>("scalersResults");
  scalersSourceToken_ = consumes<LumiScalersCollection>(edm::InputTag(scalersSource_));
       
  lumiInputTag_      = pset.getParameter<edm::InputTag>("lumiInputTag");
  lumiProducerToken_ = consumes<LumiDetails, edm::InLumi>(lumiInputTag_);

  bmtfPhInputTag_ = consumes<L1MuDTChambPhContainer>(pset.getParameter<edm::InputTag>("bmtfInputPhDigis"));
  bmtfThInputTag_ = consumes<L1MuDTChambThContainer>(pset.getParameter<edm::InputTag>("bmtfInputThDigis"));
  bmtfOutputTag_ = consumes<l1t::RegionalMuonCandBxCollection>(pset.getParameter<edm::InputTag>("bmtfOutputDigis"));
  
  OnlyBarrel_ = pset.getParameter<bool>("OnlyBarrel");
  dtExtrapolation_ = pset.getParameter<bool>("dtExtrapolation");

  localDTmuons_    = pset.getUntrackedParameter<bool>("localDTmuons","False");

  AnaTrackGlobalMu_= pset.getUntrackedParameter<bool>("AnaTrackGlobalMu","True");  // set to False when problems with tracks of global muons

  outFile_         = pset.getParameter<std::string>("outputFile");

  initialize_Tree_variables();

  //counters
  idigis       = 0;
  idtsegments  = 0;
  icscsegments = 0;
  idtltTwinMuxOut     = 0;
  idtltTwinMuxIn     = 0;
  idtltTwinMux_th  = 0;
  imuons       = 0;
  igmt         = 0;
  igtalgo      = 0; // legacy
  igttt        = 0; // legacy
  ihltFilters  = 0;
  ihltPaths    = 0;
}

void TTreeGenerator::beginLuminosityBlock(edm::LuminosityBlock const& lumiBlock,
                                          edm::EventSetup const& context)
{
   return;
}



void TTreeGenerator::analyze(const edm::Event& event, const edm::EventSetup& context)
{

   theSync->setES(context);



   edm::ESHandle<DTGeometry> dtGeomH;
   context.get<MuonGeometryRecord>().get(dtGeomH);
   const DTGeometry* dtGeom_ = dtGeomH.product();
   
   edm::ESHandle<RPCGeometry> rpcGeomH;  
   context.get<MuonGeometryRecord>().get(rpcGeomH);
   const RPCGeometry* rpcGeom_ = rpcGeomH.product();


  //retrieve the beamspot info
  if(!localDTmuons_)  
  {
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    event.getByToken(beamSpotToken_ ,recoBeamSpotHandle);
    beamspot = *recoBeamSpotHandle; 
    
    //retrieve the luminosity
    if(!runOnSimulation_) // crashes with simulation trying to get the lumiperblock   
                          // The vector<LumiScalers> but not the <LumiScalersCollection> 
    {                      
       edm::Handle<LumiScalersCollection> lumiScalers;
       event.getByToken(scalersSourceToken_, lumiScalers);
       LumiScalersCollection::const_iterator lumiIt = lumiScalers->begin();
       lumiperblock = lumiIt->instantLumi();
    }
  }

  //retrieve the collections you are interested on in the event
  edm::Handle<DTDigiCollection> dtdigis;
  edm::Handle<MuonDigiCollection<DTLayerId, DTDigiSimLink>> dtdigisSim;
  if(runOnRaw_ && !runOnSimulation_) event.getByToken(dtDigiToken_, dtdigis);
  else   // Added to cope with simulation including Digis
    if(runOnSimulation_ && runOnSimulationWithDigis_) event.getByToken(dtDigiTokenSim_, dtdigisSim);


  edm::Handle<DTRecSegment4DCollection> dtsegments4D;
  event.getByToken(dtSegmentToken_, dtsegments4D);

  context.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  edm::Handle<reco::VertexCollection> privtxs;
  if(!localDTmuons_) event.getByToken(PrimaryVertexToken_, privtxs);
  
  edm::Handle<CSCSegmentCollection> cscsegments;
  event.getByToken(cscSegmentToken_, cscsegments);

  edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxOut;
  bool hasPhiTwinMuxOut=false;
  if(runOnRaw_) hasPhiTwinMuxOut=event.getByToken(dtTrigTwinMuxOutToken_,localTriggerTwinMuxOut);

  edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxIn;
  bool hasPhiTwinMuxIn=false;
  if(runOnRaw_) hasPhiTwinMuxIn=event.getByToken(dtTrigTwinMuxInToken_,localTriggerTwinMuxIn);

  edm::Handle<L1MuDTChambThContainer> localTriggerTwinMux_Th;

  edm::Handle<reco::MuonCollection> MuList;
  if(!localDTmuons_) event.getByToken(staMuToken_,MuList);

  edm::Handle<l1t::MuonBxCollection> gmt;   // legacy
  if(!localDTmuons_) event.getByToken(gmtToken_,gmt); // legacy

  edm::Handle< L1GlobalTriggerReadoutRecord > gtrc; // legacy
  if(runOnRaw_ && !localDTmuons_) event.getByToken(gtToken_, gtrc); // legacy

  edm::Handle<edm::TriggerResults>  hltresults;
  if(!localDTmuons_) event.getByToken(triggerResultToken_, hltresults); 

  edm::Handle<trigger::TriggerEvent>  hltevent;
  if(!localDTmuons_) event.getByToken(triggerEventToken_, hltevent); 

  edm::Handle<RPCRecHitCollection> rpcHits;
  if(!localDTmuons_) event.getByToken(rpcRecHitToken_,rpcHits);

  //clear the containers
  clear_Arrays();

  // Get the propagators
  context.get<TrackingComponentsRecord>().get("SmartPropagatorAny",propagatorAlong);
  context.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", propagatorOpposite);

  //get the magnetic field
  context.get<IdealMagneticFieldRecord>().get(theBField);

  //Fill the event info block
  runnumber = event.run();
  lumiblock = event.getLuminosityBlock().luminosityBlock();
  eventNumber = event.eventAuxiliary().event();
  timestamp = event.eventAuxiliary().time().value();
  bunchXing = event.eventAuxiliary().bunchCrossing();
  orbitNum = event.eventAuxiliary().orbitNumber();

  // it's filling nothing at the moment lumiDetails->isValid() return false (30/04/2016  M.C.F) 
  if(!localDTmuons_ && !runOnSimulation_)  // Crashes when run in simulation no <LumiDetails> available         
  {
     edm::Handle<LumiDetails> lumiDetails;
     event.getLuminosityBlock().getByToken(lumiProducerToken_, lumiDetails); 
     if(lumiDetails->isValid()){
       beam1Intensity = lumiDetails->lumiBeam1Intensity(bunchXing);
       beam2Intensity = lumiDetails->lumiBeam2Intensity(bunchXing);
     }
  }

  //Primary vertex
  if(!localDTmuons_)     
  {
    if((*privtxs).size() != 0){
      PV_x = (*privtxs)[0].position().x();
      PV_y = (*privtxs)[0].position().y();
      PV_z = (*privtxs)[0].position().z();
    
      PV_xxE = (*privtxs)[0].covariance(0,0);
      PV_yyE = (*privtxs)[0].covariance(1,1);
      PV_zzE = (*privtxs)[0].covariance(2,2);
      PV_xyE = (*privtxs)[0].covariance(0,1);
      PV_xzE = (*privtxs)[0].covariance(0,2);
      PV_yzE = (*privtxs)[0].covariance(1,2);
    
      PV_normchi2 = (*privtxs)[0].chi2()/(*privtxs)[0].ndof();
    
      PV_Nvtx = (*privtxs).size();
    }
    else{
      PV_x   = -999.;
      PV_y   = -999.;
      PV_z   = -999.;
      PV_xxE = -999.;
      PV_yyE = -999.;
      PV_zzE = -999.;
      PV_xyE = -999.;
      PV_xzE = -999.;
      PV_yzE = -999.;
      PV_normchi2 = -999.;
      PV_Nvtx = -999;
    }
  }

  //DT SEGMENTS
  fill_dtsegments_variables(dtsegments4D, dtGeom_, rpcGeom_, context);

  //TwinMux
  if(runOnRaw_ && hasPhiTwinMuxIn) fill_twinmuxin_variables(localTriggerTwinMuxIn);
  if(runOnRaw_ && hasPhiTwinMuxOut) fill_twinmuxout_variables(localTriggerTwinMuxOut);

  //MUONS
  if(!localDTmuons_) fill_muons_variables(MuList);
    
  //HLT
  if(!localDTmuons_) fill_hlt_variables(event,hltresults,hltevent);
  
  analyzeUnpackingRpcRecHit(event, rpcGeom_);  
  
  tree_->Fill();
  
  return;
}

void TTreeGenerator::fill_dtsegments_variables(edm::Handle<DTRecSegment4DCollection> segments4D, 
					       const DTGeometry* dtGeom_, const RPCGeometry* rpcGeom_, 
					       const edm::EventSetup& iSetup)
{

  idtsegments = 0;
  static TVectorF dummyfloat(1); dummyfloat(0) = -999.;
  for (DTRecSegment4DCollection::id_iterator chambIt = segments4D->id_begin(); chambIt != segments4D->id_end(); ++chambIt){
    
    DTRecSegment4DCollection::range  range = segments4D->get(*chambIt);
    for (DTRecSegment4DCollection::const_iterator segment4D = range.first; segment4D!=range.second; ++segment4D){

      if(idtsegments >= dtsegmentsSize_) return;
      segm4D_wheel.push_back((*chambIt).wheel());
      segm4D_sector.push_back((*chambIt).sector());
      segm4D_station.push_back((*chambIt).station());
      const bool hasPhi = segment4D->hasPhi();
      const bool hasZed = segment4D->hasZed();
      segm4D_hasPhi.push_back(hasPhi);
      segm4D_hasZed.push_back(hasZed);
      segm4D_x_pos_loc.push_back(segment4D->localPosition().x());
      segm4D_y_pos_loc.push_back(segment4D->localPosition().y());
      segm4D_z_pos_loc.push_back(segment4D->localPosition().z());
      segm4D_x_dir_loc.push_back(segment4D->localDirection().x());
      segm4D_y_dir_loc.push_back(segment4D->localDirection().y());
      segm4D_z_dir_loc.push_back(segment4D->localDirection().z());
      
      if (hasPhi||hasZed){

	TVectorF hitExpectedPos(12);
	TVectorF hitExpectedWire(12);
	std::vector<DTWireId> wireVector;
	for (int kSL=1; kSL<4; kSL=kSL+1){
	  if ((*chambIt).station()==4 && kSL==2) continue; //FRC 21-12-2016
	  for (int kL=1; kL<5; kL++){
	    wireVector.push_back(DTWireId((*chambIt).wheel(),(*chambIt).station(),(*chambIt).sector(),kSL,kL,2));
	  }
	}

	int kkk=0;
	const DTChamber* mych = dtGeom_->chamber(*chambIt); 
	StraightLinePlaneCrossing segmentPlaneCrossing(((*mych).toGlobal(segment4D->localPosition())).basicVector(),((*mych).toGlobal(segment4D->localDirection())).basicVector(),anyDirection); 

	for(std::vector<DTWireId>::const_iterator wireIt = wireVector.begin(); wireIt!=wireVector.end(); ++wireIt){
	  const DTLayer* layer = dtGeom_->layer(*wireIt);
	  const DTChamber* chamber = dtGeom_->chamber(wireIt->layerId().chamberId());
	  pair<bool,Basic3DVector<float> > ppt = segmentPlaneCrossing.position(layer->surface());
	  bool success = ppt.first; // check for failure
	  int theExpWire=-999;
	  float theExpPos=999;

 if (success){ GlobalPoint segExrapolationToLayer(ppt.second);
	    LocalPoint  segPosAtZWireLayer = layer->toLocal(segExrapolationToLayer); 
	    LocalPoint  segPosAtZWireChamber = chamber->toLocal(segExrapolationToLayer); 
	    
	    //
	    if ((kkk<4 || kkk>7)&&hasPhi){
	      theExpPos = segPosAtZWireChamber.x();
	      theExpWire = layer->specificTopology().channel(segPosAtZWireLayer);
	    }
	    else if ((kkk>=4 && kkk<=7) &&hasZed){
	      theExpPos = segPosAtZWireChamber.y();     //theExpPos = segPosAtZWire.x();
	      //LocalPoint passPoint(-segPosAtZWire.y(),segPosAtZWire.x(),segPosAtZWire.z());
	      theExpWire = layer->specificTopology().channel(segPosAtZWireLayer);
	      
	    }
	  }
	  hitExpectedWire[kkk] = theExpWire;
	  hitExpectedPos[kkk] = theExpPos;
	  kkk++;
          if ((*chambIt).station()==4 && kkk==4) kkk+=4; //FRC 22-12-2016
	
	}// END iterator
	
	new ((*segm4D_hitsExpPos)[idtsegments]) TVectorF(hitExpectedPos);
	new ((*segm4D_hitsExpWire)[idtsegments]) TVectorF(hitExpectedWire);
      }
      
      else {
	new ((*segm4D_hitsExpPos)[idtsegments]) TVectorF(dummyfloat);
	new ((*segm4D_hitsExpWire)[idtsegments]) TVectorF(dummyfloat);
      }
      const GeomDet* geomDet = theTrackingGeometry->idToDet(segment4D->geographicalId());
      const GlobalPoint point_glb = geomDet->toGlobal(segment4D->localPosition());
      segm4D_cosx.push_back(point_glb.x());
      segm4D_cosy.push_back(point_glb.y());
      segm4D_cosz.push_back(point_glb.z());
      segm4D_phi.push_back(point_glb.phi());
      segm4D_theta.push_back(point_glb.theta());
      segm4D_eta.push_back(point_glb.eta());
      if(hasPhi) fill_dtphi_info(segment4D->phiSegment(),geomDet);
      else{
	segm4D_t0.push_back(-999.);
	segm4D_vdrift.push_back(-999.);
	segm4D_phinormchi2.push_back(-999.);
	segm4D_phinhits.push_back(-999);
	new ((*segm4D_phiHits_Pos)[idtsegments])    TVectorF(dummyfloat);
	new ((*segm4D_phiHits_PosCh)[idtsegments])  TVectorF(dummyfloat);
	new ((*segm4D_phiHits_PosErr)[idtsegments]) TVectorF(dummyfloat);
	new ((*segm4D_phiHits_Side)[idtsegments])   TVectorF(dummyfloat);
 	new ((*segm4D_phiHits_Wire)[idtsegments])   TVectorF(dummyfloat);
 	new ((*segm4D_phiHits_Layer)[idtsegments])  TVectorF(dummyfloat);
	new ((*segm4D_phiHits_SuperLayer)[idtsegments])  TVectorF(dummyfloat);
	new ((*segm4D_phiHits_Time)[idtsegments])    TVectorF(dummyfloat);
	new ((*segm4D_phiHits_TimeCali)[idtsegments])    TVectorF(dummyfloat);
      }
      if(hasZed) fill_dtz_info(segment4D->zSegment(),geomDet);
      else{
	segm4D_znormchi2.push_back(-999.);
	segm4D_znhits.push_back(-999);
	new ((*segm4D_zHits_Pos)[idtsegments])      TVectorF(dummyfloat);
	new ((*segm4D_zHits_PosCh)[idtsegments])    TVectorF(dummyfloat);
	new ((*segm4D_zHits_PosErr)[idtsegments])   TVectorF(dummyfloat);
	new ((*segm4D_zHits_Side)[idtsegments])     TVectorF(dummyfloat);
	new ((*segm4D_zHits_Wire)[idtsegments])     TVectorF(dummyfloat);
	new ((*segm4D_zHits_Layer)[idtsegments])    TVectorF(dummyfloat);
	new ((*segm4D_zHits_Time)[idtsegments])     TVectorF(dummyfloat);
	new ((*segm4D_zHits_TimeCali)[idtsegments]) TVectorF(dummyfloat);
      }
      
	if(dtExtrapolation_) {
	  
		short nDTsegmentOnRPC = 0;
		
		if(segments4D->size()>10 ||
		   segment4D->dimension() != 4 ||                    //check if the dimension of the segment is 4
		   segment4D->phiSegment()->recHits().size() < 5 ||  //check if there are 2 phi superlayers fired (at least 5 hits)
		   segment4D->zSegment()->recHits().size()   < 4     //check if there are at least 4 hits for the Z projection
		   ) {
		  
		  DT_segments_onRPC.push_back(nDTsegmentOnRPC);
		  
		  new ((*DT_extrapolated_OnRPC_BX)[idtsegments])       TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Loc_x)[idtsegments])    TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Loc_y)[idtsegments])    TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Loc_z)[idtsegments])    TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Loc_eta)[idtsegments])  TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Loc_phi)[idtsegments])  TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Glob_x)[idtsegments])   TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Glob_y)[idtsegments])   TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Glob_z)[idtsegments])   TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Glob_eta)[idtsegments]) TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Glob_phi)[idtsegments]) TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Region)[idtsegments])   TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Sector)[idtsegments])   TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Station)[idtsegments])  TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Layer)[idtsegments])    TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Roll)[idtsegments])     TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Ring)[idtsegments])     TVectorF(dummyfloat);
		  new ((*DT_extrapolated_OnRPC_Stripw)[idtsegments])   TVectorF(dummyfloat);
		}
		else {
		  
		  std::map<DTChamberId,int> DTSegmentCounter;
		  DTChamberId DTId = segment4D->chamberId();

		  TVectorF DTextrapolatedOnRPC_BX(2);
		  TVectorF DTextrapolatedOnRPC_Loc_x(2);
		  TVectorF DTextrapolatedOnRPC_Loc_y(2);
		  TVectorF DTextrapolatedOnRPC_Loc_z(2);
		  TVectorF DTextrapolatedOnRPC_Loc_eta(2);
		  TVectorF DTextrapolatedOnRPC_Loc_phi(2);
		  TVectorF DTextrapolatedOnRPC_Glob_x(2);
		  TVectorF DTextrapolatedOnRPC_Glob_y(2);
		  TVectorF DTextrapolatedOnRPC_Glob_z(2);
		  TVectorF DTextrapolatedOnRPC_Glob_eta(2);
		  TVectorF DTextrapolatedOnRPC_Glob_phi(2);
		  TVectorF DTextrapolatedOnRPC_Region(2);
		  TVectorF DTextrapolatedOnRPC_Sector(2);
		  TVectorF DTextrapolatedOnRPC_Station(2);
		  TVectorF DTextrapolatedOnRPC_Layer(2);
		  TVectorF DTextrapolatedOnRPC_Roll(2);
		  TVectorF DTextrapolatedOnRPC_Ring(2);
		  TVectorF DTextrapolatedOnRPC_Stripw(2);

		  // 			if(DTSegmentCounter[DTId]!=1 || DTId.station()==4){
		  // 				std::cout<<"DT \t \t More than one segment in this chamber, or we are in Station 4"<<std::endl;
		  // 				continue;
		  // 			}
			
		  int dtWheel = DTId.wheel(); //maybe unuseful
		  int dtStation = DTId.station();
		  int dtSector = DTId.sector();      
		  
		  LocalPoint segmentPosition= segment4D->localPosition();
		  LocalVector segmentDirection=segment4D->localDirection();
		  
		  const GeomDet* gdet=dtGeom_->idToDet(segment4D->geographicalId());
		  const BoundPlane & DTSurface = gdet->surface();	 
		  		  
		  float Xo=segmentPosition.x();
		  float Yo=segmentPosition.y();
		  float Zo=segmentPosition.z();
		  float dx=segmentDirection.x();
		  float dy=segmentDirection.y();
		  float dz=segmentDirection.z();
		  
		  ObjectMap* TheObject = ObjectMap::GetInstance(iSetup);
		  DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
		  std::set<RPCDetId> rollsForThisDT = TheObject->GetInstance(iSetup)->GetRolls(theindex);
		  assert(rollsForThisDT.size()>=1);		
		  
		  for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisDT.begin();iteraRoll != rollsForThisDT.end(); iteraRoll++){
		    const RPCRoll* rollasociated = rpcGeom_->roll(*iteraRoll);
		    RPCDetId rpcId = rollasociated->id();
		    const BoundPlane & RPCSurface = rollasociated->surface(); 
		    
		    RPCGeomServ rpcsrv(rpcId);
		    std::string nameRoll = rpcsrv.name();
		    
		    GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
		    
		    LocalPoint CenterRollinDTFrame = DTSurface.toLocal(CenterPointRollGlobal);
		    
		    float D=CenterRollinDTFrame.z();
		    
		    float X=Xo+dx*D/dz;
		    float Y=Yo+dy*D/dz;
		    float Z=D;
		    
		    const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
		    LocalPoint xmin = top_->localPosition(0.);
		    LocalPoint xmax = top_->localPosition((float)rollasociated->nstrips());
		    float rsize = fabs( xmax.x()-xmin.x() );
		    float stripl = top_->stripLength();
		    
		    float stripw = top_->pitch();
		    
		    float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));
		    
		    if(extrapolatedDistance<=MaxD){ 
		      GlobalPoint GlobalPointExtrapolated = DTSurface.toGlobal(LocalPoint(X,Y,Z));
		      LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);
		      
		      if(fabs(PointExtrapolatedRPCFrame.z()) < 1. && 
			 fabs(PointExtrapolatedRPCFrame.x()) < rsize*eyr && 
			 fabs(PointExtrapolatedRPCFrame.y()) < stripl*eyr){
			
			LocalVector segmentDirection=segment4D->localDirection();
			float dx=segmentDirection.x();
			float dz=segmentDirection.z();
			float cosal = dx/sqrt(dx*dx+dz*dz);
			float angle = acos(cosal)*180/3.1415926;
			
			RPCRecHit RPCPoint(rpcId,0,LocalPoint(PointExtrapolatedRPCFrame.x(),PointExtrapolatedRPCFrame.y(),angle));
			
			float dt_loc_x = PointExtrapolatedRPCFrame.x();
			float dt_loc_y = PointExtrapolatedRPCFrame.y();
			float dt_loc_z = PointExtrapolatedRPCFrame.z();
			float dt_loc_phi = PointExtrapolatedRPCFrame.phi();
			float dt_loc_eta = PointExtrapolatedRPCFrame.eta();
			float dt_glob_x = GlobalPointExtrapolated.x();
			float dt_glob_y = GlobalPointExtrapolated.y();
			float dt_glob_z = GlobalPointExtrapolated.z();
			float dt_glob_eta = GlobalPointExtrapolated.eta();
			float dt_glob_phi = GlobalPointExtrapolated.phi();
			
			DTextrapolatedOnRPC_BX(nDTsegmentOnRPC) = RPCPoint.BunchX();
			
			DTextrapolatedOnRPC_Loc_x(nDTsegmentOnRPC) = dt_loc_x;
			DTextrapolatedOnRPC_Loc_y(nDTsegmentOnRPC) = dt_loc_y;
			DTextrapolatedOnRPC_Loc_z(nDTsegmentOnRPC) = dt_loc_z;
			DTextrapolatedOnRPC_Loc_eta(nDTsegmentOnRPC) = dt_loc_eta;
			DTextrapolatedOnRPC_Loc_phi(nDTsegmentOnRPC) = dt_loc_phi;
			
			DTextrapolatedOnRPC_Glob_x(nDTsegmentOnRPC) = dt_glob_x;
			DTextrapolatedOnRPC_Glob_y(nDTsegmentOnRPC) = dt_glob_y;
			DTextrapolatedOnRPC_Glob_z(nDTsegmentOnRPC) = dt_glob_z;
			DTextrapolatedOnRPC_Glob_eta(nDTsegmentOnRPC) = dt_glob_eta;
			DTextrapolatedOnRPC_Glob_phi(nDTsegmentOnRPC) = dt_glob_phi;
			
			DTextrapolatedOnRPC_Region(nDTsegmentOnRPC) = rpcId.region();
			DTextrapolatedOnRPC_Sector(nDTsegmentOnRPC) = rpcId.sector();
			DTextrapolatedOnRPC_Station(nDTsegmentOnRPC) = rpcId.station();
			DTextrapolatedOnRPC_Layer(nDTsegmentOnRPC) = rpcId.layer();
			DTextrapolatedOnRPC_Roll(nDTsegmentOnRPC) = rpcId.roll();
			DTextrapolatedOnRPC_Ring(nDTsegmentOnRPC) = rpcId.ring();
			DTextrapolatedOnRPC_Stripw(nDTsegmentOnRPC) = stripw;

			nDTsegmentOnRPC++;

		      }else {
		      }//Condition for the right match
		    }else{
		    }//Distance is too big

		  }//loop over all the rolls asociated

		  DTextrapolatedOnRPC_BX.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Loc_x.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Loc_y.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Loc_z.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Loc_eta.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Loc_phi.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Glob_x.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Glob_y.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Glob_z.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Glob_eta.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Glob_phi.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Region.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Sector.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Station.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Layer.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Roll.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Ring.ResizeTo(nDTsegmentOnRPC);
		  DTextrapolatedOnRPC_Stripw.ResizeTo(nDTsegmentOnRPC);

		  DT_segments_onRPC.push_back(nDTsegmentOnRPC);
		  
		  new ((*DT_extrapolated_OnRPC_BX)[idtsegments])       TVectorF(DTextrapolatedOnRPC_BX);
		  new ((*DT_extrapolated_OnRPC_Loc_x)[idtsegments])    TVectorF(DTextrapolatedOnRPC_Loc_x);
		  new ((*DT_extrapolated_OnRPC_Loc_y)[idtsegments])    TVectorF(DTextrapolatedOnRPC_Loc_y);
		  new ((*DT_extrapolated_OnRPC_Loc_z)[idtsegments])    TVectorF(DTextrapolatedOnRPC_Loc_x);
		  new ((*DT_extrapolated_OnRPC_Loc_eta)[idtsegments])  TVectorF(DTextrapolatedOnRPC_Loc_eta);
		  new ((*DT_extrapolated_OnRPC_Loc_phi)[idtsegments])  TVectorF(DTextrapolatedOnRPC_Loc_phi);
		  new ((*DT_extrapolated_OnRPC_Glob_x)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Glob_x);
		  new ((*DT_extrapolated_OnRPC_Glob_y)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Glob_y);
		  new ((*DT_extrapolated_OnRPC_Glob_z)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Glob_z);
		  new ((*DT_extrapolated_OnRPC_Glob_eta)[idtsegments]) TVectorF(DTextrapolatedOnRPC_Glob_eta);
		  new ((*DT_extrapolated_OnRPC_Glob_phi)[idtsegments]) TVectorF(DTextrapolatedOnRPC_Glob_phi);
		  new ((*DT_extrapolated_OnRPC_Region)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Region);
		  new ((*DT_extrapolated_OnRPC_Sector)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Sector);
		  new ((*DT_extrapolated_OnRPC_Station)[idtsegments])  TVectorF(DTextrapolatedOnRPC_Station);
		  new ((*DT_extrapolated_OnRPC_Layer)[idtsegments])    TVectorF(DTextrapolatedOnRPC_Layer);
		  new ((*DT_extrapolated_OnRPC_Roll)[idtsegments])     TVectorF(DTextrapolatedOnRPC_Roll);
		  new ((*DT_extrapolated_OnRPC_Ring)[idtsegments])     TVectorF(DTextrapolatedOnRPC_Ring);
		  new ((*DT_extrapolated_OnRPC_Stripw)[idtsegments])   TVectorF(DTextrapolatedOnRPC_Stripw);

		} //else: segments4D->size()<10
	} // if: do you want extrapolation ? ?
		  
	++idtsegments;

    } // for on segment4D
  } // for on chambIt
 
  return;
  
} // close the function

void TTreeGenerator::fill_dtphi_info(const DTChamberRecSegment2D* phiSeg, const GeomDet* chamb)
{
  std::vector<DTRecHit1D> phirecHitslist = phiSeg->specificRecHits();
  //segment information
  segm4D_t0.push_back(phiSeg->t0());
  segm4D_vdrift.push_back(phiSeg->vDrift());
  segm4D_phinormchi2.push_back(phiSeg->chi2()/phiSeg->degreesOfFreedom());
  //rechits information
  const int nphirecHits = phirecHitslist.size();
  segm4D_phinhits.push_back(nphirecHits);
  TVectorF phiPosRechits(nphirecHits);
  TVectorF phiPosChRechits(nphirecHits);
  TVectorF phiPosErrRechits(nphirecHits);
  TVectorF phiSideRechits(nphirecHits);
  TVectorF phiwireRechits(nphirecHits);
  TVectorF philayerRechits(nphirecHits);
  TVectorF phisuperlayerRechits(nphirecHits);
  TVectorF phiTimeRechits(nphirecHits);
  TVectorF phiTimeCaliRechits(nphirecHits);
  int rechitscounter = 0;
  for(std::vector<DTRecHit1D>::const_iterator recHitsIt = phirecHitslist.begin(); recHitsIt!=phirecHitslist.end(); ++recHitsIt){
    const GeomDet * layer = theTrackingGeometry->idToDet(recHitsIt->wireId().layerId());
    phiPosRechits(rechitscounter)    = recHitsIt->localPosition().x();
    phiPosChRechits(rechitscounter)  = chamb->toLocal(layer->toGlobal(recHitsIt->localPosition())).x();
    phiPosErrRechits(rechitscounter) = recHitsIt->localPositionError().xx();
    phiSideRechits(rechitscounter)   = recHitsIt->lrSide();
    phiwireRechits(rechitscounter)   = recHitsIt->wireId().wire();
    philayerRechits(rechitscounter)  = recHitsIt->wireId().layerId().layer();
    phisuperlayerRechits(rechitscounter)  = recHitsIt->wireId().layerId().superlayer();
    phiTimeRechits(rechitscounter)  = recHitsIt->digiTime();
    float ttrig = theSync->offset(recHitsIt->wireId());
    phiTimeCaliRechits(rechitscounter)  = recHitsIt->digiTime()-ttrig;
    rechitscounter++;
  }
  new ((*segm4D_phiHits_Pos)[idtsegments])    TVectorF(phiPosRechits);
  new ((*segm4D_phiHits_PosCh)[idtsegments])  TVectorF(phiPosChRechits);
  new ((*segm4D_phiHits_PosErr)[idtsegments]) TVectorF(phiPosErrRechits);
  new ((*segm4D_phiHits_Side)[idtsegments])   TVectorF(phiSideRechits);
  new ((*segm4D_phiHits_Wire)[idtsegments])   TVectorF(phiwireRechits);
  new ((*segm4D_phiHits_Layer)[idtsegments])  TVectorF(philayerRechits);
  new ((*segm4D_phiHits_SuperLayer)[idtsegments])  TVectorF(phisuperlayerRechits);
  new ((*segm4D_phiHits_Time)[idtsegments])  TVectorF(phiTimeRechits);
  new ((*segm4D_phiHits_TimeCali)[idtsegments])  TVectorF(phiTimeCaliRechits);
  return;
}

void TTreeGenerator::fill_dtz_info(const DTSLRecSegment2D* zSeg,  const GeomDet* chamb)
{
  std::vector<DTRecHit1D> zrecHitslist = zSeg->specificRecHits();
  segm4D_znormchi2.push_back(zSeg->chi2()/zSeg->degreesOfFreedom());
  //rechits information
  const int nzrecHits = zrecHitslist.size();
  segm4D_znhits.push_back(nzrecHits);
  TVectorF zPosRechits(nzrecHits);
  TVectorF zPosChRechits(nzrecHits);
  TVectorF zPosErrRechits(nzrecHits);
  TVectorF zSideRechits(nzrecHits);
  TVectorF zwireRechits(nzrecHits);
  TVectorF zlayerRechits(nzrecHits);
  TVectorF zTimeRechits(nzrecHits);
  TVectorF zTimeCaliRechits(nzrecHits);
  int rechitscounter = 0;
  for(std::vector<DTRecHit1D>::const_iterator recHitsIt = zrecHitslist.begin(); recHitsIt!=zrecHitslist.end(); ++recHitsIt){
    const GeomDet * layer = theTrackingGeometry->idToDet(recHitsIt->wireId().layerId());
    zPosRechits(rechitscounter)    = recHitsIt->localPosition().y();
    zPosChRechits(rechitscounter)  = chamb->toLocal(layer->toGlobal(recHitsIt->localPosition())).y();
    zPosErrRechits(rechitscounter) = recHitsIt->localPositionError().yy();
    zSideRechits(rechitscounter)   = recHitsIt->lrSide();
    zwireRechits(rechitscounter)   = recHitsIt->wireId().wire();
    zlayerRechits(rechitscounter)  = recHitsIt->wireId().layerId().layer();
    zTimeRechits(rechitscounter)   = recHitsIt->digiTime();
    float ttrig = theSync->offset(recHitsIt->wireId());
    zTimeCaliRechits(rechitscounter)   = recHitsIt->digiTime()-ttrig;
    rechitscounter++;
  }
  new ((*segm4D_zHits_Pos)[idtsegments])    TVectorF(zPosRechits);
  new ((*segm4D_zHits_PosCh)[idtsegments])  TVectorF(zPosChRechits);
  new ((*segm4D_zHits_PosErr)[idtsegments]) TVectorF(zPosErrRechits);
  new ((*segm4D_zHits_Side)[idtsegments])   TVectorF(zSideRechits);
  new ((*segm4D_zHits_Wire)[idtsegments])   TVectorF(zwireRechits);
  new ((*segm4D_zHits_Layer)[idtsegments])  TVectorF(zlayerRechits);
  new ((*segm4D_zHits_Time)[idtsegments])   TVectorF(zTimeRechits);
  new ((*segm4D_zHits_TimeCali)[idtsegments]) TVectorF(zTimeCaliRechits);
  return;
}

void TTreeGenerator::fill_twinmuxout_variables(edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxOut)
{
  idtltTwinMuxOut = 0;
  const std::vector<L1MuDTChambPhDigi>*  phTrigs = localTriggerTwinMuxOut->getContainer();
  for(std::vector<L1MuDTChambPhDigi>::const_iterator iph = phTrigs->begin(); iph != phTrigs->end() ; ++iph){
    if(idtltTwinMuxOut >= dtltTwinMuxOutSize_) break;
    if (iph->code()!=7){
      ltTwinMuxOut_wheel.push_back(iph->whNum());
      ltTwinMuxOut_sector.push_back(iph->scNum() + 1); // DTTF[0-11] -> DT[1-12] Sector Numbering
      ltTwinMuxOut_station.push_back(iph->stNum());
      ltTwinMuxOut_quality.push_back(iph->code());
      ltTwinMuxOut_rpcbit.push_back(iph->RpcBit());
      ltTwinMuxOut_bx.push_back(iph->bxNum());
      ltTwinMuxOut_phi.push_back(iph->phi());
      ltTwinMuxOut_phiB.push_back(iph->phiB());
      ltTwinMuxOut_is2nd.push_back(iph->Ts2Tag());
      idtltTwinMuxOut++;
    }
  }
  return;
}

void TTreeGenerator::fill_twinmuxin_variables(edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxIn)
{

  idtltTwinMuxIn=0;
  const std::vector<L1MuDTChambPhDigi>*  phTrigs = localTriggerTwinMuxIn->getContainer();
  for(std::vector<L1MuDTChambPhDigi>::const_iterator iph = phTrigs->begin(); iph != phTrigs->end() ; ++iph){
    if(idtltTwinMuxIn >= dtltTwinMuxInSize_) break;
    if (iph->code()!=7){
      ltTwinMuxIn_wheel.push_back(iph->whNum());
      ltTwinMuxIn_sector.push_back(iph->scNum() + 1); // DTTF[0-11] -> DT[1-12] Sector Numbering
      ltTwinMuxIn_station.push_back(iph->stNum());
      ltTwinMuxIn_quality.push_back(iph->code());
      if (iph->Ts2Tag()==1) ltTwinMuxIn_bx.push_back(iph->bxNum()-1);
      else ltTwinMuxIn_bx.push_back(iph->bxNum());
      ltTwinMuxIn_phi.push_back(iph->phi());
      ltTwinMuxIn_phiB.push_back(iph->phiB());
      ltTwinMuxIn_is2nd.push_back(iph->Ts2Tag());
      idtltTwinMuxIn++;
    }
  }
  return;
}

void TTreeGenerator::fill_muons_variables(edm::Handle<reco::MuonCollection> MuList)
{
  imuons = 0;
  for (reco::MuonCollection::const_iterator nmuon = MuList->begin(); nmuon != MuList->end(); ++nmuon){
    if(!(nmuon->isStandAloneMuon())) continue;
    if(imuons >= STAMuSize_) break;
    const reco::TrackRef mutrackref = nmuon->outerTrack();

    STAMu_isMuGlobal.push_back(nmuon->isGlobalMuon());
    STAMu_isMuTracker.push_back(nmuon->isTrackerMuon());
    STAMu_numberOfChambers.push_back(nmuon->numberOfChambers());
    STAMu_numberOfMatches.push_back(nmuon->numberOfMatches());
    STAMu_numberOfHits.push_back(mutrackref->numberOfValidHits());

    Mu_px_mu.push_back(nmuon->px());
    Mu_py_mu.push_back(nmuon->py());
    Mu_pz_mu.push_back(nmuon->pz());
    Mu_phi_mu.push_back(nmuon->phi());
    Mu_eta_mu.push_back(nmuon->eta());

    STAMu_recHitsSize.push_back(mutrackref->recHitsSize());
    STAMu_normchi2Mu.push_back(mutrackref->chi2()/mutrackref->ndof());
    STAMu_chargeMu.push_back(mutrackref->charge());
    STAMu_dxyMu.push_back(mutrackref->dxy(beamspot.position()));
    STAMu_dzMu.push_back(mutrackref->dz(beamspot.position()));
    int segmIndex = 0;
    int segmWord = 0;

    std::vector<int> segmIndex_container;
    for (trackingRecHit_iterator recMu = mutrackref->recHitsBegin(); recMu!=mutrackref->recHitsEnd(); recMu++){
      DetId detid = (*recMu)->geographicalId(); 
      if(detid.subdetId() != MuonSubdetId::DT) continue;
      DTChamberId recChamb(detid);
      const short recWheel   = recChamb.wheel();
      const short recSector  = recChamb.sector();
      const short recStation = recChamb.station();
      //loop over the saved segments and find the position of the rechits
      //This is the quickest way to do this search: find the sector (highest number of
      //combinations), loop over the find iterator, and search for wheel and stations
      std::vector<short>::iterator sectorIt = std::find(segm4D_sector.begin(),segm4D_sector.end(),recSector);
      while(sectorIt != segm4D_sector.end()){
	segmIndex = (short) distance(segm4D_sector.begin(),sectorIt);
	if(recWheel == segm4D_wheel.at(segmIndex) && recStation == segm4D_station.at(segmIndex))
	  if(find(segmIndex_container.begin(),segmIndex_container.end(),segmIndex) == segmIndex_container.end()){
	    segmIndex_container.push_back(segmIndex);
	    segmWord |= (1 << segmIndex);
	  }
	sectorIt = std::find(sectorIt+1,segm4D_sector.end(),recSector);
      }
    }
    STAMu_segmIndex.push_back(segmWord);
    if(nmuon->isGlobalMuon() & AnaTrackGlobalMu_) { 

      const reco::TrackRef glbmutrackref = nmuon->innerTrack();
      GLBMu_normchi2Mu.push_back(glbmutrackref->chi2()/glbmutrackref->ndof());
      GLBMu_dxyMu.push_back(glbmutrackref->dxy(beamspot.position()));
      GLBMu_dzMu.push_back(glbmutrackref->dz(beamspot.position()));
      
      GLBMu_numberOfPixelHits.push_back(glbmutrackref->hitPattern().numberOfValidPixelHits());
      GLBMu_numberOfTrackerHits.push_back(glbmutrackref->hitPattern().numberOfValidTrackerHits());
      
      GLBMu_tkIsoR03.push_back(nmuon->isolationR03().sumPt);
      GLBMu_ntkIsoR03.push_back(nmuon->isolationR03().nTracks);
      GLBMu_emIsoR03.push_back(nmuon->isolationR03().emEt);
      GLBMu_hadIsoR03.push_back(nmuon->isolationR03().hadEt);

    }
    else{

      GLBMu_normchi2Mu.push_back(-999.);
      GLBMu_dxyMu.push_back(-999.);
      GLBMu_dzMu.push_back(-999.);
      GLBMu_numberOfPixelHits.push_back(-999);
      GLBMu_numberOfTrackerHits.push_back(-999);
      GLBMu_tkIsoR03.push_back(-999.);
      GLBMu_ntkIsoR03.push_back(-999.);
      GLBMu_emIsoR03.push_back(-999.);
      GLBMu_hadIsoR03.push_back(-999.);

    }

    Int_t iMatches = 0;

    int wheel[4]   = {-999, -999, -999, -999};
    int sector[4]  = {-999, -999, -999, -999};
    float x[4] = {-999., -999., -999., -999.};
    float y[4] = {-999., -999., -999., -999.};

    std::vector<int> wheels;
    std::vector<int> stations;
    std::vector<int> sectors;

    if ( nmuon->isTrackerMuon() &&
	 nmuon->isMatchesValid() )
      {
	for ( const reco::MuonChamberMatch & match : nmuon->matches() )
	  {

	    if (match.edgeX < -5. &&
		match.edgeY < -5. &&
		match.id.det() == DetId::Muon &&
		match.id.subdetId() == MuonSubdetId::DT)
	      {

		DTChamberId dtId(match.id.rawId());

		int ch  = dtId.station();
		int sec = dtId.sector();
		int wh  = dtId.wheel();

		wheels.push_back(wh);
		stations.push_back(ch);
		sectors.push_back(sec);

		x[ch -1] = match.x;
		y[ch -1] = match.y;
		sector[ch -1] = sec;
		wheel[ch -1]  = wh;

		++iMatches;

	      }
	  }
      }
	    
    TRKMu_x_MB1.push_back(x[0]);
    TRKMu_y_MB1.push_back(y[0]);
    TRKMu_sector_MB1.push_back(sector[0]);
    TRKMu_wheel_MB1.push_back(wheel[0]);

    TRKMu_x_MB2.push_back(x[1]);
    TRKMu_y_MB2.push_back(y[1]);
    TRKMu_sector_MB2.push_back(sector[1]);
    TRKMu_wheel_MB2.push_back(wheel[1]);

    TRKMu_x_MB3.push_back(x[2]);
    TRKMu_y_MB3.push_back(y[2]);
    TRKMu_sector_MB3.push_back(sector[2]);
    TRKMu_wheel_MB3.push_back(wheel[2]);

    TRKMu_x_MB4.push_back(x[3]);
    TRKMu_y_MB4.push_back(y[3]);
    TRKMu_sector_MB4.push_back(sector[3]);
    TRKMu_wheel_MB4.push_back(wheel[3]);

    // std::cout << "[TTreeGenerator::analyze] nMatches : " << iMatches << "\t";
    // for (auto & ch : stations) std::cout << ch << " "; std::cout << "\t";
    // for (auto & sec : sectors) std::cout << sec << " "; std::cout << "\t";
    // for (auto & wh : wheels) std::cout << wh << " "; std::cout << "\t";
    // std::cout << std::endl;

    if(nmuon->isCaloCompatibilityValid()) STAMu_caloCompatibility.push_back(nmuon->caloCompatibility());
    else STAMu_caloCompatibility.push_back(-999.);
    //extrapolate the muon to the MB2

    TrajectoryStateOnSurface tsos;
    tsos = cylExtrapTrkSam(mutrackref,500.);  // track at MB2 radius - extrapolation
    if (tsos.isValid()){
      static const float pig = acos(-1.);
      const double xx = tsos.globalPosition().x();
      const double yy = tsos.globalPosition().y();
      const double zz = tsos.globalPosition().z();
      const double rr       = sqrt(xx*xx + yy*yy);
      const double cosphi   = xx/rr;
      const double abspseta = -log(tan(atan(fabs(rr/zz))/2.));
      STAMu_z_mb2.push_back(zz);
      if (yy>=0) STAMu_phi_mb2.push_back(acos(cosphi));
      else       STAMu_phi_mb2.push_back(2*pig-acos(cosphi));
      if (zz>=0) STAMu_pseta_mb2.push_back(abspseta);
      else       STAMu_pseta_mb2.push_back(-abspseta);
    }
    else{
      STAMu_z_mb2.push_back(-999.);
      STAMu_phi_mb2.push_back(-999.);
      STAMu_pseta_mb2.push_back(-999.);
    }
    imuons++;
  }
  return;
}

void TTreeGenerator::fill_hlt_variables(const edm::Event &e, 
					edm::Handle<edm::TriggerResults> hltresults,
					edm::Handle<trigger::TriggerEvent> hltevent)
{
  const edm::TriggerNames TrigNames_ = e.triggerNames(*hltresults);
  const int ntrigs = hltresults->size();

  ihltPaths = 0; 

  for (int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if (hltresults->accept(itr)) {
      hlt_path.push_back(trigName);
      ++ihltPaths;      
    }
  }

  const trigger::size_type nFilters(hltevent->sizeFilters());

  ihltFilters = 0; 

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter)
    {

      if(ihltFilters >= hltSize_) break;

      std::string filterTag = hltevent->filterTag(iFilter).encode();

      if ( ( filterTag.find("Mu22") != std::string::npos ||
	     filterTag.find("Mu25") != std::string::npos ) &&
	   filterTag.find("Filtered") != std::string::npos &&
	   filterTag.find("L3")       != std::string::npos
	 )
	{

	  trigger::Keys objectKeys = hltevent->filterKeys(iFilter);
	  const trigger::TriggerObjectCollection& triggerObjects(hltevent->getObjects());

	  for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey)
	    {
	      trigger::size_type objKey = objectKeys.at(iKey);
	      const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

	      ++ihltFilters;

	      float trigObjPt = triggerObj.pt();
	      float trigObjEta = triggerObj.eta();
	      float trigObjPhi = triggerObj.phi();

	      hlt_filter.push_back(filterTag);
	      hlt_filter_pt.push_back(trigObjPt);
	      hlt_filter_eta.push_back(trigObjEta);
	      hlt_filter_phi.push_back(trigObjPhi);
	    }
	}
    }
      
  return;

}

void TTreeGenerator::analyzeUnpackingRpcRecHit(const edm::Event& event, const RPCGeometry* rpcGeom_)
{
  edm::Handle<RPCRecHitCollection> UnpackingRpcHits;
  event.getByToken(UnpackingRpcRecHitToken_, UnpackingRpcHits);
  RPCRecHitCollection::const_iterator recHit;
  irpcrechits_TwinMux=0;
  for(recHit = UnpackingRpcHits->begin(); recHit != UnpackingRpcHits->end(); recHit++){ 
    int cls = recHit->clusterSize();
    int firststrip = recHit->firstClusterStrip();
    int bx = recHit->BunchX();
    RPCDetId rpcId = recHit->rpcId();
    int region = rpcId.region();
    int stat = rpcId.station();
    int sect = rpcId.sector();
    int layer = rpcId.layer();
    int subsector = rpcId.subsector();
    int roll = rpcId.roll();
    int ring = rpcId.ring();

    LocalPoint recHitPos=recHit->localPosition();
    float loc_x = recHitPos.x();
    float loc_y = recHitPos.y();
    float loc_z = recHitPos.z();
	float loc_phi = recHitPos.phi();
	float loc_eta = recHitPos.eta();

    const RPCRoll * rollasociated = rpcGeom_->roll(rpcId);
    const BoundPlane & RPCSurface = rollasociated->surface(); 
    GlobalPoint RPCRecHitInGlobal = RPCSurface.toGlobal(recHitPos);
    float glob_x = RPCRecHitInGlobal.x();
    float glob_y = RPCRecHitInGlobal.y();
    float glob_z = RPCRecHitInGlobal.z();
    float glob_eta = RPCRecHitInGlobal.eta();
    float glob_phi = RPCRecHitInGlobal.phi();
    
	if(OnlyBarrel_ && region != 0) continue;
	    
    RpcRechit_TwinMux_region.push_back(region);
    RpcRechit_TwinMux_clusterSize.push_back(cls);
    RpcRechit_TwinMux_strip.push_back(firststrip);
    RpcRechit_TwinMux_bx.push_back(bx);
    RpcRechit_TwinMux_station.push_back(stat);
    RpcRechit_TwinMux_sector.push_back(sect);
    RpcRechit_TwinMux_layer.push_back(layer);
    RpcRechit_TwinMux_subsector.push_back(subsector);
    RpcRechit_TwinMux_roll.push_back(roll);
    RpcRechit_TwinMux_ring.push_back(ring);
    RpcRechit_TwinMux_Loc_x.push_back(loc_x);
    RpcRechit_TwinMux_Loc_y.push_back(loc_y);
    RpcRechit_TwinMux_Loc_z.push_back(loc_z);
	RpcRechit_TwinMux_Loc_eta.push_back(loc_eta);
	RpcRechit_TwinMux_Loc_phi.push_back(loc_phi);    


    RpcRechit_TwinMux_Glob_x.push_back(glob_x);
    RpcRechit_TwinMux_Glob_y.push_back(glob_y);
    RpcRechit_TwinMux_Glob_z.push_back(glob_z);
	RpcRechit_TwinMux_Glob_eta.push_back(glob_eta);
	RpcRechit_TwinMux_Glob_phi.push_back(glob_phi);
    irpcrechits_TwinMux++;
  }
  return;
}

void TTreeGenerator::beginJob()
{
  outFile = new TFile(outFile_.c_str(), "RECREATE", "");
  outFile->cd();

  tree_ = new TTree ("DTTree", "CMSSW DT tree");

  //Event info
  tree_->Branch("runnumber",&runnumber,"runnumber/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
  tree_->Branch("eventNumber",&eventNumber,"eventNumber/I");
  tree_->Branch("timestamp",&timestamp,"timestamp/F");
  tree_->Branch("bunchXing",&bunchXing,"bunchXing/I");
  tree_->Branch("orbitNum",&orbitNum,"orbitNum/I");

  //Primary vertex
  tree_->Branch("PV_x",&PV_x,"PV_x/F");
  tree_->Branch("PV_y",&PV_y,"PV_y/F");
  tree_->Branch("PV_z",&PV_z,"PV_z/F");

  tree_->Branch("PV_xxE",&PV_xxE,"PV_xxE/F");
  tree_->Branch("PV_yyE",&PV_yyE,"PV_yyE/F");
  tree_->Branch("PV_zzE",&PV_zzE,"PV_zzE/F");
  tree_->Branch("PV_xyE",&PV_xyE,"PV_xyE/F");
  tree_->Branch("PV_xzE",&PV_xzE,"PV_xzE/F");
  tree_->Branch("PV_yzE",&PV_yzE,"PV_yzE/F");

  tree_->Branch("PV_normchi2",&PV_normchi2,"PV_normch2/F");
  tree_->Branch("PV_Nvtx",&PV_Nvtx,"PV_Nvtx/F");

  //luminosity
  tree_->Branch("lumiperblock",&lumiperblock,"lumiperblock/F");
  tree_->Branch("beam1Intensity",&beam1Intensity,"beam1Intensity/F");
  tree_->Branch("beam2Intensity",&beam2Intensity,"beam2Intensity/F");

  //HLT
  tree_->Branch("hlt_path",&hlt_path,32000,-1);
  tree_->Branch("hlt_filter",&hlt_filter,32000,-1);

  tree_->Branch("hlt_filter_phi",&hlt_filter_phi);
  tree_->Branch("hlt_filter_eta",&hlt_filter_eta);
  tree_->Branch("hlt_filter_pt", &hlt_filter_pt);

  //DT segment variables
  tree_->Branch("dtsegm4D_wheel",&segm4D_wheel);
  tree_->Branch("dtsegm4D_sector",&segm4D_sector);
  tree_->Branch("dtsegm4D_station",&segm4D_station);

  tree_->Branch("dtsegm4D_hasPhi",&segm4D_hasPhi);
  tree_->Branch("dtsegm4D_hasZed",&segm4D_hasZed);
  tree_->Branch("dtsegm4D_x_pos_loc",&segm4D_x_pos_loc);
  tree_->Branch("dtsegm4D_y_pos_loc",&segm4D_y_pos_loc);
  tree_->Branch("dtsegm4D_z_pos_loc",&segm4D_z_pos_loc);
  tree_->Branch("dtsegm4D_x_dir_loc",&segm4D_x_dir_loc);
  tree_->Branch("dtsegm4D_y_dir_loc",&segm4D_y_dir_loc);
  tree_->Branch("dtsegm4D_z_dir_loc",&segm4D_z_dir_loc);
  tree_->Branch("dtsegm4D_cosx",&segm4D_cosx);
  tree_->Branch("dtsegm4D_cosy",&segm4D_cosy);
  tree_->Branch("dtsegm4D_cosz",&segm4D_cosz);
  tree_->Branch("dtsegm4D_phi",&segm4D_phi);
  tree_->Branch("dtsegm4D_theta",&segm4D_theta);
  tree_->Branch("dtsegm4D_eta",&segm4D_eta);
  tree_->Branch("dtsegm4D_t0",&segm4D_t0);
  tree_->Branch("dtsegm4D_vdrift",&segm4D_vdrift);
  tree_->Branch("dtsegm4D_phinormchisq",&segm4D_phinormchi2);
  tree_->Branch("dtsegm4D_phinhits",&segm4D_phinhits);
  tree_->Branch("dtsegm4D_znormchisq",&segm4D_znormchi2);
  tree_->Branch("dtsegm4D_znhits",&segm4D_znhits);

  tree_->Branch("dtsegm4D_hitsExpPos", &segm4D_hitsExpPos, 2048000,0);
  tree_->Branch("dtsegm4D_hitsExpWire", &segm4D_hitsExpWire, 2048000,0);

  //rechits info
  tree_->Branch("dtsegm4D_phi_hitsPos",&segm4D_phiHits_Pos,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsPosCh",&segm4D_phiHits_PosCh,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsPosErr",&segm4D_phiHits_PosErr,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsSide",&segm4D_phiHits_Side,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsWire",&segm4D_phiHits_Wire,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsLayer",&segm4D_phiHits_Layer,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsSuperLayer",&segm4D_phiHits_SuperLayer,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsTime",&segm4D_phiHits_Time,2048000,0);
  tree_->Branch("dtsegm4D_phi_hitsTimeCali",&segm4D_phiHits_TimeCali,2048000,0);

  tree_->Branch("dtsegm4D_z_hitsPos",&segm4D_zHits_Pos,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsPosCh",&segm4D_zHits_PosCh,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsPosErr",&segm4D_zHits_PosErr,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsSide",&segm4D_zHits_Side,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsWire",&segm4D_zHits_Wire,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsLayer",&segm4D_zHits_Layer,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsTime",&segm4D_zHits_Time,2048000,0);
  tree_->Branch("dtsegm4D_z_hitsTimeCali",&segm4D_zHits_TimeCali,2048000,0);
  
  tree_->Branch("NDTsegmentonRPC", &DT_segments_onRPC);

  tree_->Branch("DTextrapolatedOnRPCBX", &DT_extrapolated_OnRPC_BX,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLocX", &DT_extrapolated_OnRPC_Loc_x,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLocY", &DT_extrapolated_OnRPC_Loc_y,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLocZ", &DT_extrapolated_OnRPC_Loc_z,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLocEta", &DT_extrapolated_OnRPC_Loc_eta,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLocPhi", &DT_extrapolated_OnRPC_Loc_phi,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCGlobX", &DT_extrapolated_OnRPC_Glob_x,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCGlobY", &DT_extrapolated_OnRPC_Glob_y,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCGlobZ", &DT_extrapolated_OnRPC_Glob_z,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCGlobEta", &DT_extrapolated_OnRPC_Glob_eta,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCGlobPhi", &DT_extrapolated_OnRPC_Glob_phi,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCRegion", &DT_extrapolated_OnRPC_Region,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCSector", &DT_extrapolated_OnRPC_Sector,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCStation", &DT_extrapolated_OnRPC_Station,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCLayer", &DT_extrapolated_OnRPC_Layer,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCRoll", &DT_extrapolated_OnRPC_Roll,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCRing", &DT_extrapolated_OnRPC_Ring,2048000,0);
  tree_->Branch("DTextrapolatedOnRPCStripw", &DT_extrapolated_OnRPC_Stripw,2048000,0);

  //Twinmux variables
  tree_->Branch("ltTwinMuxIn_wheel",&ltTwinMuxIn_wheel);
  tree_->Branch("ltTwinMuxIn_sector",&ltTwinMuxIn_sector);
  tree_->Branch("ltTwinMuxIn_station",&ltTwinMuxIn_station);
  tree_->Branch("ltTwinMuxIn_quality",&ltTwinMuxIn_quality);
  tree_->Branch("ltTwinMuxIn_bx",&ltTwinMuxIn_bx);
  tree_->Branch("ltTwinMuxIn_phi",&ltTwinMuxIn_phi);
  tree_->Branch("ltTwinMuxIn_phiB",&ltTwinMuxIn_phiB);
  tree_->Branch("ltTwinMuxIn_is2nd",&ltTwinMuxIn_is2nd);

  tree_->Branch("ltTwinMuxOut_wheel",&ltTwinMuxOut_wheel);
  tree_->Branch("ltTwinMuxOut_sector",&ltTwinMuxOut_sector);
  tree_->Branch("ltTwinMuxOut_station",&ltTwinMuxOut_station);
  tree_->Branch("ltTwinMuxOut_quality",&ltTwinMuxOut_quality);
  tree_->Branch("ltTwinMuxOut_rpcbit",&ltTwinMuxOut_rpcbit);
  tree_->Branch("ltTwinMuxOut_bx",&ltTwinMuxOut_bx);
  tree_->Branch("ltTwinMuxOut_phi",&ltTwinMuxOut_phi);
  tree_->Branch("ltTwinMuxOut_phiB",&ltTwinMuxOut_phiB);
  tree_->Branch("ltTwinMuxOut_is2nd",&ltTwinMuxOut_is2nd);

  //muon variables
  tree_->Branch("Mu_isMuGlobal",&STAMu_isMuGlobal);
  tree_->Branch("Mu_isMuTracker",&STAMu_isMuTracker);
  tree_->Branch("Mu_numberOfChambers_sta",&STAMu_numberOfChambers);
  tree_->Branch("Mu_numberOfMatches_sta",&STAMu_numberOfMatches);
  tree_->Branch("Mu_numberOfHits_sta",&STAMu_numberOfHits);
  tree_->Branch("Mu_segmentIndex_sta",&STAMu_segmIndex);
  tree_->Branch("Mu_px",&Mu_px_mu);
  tree_->Branch("Mu_py",&Mu_py_mu);
  tree_->Branch("Mu_pz",&Mu_pz_mu);
  tree_->Branch("Mu_phi",&Mu_phi_mu);
  tree_->Branch("Mu_eta",&Mu_eta_mu);
  tree_->Branch("Mu_recHitsSize",&STAMu_recHitsSize);
  tree_->Branch("Mu_normchi2_sta",&STAMu_normchi2Mu);
  tree_->Branch("Mu_charge",&STAMu_chargeMu);
  tree_->Branch("Mu_dxy_sta",&STAMu_dxyMu);
  tree_->Branch("Mu_dz_sta",&STAMu_dzMu);
  tree_->Branch("Mu_normchi2_glb",&GLBMu_normchi2Mu);
  tree_->Branch("Mu_dxy_glb",&GLBMu_dxyMu);
  tree_->Branch("Mu_dz_glb",&GLBMu_dzMu);
  tree_->Branch("Mu_numberOfPixelHits_glb",&GLBMu_numberOfPixelHits);
  tree_->Branch("Mu_numberOfTrackerHits_glb",&GLBMu_numberOfTrackerHits);
  tree_->Branch("Mu_tkIsoR03_glb",&GLBMu_tkIsoR03);
  tree_->Branch("Mu_ntkIsoR03_glb",&GLBMu_ntkIsoR03);
  tree_->Branch("Mu_emIsoR03_glb",&GLBMu_emIsoR03);
  tree_->Branch("Mu_hadIsoR03_glb",&GLBMu_hadIsoR03);
  tree_->Branch("STAMu_caloCompatibility",&STAMu_caloCompatibility);
  tree_->Branch("Mu_z_mb2_mu",&STAMu_z_mb2);
  tree_->Branch("Mu_phi_mb2_mu",&STAMu_phi_mb2);
  tree_->Branch("Mu_pseta_mb2_mu",&STAMu_pseta_mb2);

  tree_->Branch("TRKMu_x_MB1",&TRKMu_x_MB1);
  tree_->Branch("TRKMu_y_MB1",&TRKMu_y_MB1);
  tree_->Branch("TRKMu_sector_MB1",&TRKMu_sector_MB1);
  tree_->Branch("TRKMu_wheel_MB1",&TRKMu_wheel_MB1);

  tree_->Branch("TRKMu_x_MB2",&TRKMu_x_MB2);
  tree_->Branch("TRKMu_y_MB2",&TRKMu_y_MB2);
  tree_->Branch("TRKMu_sector_MB2",&TRKMu_sector_MB2);
  tree_->Branch("TRKMu_wheel_MB2",&TRKMu_wheel_MB2);

  tree_->Branch("TRKMu_x_MB3",&TRKMu_x_MB3);
  tree_->Branch("TRKMu_y_MB3",&TRKMu_y_MB3);
  tree_->Branch("TRKMu_sector_MB3",&TRKMu_sector_MB3);
  tree_->Branch("TRKMu_wheel_MB3",&TRKMu_wheel_MB3);

  tree_->Branch("TRKMu_x_MB4",&TRKMu_x_MB4);
  tree_->Branch("TRKMu_y_MB4",&TRKMu_y_MB4);
  tree_->Branch("TRKMu_sector_MB4",&TRKMu_sector_MB4);
  tree_->Branch("TRKMu_wheel_MB4",&TRKMu_wheel_MB4);
  
  //counters
  tree_->Branch("Ndtsegments",&idtsegments,"Ndtsegments/S");
  tree_->Branch("NdtltTwinMuxOut",&idtltTwinMuxOut,"NdtltTwinMuxOut/S");
  tree_->Branch("NdtltTwinMuxIn",&idtltTwinMuxIn,"NdtltTwinMuxIn/S");
  tree_->Branch("Nmuons",&imuons,"Nmuons/S");
  tree_->Branch("NhltPaths",&ihltPaths,"Nhlt/S");
  tree_->Branch("NhltFilters",&ihltFilters,"Nhlt/S");
  tree_->Branch("NrpcRecHits",&irpcrechits,"NrpcRecHits/S");
  
  //Unpacking RPC RecHit
  tree_->Branch("NirpcrechitsTwinMux", &irpcrechits_TwinMux);
  tree_->Branch("RpcRecHitTwinMuxRegion", &RpcRechit_TwinMux_region);
  tree_->Branch("RpcRecHitTwinMuxClusterSize", &RpcRechit_TwinMux_clusterSize);
  tree_->Branch("RpcRecHitTwinMuxStrip", &RpcRechit_TwinMux_strip);
  tree_->Branch("RpcRecHitTwinMuxBx", &RpcRechit_TwinMux_bx);
  tree_->Branch("RpcRecHitTwinMuxStation", &RpcRechit_TwinMux_station);
  tree_->Branch("RpcRecHitTwinMuxSector", &RpcRechit_TwinMux_sector);
  tree_->Branch("RpcRecHitTwinMuxLayer", &RpcRechit_TwinMux_layer);
  tree_->Branch("RpcRecHitTwinMuxSubsector", &RpcRechit_TwinMux_subsector);
  tree_->Branch("RpcRecHitTwinMuxRoll", &RpcRechit_TwinMux_roll);
  tree_->Branch("RpcRecHitTwinMuxRing", &RpcRechit_TwinMux_ring);
  tree_->Branch("RpcRechitTwinMuxLocX", &RpcRechit_TwinMux_Loc_x);
  tree_->Branch("RpcRechitTwinMuxLocY", &RpcRechit_TwinMux_Loc_y);
  tree_->Branch("RpcRechitTwinMuxLocZ", &RpcRechit_TwinMux_Loc_z);
  tree_->Branch("RpcRechitTwinMuxLocEta", &RpcRechit_TwinMux_Loc_eta);
  tree_->Branch("RpcRechitTwinMuxLocPhi", &RpcRechit_TwinMux_Loc_phi);
  tree_->Branch("RpcRechitTwinMuxGlobX", &RpcRechit_TwinMux_Glob_x);
  tree_->Branch("RpcRechitTwinMuxGlobY", &RpcRechit_TwinMux_Glob_y);
  tree_->Branch("RpcRechitTwinMuxGlobZ", &RpcRechit_TwinMux_Glob_z);
  tree_->Branch("RpcRechitTwinMuxGlobEta", &RpcRechit_TwinMux_Glob_eta);
  tree_->Branch("RpcRechitTwinMuxGlobPhi", &RpcRechit_TwinMux_Glob_phi);
  
  return;
}

void TTreeGenerator::endJob()
{
  outFile->cd();
  tree_->Write();
  outFile->Close();

  return;
}

inline void TTreeGenerator::clear_Arrays()
{
  //digi variables
  // digi_wheel.clear();
  // digi_sector.clear();
  // digi_station.clear();
  // digi_sl.clear();
  // digi_layer.clear();
  // digi_wire.clear();
  // digi_time.clear();

  //DT segment variables 
  segm4D_wheel.clear();
  segm4D_sector.clear();
  segm4D_station.clear();
  segm4D_hasPhi.clear();
  segm4D_hasZed.clear();
  segm4D_x_pos_loc.clear();
  segm4D_y_pos_loc.clear();
  segm4D_z_pos_loc.clear();
  segm4D_x_dir_loc.clear();
  segm4D_y_dir_loc.clear();
  segm4D_z_dir_loc.clear();
  segm4D_cosx.clear();
  segm4D_cosy.clear();
  segm4D_cosz.clear();
  segm4D_phi.clear();
  segm4D_theta.clear();
  segm4D_eta.clear();
  segm4D_t0.clear();
  segm4D_vdrift.clear();
  segm4D_phinormchi2.clear();
  segm4D_phinhits.clear();
  segm4D_znormchi2.clear();
  segm4D_znhits.clear();

  segm4D_hitsExpPos->Clear();
  segm4D_hitsExpWire->Clear();

  segm4D_phiHits_Pos->Clear();
  segm4D_phiHits_PosCh->Clear();
  segm4D_phiHits_PosErr->Clear();
  segm4D_phiHits_Side->Clear();
  segm4D_phiHits_Wire->Clear();
  segm4D_phiHits_Layer->Clear();
  segm4D_phiHits_SuperLayer->Clear();
  segm4D_phiHits_Time->Clear();
  segm4D_phiHits_TimeCali->Clear();
  segm4D_hitsExpPos->Clear();
  segm4D_hitsExpWire->Clear();

  segm4D_zHits_Pos->Clear();
  segm4D_zHits_PosCh->Clear();
  segm4D_zHits_PosErr->Clear();
  segm4D_zHits_Side->Clear();
  segm4D_zHits_Wire->Clear();
  segm4D_zHits_Layer->Clear();
  segm4D_zHits_Time->Clear();
  segm4D_zHits_TimeCali->Clear();

  DT_segments_onRPC.clear();
  DT_extrapolated_OnRPC_BX->Clear();
  DT_extrapolated_OnRPC_Loc_x->Clear();
  DT_extrapolated_OnRPC_Loc_y->Clear();
  DT_extrapolated_OnRPC_Loc_z->Clear();
  DT_extrapolated_OnRPC_Loc_eta->Clear();
  DT_extrapolated_OnRPC_Loc_phi->Clear();
  DT_extrapolated_OnRPC_Glob_x->Clear();
  DT_extrapolated_OnRPC_Glob_y->Clear();
  DT_extrapolated_OnRPC_Glob_z->Clear();
  DT_extrapolated_OnRPC_Glob_eta->Clear();
  DT_extrapolated_OnRPC_Glob_phi->Clear();
  DT_extrapolated_OnRPC_Region->Clear();
  DT_extrapolated_OnRPC_Sector->Clear();
  DT_extrapolated_OnRPC_Station->Clear();
  DT_extrapolated_OnRPC_Layer->Clear();
  DT_extrapolated_OnRPC_Roll->Clear();
  DT_extrapolated_OnRPC_Ring->Clear();
  DT_extrapolated_OnRPC_Stripw->Clear();

  //TM Variables
  ltTwinMuxIn_wheel.clear();
  ltTwinMuxIn_sector.clear();
  ltTwinMuxIn_station.clear();
  ltTwinMuxIn_quality.clear();
  ltTwinMuxIn_bx.clear();
  ltTwinMuxIn_phi.clear();
  ltTwinMuxIn_phiB.clear();
  ltTwinMuxIn_is2nd.clear();

  ltTwinMuxOut_wheel.clear();
  ltTwinMuxOut_sector.clear();
  ltTwinMuxOut_station.clear();
  ltTwinMuxOut_quality.clear();
  ltTwinMuxOut_rpcbit.clear();
  ltTwinMuxOut_bx.clear();
  ltTwinMuxOut_phi.clear();
  ltTwinMuxOut_phiB.clear();
  ltTwinMuxOut_is2nd.clear();

  //muon variables
  STAMu_isMuGlobal.clear();
  STAMu_isMuTracker.clear();
  STAMu_numberOfChambers.clear();
  STAMu_numberOfMatches.clear();
  STAMu_numberOfHits.clear();
  STAMu_segmIndex.clear();

  Mu_px_mu.clear();
  Mu_py_mu.clear();
  Mu_pz_mu.clear();
  Mu_phi_mu.clear();
  Mu_eta_mu.clear();
  STAMu_recHitsSize.clear();
  STAMu_normchi2Mu.clear();
  STAMu_chargeMu.clear();
  STAMu_dxyMu.clear();
  STAMu_dzMu.clear();

  GLBMu_normchi2Mu.clear();
  GLBMu_dxyMu.clear();
  GLBMu_dzMu.clear();

  GLBMu_numberOfPixelHits.clear();
  GLBMu_numberOfTrackerHits.clear();

  GLBMu_tkIsoR03.clear();
  GLBMu_ntkIsoR03.clear();
  GLBMu_emIsoR03.clear();
  GLBMu_hadIsoR03.clear();

  STAMu_caloCompatibility.clear();

  STAMu_z_mb2.clear();
  STAMu_phi_mb2.clear();
  STAMu_pseta_mb2.clear();

  TRKMu_x_MB1.clear();
  TRKMu_y_MB1.clear();
  TRKMu_sector_MB1.clear();
  TRKMu_wheel_MB1.clear();

  TRKMu_x_MB2.clear();
  TRKMu_y_MB2.clear();
  TRKMu_sector_MB2.clear();
  TRKMu_wheel_MB2.clear();

  TRKMu_x_MB3.clear();
  TRKMu_y_MB3.clear();
  TRKMu_sector_MB3.clear();
  TRKMu_wheel_MB3.clear();

  TRKMu_x_MB4.clear();
  TRKMu_y_MB4.clear();
  TRKMu_sector_MB4.clear();
  TRKMu_wheel_MB4.clear();

  //HLT
  hlt_path.clear();
  hlt_filter.clear();

  hlt_filter_phi.clear();
  hlt_filter_eta.clear();
  hlt_filter_pt.clear();
  
  // Unpacking RPC rec hits
  RpcRechit_TwinMux_region.clear();
  RpcRechit_TwinMux_clusterSize.clear();
  RpcRechit_TwinMux_strip.clear();
  RpcRechit_TwinMux_bx.clear();
  RpcRechit_TwinMux_station.clear();
  RpcRechit_TwinMux_sector.clear();
  RpcRechit_TwinMux_layer.clear();
  RpcRechit_TwinMux_subsector.clear();
  RpcRechit_TwinMux_roll.clear();
  RpcRechit_TwinMux_ring.clear();
  RpcRechit_TwinMux_Loc_x.clear();
  RpcRechit_TwinMux_Loc_y.clear();
  RpcRechit_TwinMux_Loc_z.clear(); 
  RpcRechit_TwinMux_Loc_eta.clear(); 
  RpcRechit_TwinMux_Loc_phi.clear(); 
  RpcRechit_TwinMux_Glob_x.clear();
  RpcRechit_TwinMux_Glob_y.clear();
  RpcRechit_TwinMux_Glob_z.clear(); 
  RpcRechit_TwinMux_Glob_eta.clear(); 
  RpcRechit_TwinMux_Glob_phi.clear(); 
    
  return;
}

void TTreeGenerator::initialize_Tree_variables()
{
  //Event variables
  runnumber   = 0;
  lumiblock   = 0;
  eventNumber = 0;
  timestamp   = 0.;
  bunchXing   = 0;
  orbitNum    = 0;

  PV_x = 0.;
  PV_y = 0.;
  PV_z = 0.;
  PV_xxE = 0.;
  PV_yyE = 0.;
  PV_zzE = 0.;
  PV_xyE = 0.;
  PV_xzE = 0.;
  PV_yzE = 0.;
  PV_normchi2 = 0.;
  PV_Nvtx = 0.;

  lumiperblock = 0.;
  beam1Intensity = -1.;
  beam2Intensity = -1.;

  segm4D_phiHits_Pos    = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_PosCh  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_PosErr = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Side   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Wire   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Layer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_SuperLayer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_Time   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_phiHits_TimeCali   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_hitsExpPos     = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_hitsExpWire    = new TClonesArray("TVectorF",dtsegmentsSize_);

  segm4D_zHits_Pos    = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_PosCh  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_PosErr = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Side   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Wire   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Layer  = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_Time   = new TClonesArray("TVectorF",dtsegmentsSize_);
  segm4D_zHits_TimeCali   = new TClonesArray("TVectorF",dtsegmentsSize_);

  DT_extrapolated_OnRPC_BX       = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Loc_x    = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Loc_y    = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Loc_z    = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Loc_eta  = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Loc_phi  = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Glob_x   = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Glob_y   = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Glob_z   = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Glob_eta = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Glob_phi = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Region   = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Sector   = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Station  = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Layer    = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Roll     = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Ring     = new TClonesArray("TVectorF",dtsegmentsSize_);
  DT_extrapolated_OnRPC_Stripw   = new TClonesArray("TVectorF",dtsegmentsSize_);

  return;
}

TrajectoryStateOnSurface TTreeGenerator::cylExtrapTrkSam(reco::TrackRef track, const float rho) const
{
  Cylinder::PositionType pos(0.,0.,0.);
  Cylinder::RotationType rot;
  Cylinder::CylinderPointer myCylinder = Cylinder::build(pos, rot, rho);

  FreeTrajectoryState recoStart = freeTrajStateMuon(track);
  TrajectoryStateOnSurface recoProp;
  recoProp = propagatorAlong->propagate(recoStart, *myCylinder);
  if (!recoProp.isValid()) {
    recoProp = propagatorOpposite->propagate(recoStart, *myCylinder);
  }
  return recoProp;
}

FreeTrajectoryState TTreeGenerator::freeTrajStateMuon(const reco::TrackRef track) const
{
  const GlobalPoint  innerPoint(track->innerPosition().x(),track->innerPosition().y(),track->innerPosition().z());
  const GlobalVector innerVec  (track->innerMomentum().x(),track->innerMomentum().y(),track->innerMomentum().z());  
  
  FreeTrajectoryState recoStart(innerPoint, innerVec, track->charge(), &*theBField);
  
  return recoStart;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTreeGenerator);