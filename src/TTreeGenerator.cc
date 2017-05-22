//
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

#include "CalibMuon/DTDigiSync/interface/DTTTrigBaseSync.h"
#include "CalibMuon/DTDigiSync/interface/DTTTrigSyncFactory.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
#include "UserCode/DTDPGAnalysis/interface/TTreeGenerator.h"
#include "UserCode/DTDPGAnalysis/interface/DTSegmentToRPC.h"

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

  triggerTag_      = pset.getParameter<edm::InputTag>("TriggerTag");
  triggerToken_    = consumes<edm::TriggerResults>(edm::InputTag(triggerTag_));

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

// unpacking Rpc RecHit
// UnpackingRpcRecHitLabel_  = pset.getParameter<edm::InputTag>("UnpackingRpcRecHitLabel");
//   UnpackingRpcRecHitToken_  = consumes<RPCRecHitCollection>(edm::InputTag(UnpackingRpcRecHitLabel_));

//***  Moved to begin of function, needed for stablishing the digi collection to be used
//  runOnRaw_        = pset.getParameter<bool>("runOnRaw");
//  runOnSimulation_ = pset.getParameter<bool>("runOnSimulation");
//  runOnSimulationWithDigis_ = pset.getParameter<bool>("runOnSimulationWithDigis");
//  ****

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
  ihlt         = 0;
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
    //event.getByLabel(beamSpotTag_ ,recoBeamSpotHandle);  // Doesn't work after 75X
    event.getByToken(beamSpotToken_ ,recoBeamSpotHandle);
    beamspot = *recoBeamSpotHandle; 
    
    //retrieve the luminosity
    if(!runOnSimulation_) // crashes with simulation trying to get the lumiperblock   
                          // The vector<LumiScalers> but not the <LumiScalersCollection> 
    {                      
       edm::Handle<LumiScalersCollection> lumiScalers;
       //event.getByLabel(scalersSource_, lumiScalers);  // Doesn't work after 75X
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

//   edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxOut;
//   bool hasPhiTwinMuxOut=false;
//   if(runOnRaw_) hasPhiTwinMuxOut=event.getByToken(dtTrigTwinMuxOutToken_,localTriggerTwinMuxOut);
// 
//   edm::Handle<L1MuDTChambPhContainer> localTriggerTwinMuxIn;
//   bool hasPhiTwinMuxIn=false;
//   if(runOnRaw_) hasPhiTwinMuxIn=event.getByToken(dtTrigTwinMuxInToken_,localTriggerTwinMuxIn);
// 
//   edm::Handle<L1MuDTChambThContainer> localTriggerTwinMux_Th;
//   bool hasThetaTwinMux=false;
//   if(runOnRaw_) hasThetaTwinMux=event.getByToken(dtTrigTwinMux_ThToken_,localTriggerTwinMux_Th);

  edm::Handle<reco::MuonCollection> MuList;
  if(!localDTmuons_) event.getByToken(staMuToken_,MuList);

  edm::Handle<l1t::MuonBxCollection> gmt;   // legacy
  if(!localDTmuons_) event.getByToken(gmtToken_,gmt); // legacy

  edm::Handle< L1GlobalTriggerReadoutRecord > gtrc; // legacy
  if(runOnRaw_ && !localDTmuons_) event.getByToken(gtToken_, gtrc); // legacy

  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  context.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  edm::Handle<edm::TriggerResults>  hltresults;
  if(!localDTmuons_) event.getByToken(triggerToken_, hltresults); 

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

  //DIGIS
  if(runOnRaw_ && !runOnSimulation_) fill_digi_variables(dtdigis);
  else   // Added to cope with simulation including Digis
    if(runOnSimulation_ && runOnSimulationWithDigis_) fill_digi_variablesSim(dtdigisSim);


  //DT SEGMENTS
//   fill_dtsegments_variables(dtsegments4D, dtGeom_);
  fill_dtsegments_variables(dtsegments4D, dtGeom_, rpcGeom_, context);

  //CSC SEGMENTS
//   if(!localDTmuons_) fill_cscsegments_variables(cscsegments);

  //TwinMux
//   if(runOnRaw_ && hasPhiTwinMuxIn) fill_twinmuxin_variables(localTriggerTwinMuxIn);
//   if(runOnRaw_ && hasPhiTwinMuxOut) fill_twinmuxout_variables(localTriggerTwinMuxOut);
//   if(runOnRaw_ && hasThetaTwinMux) fill_twinmuxth_variables(localTriggerTwinMux_Th);
// 
  //MUONS
  if(!localDTmuons_) fill_muons_variables(MuList);

  //GMT
//   if(!localDTmuons_) fill_gmt_variables(gmt); // legacy

  //GT
//   if(runOnRaw_ && !localDTmuons_) fill_gt_variables(gtrc,menu); // legacy
    
  //HLT
//   if(!localDTmuons_) fill_hlt_variables(event,hltresults);

  // RPC
  if(!localDTmuons_) fill_rpc_variables(event,rpcHits);
  
//    analyzeBMTF(event);

//    analyzeRPCunpacking(event);

   analyzeUnpackingRpcRecHit(event, rpcGeom_);
   
//    extrapolate_DTsegment_onRPC(dtsegments4D, dtGeom_, rpcGeom_, context);
  
  tree_->Fill();

  return;
}

void TTreeGenerator::fill_digi_variables(edm::Handle<DTDigiCollection> dtdigis)
{

  idigis = 0;
  for (DTDigiCollection::DigiRangeIterator dtLayerIdIt = dtdigis->begin(); dtLayerIdIt!=dtdigis->end(); dtLayerIdIt++){
    for (DTDigiCollection::const_iterator digiIt = ((*dtLayerIdIt).second).first;digiIt!=((*dtLayerIdIt).second).second; ++digiIt){
      if(idigis >= digisSize_) return;
      digi_wheel.push_back((*dtLayerIdIt).first.wheel());
      digi_sector.push_back((*dtLayerIdIt).first.sector());
      digi_station.push_back((*dtLayerIdIt).first.station());
      digi_sl.push_back((*dtLayerIdIt).first.superLayer());
      digi_layer.push_back((*dtLayerIdIt).first.layer());
      digi_wire.push_back((*digiIt).wire());
      digi_time.push_back((*digiIt).time());
      idigis++;
    }
  }
  return;
}


void TTreeGenerator::fill_digi_variablesSim(edm::Handle<DTDigiSimLinkCollection> dtdigisSim)
{

  idigis = 0;
  for (DTDigiSimLinkCollection::DigiRangeIterator dtLayerIdIt = dtdigisSim->begin(); dtLayerIdIt!=dtdigisSim->end(); dtLayerIdIt++){
    for (DTDigiSimLinkCollection::const_iterator digiIt = ((*dtLayerIdIt).second).first;digiIt!=((*dtLayerIdIt).second).second; ++digiIt){
      if(idigis >= digisSize_) return;
      digi_wheel.push_back((*dtLayerIdIt).first.wheel());
      digi_sector.push_back((*dtLayerIdIt).first.sector());
      digi_station.push_back((*dtLayerIdIt).first.station());
      digi_sl.push_back((*dtLayerIdIt).first.superLayer());
      digi_layer.push_back((*dtLayerIdIt).first.layer());
      digi_wire.push_back((*digiIt).wire());
      digi_time.push_back((*digiIt).time());
      idigis++;
    }
  }
  return;
}
// void TTreeGenerator::fill_dtsegments_variables(edm::Handle<DTRecSegment4DCollection> segments4D, const DTGeometry* dtGeom_)
void TTreeGenerator::fill_dtsegments_variables(edm::Handle<DTRecSegment4DCollection> segments4D, const DTGeometry* dtGeom_, const RPCGeometry* rpcGeom_, const edm::EventSetup& iSetup)
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
      const GlobalVector point_glb = geomDet->toGlobal(segment4D->localDirection());
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
      idtsegments++;

  
  if(dtExtrapolation_){
  
  	DT_segment_onRPC = 0;
	
	if(segments4D->size()>10){ //Too many segments in this event we are not doing the extrapolation
		std::cout<<"-----------------------------Too many segments in this event we are not doing the extrapolation"<<std::endl;
	}
	else{   
		
	    std::map<DTChamberId,int> DTSegmentCounter;
		
		DTChamberId DTId = segment4D->chamberId();
		
// 		if(DTSegmentCounter[DTId]!=1 || DTId.station()==4){
// 			std::cout<<"DT \t \t More than one segment in this chamber, or we are in Station 4"<<std::endl;
// 			continue;
// 		}

		std::cout<<" ok sto estrapolando "<<std::endl;
				
		int dtWheel = DTId.wheel();
		int dtStation = DTId.station();
		int dtSector = DTId.sector();      
		
		LocalPoint segmentPosition= segment4D->localPosition();
		LocalVector segmentDirection=segment4D->localDirection();
		
		const GeomDet* gdet=dtGeom_->idToDet(segment4D->geographicalId());
		const BoundPlane & DTSurface = gdet->surface();	 
		
		if(segment4D->dimension()!=4) continue;  //check if the dimension of the segment is 4
		if((segment4D->phiSegment()->recHits()).size()<7) continue; //check if there are 2 phi superlayers fired (at least 7hits)
		if((segment4D->zSegment()->recHits()).size()<4) continue; //check if there are at least 4 hits for the Z projection
   
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
				float dt_glob_x = GlobalPointExtrapolated.x();
				float dt_glob_y = GlobalPointExtrapolated.y();
				float dt_glob_z = GlobalPointExtrapolated.z();
				float dt_glob_eta = GlobalPointExtrapolated.eta();
				float dt_glob_phi = GlobalPointExtrapolated.phi();
				
				DT_extrapolated_OnRPC_BX.push_back(RPCPoint.BunchX());
				
				DT_extrapolated_OnRPC_Loc_x.push_back(dt_loc_x);
				DT_extrapolated_OnRPC_Loc_y.push_back(dt_loc_y);
				DT_extrapolated_OnRPC_Loc_z.push_back(dt_loc_z);
				
				DT_extrapolated_OnRPC_Glob_x.push_back(dt_glob_x);
				DT_extrapolated_OnRPC_Glob_y.push_back(dt_glob_y);
				DT_extrapolated_OnRPC_Glob_z.push_back(dt_glob_z);
				DT_extrapolated_OnRPC_Glob_eta.push_back(dt_glob_eta);
				DT_extrapolated_OnRPC_Glob_phi.push_back(dt_glob_phi);
				
				DT_extrapolated_OnRPC_Region.push_back(rpcId.region());
				DT_extrapolated_OnRPC_Sector.push_back(rpcId.sector());
				DT_extrapolated_OnRPC_Station.push_back(rpcId.station());
				DT_extrapolated_OnRPC_Layer.push_back(rpcId.layer());
				DT_extrapolated_OnRPC_Roll.push_back(rpcId.roll());
				DT_extrapolated_OnRPC_Ring.push_back(rpcId.ring());
				DT_extrapolated_OnRPC_Stripw.push_back(stripw);
						
				DT_segment_onRPC++;
		  
			  }else {
/*			    if(debug)*/ std::cout<<"DT \t \t \t \t No the prediction is outside of this roll"<<std::endl;
			  }//Condition for the right match
			}else{
/*			  if(debug) */ std::cout<<"DT \t \t \t No, Exrtrapolation too long!, canceled"<<std::endl;
			}//Distance is too big

		}//loop over all the rolls asociated  
	} //else: segments4D->size()<10
  } // if: do you want extrapolation ?
  
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

void TTreeGenerator::fill_cscsegments_variables(edm::Handle<CSCSegmentCollection> cscsegments)
{
  icscsegments = 0;
  for (CSCSegmentCollection::id_iterator chambIt = cscsegments->id_begin(); chambIt != cscsegments->id_end(); ++chambIt){
    CSCSegmentCollection::range  range = cscsegments->get(*chambIt);
    for (CSCSegmentCollection::const_iterator cscsegment = range.first; cscsegment!=range.second; ++cscsegment){
      if(icscsegments >= cscsegmentsSize_) return;
      cscsegm_ring.push_back((*chambIt).ring());
      cscsegm_chamber.push_back((*chambIt).chamber());
      cscsegm_station.push_back((*chambIt).station());
      const GeomDet* geomDet = theTrackingGeometry->idToDet(cscsegment->geographicalId());
      const GlobalVector point_glb = geomDet->toGlobal(cscsegment->localDirection());
      cscsegm_cosx.push_back(point_glb.x());
      cscsegm_cosy.push_back(point_glb.y());
      cscsegm_cosz.push_back(point_glb.y());
      cscsegm_normchi2.push_back(cscsegment->chi2()/cscsegment->degreesOfFreedom());
      cscsegm_nRecHits.push_back(cscsegment->nRecHits());
      icscsegments++;
    }
  }
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

void TTreeGenerator::fill_twinmuxth_variables(edm::Handle<L1MuDTChambThContainer> localTriggerTwinMux_Th)
{
  idtltTwinMux_th = 0;
  const std::vector<L1MuDTChambThDigi>*  thTrigs = localTriggerTwinMux_Th->getContainer();
  for(std::vector<L1MuDTChambThDigi>::const_iterator ith = thTrigs->begin(); ith != thTrigs->end() ; ++ith){
    if(idtltTwinMux_th >= dtltTwinMuxThSize_) break;
    ltTwinMux_thWheel.push_back(ith->whNum());
    ltTwinMux_thSector.push_back(ith->scNum() + 1); // DTTF[0-11] -> DT[1-12] Sector Numbering
    ltTwinMux_thStation.push_back(ith->stNum());
    ltTwinMux_thBx.push_back(ith->bxNum());
    unsigned short int thcode=0;
    for (int pos=0; pos<7; pos++)
      if (ith->code(pos))
      thcode=thcode | (0x1<<pos);
    ltTwinMux_thHits.push_back(thcode);
    idtltTwinMux_th++;
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
      // This part  gives problems in MWGR16. All the calls glbmutrackref->... give similar error  (different Product ID number)
      // RefCore: A request to resolve a reference to a product of type 'std::vector<reco::Track>' with ProductID '2:533'
      // can not be satisfied because the product cannot be found.
      // The problem is due to the fact that the Muon collections exists, 
      // vector<reco::Muon>                    "muons"                     ""                "RECO"
      // and also the track collecions for STA, 
      // but not for GLOBAL!!!
      //  vector<reco::Track>                   "globalCosmicMuons"         ""                "RECO"    
      // AnalyzeTracksFromGlobalMuons_
      // An extra variable is add  to control from python the possibility of skiping this part of the code
      //
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

void TTreeGenerator::fill_gmt_variables(const edm::Handle<l1t::MuonBxCollection> & gmt)
{
  igmt = 0;

  for (int ibx = gmt->getFirstBX(); ibx <= gmt->getLastBX(); ++ibx) 
    {
      for (auto l1MuIt = gmt->begin(ibx); l1MuIt != gmt->end(ibx); ++l1MuIt)
	{

	  gmt_bx.push_back(ibx);
	  gmt_phi.push_back(l1MuIt->phi());
	  gmt_eta.push_back(l1MuIt->eta());
	  gmt_pt.push_back(l1MuIt->pt());
	  gmt_qual.push_back(l1MuIt->hwQual());
	  gmt_tf_idx.push_back(l1MuIt->tfMuonIndex());
	  gmt_charge.push_back(l1MuIt->hwChargeValid() ? l1MuIt->charge() : 0);
	  igmt++;
	}
    }
  
  return;
}

void TTreeGenerator::fill_gt_variables(edm::Handle<L1GlobalTriggerReadoutRecord> gtrr, const L1GtTriggerMenu* menu) // CB FIXME speedup // legacy
{
  igtalgo = 0;
  igttt   = 0;
  const AlgorithmMap algoMap = menu->gtAlgorithmMap();
  const AlgorithmMap ttMap   = menu->gtTechnicalTriggerMap();
  for(int ibx=0; ibx<5; ibx++) {
    DecisionWord gtDecisionWord = gtrr->decisionWord(ibx-2);
    if (!gtDecisionWord.empty()) {
      CItAlgo algoMapIt  = algoMap.begin();
      CItAlgo algoMapEnd = algoMap.end();
      for(; algoMapIt!=algoMapEnd; ++algoMapIt) {
	if( menu->gtAlgorithmResult((*algoMapIt).first, gtDecisionWord) ) {
	  // gt_algo_bit.push_back(TString((*algoMapIt).first));
	  gt_algo_bit.push_back((algoMapIt->second).algoBitNumber());
	  gt_algo_bx.push_back(ibx-2);
	  igtalgo++;
	}
      }
    }
    TechnicalTriggerWord gtTTWord = gtrr->technicalTriggerWord(ibx-2);
    if (!gtTTWord.empty()) {
      CItAlgo ttMapIt  = ttMap.begin();
      CItAlgo ttMapEnd = ttMap.end();
      for(; ttMapIt!=ttMapEnd; ++ttMapIt) {
	int bitNumber = (ttMapIt->second).algoBitNumber();
	if (gtTTWord.at(bitNumber)) {
	  //	  gt_tt_bit.push_back(TString(ttMapIt->first));
	  gt_tt_bit.push_back((ttMapIt->second).algoBitNumber());
	  gt_tt_bx.push_back(ibx-2);
	  igttt++;
	}
      }
    }
  }
  return;
}

void TTreeGenerator::fill_hlt_variables(const edm::Event &e, edm::Handle<edm::TriggerResults> hltresults)
{
  ihlt = 0; 
  const edm::TriggerNames TrigNames_ = e.triggerNames(*hltresults);
  const int ntrigs = hltresults->size();

  for (int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if (hltresults->accept(itr)) {
      hlt_path.push_back(trigName);
      ihlt++;
    }
  }
  return;
}

void TTreeGenerator::fill_rpc_variables(const edm::Event &e, edm::Handle<RPCRecHitCollection> rpcrechits){
  RPCRecHitCollection::const_iterator recHit;
  irpcrechits=0;
  for(recHit = rpcrechits->begin(); recHit != rpcrechits->end(); recHit++){ 
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
    rpc_region.push_back(region);
    rpc_clusterSize.push_back(cls);
    rpc_strip.push_back(firststrip);
    rpc_bx.push_back(bx);
    rpc_station.push_back(stat);
    rpc_sector.push_back(sect);
    rpc_layer.push_back(layer);
    rpc_subsector.push_back(subsector);
    rpc_roll.push_back(roll);
    rpc_ring.push_back(ring);
    irpcrechits++;
    
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
    


    const RPCRoll * rollasociated = rpcGeom_->roll(rpcId);
    const BoundPlane & RPCSurface = rollasociated->surface(); 
    GlobalPoint RPCRecHitInGlobal = RPCSurface.toGlobal(recHitPos);
    float glob_x = RPCRecHitInGlobal.x();
    float glob_y = RPCRecHitInGlobal.y();
    float glob_z = RPCRecHitInGlobal.z();
    float glob_eta = RPCRecHitInGlobal.eta();
    float glob_phi = RPCRecHitInGlobal.phi();
    RpcRechit_TwinMux_Glob_x.push_back(glob_x);
    RpcRechit_TwinMux_Glob_y.push_back(glob_y);
    RpcRechit_TwinMux_Glob_z.push_back(glob_z);
	RpcRechit_TwinMux_Glob_eta.push_back(glob_eta);
	RpcRechit_TwinMux_Glob_phi.push_back(glob_phi);
    
    
    irpcrechits_TwinMux++;
  }
  return;
}

void TTreeGenerator::analyzeRPCunpacking(const edm::Event& event)
{
	edm::Handle<MuonDigiCollection<RPCDetId,RPCDigi> > rpc;
	event.getByToken(rpcToken_, rpc);
	auto chamber = rpc->begin();
	auto chend  = rpc->end();
	irpcdigi_TwinMux=0;
	for( ; chamber != chend; ++chamber ) {
	
		if(OnlyBarrel_ && (*chamber).first.region() != 0) continue;

   		auto digi = (*chamber).second.first;
	    auto dend = (*chamber).second.second;
		for( ; digi != dend; ++digi ) {
			RpcDigi_TwinMux_bx.push_back(digi->bx());
			RpcDigi_TwinMux_strip.push_back(digi->strip());
		}
		RpcDigi_TwinMux_region.push_back((*chamber).first.region()); 
		RpcDigi_TwinMux_ring.push_back((*chamber).first.ring());
		RpcDigi_TwinMux_station.push_back((*chamber).first.station());
		RpcDigi_TwinMux_layer.push_back((*chamber).first.layer());
		RpcDigi_TwinMux_sector.push_back((*chamber).first.sector());
		RpcDigi_TwinMux_subsector.push_back((*chamber).first.subsector()); 
		RpcDigi_TwinMux_roll.push_back((*chamber).first.roll());
		RpcDigi_TwinMux_trIndex.push_back((*chamber).first.trIndex());
		RpcDigi_TwinMux_det.push_back((*chamber).first.det());
		RpcDigi_TwinMux_subdetId.push_back((*chamber).first.subdetId()); 
		RpcDigi_TwinMux_rawId.push_back((*chamber).first.rawId()); 
		irpcdigi_TwinMux++;
	}  

}

void TTreeGenerator::extrapolate_DTsegment_onRPC(edm::Handle<DTRecSegment4DCollection> segments4D, const DTGeometry* dtGeom_, const RPCGeometry* rpcGeom_, const edm::EventSetup& iSetup)
{
// 	DT_segment_onRPC = 0;
// 	
// 	if(segments4D->size()>10){ //Too many segments in this event we are not doing the extrapolation
// 		std::cout<<"-----------------------------Too many segments in this event we are not doing the extrapolation"<<std::endl;
// 	}
// 	else{
// 		    std::map<DTChamberId,int> DTSegmentCounter;
// 		    DTRecSegment4DCollection::const_iterator segment;  
//   
//     		for (segment = segments4D->begin();segment!=segments4D->end(); ++segment){
// 			      DTSegmentCounter[segment->chamberId()]++;
// 		    }
// 		    
// 		    for (segment = segments4D->begin(); segment != segments4D->end(); ++segment){
// 		    
// 		    	DTChamberId DTId = segment->chamberId();
//       
// // 				if(debug) std::cout<<"DT  \t \t This Segment is in Chamber id: "<<DTId<<std::endl;
// // 				if(debug) std::cout<<"DT  \t \t Number of segments in this DT = "<<DTSegmentCounter[DTId]<<std::endl;
// // 				if(debug) std::cout<<"DT  \t \t Is the only one in this DT? and is not in the 4th Station?"<<std::endl;
// 
//       			if(DTSegmentCounter[DTId]!=1 || DTId.station()==4){
// // 					std::cout<<"DT \t \t More than one segment in this chamber, or we are in Station 4"<<std::endl;
// 				    continue;
// 				}
// 				
// 				int dtWheel = DTId.wheel();
// 				int dtStation = DTId.station();
// 				int dtSector = DTId.sector();      
// 				
// 				LocalPoint segmentPosition= segment->localPosition();
// 				LocalVector segmentDirection=segment->localDirection();
// 				
// 				const GeomDet* gdet=dtGeom_->idToDet(segment->geographicalId());
// 				const BoundPlane & DTSurface = gdet->surface();	 
// 				
// 				if(segment->dimension()!=4) continue;  //check if the dimension of the segment is 4
// 				if((segment->phiSegment()->recHits()).size()<7) continue; //check if there are 2 phi superlayers fired (at least 7hits)
// 				if((segment->zSegment()->recHits()).size()<4) continue; //check if there are at least 4 hits for the Z projection
// 				
// // 				if(debug) std::cout<<"DT  \t \t yes"<<std::endl;
// // 				if(debug) std::cout<<"DT  \t \t DT Segment Dimension "<<segment->dimension()<<std::endl; 
//     
// 				float Xo=segmentPosition.x();
// 				float Yo=segmentPosition.y();
// 				float Zo=segmentPosition.z();
// 				float dx=segmentDirection.x();
// 				float dy=segmentDirection.y();
// 				float dz=segmentDirection.z();
//      
// // 				if(debug)  std::cout<<"Calling to Object Map class"<<std::endl;
// 				ObjectMap* TheObject = ObjectMap::GetInstance(iSetup);
// // 				if(debug) std::cout<<"Creating the DTIndex"<<std::endl;
// 				DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
// // 				if(debug) std::cout<<"Getting the Rolls for the given index"<<std::endl;
// 				std::set<RPCDetId> rollsForThisDT = TheObject->GetInstance(iSetup)->GetRolls(theindex);
// // 				if(debug) std::cout<<"DT  \t \t Number of rolls for this DT = "<<rollsForThisDT.size()<<std::endl;
// 				assert(rollsForThisDT.size()>=1);				
// 
// // 				if(debug) std::cout<<"DT  \t \t Loop over all the rolls asociated to this DT"<<std::endl;
// 				for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisDT.begin();iteraRoll != rollsForThisDT.end(); iteraRoll++){
// 					const RPCRoll* rollasociated = rpcGeom_->roll(*iteraRoll);
// 					RPCDetId rpcId = rollasociated->id();
// 					const BoundPlane & RPCSurface = rollasociated->surface(); 
// 	
// 					RPCGeomServ rpcsrv(rpcId);
// 					std::string nameRoll = rpcsrv.name();
// 	
// // 					if(debug) std::cout<<"DT  \t \t \t RollName: "<<nameRoll<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t Doing the extrapolation to this roll"<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t DT Segment Direction in DTLocal "<<segmentDirection<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t DT Segment Point in DTLocal "<<segmentPosition<<std::endl;
// // 	
// 					GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
// 
// 					LocalPoint CenterRollinDTFrame = DTSurface.toLocal(CenterPointRollGlobal);
// 	
// // 					if(debug) std::cout<<"DT  \t \t \t Center (0,0,0) Roll In DTLocal"<<CenterRollinDTFrame<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t Center (0,0,0) of the Roll in Global"<<CenterPointRollGlobal<<std::endl;
// // 	
// 					float D=CenterRollinDTFrame.z();
// 	
// 					float X=Xo+dx*D/dz;
// 					float Y=Yo+dy*D/dz;
// 					float Z=D;
// 	
// 					const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
// // 					std::cout<<"----------------------------------------------------------------------------------------------------stripw = "<<top_->pitch()<<std::endl;
// 					LocalPoint xmin = top_->localPosition(0.);
// // 					if(debug) std::cout<<"DT  \t \t \t xmin of this  Roll "<<xmin<<"cm"<<std::endl;
// 					LocalPoint xmax = top_->localPosition((float)rollasociated->nstrips());
// // 					if(debug) std::cout<<"DT  \t \t \t xmax of this  Roll "<<xmax<<"cm"<<std::endl;
// 					float rsize = fabs( xmax.x()-xmin.x() );
// // 					if(debug) std::cout<<"DT  \t \t \t Roll Size "<<rsize<<"cm"<<std::endl;
// 					float stripl = top_->stripLength();
// 	
// 					float stripw = top_->pitch();
// 	
// // 					if(debug) std::cout<<"DT  \t \t \t Strip Lenght "<<stripl<<"cm"<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t Strip Width "<<stripw<<"cm"<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t X Predicted in DTLocal= "<<X<<"cm"<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t Y Predicted in DTLocal= "<<Y<<"cm"<<std::endl;
// // 					if(debug) std::cout<<"DT  \t \t \t Z Predicted in DTLocal= "<<Z<<"cm"<<std::endl;
// 	
// 					float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));
// 	
// // 					if(debug) std::cout<<"DT  \t \t \t Is the distance of extrapolation less than MaxD? ="<<extrapolatedDistance<<"cm"<<"MaxD="<<MaxD<<"cm"<<std::endl;
// 	
// 					if(extrapolatedDistance<=MaxD){ 
// // 					  if(debug) std::cout<<"DT  \t \t \t yes"<<std::endl;   
// 					  GlobalPoint GlobalPointExtrapolated = DTSurface.toGlobal(LocalPoint(X,Y,Z));
// // 					  if(debug) std::cout<<"DT  \t \t \t Point ExtraPolated in Global"<<GlobalPointExtrapolated<< std::endl;
// 					  LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);
// 					
// // 					  if(debug) std::cout<<"DT  \t \t \t Point Extrapolated in RPCLocal"<<PointExtrapolatedRPCFrame<< std::endl;
// // 					  if(debug) std::cout<<"DT  \t \t \t Corner of the Roll = ("<<rsize*eyr<<","<<stripl*eyr<<")"<<std::endl;
// // 					  if(debug) std::cout<<"DT \t \t \t Info About the Point Extrapolated in X Abs ("<<fabs(PointExtrapolatedRPCFrame.x())<<","
// // 							     <<fabs(PointExtrapolatedRPCFrame.y())<<","<<fabs(PointExtrapolatedRPCFrame.z())<<")"<<std::endl;
// // 					  if(debug) std::cout<<"DT  \t \t \t Does the extrapolation go inside this roll?"<<std::endl;
// 					
// 					  if(fabs(PointExtrapolatedRPCFrame.z()) < 1. && 
// 					     fabs(PointExtrapolatedRPCFrame.x()) < rsize*eyr && 
// 					     fabs(PointExtrapolatedRPCFrame.y()) < stripl*eyr){
// // 					    if(debug) std::cout<<"DT  \t \t \t \t yes"<<std::endl;	
// // 					    if(debug) std::cout<<"DT  \t \t \t \t Creating the RecHit"<<std::endl;
// 
// 					    LocalVector segmentDirection=segment->localDirection();
// 					    float dx=segmentDirection.x();
// 					    float dz=segmentDirection.z();
// 					    float cosal = dx/sqrt(dx*dx+dz*dz);
// 					    float angle = acos(cosal)*180/3.1415926;
// 					    					    
// 					    RPCRecHit RPCPoint(rpcId,0,LocalPoint(PointExtrapolatedRPCFrame.x(),PointExtrapolatedRPCFrame.y(),angle));
// 
// 						float dt_loc_x = PointExtrapolatedRPCFrame.x();
// 						float dt_loc_y = PointExtrapolatedRPCFrame.y();
// 						float dt_loc_z = PointExtrapolatedRPCFrame.z();
// 						float dt_glob_x = GlobalPointExtrapolated.x();
// 						float dt_glob_y = GlobalPointExtrapolated.y();
// 						float dt_glob_z = GlobalPointExtrapolated.z();
// 						float dt_glob_eta = GlobalPointExtrapolated.eta();
// 						float dt_glob_phi = GlobalPointExtrapolated.phi();
// 						
// 						DT_extrapolated_OnRPC_BX.push_back(RPCPoint.BunchX());
// 						
// 						DT_extrapolated_OnRPC_Loc_x.push_back(dt_loc_x);
// 						DT_extrapolated_OnRPC_Loc_y.push_back(dt_loc_y);
// 						DT_extrapolated_OnRPC_Loc_z.push_back(dt_loc_z);
// 						
// 						DT_extrapolated_OnRPC_Glob_x.push_back(dt_glob_x);
// 						DT_extrapolated_OnRPC_Glob_y.push_back(dt_glob_y);
// 						DT_extrapolated_OnRPC_Glob_z.push_back(dt_glob_z);
// 						DT_extrapolated_OnRPC_Glob_eta.push_back(dt_glob_eta);
// 						DT_extrapolated_OnRPC_Glob_phi.push_back(dt_glob_phi);
// 						
// 						DT_extrapolated_OnRPC_Region.push_back(rpcId.region());
// 						DT_extrapolated_OnRPC_Sector.push_back(rpcId.sector());
// 						DT_extrapolated_OnRPC_Station.push_back(rpcId.station());
// 						DT_extrapolated_OnRPC_Layer.push_back(rpcId.layer());
// 						DT_extrapolated_OnRPC_Roll.push_back(rpcId.roll());
// 						DT_extrapolated_OnRPC_Ring.push_back(rpcId.ring());
// 						DT_extrapolated_OnRPC_Stripw.push_back(stripw);
// 												
// 						DT_segment_onRPC++;
// 						
// // 					    if(debug) std::cout<<"DT  \t \t \t \t Clearing the vector"<<std::endl; 
// // 					    RPCPointVector.clear();
// // 					    if(debug) std::cout<<"DT  \t \t \t \t Pushing back"<<std::endl;	
// // 					    RPCPointVector.push_back(RPCPoint); 
// // 					    if(debug) std::cout<<"DT  \t \t \t \t Putting the vector"<<std::endl;	
// // 					    _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());
// 
// // 					    if(debug) std::cout<<"DT \t \t \t \t Filling container with "<<nameRoll
// // 				 			<<" Point.x="<<PointExtrapolatedRPCFrame.x()<<" Point.y="<<PointExtrapolatedRPCFrame.y()<<" size="<<RPCPointVector.size()<<std::endl;
// 					  
// 					  }else {
// /*					    if(debug)*/ std::cout<<"DT \t \t \t \t No the prediction is outside of this roll"<<std::endl;
// 					  }//Condition for the right match
// 					}else{
// /*					  if(debug) */ std::cout<<"DT \t \t \t No, Exrtrapolation too long!, canceled"<<std::endl;
// 					}//D so big
// 
// 
// 				}//loop over all the rolls asociated
// 
// 		    
// 		    }
//       


//       }
//     }
//   
//     if(incldtMB4){
//       if(all4DSegments->size()>0){     
// 	if(debug) std::cout<<"MB4 \t \t Loop Over all4DSegments "<<all4DSegments->size()<<std::endl;
// 	extrapolatedRolls.clear();
// 	for (segment = all4DSegments->begin(); segment != all4DSegments->end(); ++segment){ 
//     
// 	  DTChamberId DTId = segment->chamberId();
// 	
// 	  if(debug) std::cout<<"MB4 \t \t This Segment is in Chamber id: "<<DTId<<std::endl;
// 	  if(debug) std::cout<<"MB4 \t \t Number of segments in this DT = "<<DTSegmentCounter[DTId]<<std::endl;
// 	  if(debug) std::cout<<"MB4 \t \t \t Is the only one in this DT? and is in the Station 4?"<<std::endl;
// 
// 	  if(DTSegmentCounter[DTId] == 1 && DTId.station()==4){
// 
// 	    if(debug) std::cout<<"MB4 \t \t \t yes"<<std::endl;
// 	    int dtWheel = DTId.wheel();
// 	    int dtStation = DTId.station();
// 	    int dtSector = DTId.sector();
//       
// 	    LocalPoint segmentPosition= segment->localPosition();
// 	    LocalVector segmentDirection=segment->localDirection();
// 	    
// 	    if(debug) std::cout<<"MB4 \t \t \t \t The Segment in MB4 is 2D?"<<std::endl;
// 	    if(segment->dimension()==2){
// 	      if(debug) std::cout<<"MB4 \t \t \t \t yes"<<std::endl;
// 	      LocalVector segmentDirectionMB4=segmentDirection;
// 	      LocalPoint segmentPositionMB4=segmentPosition;
// 	    
// 	      const BoundPlane& DTSurface4 = dtGeo->idToDet(DTId)->surface();
// 	        
// 	      DTRecSegment4DCollection::const_iterator segMB3;  
// 	        
// 	      if(debug) std::cout<<"MB4 \t \t \t \t Loop on segments in =sector && MB3 && adjacent sectors && y dim=4"<<std::endl;
// 	      for(segMB3=all4DSegments->begin();segMB3!=all4DSegments->end();++segMB3){
// 	            
// 		DTChamberId dtid3 = segMB3->chamberId();  
// 		
// 		if(debug) std::cout<<"MB4  \t \t \t \t Segment in Chamber ="<<dtid3<<std::endl;
// 	            
// 		if(distsector(dtid3.sector(),DTId.sector())<=1 //The DT sector could be 13 or 14 and because is corrected in the calculation of the distance.
// 		   && distwheel(dtid3.wheel(),DTId.wheel())<=1 //The we could have segments in neighbohr wheels in pp collisions 
// 		   && dtid3.station()==3
// 		   && DTSegmentCounter[dtid3] == 1
// 		   && segMB3->dimension()==4){
// 		  
// 		  if(debug) std::cout<<"MB4  \t \t \t \t distsector ="<<distsector(dtid3.sector(),DTId.sector())<<std::endl;
// 		  if(debug) std::cout<<"MB4  \t \t \t \t distwheel ="<<distwheel(dtid3.wheel(),DTId.wheel())<<std::endl;
// 
// 		  const GeomDet* gdet3=dtGeo->idToDet(segMB3->geographicalId());
// 		  const BoundPlane & DTSurface3 = gdet3->surface();
// 
// 		  LocalVector segmentDirectionMB3 =  segMB3->localDirection();
// 		  GlobalPoint segmentPositionMB3inGlobal = DTSurface3.toGlobal(segMB3->localPosition());
// 		  GlobalPoint segmentPositionMB4inGlobal = DTSurface4.toGlobal(segmentPosition);
// 		
// 		  LocalVector segDirMB3inMB4Frame=DTSurface4.toLocal(DTSurface3.toGlobal(segmentDirectionMB3));
// 		
// 		  GlobalVector segDirMB4inGlobalFrame=DTSurface4.toGlobal(segmentDirectionMB4);
// 		  GlobalVector segDirMB3inGlobalFrame=DTSurface3.toGlobal(segmentDirectionMB3);
// 		
// 		  float dx=segDirMB4inGlobalFrame.x();
// 		  float dy=segDirMB4inGlobalFrame.y();
// 				
// 		  float dx3=segDirMB3inGlobalFrame.x();
// 		  float dy3=segDirMB3inGlobalFrame.y();
// 				
// 		  double cosAng=fabs(dx*dx3+dy*dy3/sqrt((dx3*dx3+dy3*dy3)*(dx*dx+dy*dy)));
// 		  
// 		  if(debug) std::cout<<"MB4 \t \t \t \t cosAng"<<cosAng<<"Beetween "<<dtid3<<" and "<<DTId<<std::endl;
// 		  
// 		  if(debug){
// 		    std::cout<<"MB4 \t \t \t \t dx="<<dx<<" dy="<<dy<<std::endl;
// 		    std::cout<<"MB4 \t \t \t \t dx3="<<dx3<<" dy3="<<dy<<std::endl;
// 		    std::cout<<"MB4 \t \t \t \t cosAng="<<cosAng<<std::endl;
// 		  }
// 
// 		  float DistanceBetweenSegments = ((segmentPositionMB3inGlobal) - (segmentPositionMB4inGlobal)).mag();
// 		  
// 		  if(cosAng>MinCosAng && DistanceBetweenSegments < MaxDistanceBetweenSegments){
// 		    
// 		    if(debug) std::cout<<"MB4 \t \t \t \t Distance between segments="<<DistanceBetweenSegments<<std::endl;
// 		  
// 		    if(debug) std::cout<<"MB4 \t \t We found compatible Segments (similar direction and close enough) in "<<dtid3<<" and "<<DTId<<std::endl;
// 
// 		    if(dtSector==13){
// 		      dtSector=4;
// 		    }
// 		    if(dtSector==14){
// 		      dtSector=10;
// 		    }
// 		   
// 		    if(debug)  std::cout<<"Calling to Object Map class"<<std::endl;
// 		    ObjectMap* TheObject = ObjectMap::GetInstance(iSetup);
// 		    if(debug) std::cout<<"Creating the DTIndex"<<std::endl;
// 		    DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
// 		    if(debug) std::cout<<"Getting the Rolls for the given index"<<std::endl;
// 
// 		    std::set<RPCDetId> rollsForThisDT = TheObject->GetInstance(iSetup)->GetRolls(theindex);
// 
// 		    if(debug) std::cout<<"MB4 \t \t Number of rolls for this DT = "<<rollsForThisDT.size()<<std::endl;
// 		    
// 		    assert(rollsForThisDT.size()>=1);
// 		        
// 		    if(debug) std::cout<<"MB4  \t \t Loop over all the rolls asociated to this DT"<<std::endl;
// 		    for (std::set<RPCDetId>::iterator iteraRoll=rollsForThisDT.begin();iteraRoll != rollsForThisDT.end(); iteraRoll++){
// 		      const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll); //roll asociated to MB4
// 		      RPCDetId rpcId = rollasociated->id();
// 		      const BoundPlane & RPCSurfaceRB4 = rollasociated->surface(); //surface MB4
// 
// 		      RPCGeomServ rpcsrv(rpcId);
// 		      std::string nameRoll = rpcsrv.name();
// 
// 		      if(debug) std::cout<<"MB4  \t \t \t RollName: "<<nameRoll<<std::endl;
// 		      if(debug) std::cout<<"MB4  \t \t \t Doing the extrapolation to this roll"<<std::endl;
// 		        
// 		      GlobalPoint CenterPointRollGlobal=RPCSurfaceRB4.toGlobal(LocalPoint(0,0,0));
// 		      LocalPoint CenterRollinMB4Frame = DTSurface4.toLocal(CenterPointRollGlobal); //In MB4
// 		      LocalPoint segmentPositionMB3inMB4Frame = DTSurface4.toLocal(segmentPositionMB3inGlobal); //In MB4
// 		      //LocalPoint segmentPositionMB3inRB4Frame = RPCSurfaceRB4.toLocal(segmentPositionMB3inGlobal); //In MB4
// 		      LocalVector segmentDirectionMB3inMB4Frame = DTSurface4.toLocal(segDirMB3inGlobalFrame); //In MB4
// 		        
// 		      //The exptrapolation is done in MB4 frame. for local x and z is done from MB4,
// 		      float Dxz=CenterRollinMB4Frame.z();
// 		      float Xo4=segmentPositionMB4.x();
// 		      float dxl=segmentDirectionMB4.x(); //dx local for MB4 segment in MB4 Frame
// 		      float dzl=segmentDirectionMB4.z(); //dx local for MB4 segment in MB4 Frame
// 		        
// 		      float X=Xo4+dxl*Dxz/dzl; //In MB4 frame
// 		      float Z=Dxz;//In MB4 frame
// 		        
// 		      //for local y is done from MB3 
// 		      float Yo34 = segmentPositionMB3inMB4Frame.y();
// 		      float dy34 = segmentDirectionMB3inMB4Frame.y();
// 		      float dz34 = segmentDirectionMB3inMB4Frame.z();
// 		      float Dy=Dxz-(segmentPositionMB3inMB4Frame.z()); //Distance beetween the segment in MB3 and the RB4 surface
// 
// 		      if(debug) std::cout<<"MB4 \t \t \t The distance to extrapolate in Y from MB3 is "<<Dy<<"cm"<<std::endl;
// 		        
// 		      float Y=Yo34+dy34*Dy/dz34;//In MB4 Frame
// 		          
// 		      const RectangularStripTopology* top_
// 			=dynamic_cast<const RectangularStripTopology*>(&(rollasociated->topology())); //Topology roll asociated MB4
// 		      LocalPoint xmin = top_->localPosition(0.);
// 		      LocalPoint xmax = top_->localPosition((float)rollasociated->nstrips());
// 		      float rsize = fabs( xmax.x()-xmin.x() );
// 		      float stripl = top_->stripLength();
// 		      float stripw = top_->pitch();
// 
// 			    
// 		      if(debug) std::cout<<"MB4 \t \t \t Strip Lenght "<<stripl<<"cm"<<std::endl;
// 		      if(debug) std::cout<<"MB4 \t \t \t Strip Width "<<stripw<<"cm"<<std::endl;
// 
// 		      if(debug) std::cout<<"MB4 \t \t \t X Predicted in MB4DTLocal= "<<X<<"cm"<<std::endl;
// 		      if(debug) std::cout<<"MB4 \t \t \t Y Predicted in MB4DTLocal= "<<Y<<"cm"<<std::endl;
// 		      if(debug) std::cout<<"MB4 \t \t \t Z Predicted in MB4DTLocal= "<<Z<<"cm"<<std::endl;
// 
// 		      float extrapolatedDistance = sqrt((Y-Yo34)*(Y-Yo34)+Dy*Dy);
// 
// 		      if(debug) std::cout<<"MB4 \t \t \t segmentPositionMB3inMB4Frame"<<segmentPositionMB3inMB4Frame<<std::endl;
// 		      if(debug) std::cout<<"MB4 \t \t \t segmentPositionMB4inMB4Frame"<<segmentPosition<<std::endl;
// 
// 		      if(debug) std::cout<<"MB4 \t \t \t segmentDirMB3inMB4Frame"<<segDirMB3inMB4Frame<<std::endl;
// 		      if(debug) std::cout<<"MB4 \t \t \t segmentDirMB4inMB4Frame"<<segmentDirectionMB4<<std::endl;
// 			    
// 		      if(debug) std::cout<<"MB4 \t \t \t CenterRB4PositioninMB4Frame"<<CenterRollinMB4Frame<<std::endl;
// 			    
// 		      if(debug) std::cout<<"MB4 \t \t \t Is the extrapolation distance ="<<extrapolatedDistance<<"less than "<<MaxDrb4<<std::endl;
//     
// 
// 		      if(extrapolatedDistance<=MaxDrb4){ 
// 			if(debug) std::cout<<"MB4 \t \t \t yes"<<std::endl;
// 
// 			GlobalPoint GlobalPointExtrapolated = DTSurface4.toGlobal(LocalPoint(X,Y,Z));
// 			        
// 			if(debug) std::cout<<"MB4 \t \t \t Point ExtraPolated in Global"<<GlobalPointExtrapolated<< std::endl;
// 			        
// 			LocalPoint PointExtrapolatedRPCFrame = RPCSurfaceRB4.toLocal(GlobalPointExtrapolated);
// 
// 			if(debug) std::cout<<"MB4 \t \t \t Point Extrapolated in RPCLocal"<<PointExtrapolatedRPCFrame<< std::endl;
// 			if(debug) std::cout<<"MB4 \t \t \t Corner of the Roll = ("<<rsize*eyr<<","<<stripl*eyr<<")"<<std::endl;
// 			if(debug) std::cout<<"MB4 \t \t \t Info About the Point Extrapolated in X Abs ("<<fabs(PointExtrapolatedRPCFrame.x())<<","
// 					   <<fabs(PointExtrapolatedRPCFrame.y())<<","<<fabs(PointExtrapolatedRPCFrame.z())<<")"<<std::endl;
// 			  
// 			if(debug) std::cout<<"MB4 \t \t \t Does the extrapolation go inside this roll?"<<std::endl;
// 			  
// 			if(fabs(PointExtrapolatedRPCFrame.z()) < 5.  &&
// 			   fabs(PointExtrapolatedRPCFrame.x()) < rsize*eyr &&
// 			   fabs(PointExtrapolatedRPCFrame.y()) < stripl*eyr){
// 			  if(debug) std::cout<<"MB4 \t \t \t \t yes"<<std::endl;
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Creating the RecHit"<<std::endl;
// 
// 			  LocalVector segmentDirection=segment->localDirection();
// 			  float dx=segmentDirection.x();
// 			  float dz=segmentDirection.z();
// 			  float cosal = dx/sqrt(dx*dx+dz*dz);
// 			  float angle = acos(cosal)*180/3.1415926;
// 
// 			  RPCRecHit RPCPointMB4(rpcId,0,LocalPoint(PointExtrapolatedRPCFrame.x(),PointExtrapolatedRPCFrame.y(),angle));
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Clearing the RPCPointVector"<<std::endl;
// 			  RPCPointVector.clear();
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Pushing Back"<<std::endl;
// 			  RPCPointVector.push_back(RPCPointMB4);
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Putting for "<<rpcId<<std::endl;
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Filling container with "<<nameRoll
// 					     <<" Point.x="<<PointExtrapolatedRPCFrame.x()<<" Point.y="<<PointExtrapolatedRPCFrame.y()<<" size="<<RPCPointVector.size()<<std::endl;
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Number of rolls already extrapolated in RB4 = "<<extrapolatedRolls.size()<<std::endl;
// 			  if(find (extrapolatedRolls.begin(),extrapolatedRolls.end(),rpcId.rawId()) == extrapolatedRolls.end()){
// 			    extrapolatedRolls.push_back(rpcId.rawId());
// 			    _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());
// 			  }else{
// 			    if(debug) std::cout<<"MB4 \t \t \t \t roll already extrapolated "<<rpcId<<std::endl;
// 			  }
// 			  if(debug) std::cout<<"MB4 \t \t \t \t Extrapolations done after this point = "<<extrapolatedRolls.size()<<std::endl;
// 			  if(debug) for(uint32_t m=0;m<extrapolatedRolls.size();m++) std::cout<<"MB4 \t \t \t \t"<< extrapolatedRolls.at(m)<<std::endl;
// 			}else{
// 			  if(debug) std::cout<<"MB4 \t \t \t \t No the prediction is outside of this roll"<<std::endl;
// 			}
// 		      }//Condition for the right match
// 		      else{
// 			if(debug) std::cout<<"MB4 \t \t \t No, Exrtrapolation too long!, canceled"<<std::endl;
// 		      }
// 		    }//loop over all the rollsasociated
// 		  }else{
// 		    if(debug) std::cout<<"MB4 \t \t \t \t I found segments in MB4 and MB3 adjacent wheel and/or sector but not compatibles, Diferent Directions"<<std::endl;
// 		  }
// 		}else{//if dtid3.station()==3&&dtid3.sector()==DTId.sector()&&dtid3.wheel()==DTId.wheel()&&segMB3->dim()==4
// 		  if(debug) std::cout<<"MB4 \t \t \t No the same station or same wheel or segment dim in mb3 not 4D"<<std::endl;
// 		}
// 	      }//loop over all the segments looking for one in MB3 
// 	    }else{
// 	      if(debug) std::cout<<"MB4 \t \t \t Is NOT a 2D Segment"<<std::endl;
// 	    }
// 	  }else{
// 	    if(debug) std::cout<<"MB4 \t \t \t \t There is not just one segment or is not in station 4"<<std::endl;
// 	  }//De aca para abajo esta en dtpart.inl
// 	}
//       }else{
// 	if(debug) std::cout<<"MB4 \t This event doesn't have 4D Segment"<<std::endl;
//       }
//     }
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
// 	}


	
	
// 	
// 	for (DTRecSegment4DCollection::const_iterator segment = all4DSegments->begin(); segment != all4DSegments->end(); ++segment){
// 		std::cout<<"----------------------------------------------------------------------ciao"<<std::endl;
// 	}

}

void TTreeGenerator::analyzeBMTF(const edm::Event& event)
{
///Output of BMTF
  int ctr = 0;
  edm::Handle<l1t::RegionalMuonCandBxCollection> mycoll;
  event.getByToken(bmtfOutputTag_,mycoll);

  if (mycoll.isValid()) {
        int firstbx = (*mycoll).getFirstBX();
        int lastbx  = (*mycoll).getLastBX() + 1;
        for(int i=firstbx; i<lastbx; i++){
  				for (auto mu = (*mycoll).begin(i); mu != (*mycoll).end(i); ++mu) {
		    	  ctr++;
			      Bmtf_Pt.push_back(mu->hwPt()*0.5);
			      Bmtf_Eta.push_back(mu->hwEta()*0.010875);
			      Bmtf_FineBit.push_back(mu->hwHF());
			      Bmtf_Phi.push_back(mu->hwPhi());
			      //2*mu->hwPhi()*TMath::Pi()/576;
			      int phi = 0;
			      phi = mu->processor()*48 + mu->hwPhi();
			      phi += 576 - 24;
			      phi = phi % 576;
			      Bmtf_GlobalPhi.push_back(phi);
		    	  Bmtf_qual.push_back(mu->hwQual());
	    		  Bmtf_ch.push_back(mu->hwSign());
    			  Bmtf_bx.push_back(i);
			      Bmtf_processor.push_back(mu->processor());
			      std::map<int, int>  trAdd;
	    		  trAdd = mu->trackAddress();
		    	  int wheel = pow(-1,trAdd[0]) * trAdd[1];
    			  Bmtf_wh.push_back(wheel);
			      Bmtf_trAddress.push_back(trAdd[2]);
    			  Bmtf_trAddress.push_back(trAdd[3]);
	    		  Bmtf_trAddress.push_back(trAdd[4]);
			      Bmtf_trAddress.push_back(trAdd[5]);
    			} // for mu
		        Bmtf_Size = ctr;
 	 } // for i
  } //if
  else 
      edm::LogInfo("L1Prompt") << "can't find L1MuMBTrackContainer";
      

  edm::Handle<L1MuDTChambPhContainer> bmtfPhInputs;  
  event.getByToken(bmtfPhInputTag_, bmtfPhInputs);
  
  if(!(bmtfPhInputs.isValid())) std::cout<<"no  ok"<<std::endl;
  
   L1MuDTChambPhContainer::Phi_Container const *PhContainer = bmtfPhInputs->getContainer();

   Bmtf_phSize = PhContainer->size();
   int iphtr=0;
   for( L1MuDTChambPhContainer::Phi_Container::const_iterator DTPhDigiItr =  PhContainer->begin() ;
       DTPhDigiItr != PhContainer->end(); ++DTPhDigiItr ){
		      Bmtf_phBx.push_back     (  DTPhDigiItr->bxNum() );
    		  Bmtf_phTs2Tag.push_back     ( DTPhDigiItr->Ts2Tag() );
		      Bmtf_phWh.push_back     (  DTPhDigiItr->whNum() );
	    	  Bmtf_phSe.push_back     (  DTPhDigiItr->scNum() );
	    	  Bmtf_phSt.push_back     (  DTPhDigiItr->stNum() );
    		  Bmtf_phAng.push_back    (  DTPhDigiItr->phi()   );
		      Bmtf_phBandAng.push_back(  DTPhDigiItr->phiB()  );
    		  Bmtf_phCode.push_back   (  DTPhDigiItr->code()  );
		      iphtr++;
    }



  edm::Handle<L1MuDTChambThContainer > bmtfThInputs;
  event.getByToken(bmtfThInputTag_, bmtfThInputs); 
  
   L1MuDTChambThContainer::The_Container const *ThContainer = bmtfThInputs->getContainer();
   
   if(!(bmtfThInputs.isValid())) std::cout<<"no  ok"<<std::endl;

   int ithtr=0;
   Bmtf_thSize = ThContainer->size();

   for( L1MuDTChambThContainer::The_Container::const_iterator
	 DTThDigiItr =  ThContainer->begin() ;
       DTThDigiItr != ThContainer->end() ;
       ++DTThDigiItr )
     {

      Bmtf_thBx.push_back( DTThDigiItr->bxNum()  );
      Bmtf_thWh.push_back( DTThDigiItr->whNum() );
      Bmtf_thSe.push_back( DTThDigiItr->scNum() );
      Bmtf_thSt.push_back( DTThDigiItr->stNum() );

      ostringstream  ss1, ss2; 
      ss1.clear(); ss2.clear();
      ss1<<"9"; ss2<<"9";

      for(int j=0; j<7; j++){
        ss1<<DTThDigiItr->position(j);
        ss2<<DTThDigiItr->code(j) ;
      }
      Bmtf_thTheta.push_back(stoi(ss1.str())) ;
      Bmtf_thCode.push_back(stoi(ss2.str()));

      ithtr++;

    }     

  
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

  //digi variables
  tree_->Branch("digi_wheel",&digi_wheel);
  tree_->Branch("digi_sector",&digi_sector);
  tree_->Branch("digi_station",&digi_station);
  tree_->Branch("digi_sl",&digi_sl);
  tree_->Branch("digi_layer",&digi_layer);
  tree_->Branch("digi_wire",&digi_wire);
  tree_->Branch("digi_time",&digi_time);

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

  //CSC segment variables
  tree_->Branch("cscsegm_ring",&cscsegm_ring);
  tree_->Branch("cscsegm_chamber",&cscsegm_chamber);
  tree_->Branch("cscsegm_station",&cscsegm_station);
  tree_->Branch("cscsegm_cosx",&cscsegm_cosx);
  tree_->Branch("cscsegm_cosy",&cscsegm_cosy);
  tree_->Branch("cscsegm_cosz",&cscsegm_cosz);
  tree_->Branch("cscsegm_phi",&cscsegm_phi);
  tree_->Branch("cscsegm_eta",&cscsegm_eta);
  tree_->Branch("cscsegm_normchisq",&cscsegm_normchi2);
  tree_->Branch("cscsegm_nRecHits",&cscsegm_nRecHits);
 
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

  tree_->Branch("ltTwinMux_thWheel",&ltTwinMux_thWheel);
  tree_->Branch("ltTwinMux_thSector",&ltTwinMux_thSector);
  tree_->Branch("ltTwinMux_thStation",&ltTwinMux_thStation);
  tree_->Branch("ltTwinMux_thBx",&ltTwinMux_thBx);
  tree_->Branch("ltTwinMux_thHits",&ltTwinMux_thHits);

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

  //GMT
  tree_->Branch("gmt_bx",&gmt_bx);
  tree_->Branch("gmt_phi",&gmt_phi);
  tree_->Branch("gmt_eta",&gmt_eta);
  tree_->Branch("gmt_pt",&gmt_pt);
  tree_->Branch("gmt_charge",&gmt_charge);
  tree_->Branch("gmt_qual",&gmt_qual);
  tree_->Branch("gmt_tf_idx",&gmt_tf_idx);

  //GT  // legacy
  tree_->Branch("gt_algo_bit",&gt_algo_bit);
  tree_->Branch("gt_algo_bx",&gt_algo_bx);
  tree_->Branch("gt_tt_bit",&gt_tt_bit);
  tree_->Branch("gt_tt_bx",&gt_tt_bx);

  //RPC
  tree_->Branch(  "rpc_region",&rpc_region);			   
  tree_->Branch(  "rpc_clusterSize",&rpc_clusterSize);		   
  tree_->Branch(  "rpc_strip",&rpc_strip);			   
  tree_->Branch(  "rpc_bx",&rpc_bx);			   
  tree_->Branch(  "rpc_station",&rpc_station);		   
  tree_->Branch(  "rpc_sector",&rpc_sector);		   
  tree_->Branch(  "rpc_layer",&rpc_layer);			   
  tree_->Branch(  "rpc_subsector",&rpc_subsector);		   
  tree_->Branch(  "rpc_roll",&rpc_roll);			   
  tree_->Branch(  "rpc_ring",&rpc_ring);                     
  
  //counters
  tree_->Branch("Ndigis",&idigis,"Ndigis/S");
  tree_->Branch("Ndtsegments",&idtsegments,"Ndtsegments/S");
  tree_->Branch("Ncscsegments",&icscsegments,"Ncscsegments/S");
  tree_->Branch("NdtltTwinMuxOut",&idtltTwinMuxOut,"NdtltTwinMuxOut/S");
  tree_->Branch("NdtltTwinMux_th",&idtltTwinMux_th,"NdtltTwinMux_th/S");
  tree_->Branch("NdtltTwinMuxIn",&idtltTwinMuxIn,"NdtltTwinMuxIn/S");
  tree_->Branch("Nmuons",&imuons,"Nmuons/S");
  tree_->Branch("Ngmt",&igmt,"Ngmt/S");
  tree_->Branch("Ngtalgo",&igtalgo,"Ngtalgo/S");
  tree_->Branch("Ngttechtrig",&igttt,"Ngttt/S");
  tree_->Branch("Nhlt",&ihlt,"Nhlt/S");
  tree_->Branch("NrpcRecHits",&irpcrechits,"NrpcRecHits/S");
  
   // Bmtf
   tree_->Branch("bmtfPt", &Bmtf_Pt);
   tree_->Branch("bmtfEta", &Bmtf_Eta);
   tree_->Branch("bmtfPhi", &Bmtf_Phi);
   tree_->Branch("bmtfGlobalPhi", &Bmtf_GlobalPhi);
   tree_->Branch("bmftFineBit", &Bmtf_FineBit);
   tree_->Branch("bmtfqual", &Bmtf_qual);
   tree_->Branch("bmtfch", &Bmtf_ch);
   tree_->Branch("bmtfbx", &Bmtf_bx);
   tree_->Branch("bmtfprocessor", &Bmtf_processor);
   tree_->Branch("bmtfwh", &Bmtf_wh);
   tree_->Branch("bmtftrAddress", &Bmtf_trAddress);
   tree_->Branch("bmtfSize", &Bmtf_Size);
   
   tree_->Branch("bmtfPhSize", &Bmtf_phSize);
   tree_->Branch("bmtfPhBx", &Bmtf_phBx);
   tree_->Branch("bmtfPhWh", &Bmtf_phWh);
   tree_->Branch("bmtfPhSe", &Bmtf_phSe);
   tree_->Branch("bmtfPhSt", &Bmtf_phSt);
   tree_->Branch("bmtfPhAng", &Bmtf_phAng);
   tree_->Branch("bmtfPhBandAng", &Bmtf_phBandAng);
   tree_->Branch("bmtfPhCode", &Bmtf_phCode);
   tree_->Branch("bmtfPhTs2Tag", &Bmtf_phTs2Tag);
   
   tree_->Branch("bmtfThSize", &Bmtf_thSize);
   tree_->Branch("bmtfThBx", &Bmtf_thBx);
   tree_->Branch("bmtfThWh", &Bmtf_thWh);
   tree_->Branch("bmtfThSe", &Bmtf_thSe);
   tree_->Branch("bmtfThSt", &Bmtf_thSt);
   tree_->Branch("bmtfThTheta", &Bmtf_thTheta);
   tree_->Branch("bmtfThCode", &Bmtf_thCode);
  

   // Unpacking RPC
   tree_->Branch("NirpcdigiTwinMux", &irpcdigi_TwinMux);
   tree_->Branch("RpcDigiTwinMuxBx", &RpcDigi_TwinMux_bx);
   tree_->Branch("RpcDigiTwinMuxStrip", &RpcDigi_TwinMux_strip);
   tree_->Branch("RpcDigiTwinMuxRegion", &RpcDigi_TwinMux_region);
   tree_->Branch("RpcDigiTwinMuxRing", &RpcDigi_TwinMux_ring);
   tree_->Branch("RpcDigiTwinMuxStation", &RpcDigi_TwinMux_station);
   tree_->Branch("RpcDigiTwinMuxLayer", &RpcDigi_TwinMux_layer);
   tree_->Branch("RpcDigiTwinMuxSector", &RpcDigi_TwinMux_sector);
   tree_->Branch("RpcDigiTwinMuxSubSector", &RpcDigi_TwinMux_subsector);
   tree_->Branch("RpcDigiTwinMuxRoll", &RpcDigi_TwinMux_roll);
   tree_->Branch("RpcDigiTwinMuxTrIndex", &RpcDigi_TwinMux_trIndex);
   tree_->Branch("RpcDigiTwinMuxDet", &RpcDigi_TwinMux_det);
   tree_->Branch("RpcDigiTwinMuxSubdetId", &RpcDigi_TwinMux_subdetId);
   tree_->Branch("RpcDigiTwinMuxRawId", &RpcDigi_TwinMux_rawId);

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
   tree_->Branch("RpcRechitTwinMuxGlobX", &RpcRechit_TwinMux_Glob_x);
   tree_->Branch("RpcRechitTwinMuxGlobY", &RpcRechit_TwinMux_Glob_y);
   tree_->Branch("RpcRechitTwinMuxGlobZ", &RpcRechit_TwinMux_Glob_z);
   tree_->Branch("RpcRechitTwinMuxGlobEta", &RpcRechit_TwinMux_Glob_eta);
   tree_->Branch("RpcRechitTwinMuxGlobPhi", &RpcRechit_TwinMux_Glob_phi);
      
   tree_->Branch("DTextrapolatedOnRPCBX", &DT_extrapolated_OnRPC_BX);
   tree_->Branch("DTextrapolatedOnRPCLocX", &DT_extrapolated_OnRPC_Loc_x);
   tree_->Branch("DTextrapolatedOnRPCLocY", &DT_extrapolated_OnRPC_Loc_y);
   tree_->Branch("DTextrapolatedOnRPCLocZ", &DT_extrapolated_OnRPC_Loc_z);
   tree_->Branch("DTextrapolatedOnRPCGlobX", &DT_extrapolated_OnRPC_Glob_x);
   tree_->Branch("DTextrapolatedOnRPCGlobY", &DT_extrapolated_OnRPC_Glob_y);
   tree_->Branch("DTextrapolatedOnRPCGlobZ", &DT_extrapolated_OnRPC_Glob_z);
   tree_->Branch("DTextrapolatedOnRPCGlobEta", &DT_extrapolated_OnRPC_Glob_eta);
   tree_->Branch("DTextrapolatedOnRPCGlobPhi", &DT_extrapolated_OnRPC_Glob_phi);
   tree_->Branch("DTextrapolatedOnRPCRegion", &DT_extrapolated_OnRPC_Region);
   tree_->Branch("DTextrapolatedOnRPCSector", &DT_extrapolated_OnRPC_Sector);
   tree_->Branch("DTextrapolatedOnRPCStation", &DT_extrapolated_OnRPC_Station);
   tree_->Branch("DTextrapolatedOnRPCLayer", &DT_extrapolated_OnRPC_Layer);
   tree_->Branch("DTextrapolatedOnRPCRoll", &DT_extrapolated_OnRPC_Roll);
   tree_->Branch("DTextrapolatedOnRPCRing", &DT_extrapolated_OnRPC_Ring);
   tree_->Branch("DTextrapolatedOnRPCStripw", &DT_extrapolated_OnRPC_Stripw);
   tree_->Branch("NDTsegmentonRPC", &DT_segment_onRPC);
      

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
  digi_wheel.clear();
  digi_sector.clear();
  digi_station.clear();
  digi_sl.clear();
  digi_layer.clear();
  digi_wire.clear();
  digi_time.clear();

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

  //CSC segment variables 
  cscsegm_ring.clear();
  cscsegm_chamber.clear();
  cscsegm_station.clear();
  cscsegm_cosx.clear();
  cscsegm_cosy.clear();
  cscsegm_cosz.clear();
  cscsegm_phi.clear();
  cscsegm_eta.clear();
  cscsegm_normchi2.clear();
  cscsegm_nRecHits.clear();

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

  ltTwinMux_thWheel.clear();
  ltTwinMux_thSector.clear();
  ltTwinMux_thStation.clear();
  ltTwinMux_thBx.clear();
  ltTwinMux_thHits.clear();

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

  //GMT
  gmt_bx.clear();
  gmt_phi.clear();
  gmt_eta.clear();
  gmt_pt.clear();
  gmt_qual.clear();
  gmt_charge.clear();
  gmt_tf_idx.clear();

  //GT
  gt_algo_bit.clear();
  gt_algo_bx.clear();
  gt_tt_bit.clear();
  gt_tt_bx.clear();

  //HLT
  hlt_path.clear();

  // RPC rec hits
  rpc_region.clear();
  rpc_clusterSize.clear();
  rpc_strip.clear();
  rpc_bx.clear();
  rpc_station.clear();
  rpc_sector.clear();
  rpc_layer.clear();
  rpc_subsector.clear();
  rpc_roll.clear();
  rpc_ring.clear();
  
        //Bmtf_Size.clear();
   Bmtf_Pt.clear();
   Bmtf_Eta.clear();
   Bmtf_Phi.clear();
   Bmtf_GlobalPhi.clear();
   Bmtf_qual.clear();
   Bmtf_ch.clear();
   Bmtf_bx.clear();
   Bmtf_processor.clear();
   Bmtf_trAddress.clear();
   Bmtf_wh.clear();
   Bmtf_FineBit.clear();

	Bmtf_phBx.clear();
	Bmtf_phWh.clear();
	Bmtf_phSe.clear();
	Bmtf_phSt.clear();
	Bmtf_phAng.clear();
	Bmtf_phBandAng.clear();
	Bmtf_phCode.clear();
	Bmtf_phTs2Tag.clear();

	Bmtf_thBx.clear();
	Bmtf_thWh.clear();
	Bmtf_thSe.clear();
  	Bmtf_thSt.clear();
	Bmtf_thTheta.clear();
	Bmtf_thCode.clear();

   RpcDigi_TwinMux_bx.clear();
   RpcDigi_TwinMux_strip.clear();
   RpcDigi_TwinMux_region.clear();
   RpcDigi_TwinMux_ring.clear();
   RpcDigi_TwinMux_station.clear();
   RpcDigi_TwinMux_layer.clear();
   RpcDigi_TwinMux_sector.clear();
   RpcDigi_TwinMux_subsector.clear();
   RpcDigi_TwinMux_roll.clear();
   RpcDigi_TwinMux_trIndex.clear();
   RpcDigi_TwinMux_det.clear();
   RpcDigi_TwinMux_subdetId.clear();
   RpcDigi_TwinMux_rawId.clear();

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
   RpcRechit_TwinMux_Glob_x.clear();
   RpcRechit_TwinMux_Glob_y.clear();
   RpcRechit_TwinMux_Glob_z.clear(); 
   RpcRechit_TwinMux_Glob_eta.clear(); 
   RpcRechit_TwinMux_Glob_phi.clear(); 

   DT_extrapolated_OnRPC_BX.clear();
   DT_extrapolated_OnRPC_Loc_x.clear();
   DT_extrapolated_OnRPC_Loc_y.clear();
   DT_extrapolated_OnRPC_Loc_z.clear();
   DT_extrapolated_OnRPC_Glob_x.clear();
   DT_extrapolated_OnRPC_Glob_y.clear();
   DT_extrapolated_OnRPC_Glob_z.clear();
   DT_extrapolated_OnRPC_Glob_eta.clear();
   DT_extrapolated_OnRPC_Glob_phi.clear();
   DT_extrapolated_OnRPC_Region.clear();
   DT_extrapolated_OnRPC_Sector.clear();
   DT_extrapolated_OnRPC_Station.clear();
   DT_extrapolated_OnRPC_Layer.clear();
   DT_extrapolated_OnRPC_Roll.clear();
   DT_extrapolated_OnRPC_Ring.clear();
   DT_extrapolated_OnRPC_Stripw.clear();
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
 
