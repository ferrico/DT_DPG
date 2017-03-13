//
// Original Author:  Mario Pelliccioni
//         Created:  Tue May 4 15:56:24 CEST 2010

// user include files
#include "UserCode/DTDPGAnalysis/interface/DTMuonSelection.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include <iostream>
#include <stdio.h>

using namespace std;

DTMuonSelection::DTMuonSelection(const edm::ParameterSet& iConfig)
{

  etaMin               = iConfig.getParameter<double>("etaMin");
  etaMax               = iConfig.getParameter<double>("etaMax");
  ptMin                = iConfig.getParameter<double>("ptMin");
  tightness            = iConfig.getParameter<int>("tightness");

  dtSegmentLabel       = iConfig.getParameter<edm::InputTag>("dtSegmentLabel");
  dtSegmentToken_      = consumes<DTRecSegment4DCollection>(edm::InputTag(dtSegmentLabel)); 
  //muonList           = iConfig.getParameter<edm::InputTag>("Muons");
  muonList             = iConfig.getParameter<edm::InputTag>("src");
  muonListToken_       = consumes<reco::MuonCollection>(edm::InputTag(muonList));
  SAmuonList           = iConfig.getParameter<edm::InputTag>("SAMuons");
  SAmuonToken_         = consumes<reco::TrackCollection>(edm::InputTag(SAmuonList));

  rng = new TRandom3(0);

  char sel[100];

       if (tightness==0) sprintf (sel,"whatever DT segment (loose Minimum Bias selection)");
  else if (tightness==1) sprintf (sel,"whatever muon (loose cosmic selection)");
  else if (tightness==2) sprintf (sel,"a Global muon (standard collision selection)");
  else if (tightness>2)  sprintf (sel," event number %u",tightness);

  cout<<endl<<"********************************************************"<<endl
            <<"* WARNING applying DTMuon filter                       *"<<endl
            <<"* requiring "<<       sel                                <<endl;
  if (tightness>0 && tightness<3)  {
      cout  <<"* with pt > "<<ptMin<<" and "                            <<endl
            <<"* "<<etaMin<<" < eta < "<<etaMax                         <<endl;
  }
  cout      <<"********************************************************"<<endl;
}


DTMuonSelection::~DTMuonSelection() { }


bool DTMuonSelection::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool result = false;

  if (tightness==0) {

   //Retrieve the DT segment list
   edm::Handle<DTRecSegment4DCollection> dtsegments4D;
   //iEvent.getByLabel(dtSegmentLabel, dtsegments4D);;// Doesn't work after 75
   iEvent.getByToken(dtSegmentToken_,dtsegments4D);

  for (DTRecSegment4DCollection::id_iterator chambIt = dtsegments4D->id_begin(); chambIt != dtsegments4D->id_end(); ++chambIt){
    result=true;
    break;
  }
  return result;
  }

  else if  (tightness==1) {

   // Get the RecTrack collection from the event
    edm::Handle<reco::TrackCollection> staTracks;
    //iEvent.getByLabel(SAmuonList, staTracks);
    iEvent.getByToken(SAmuonToken_,staTracks);

    reco::TrackCollection::const_iterator staTrack;

    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack){
      if (staTrack->pt() > ptMin && staTrack->eta() > etaMin && staTrack->eta() < etaMax){
        result = true;
        break;
      }
    }
  }

  else if  (tightness==2) {

    //Retrieve the muons list
    edm::Handle<reco::MuonCollection> MuHandle;
    //iEvent.getByLabel(muonList,MuHandle);// Doesn't work after 75
    iEvent.getByToken(muonListToken_,MuHandle);

    for (reco::MuonCollection::const_iterator nmuon = MuHandle->begin(); nmuon != MuHandle->end(); ++nmuon){

      if(nmuon->pt() > ptMin && nmuon->eta() > etaMin && nmuon->eta() < etaMax && nmuon->isGlobalMuon()){
        result = true;
        break;
      }
    }
  }
  else if (tightness>2) {
    if (iEvent.eventAuxiliary().event() == tightness) result = true;

  }
  //too many events in the dataset, throw away half of them
  //const Double_t rng_choice = rng->Uniform();
  //if(rng_choice < 0.5) return false;

  return result;
}



// ------------ method called once each job just before starting event loop  ------------
void  DTMuonSelection::beginJob() {
}



// ------------ method called once each job just after ending the event loop  ------------
void  DTMuonSelection::endJob() {
}
