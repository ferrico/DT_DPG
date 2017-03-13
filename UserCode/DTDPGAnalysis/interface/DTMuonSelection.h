//
// Original Author:  Mario Pelliccioni, Gianluca Cerminara
//         Created:  Tue Sep  9 15:56:24 CEST 2008

#include <TRandom3.h>

// user include files
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class DTMuonSelection : public edm::EDFilter {
public:

  explicit DTMuonSelection(const edm::ParameterSet&);

  ~DTMuonSelection();
  
private:
  virtual void beginJob() ;

  virtual bool filter(edm::Event&, const edm::EventSetup&);

  virtual void endJob() ;
  
  edm::InputTag muonList;
  edm::InputTag SAmuonList;
  edm::InputTag dtSegmentLabel;

  edm::EDGetTokenT<reco::MuonCollection>muonListToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection>dtSegmentToken_;
  edm::EDGetTokenT<reco::TrackCollection>SAmuonToken_;

  unsigned int tightness;
  double etaMin;
  double etaMax;
  double ptMin;

  TRandom3 *rng;
};
