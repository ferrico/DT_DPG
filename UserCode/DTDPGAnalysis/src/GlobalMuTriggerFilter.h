//
// Original Author:  Maria Fouz Iglesias
//         Created:  Thu Oct 30 12:08:09 CEST 2010
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutRecord.h"

#include <string>
#include <iostream>

using namespace edm;
using namespace std;

//
// class decleration
//

class GlobalMuTriggerFilter : public edm::EDFilter {
   public:
      explicit GlobalMuTriggerFilter(const edm::ParameterSet&);
      ~GlobalMuTriggerFilter();


   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag theL1MuonGlobalTriggerLabel;
      edm::EDGetTokenT<L1MuGMTReadoutCollection> theL1MuonGlobalTriggerToken_ ;

};

