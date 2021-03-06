/*
 *  See header file for a description of this class.
 *
 *  Update for reading TwinMux instead of DCC - February 2016
 *
 *  $Date: 2013/07/04 17:16:25 $
 *  \author M.C Fouz   
 */


#include <UserCode/DTDPGAnalysis/src/DTDPGCreateWheelSummary.h>

// Framework
#include <FWCore/Framework/interface/Event.h>
#include "DataFormats/Common/interface/Handle.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>


// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "boost/filesystem.hpp"
#include "TGraph.h"
#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "TPaletteAxis.h"

using namespace boost::filesystem;
using namespace edm;
using namespace std;


DTDPGCreateWheelSummary::DTDPGCreateWheelSummary(const edm::ParameterSet& ps){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Constructor";

  myParameters = ps;
  myRunNumber  = ps.getUntrackedParameter<int>("runNumber",0);
  ProcessDDUTrigger = ps.getUntrackedParameter<bool>("IncludeDDUTrigger","False");

  iHwMax=1;
  if(ProcessDDUTrigger)iHwMax=2;

}

DTDPGCreateWheelSummary::~DTDPGCreateWheelSummary(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Destructor called";

}


void DTDPGCreateWheelSummary::beginJob(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateWheelSummary]: beginJob";

}


void DTDPGCreateWheelSummary::beginRun(const edm::Run& run,const edm::EventSetup& context){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateWheelSummary]: beginRun";

  context.get<MuonGeometryRecord>().get(myMuonGeom);

}

void DTDPGCreateWheelSummary::analyze(const edm::Event& e, const edm::EventSetup& context){
  
}

void DTDPGCreateWheelSummary::createGifFile(string fileName, TCanvas *canvas, bool isExtraFile) {

  stringstream gifBase, gifTag;
  gifBase << "Run" << myRunNumber ;
  if (isExtraFile) { gifBase << "/ExtraPlots"; }
  try {
    create_directories(gifBase.str());
  } catch(const std::exception & ex) {
    throw cms::Exception("DTDPGCreateSummaryError")<< "[DTDPGCreateWheelSummary]: Excepiton " << ex.what() << " thrown creating " << gifBase.str() << " directory" << endl;
  }
  
  gifTag << "_r" << myRunNumber;
  string gifFile = gifBase.str() + "/" + fileName + gifTag.str() + ".gif";
  canvas->Print(gifFile.c_str());

}

void DTDPGCreateWheelSummary::createGifFile(string fileName, TCanvas *canvas, int wh, bool isExtraFile) {

  stringstream gifBase, gifTag;
  gifBase << "Run" << myRunNumber << "/Wheel" << showpos << wh;
  if (isExtraFile) { gifBase << "/ExtraPlots"; }

  try {
    create_directories(gifBase.str());
  } catch(const std::exception & ex) {
    throw cms::Exception("DTDPGCreateSummaryError")<< "[DTDPGCreateWheelSummary]: Excepiton " << ex.what() << " thrown creating " << gifBase.str() << " directory" << endl;
  }
  
  gifTag << "_r" << myRunNumber <<"_W" << wh;
  string gifFile = gifBase.str() + "/" + fileName + gifTag.str() + ".gif";
  canvas->Print(gifFile.c_str());

}

void DTDPGCreateWheelSummary::endJob(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateWheelSummary] endJob called!";

  // The root file which contain the histos
  myFile = new TFile(myParameters.getUntrackedParameter<string>("rootFileName", "DQM_DT.root").c_str(), "READ");


  //FRC: prepare for automatc writing out of dead channels!

  char deadname[100];
  sprintf (deadname, "Run%i/DeadChannelList_r%i.txt",myRunNumber,myRunNumber);
  edm::LogVerbatim ("DTDPGSummary") << " opening txt dead channel file ";
  DeadChannelList = new ofstream(deadname);

  char cMainFolder[30];sprintf(cMainFolder,"DQMData/Run %d", myRunNumber);
  TFolder *mfolder=(TFolder*)myFile->Get(cMainFolder);
  if(!mfolder)
  {
    throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: Folder = " << cMainFolder << " does not exist!!! Check the run number" << endl;
  }

  myMainFolder.append(cMainFolder);
  myMainFolder.append("/DT/Run summary/");

  myL1TFolder.append(cMainFolder);
  myL1TFolder.append("/L1T/Run summary/");


  char SLname[3][20]={"Phi1","Theta","Phi2"};

  //Change standard palette colors
  const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  const Int_t NCont = 155;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.99, 0.40, 0.00, 0.80, 0.95 };
  Double_t green[NRGBs] = { 0.99, 0.40, 0.40, 0.10, 0.00 };
  Double_t blue[NRGBs]  = { 0.20, 0.00, 0.90, 0.10, 0.00 };
  // It needs to be initialize if not root gives a warning and cmssw an error and
  // the program stops
  TColor::InitializeColors();
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetCanvasColor(0);
  //gStyle->SetPadColor(0);


  // >>>>>>>>>>>>>>>>   RECO HITS PLOTS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary histos with number of reco hits      
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  // CREATE HISTOS...................................
  TH2F *  RateHitsSeg[5];
  TH2F *  TotHitsSeg[5];

  for(int iw=-2;iw<=2;iw++)
    {
      char thehtit[240]; char thehistoname[240];

      sprintf(thehtit,"RateHitsSegW%d",iw);
      sprintf(thehistoname,"Number of hits in Reco Segment - Percentage W%d",iw);
      RateHitsSeg[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,42,-1,83);
      sprintf(thehtit,"TotHitsSegW%d",iw);
      sprintf(thehistoname,"Number of hits in Reco Segment  W%d",iw);
      TotHitsSeg[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,42,-1,83);

      // Put Labels
      for(int isec=1;isec<15;isec++)
	{
	  char  EWlabel[40];
	  sprintf(EWlabel,"Sect-%d",isec);
	  RateHitsSeg[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	}
      for(int inb=1;inb<42+1;inb++)RateHitsSeg[iw+2]->GetYaxis()->SetBinLabel(inb,"");
      for(int imb=1;imb<5;imb++)
	{
	  RateHitsSeg[iw+2]->GetYaxis()->SetBinLabel(3+(imb-1)*10,"3-4");
	  RateHitsSeg[iw+2]->GetYaxis()->SetBinLabel(5+(imb-1)*10,"7-8");
	  RateHitsSeg[iw+2]->GetYaxis()->SetBinLabel(7+(imb-1)*10,"11-12");
	  RateHitsSeg[iw+2]->GetYaxis()->SetBinLabel(9+(imb-1)*10,"15-16");

	  TotHitsSeg[iw+2]->GetYaxis()->SetBinLabel(3+(imb-1)*10,"3-4");
	  TotHitsSeg[iw+2]->GetYaxis()->SetBinLabel(5+(imb-1)*10,"7-8");
	  TotHitsSeg[iw+2]->GetYaxis()->SetBinLabel(7+(imb-1)*10,"11-12");
	  TotHitsSeg[iw+2]->GetYaxis()->SetBinLabel(9+(imb-1)*10,"15-16");
	}
    }// End loop on wheels


  // FILL   HISTOS...................................
  //vector<DTChamber*>::const_iterator chIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chEnd = myMuonGeom->chambers().end();
  for (; chIt != chEnd; ++chIt) {
    DTChamberId ch = (*chIt)->id();
    stringstream wheel; wheel << ch.wheel();
    stringstream station; station << ch.station();
    stringstream sector; sector << ch.sector();
    int iw  = ch.wheel();
    int isec= ch.sector();
    int imb = ch.station();
   
    string recoFolder = myMainFolder + "02-Segments/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() ; 
    string histoName = recoFolder + "/h4DSegmNHits_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();  

    //TH2F * HSegHits = (TH2F*) myFile -> FindObjectAny(histo_name.c_str());
    TH2F * HSegHits = (TH2F*) myFile -> Get(histoName.c_str());
    if(HSegHits==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoName.c_str() << " not found";
    else
      {
	int totbin=HSegHits->GetXaxis()->GetNbins();
	for(int inb=1;inb<totbin+1;inb++)
	  {
	    float xcenter=HSegHits->GetBinCenter(inb);
	    float xval=HSegHits->GetBinContent(inb);
	    if(xcenter<19 && xval>0)
	      TotHitsSeg[iw+2]->Fill(float(isec),xcenter+(imb-1)*20,xval);
	  }
      }

  }// End loop on DT Chambers

  // COMPUTE RATES ..................................
  for(int iw=0;iw<5;iw++)
    for(int isec=1;isec<16;isec++)
      for(int imb=0;imb<15;imb++)
	{
	  float summAll=0;
	  for(int ih=1;ih<11;ih++)
	    summAll +=TotHitsSeg[iw]->GetBinContent(isec,imb*10+ih);

	  if(summAll>0)
	    for(int ih=1;ih<11;ih++)
	      {
		float xval=TotHitsSeg[iw]->GetBinContent(isec,imb*10+ih);
		RateHitsSeg[iw]->SetBinContent(isec,imb*10+ih,100*xval/summAll);
	      }
	}// end loop on wheels && sectors && chambers



  // DRAW   HISTOS...................................
  TCanvas * c0=new TCanvas("c0","Percentage of Hits on Segments",500,0,900,950);
  c0->SetTopMargin(0.36);
  c0->SetLeftMargin(0.05);
  c0->SetRightMargin(0.05);
  c0->SetBottomMargin(0.20);
  c0->SetFillColor(0);
  c0->Divide(2,3,0.001,0.001);


  TCanvas  * c1[5];
  c1[0]=new TCanvas("c1_1","Percentage of Hits on Segments",500,0,900,950);
  c1[1]=new TCanvas("c1_2","Percentage of Hits on Segments",500,0,900,950);
  c1[2]=new TCanvas("c1_3","Percentage of Hits on Segments",500,0,900,950);
  c1[3]=new TCanvas("c1_4","Percentage of Hits on Segments",500,0,900,950);
  c1[4]=new TCanvas("c1_5","Percentage of Hits on Segments",500,0,900,950);


  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
    {
      c1[iw+2]->SetTopMargin(0.10);
      c1[iw+2]->cd()->SetRightMargin(0.12);
      c1[iw+2]->cd()->SetFillColor(0);
      c1[iw+2]->cd()->SetBorderMode(0);
      c1[iw+2]->cd()->SetLineColor(4);

      c1[iw+2]->cd();
      RateHitsSeg[iw+2]->SetMaximum(90.);
      RateHitsSeg[iw+2]->SetMinimum(5.2);
      RateHitsSeg[iw+2]->GetXaxis()->SetLabelSize(0.03);
      RateHitsSeg[iw+2]->GetYaxis()->SetLabelSize(0.03);
      RateHitsSeg[iw+2]->Draw("colz");

      char titlename[300];
      sprintf(titlename,"Hits percentage on RecoSegments W%d",iw);
      TPaveLabel * title = new TPaveLabel(2.1,82,13.5,88,titlename);
      title->SetFillColor(0);
      title->SetTextColor(4);
      title->Draw();

      char Tit[100];
      int step=20;
      for(int imb=1;imb<5;imb++)
	{
	  sprintf(Tit,"MB-%d",imb);
	  TText * tTit  = new TText(13.1,10+(imb-1)*step,Tit) ;
	  tTit->SetTextSize(0.04);
	  tTit->Draw();
	}


      createGifFile("RateHitsSeg",c1[iw+2],iw);

      c0->cd(iw+3)->cd();
      RateHitsSeg[iw+2]->GetXaxis()->SetLabelSize(0.05);
      RateHitsSeg[iw+2]->GetYaxis()->SetLabelSize(0.06);
      RateHitsSeg[iw+2]->Draw("colz");
      title->Draw();
      for(int imb=1;imb<5;imb++)
	{
	  sprintf(Tit,"MB-%d",imb);
	  TText * tTit  = new TText(13.1,10+(imb-1)*step,Tit) ;
	  tTit->SetTextSize(0.04);
	  tTit->Draw();
	}



    }// End loop on wheels
  createGifFile("RateHitsSeg",c0);


  // >>>>>>>>>>>>>>>>   CHAMBER EFFICIENCY PLOTS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary histos with efficiency per chamber/superlayer/layer     
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  TH2F * EffForWheels[5];
  TH2F * SegForWheels[5];
  TH2F * EffForWheelsLayers[5];
  TH2F * SegForWheelsLayers[5];
  TH2F * EffForSectorsCells[5][14];
  TH2F * SegForSectorsCells[5][14];

  TH2F * EffSegReco[5];
  TH2F * SegReco[5];

  TH1F * MeanEffCells= new TH1F("MeanEffCells","Mean Efficiency per Cell",101, 0.,101.);
  TH1F * MeanEffLayers= new TH1F("MeanEffLayers","Mean Efficiency per Layer",101, 0.,101.);
  TH1F * MeanEffSuperLayers= new TH1F("MeanEffSuperLayers","Mean Efficiency per SuperLayer",101, 0.,101.);

  // CREATE HISTOS...................................
  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
    {
      char thehtit[240]; char thehistoname[240];
      sprintf(thehtit,"EfW%d",iw);
      sprintf(thehistoname,"Efficiencies per SuperLayer W%d",iw);
      EffForWheels[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,17,0.,17.);
      sprintf(thehtit,"SegW%d",iw);
      sprintf(thehistoname,"Segments per SuperLayer W%d",iw);
      SegForWheels[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,17,0.,17.);

      sprintf(thehtit,"EfW%dLay",iw);
      sprintf(thehistoname,"Efficiencies per Layer W%d",iw);
      EffForWheelsLayers[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,75,0.,75.);
      sprintf(thehtit,"SegW%dLay",iw);
      sprintf(thehistoname,"Segments per Layer W%d",iw);
      SegForWheelsLayers[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,75,0.,75.);


      sprintf(thehtit,"EffSegRecoW%d",iw);
      sprintf(thehistoname,"Extrapolated Segment Reco Efficiencies W%d",iw);
      EffSegReco[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,4,1.,5.);
      sprintf(thehtit,"SegRecoW%d",iw);
      sprintf(thehistoname,"Extrapolated Segment RecoW%d",iw);
      SegReco[iw+2]= new TH2F(thehtit,thehistoname,14, 1.,15.,4,1.,5.);

      for(int isec=1;isec<15;isec++)
	{
	  sprintf(thehtit,"EfW%dSector%d",iw,isec);
	  sprintf(thehistoname,"Efficiencies per Cell W%dS%d",iw,isec);
	  EffForSectorsCells[iw+2][isec-1]= new TH2F(thehtit,thehistoname,100, 0.,100.,75,0.,75.);
	  sprintf(thehtit,"SegW%dSector%d",iw,isec);
	  sprintf(thehistoname,"Segments per Cell W%dS%d",iw,isec);
	  SegForSectorsCells[iw+2][isec-1]= new TH2F(thehtit,thehistoname,100, 0.,100.,75,0.,75.);
	}

      // Put Labels
      for(int isec=1;isec<15;isec++)
	{
	  char  EWlabel[40];
	  sprintf(EWlabel,"Sect-%d",isec);
	  EffForWheels[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	  SegForWheels[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	  EffForWheelsLayers[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	  SegForWheelsLayers[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	  EffSegReco[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	  SegReco[iw+2]->GetXaxis()->SetBinLabel(isec,EWlabel);
	}

      for(int imb=1;imb<5;imb++)
	{
	  char  EWlabel[40];
	  sprintf(EWlabel,"MB%d",imb);
	  EffSegReco[iw+2]->GetYaxis()->SetBinLabel(imb,EWlabel);
	  SegReco[iw+2]->GetYaxis()->SetBinLabel(imb,EWlabel);

	  for(int isl=1;isl<4;isl++)
	    {
	      sprintf(EWlabel,"MB%d_SL%s",imb,SLname[isl-1]);
	      EffForWheels[iw+2]->GetYaxis()->SetBinLabel((imb-1)*4+isl+1,EWlabel);
	      SegForWheels[iw+2]->GetYaxis()->SetBinLabel((imb-1)*4+isl+1,EWlabel);
	      for(int il=1;il<4;il++)
		{
		  EffForWheelsLayers[iw+2]->GetYaxis()->SetBinLabel((imb-1)*16+isl*4+il,"");
		  SegForWheelsLayers[iw+2]->GetYaxis()->SetBinLabel((imb-1)*16+isl*4+il,"");

		  for(int isec=1;isec<15;isec++)
		    {
		      EffForSectorsCells[iw+2][isec-1]->GetYaxis()->SetBinLabel((imb-1)*16+isl*4+il,"");
		      SegForSectorsCells[iw+2][isec-1]->GetYaxis()->SetBinLabel((imb-1)*16+isl*4+il,"");
		    }
		}
	    }
	}

    }// End Looop on wheels

  // FILL   HISTOS...................................
  chIt = myMuonGeom->chambers().begin();
  chEnd = myMuonGeom->chambers().end();
  for (; chIt != chEnd; ++chIt) {
    DTChamberId ch = (*chIt)->id();
    stringstream wheel; wheel << ch.wheel();
    stringstream station; station << ch.station();
    stringstream sector; sector << ch.sector();
    int iw  = ch.wheel();
    int isec= ch.sector();
    int imb = ch.station();
   
    // Chamber Efficiency Task
    string recoFolder = myMainFolder + "01-DTChamberEfficiency/Task/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str()  ;
    string histoNameReco  = recoFolder + "/hEffGoodSegVsPosDen_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();  
    string histoNameFound = recoFolder + "/hEffGoodCloseSegVsPosNum_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();  
    TH2F * HSegReco  = (TH2F*) myFile -> Get(histoNameReco.c_str());
    TH2F * HSegFound = (TH2F*) myFile -> Get(histoNameFound.c_str());

    if(HSegReco==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameReco.c_str() << " not found" << endl;
    else
      SegReco[iw+2]->Fill(isec,imb,HSegReco->GetEntries());

    if(HSegFound==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameFound.c_str() << " not found" << endl;
    else
      EffSegReco[iw+2]->Fill(isec,imb,HSegFound->GetEntries());

    // Efficiency Task
    for(int isl=1;isl<4;isl++) // Loop on SLs 
      {
	if(ch.station()==4 && isl==2) isl++;  // skip theta on MB4
	float all=0;
	float found=0;
	stringstream cisl; cisl<<isl;

	recoFolder = myMainFolder + "DTEfficiencyTask/Wheel" + wheel.str() + "/Station" + station.str() + "/Sector" + sector.str() + "/SuperLayer" + cisl.str() ;

	for(int il=1;il<5;il++) // Loop on Layers 
	  {
	    float lfound=0; float lall=0;
	    stringstream cil; cil<<il;
	    histoNameReco = recoFolder + "/hRecSegmOccupancy_W" + wheel.str() + "_St" + station.str() 
	      + "_Sec" + sector.str() + "_SL" + cisl.str()+ "_L" + cil.str();  
	    TH1F * HSeg = (TH1F*) myFile -> Get(histoNameReco.c_str());
	    if(HSeg==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = "  << histoNameReco.c_str() << " not found" << endl;
	    else
	      {
		all+=HSeg->GetEntries();
		lall=HSeg->GetEntries();
	      }

	    histoNameFound = recoFolder + "/hEffUnassOccupancy_W" + wheel.str() + "_St" + station.str() 
	      + "_Sec" + sector.str() + "_SL" + cisl.str()+ "_L" + cil.str();  
	    TH1F * HFound = (TH1F*) myFile -> Get(histoNameFound.c_str());
	    if(HFound==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameFound.c_str() << " not found" << endl;
	    else
	      {
		found+=HFound->GetEntries();
		lfound=HFound->GetEntries();
	      }

	    if(lall>0)
	      { 
		EffForWheelsLayers[iw+2]->Fill(isec,(imb-1)*17+isl*5+il,lfound);
		SegForWheelsLayers[iw+2]->Fill(isec,(imb-1)*17+isl*5+il,lall);
	      }


	    if(HSeg!=NULL && HFound!=NULL)
	      {
		int nb1=HSeg->GetNbinsX();
		int nb2=HFound->GetNbinsX();
		if(nb1 != nb2){ throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR : mistmach on number of bins" << endl;}
		for(int ibin=1;ibin<nb1;ibin++)
		  {
		    float xpos1=HSeg->GetBinCenter(ibin);
		    float xpos2=HFound->GetBinCenter(ibin);
		    if(xpos1 != xpos2){ throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR :  mistmach on xpos" << endl;}
		    float xval1=HSeg->GetBinContent(ibin);
		    float xval2=HFound->GetBinContent(ibin);
                    char printline[50];
		    if(xval1>0)
		      {
			EffForSectorsCells[iw+2][isec-1]->Fill(xpos1,(imb-1)*17+isl*5+il,xval2);
			SegForSectorsCells[iw+2][isec-1]->Fill(xpos1,(imb-1)*17+isl*5+il,xval1);

			//FRC: write out dead channels
                        // a channel is defined as dead if epsilon<0.5 and there are at least 6 crossing segments:
                        // this corresponds to a maximum Depsilon of 20.4% (at epsilon=0.5)
                        if ( xval1>5 && xval2/xval1<0.5 ) {

                              double effi=xval2/xval1;
                              //double deffi=sqrt(xval2*(xval1-xval2)/xval1)/xval1; // this is the canonical binomial error
                              double deffi = 1/2./sqrt(xval1); // this is the maximum error (at effi=0.5) for a given xval1

                              sprintf (printline,"%3i%3u%4u%3u%3u%4u%7.2f%7.2f",iw,imb,isec,isl,il,int(xpos1),effi*100,deffi*100);
                              *DeadChannelList<<printline<<endl;
			      //*DeadChannelList<<iw<<" "<<isec<<" "<<imb<<" "<<isl<<" "<<il<<" "<<xpos1<<" "<<effi*100<<" "<<deffi*100<<endl;
			}
		        else if (xval1<6) {
                          sprintf (printline,"%3i%2u%3u%2u%2u%4u -1 %4u",iw,imb,isec,isl,il,int(xpos1),int(xval1));
		          *DeadChannelList<<printline<<endl;  // not enough segments crossing this cell!
		        }		       
		      }
		    else {
                      sprintf (printline,"%3i%2u%3u%2u%2u%4u -2 %4u",iw,imb,isec,isl,il,int(xpos1),int(xval1));
		      *DeadChannelList<<printline<<endl;  // no segments crossing this cell!
		    }
		  }
	      }
	  }// end loop layers
	if(all>0)
	  {
	    EffForWheels[iw+2]->Fill(isec,(imb-1)*4+isl,found);
	    SegForWheels[iw+2]->Fill(isec,(imb-1)*4+isl,all);
	  }
      }//End loop on SLs
  }// End loop on DT Chambers

  // COMPUTE EFFICIENCIES ...........................
  for(int iw=0;iw<5;iw++)
    {
      int nb11=EffForWheels[iw]->GetNbinsX();
      int nb21=EffForWheels[iw]->GetNbinsY();
      int nb12=SegForWheels[iw]->GetNbinsX();
      int nb22=SegForWheels[iw]->GetNbinsY();
      if(nb11 != nb12 || nb21 !=nb22 )
	{throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR : mistmach on number of bins" << endl;}

      for(int ibnx=1; ibnx<nb11+1; ibnx++)
	for(int ibny=1; ibny<nb21+1; ibny++)
	  {
	    float xval1=EffForWheels[iw]->GetBinContent(ibnx,ibny);
	    float xval2=SegForWheels[iw]->GetBinContent(ibnx,ibny);
	    if(xval2>0)
          {
	      EffForWheels[iw]->SetBinContent(ibnx,ibny,100*xval1/xval2);
            MeanEffSuperLayers->Fill(100*xval1/xval2);
          }
	  }

      nb11=EffForWheelsLayers[iw]->GetNbinsX();
      nb21=EffForWheelsLayers[iw]->GetNbinsY();
      nb12=SegForWheelsLayers[iw]->GetNbinsX();
      nb22=SegForWheelsLayers[iw]->GetNbinsY();
      if(nb11 != nb12 || nb21 !=nb22 )
	{throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR : mistmach on number of bins" << endl;}

      for(int ibnx=1; ibnx<nb11+1; ibnx++)
	for(int ibny=1; ibny<nb21+1; ibny++)
	  {
	    float xval1=EffForWheelsLayers[iw]->GetBinContent(ibnx,ibny);
	    float xval2=SegForWheelsLayers[iw]->GetBinContent(ibnx,ibny);
	    if(xval2>0)
          {
	      EffForWheelsLayers[iw]->SetBinContent(ibnx,ibny,100*xval1/xval2);
            MeanEffLayers->Fill(100*xval1/xval2);
          }
	  }

      nb11=EffSegReco[iw]->GetNbinsX();
      nb21=EffSegReco[iw]->GetNbinsY();
      nb12=SegReco[iw]->GetNbinsX();
      nb22=SegReco[iw]->GetNbinsY();
      if(nb11 != nb12 || nb21 !=nb22 )
	{throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR : mistmach on number of bins" << endl;}

      for(int ibnx=1; ibnx<nb11+1; ibnx++)
	for(int ibny=1; ibny<nb21+1; ibny++)
	  {
	    float xval1=EffSegReco[iw]->GetBinContent(ibnx,ibny);
	    float xval2=SegReco[iw]->GetBinContent(ibnx,ibny);
	    if(xval2>0)
	      EffSegReco[iw]->SetBinContent(ibnx,ibny,100*xval1/xval2);
	  }

      for(int isec=0;isec<14;isec++)
	{
	  nb11=EffForSectorsCells[iw][isec]->GetNbinsX();
	  nb21=EffForSectorsCells[iw][isec]->GetNbinsY();
	  nb12=SegForSectorsCells[iw][isec]->GetNbinsX();
	  nb22=SegForSectorsCells[iw][isec]->GetNbinsY();
	  if(nb11 != nb12 || nb21 !=nb22 )
	    {throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateWheelSummary]: ERROR : mistmach on number of bins" << endl;}

	  for(int ibnx=1; ibnx<nb11+1; ibnx++)
	    for(int ibny=1; ibny<nb21+1; ibny++)
	      {
		float xval1=EffForSectorsCells[iw][isec]->GetBinContent(ibnx,ibny);
		float xval2=SegForSectorsCells[iw][isec]->GetBinContent(ibnx,ibny);
		if(xval2>0)
            {
		  EffForSectorsCells[iw][isec]->SetBinContent(ibnx,ibny,100*xval1/xval2);
              MeanEffCells->Fill(100*xval1/xval2);
            }

	      }
	}
    } // END Compute Efficiencies  

  // DRAW   HISTOS...................................

  // color palette to have green when efficiency=1 
  int col[20]={100,100,98,98,97,97,96,96,95,94,93,92,91,90,89,88,87,86,85,78};
  gStyle->SetPalette(20,col);


  // SuperLayer Efficiencies
  TCanvas * ct2=new TCanvas("ct2","Efficiencies",500,0,900,950);

  ct2->SetTopMargin(0.26);
  ct2->SetLeftMargin(0.12);
  ct2->SetRightMargin(0.05);
  ct2->SetBottomMargin(0.20);
  ct2->SetFillColor(0);
  ct2->Divide(1,5,0.,0.05);

  char titlename[300];
  sprintf(titlename,"SuperLayer Efficiencies");
  TPaveLabel* title = new TPaveLabel(0.1,0.945,0.9,0.98,titlename);
  title->SetFillColor(0);
  title->SetTextColor(4);
  title->Draw();

  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
    {
      ct2->cd(iw+3)->SetRightMargin(0.1);
      ct2->cd(iw+3)->SetFillColor(0);
      ct2->cd(iw+3)->SetBorderMode(0);
      ct2->cd(iw+3)->SetLineColor(4);

      EffForWheels[iw+2]->SetMaximum(100.);
      EffForWheels[iw+2]->SetMinimum(70.2);
      EffForWheels[iw+2]->GetXaxis()->SetLabelSize(0.1);
      EffForWheels[iw+2]->GetYaxis()->SetLabelSize(0.11);
      EffForWheels[iw+2]->GetZaxis()->SetLabelSize(0.11); // Serves to change palette text size 
      if(iw==2) EffForWheels[iw+2]->GetYaxis()->SetLabelSize(0.09);
      EffForWheels[iw+2]->Draw("colz");

      char WheelTit[100]; sprintf(WheelTit,"Wheel %d",iw);
      TPaveLabel * box;
      if(iw==2)box = new TPaveLabel(0.80,0.30,0.92,0.45,WheelTit,"NDC");
      else
	box = new TPaveLabel(0.80,0.10,0.92,0.25,WheelTit,"NDC");
      box->SetFillColor(0);
      box->Draw();

    }// End loop on wheels
  createGifFile("EfficiencyperSuperLayer",ct2);

  // Layer Efficiencies
  TCanvas * ct3=new TCanvas("ct3","Efficiencies",500,0,900,950);

  ct3->SetTopMargin(0.26);
  ct3->SetLeftMargin(0.12);
  ct3->SetRightMargin(0.05);
  ct3->SetBottomMargin(0.20);
  ct3->SetFillColor(0);
  ct3->Divide(1,5,0.,0.05);

  sprintf(titlename,"Layer Efficiencies");
  title = new TPaveLabel(0.1,0.945,0.9,0.98,titlename);
  title->SetFillColor(0);
  title->SetTextColor(4);
  title->Draw();

  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
    {
      ct3->cd(iw+3)->SetRightMargin(0.1);
      ct3->cd(iw+3)->SetFillColor(0);
      ct3->cd(iw+3)->SetBorderMode(0);
      ct3->cd(iw+3)->SetLineColor(4);

      EffForWheelsLayers[iw+2]->SetMaximum(101.);
      EffForWheelsLayers[iw+2]->SetMinimum(70.2);
      EffForWheelsLayers[iw+2]->GetXaxis()->SetLabelSize(0.1);
      EffForWheelsLayers[iw+2]->GetYaxis()->SetLabelSize(0.11);
      EffForWheelsLayers[iw+2]->GetZaxis()->SetLabelSize(0.11); // Serves to change palette text size 
      if(iw==2) EffForWheelsLayers[iw+2]->GetYaxis()->SetLabelSize(0.09);
      EffForWheelsLayers[iw+2]->Draw("colz");


      char WheelTit[100]; sprintf(WheelTit,"Wheel %d",iw);
      TPaveLabel * box;
      if(iw==2)box = new TPaveLabel(0.80,0.30,0.92,0.45,WheelTit,"NDC");
      else
	box = new TPaveLabel(0.80,0.10,0.92,0.25,WheelTit,"NDC");
      box->SetFillColor(0);
      box->Draw();

      for(int imb=0;imb<4;imb++)
	{
          int step=imb*18;
          char MBName[10];sprintf(MBName,"MB%d",imb+1);
          TText * MB   = new TText(-0.65,10+step,MBName);
          MB->SetTextSize(0.12);
          if(iw==2)MB->SetTextSize(0.10);
          MB->SetTextColor(imb+1);
          MB->Draw();

          TText * phi1   = new TText(0.3,5+step,"Phi1");
          TText * theta  = new TText(0.3,10+step,"Theta");
          TText * phi2   = new TText(0.3,15+step,"Phi2") ;
          phi1->SetTextSize(0.075); theta->SetTextSize(0.075); phi2->SetTextSize(0.075);

          if(iw==2){phi1->SetTextSize(0.06); theta->SetTextSize(0.06); phi2->SetTextSize(0.06);}

          phi1->SetTextColor(imb+1); theta->SetTextColor(imb+1); phi2->SetTextColor(imb+1);

          phi1->Draw(); phi2->Draw();
          if(imb<3)theta->Draw();
	}

    }// End loop on wheels
  createGifFile("EfficiencyperLayer",ct3);

  // Cell Efficiencies
  TCanvas * ct4[5];
  ct4[0]=new TCanvas("ct4_1","Efficiencies W-2",500,0,900,950);
  ct4[1]=new TCanvas("ct4_2","Efficiencies W-1",500,0,900,950);
  ct4[2]=new TCanvas("ct4_3","Efficiencies W0 ",500,0,900,950);
  ct4[3]=new TCanvas("ct4_4","Efficiencies W1 ",500,0,900,950);
  ct4[4]=new TCanvas("ct4_5","Efficiencies W2 ",500,0,900,950);

  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
    {
      ct4[iw+2]->cd()->SetTopMargin(0.36);
      ct4[iw+2]->cd()->SetLeftMargin(0.15);
      ct4[iw+2]->cd()->SetRightMargin(0.05);
      ct4[iw+2]->cd()->SetBottomMargin(0.20);
      ct4[iw+2]->cd()->SetFillColor(0);
      ct4[iw+2]->cd()->Divide(2,6,0.,0.05);

      char titlename[300];
      sprintf(titlename,"Single Cell Efficiencies Wheel %d",iw);
      //title = new TPaveLabel(0.1,0.945,0.9,0.98,titlename);
      title = new TPaveLabel(0.1,0.955,0.9,0.99,titlename);
      title->SetFillColor(0);
      title->SetTextColor(4);
      title->Draw();

      for(int isec=1;isec<13;isec++) // Loop on sectors 
	{
	  int ctid=1+(isec-1)*2;
	  if(isec>6) ctid=(isec-6)*2;
	  ct4[iw+2]->cd(ctid)->SetRightMargin(0.1);
	  ct4[iw+2]->cd(ctid)->SetFillColor(0);
	  ct4[iw+2]->cd(ctid)->SetBorderMode(0);
	  ct4[iw+2]->cd(ctid)->SetLineColor(4);
      
	  EffForSectorsCells[iw+2][isec-1]->SetMaximum(101.);
	  EffForSectorsCells[iw+2][isec-1]->SetMinimum(70.2);
	  EffForSectorsCells[iw+2][isec-1]->GetXaxis()->SetLabelSize(0.1);
	  EffForSectorsCells[iw+2][isec-1]->GetYaxis()->SetLabelSize(0.11);
	  EffForSectorsCells[iw+2][isec-1]->GetZaxis()->SetLabelSize(0.11);// Serves to change palette text size 
	  if(iw==2) EffForSectorsCells[iw+2][isec-1]->GetYaxis()->SetLabelSize(0.09);
	  if(isec==6 || isec==12){
	    EffForSectorsCells[iw+2][isec-1]->GetXaxis()->SetTitleSize(0.1);
	    EffForSectorsCells[iw+2][isec-1]->GetXaxis()->SetTitleOffset(0.85);
	    EffForSectorsCells[iw+2][isec-1]->GetXaxis()->SetTitle("# cell");
	  }
	  EffForSectorsCells[iw+2][isec-1]->Draw("colz");
      
	  char WheelTit[100]; sprintf(WheelTit,"Sector %d",isec);
	  TPaveLabel * box;
	  if(isec==6 || isec==12)box = new TPaveLabel(0.70,0.40,0.89,0.25,WheelTit,"NDC");
	  else
	    box = new TPaveLabel(0.70,0.05,0.89,0.20,WheelTit,"NDC");
	  box->SetFillColor(0);
	  box->Draw();
      
	  if(isec<7)
	    for(int imb=0;imb<4;imb++)
	      {
		int step=imb*18;
		char MBName[10];sprintf(MBName,"MB%d",imb+1);
		TText * MB   = new TText(-20.3,10+step,MBName);
		MB->SetTextSize(0.12);
		if(isec==6 || isec==12)MB->SetTextSize(0.10);
		MB->SetTextColor(imb+1);
		MB->Draw();
      
      
		TText * phi1   = new TText(-8.5,5+step,"Phi1");
		TText * theta  = new TText(-8.5,10+step,"Theta");
		TText * phi2   = new TText(-8.5,15+step,"Phi2") ;
		phi1->SetTextSize(0.075); theta->SetTextSize(0.075); phi2->SetTextSize(0.075);
      
		if(isec==6 || isec==12){phi1->SetTextSize(0.06); theta->SetTextSize(0.06); phi2->SetTextSize(0.06);}
      
		phi1->SetTextColor(imb+1); theta->SetTextColor(imb+1); phi2->SetTextColor(imb+1);
      
		phi1->Draw(); phi2->Draw();
		if(imb<3)theta->Draw();
	      }
      
	}// End loop on sectors

      createGifFile("EfficiencyperCell",ct4[iw+2],iw);

    }// End loop on wheels

  gStyle->SetOptStat(1);
  gStyle->SetStatX(0.45);
  gStyle->SetStatW(0.18);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.18);
  
  gStyle->SetOptFit(1);
  //Average of efficiencies 
  TCanvas * ct201=new TCanvas("ct201","Efficiencies per Cell",500,0,900,950);
  ct201->cd()->SetTopMargin(0.10);
  ct201->cd()->SetRightMargin(0.12);
  ct201->cd()->SetFillColor(0);
  ct201->cd()->SetBorderMode(0);
  ct201->cd()->SetLineColor(4);

  ct201->cd();

  MeanEffCells->SetStats(1);
  MeanEffCells->SetLineColor(2);
  MeanEffCells->SetFillColor(5);
  MeanEffCells->Draw();
  try { MeanEffCells->Fit("gaus","Q"); } 
  catch(const cms::Exception&) {
     edm::LogError("DTDPGSummary") << "[DTDPGCreateWheelSummary]:  Error fitting Mean Efficiency Cells";
  }
  MeanEffCells->Draw("same");
  createGifFile("MeanEffCell",ct201,true);

  TCanvas * ct202=new TCanvas("ct202","Efficiencies per Layer",500,0,900,950);
  ct202->cd()->SetTopMargin(0.10);
  ct202->cd()->SetRightMargin(0.12);
  ct202->cd()->SetFillColor(0);
  ct202->cd()->SetBorderMode(0);
  ct202->cd()->SetLineColor(4);

  ct202->cd();

  MeanEffLayers->SetStats(1);
  MeanEffLayers->SetLineColor(2);
  MeanEffLayers->SetFillColor(5);
  MeanEffLayers->Draw();
  try { MeanEffLayers->Fit("gaus","Q"); } 
  catch(const cms::Exception&) {
     edm::LogError("DTDPGSummary") << "[DTDPGCreateWheelSummary]:  Error fitting Mean Efficiency Layers";
  }
  MeanEffLayers->Draw("same");
  createGifFile("MeanEffLayer",ct202,true);


  TCanvas * ct203=new TCanvas("ct203","Efficiencies per SuperLayer",500,0,900,950);
  ct203->cd()->SetTopMargin(0.10);
  ct203->cd()->SetRightMargin(0.12);
  ct203->cd()->SetFillColor(0);
  ct203->cd()->SetBorderMode(0);
  ct203->cd()->SetLineColor(4);

  ct203->cd();

  MeanEffSuperLayers->SetStats(1);
  MeanEffSuperLayers->SetLineColor(2);
  MeanEffSuperLayers->SetFillColor(5);
  MeanEffSuperLayers->Draw();
  try { MeanEffSuperLayers->Fit("gaus","Q"); } 
  catch(const cms::Exception&) {
     edm::LogError("DTDPGSummary") << "[DTDPGCreateWheelSummary]:  Error fitting Mean Efficiency SuperLayers";
  }
  MeanEffSuperLayers->Draw("same");
  createGifFile("MeanEffSuperLayer",ct203,true);


  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //ROB Errors
  TCanvas * ct101[5];
  ct101[0]=new TCanvas("ct101_1","ROS Errors W-2",500,0,1150,950);
  ct101[1]=new TCanvas("ct101_2","ROS Errors W-1",500,0,1150,950);
  ct101[2]=new TCanvas("ct101_3","ROS Errors W0", 500,0,1150,950);
  ct101[3]=new TCanvas("ct101_4","ROS Errors W1", 500,0,1150,950);
  ct101[4]=new TCanvas("ct101_5","ROS Errors W2", 500,0,1150,950);

  //FRC divide stations
  TLine* lMB1 = new TLine (0., 6.,17., 6.);
  TLine* lMB2 = new TLine (0.,12.,17.,12.);
  TLine* lMB3 = new TLine (0.,18.,17.,18.);
  TLine* lMB4 = new TLine (0.,24.,17.,24.);
  lMB1->SetLineWidth(2.);
  lMB2->SetLineWidth(2.);
  lMB3->SetLineWidth(2.);
  lMB4->SetLineWidth(2.);
  TText* tMB1 = new  TText ( 7., 2.,"MB1");
  TText* tMB2 = new  TText ( 8., 8.,"MB2");
  TText* tMB3 = new  TText ( 9.,14.,"MB3");
  TText* tMB4 = new  TText (10.,20.,"MB4");
  tMB1->SetTextSize(0.06);
  tMB2->SetTextSize(0.06);
  tMB3->SetTextSize(0.06);
  tMB4->SetTextSize(0.06);


  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
  {
   ct101[iw+2]->cd()->SetTopMargin(0.36);
   ct101[iw+2]->cd()->SetLeftMargin(0.05);
   ct101[iw+2]->cd()->SetRightMargin(0.05);
   ct101[iw+2]->cd()->SetBottomMargin(0.20);
   ct101[iw+2]->cd()->Divide(3,4,0.0001,0.0001);

   gStyle->SetOptTitle(1);
   gStyle->SetPalette(1);
  
   gStyle->SetTitleH(0.1);
   gStyle->SetTitleW(0.97);
   gStyle->SetTitleH(0.1);

   gStyle->SetTitleBorderSize(2);

   for(int isec=1;isec<13;isec++) // Loop on sectors 
   {
    ct101[iw+2]->cd(isec);
    ct101[iw+2]->cd(isec)->SetBottomMargin(0.15);
    ct101[iw+2]->cd(isec)->SetGrid();
 
    stringstream dduID;
    //if(isec<7)dduID << 770+2+iw;  // for 10 DDUs
    //else dduID <<  775+2+iw;      // for 10 DDUs
    dduID << 770+2+iw; // for 5 DDUs
    stringstream ROSID; ROSID << isec;

    string histoName = myMainFolder + "00-DataIntegrity/FED" + dduID.str() + "/ROS" + ROSID.str() + "/FED" + dduID.str() + "_ROS" + ROSID.str() + "_ROSError";
    TH1F *histoROSError = (TH1F*) myFile -> Get(histoName.c_str());
    if(histoROSError) {
     histoROSError->SetStats( 0 );
     histoROSError->GetXaxis()->SetBinLabel(1,"#splitline{Link}{TimeOut}");
     histoROSError->GetXaxis()->SetBinLabel(2,"#splitline{Ev.Id.}{Mis.}");
     histoROSError->GetXaxis()->SetBinLabel(3,"#splitline{FIFO}{#splitline{almost}{full}}");
     histoROSError->GetXaxis()->SetBinLabel(4,"#splitline{FIFO}{full}");
     histoROSError->GetXaxis()->SetBinLabel(5,"#splitline{CEROS}{timeout}");
     histoROSError->GetXaxis()->SetBinLabel(6,"#splitline{Max.}{wds}");
     histoROSError->GetXaxis()->SetBinLabel(7,"#splitline{TDC}{#splitline{parity}{err}}");
     histoROSError->GetXaxis()->SetBinLabel(8,"#splitline{BX ID}{Mis.}");
     histoROSError->GetXaxis()->SetBinLabel(9,"#splitline{Ch}{blocked}");
     histoROSError->GetXaxis()->SetBinLabel(10,"#splitline{Ev ID}{Mis.}");
     histoROSError->GetXaxis()->SetBinLabel(11,"#splitline{CEROS}{blocked}");
     histoROSError->GetXaxis()->SetBinLabel(12,"#splitline{TDC}{Fatal}");
     histoROSError->GetXaxis()->SetBinLabel(13,"#splitline{TDC}{#splitline{FIFO}{Ov.}}");
     histoROSError->GetXaxis()->SetBinLabel(14,"#splitline{L1}{#splitline{Buffer}{Ov.}}");
     histoROSError->GetXaxis()->SetBinLabel(15,"#splitline{TDCL1A}{#splitline{FIFO}{Ov.}}");
     histoROSError->GetXaxis()->SetBinLabel(16,"#splitline{TDC}{#splitline{hit}{err.}}");
     histoROSError->GetXaxis()->SetBinLabel(17,"#splitline{TDC}{#splitline{hit}{rej.}}");
     histoROSError->GetYaxis()->SetLabelSize(0.045);
     histoROSError->GetXaxis()->SetLabelSize(0.04);
     histoROSError->LabelsOption("h","X");
     histoROSError->Draw("colz");
     lMB1->Draw();
     lMB2->Draw();
     lMB3->Draw();
     lMB4->Draw();
     tMB1->Draw();
     tMB2->Draw();
     tMB3->Draw();
     tMB4->Draw();
    }
   }// end loop on sectors

   createGifFile("DataIntegrity_ROBErrors",ct101[iw+2],iw);
  }// end loop on wheels


  TCanvas * ct102[5];
  ct102[0]=new TCanvas("ct102_1","ROS EventLenght W-2",500,0,1150,950);
  ct102[1]=new TCanvas("ct102_2","ROS EventLenght W-1",500,0,1150,950);
  ct102[2]=new TCanvas("ct102_3","ROS EventLenght W0", 500,0,1150,950);
  ct102[3]=new TCanvas("ct102_4","ROS EventLenght W1", 500,0,1150,950);
  ct102[4]=new TCanvas("ct102_5","ROS EventLenght W2", 500,0,1150,950);

  for(int iw=-2;iw<2+1;iw++) // Loop on wheels
  {
   ct102[iw+2]->cd()->SetTopMargin(0.36);
   ct102[iw+2]->cd()->SetLeftMargin(0.05);
   ct102[iw+2]->cd()->SetRightMargin(0.05);
   ct102[iw+2]->cd()->SetBottomMargin(0.20);
   ct102[iw+2]->cd()->Divide(3,4,0.0001,0.0001);

   for(int isec=1;isec<13;isec++) // Loop on sectors 
   {
    ct102[iw+2]->cd(isec);
    ct102[iw+2]->cd(isec)->SetLogy(1);

    gStyle->SetOptStat(111111);
    gStyle->SetStatX(0.88);
    gStyle->SetStatW(0.32);
    //gStyle->SetStatY(0.98);
    gStyle->SetStatY(0.88);
    gStyle->SetStatH(0.40);

    stringstream dduID;
    //if(isec<7)dduID << 770+2+iw; // for 10 DDUs
    //else dduID <<  775+2+iw;     // for 10 DDUs
    dduID << 770+2+iw;   // for 5 DDUs
    stringstream ROSID; ROSID << isec;

    string histoName = myMainFolder + "00-DataIntegrity/FED" + dduID.str() + "/ROS" + ROSID.str() + "/FED" + dduID.str() + "_ROS" + ROSID.str() + "_ROSEventLenght";

    TH1F *histoROSEL = (TH1F*) myFile -> Get(histoName.c_str());
    if(histoROSEL) 
    {
     histoROSEL->GetYaxis()->SetLabelSize(0.07);
     histoROSEL->GetXaxis()->SetLabelSize(0.07);
     histoROSEL->Draw();
    }

   }// end loop on sectors

   createGifFile("ROBEventLength",ct102[iw+2],iw);
   delete ct102[iw+2];
  }// end loop on wheels



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // >>>>>>>>>>>>>>>>   TRIGGER PLOTS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary histos with trigger efficiency (All/HHHL),  fracion of Correlated/H triggers,
  // LUTs residuals mean and RMS, BX distribution mean and TwinMux BX and quality distribution
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // CREATE HISTOS...................................

  char thehtit[240]; char thehistoname[240];

  TH2F * EffTrigPhi[2];
  TH1F * EffTrigPhi1d[2];
  TH1F * EffTrigPhi1dMB4[2];
  TH2F * EffTrigPhiHHHL[2];
  TH1F * EffTrigPhiHHHL1d[2];


  for (int iHw=0;iHw<iHwMax;++iHw) {
    //string hwSource       = iHw ? "DDU" : "TM";
    string hwSource;
    if(iHw==0) hwSource   = "TM";
    if(iHw==1) hwSource   =  "DDU";
    string hTitle         = "TrigPhiEff_" + hwSource;
    string hName          = "Trigger Efficiency - Phi (" + hwSource + " any Quality & Slot)";
    EffTrigPhi[iHw]       = new TH2F(hTitle.c_str(),hName.c_str(),14, 1.,15.,26,0.,26.);
    hTitle                = "TrigPhiEff1d_" + hwSource;
    EffTrigPhi1d[iHw]     = new TH1F(hTitle.c_str(),hName.c_str(),100,0,100);
    hTitle                = "TrigPhiEff1dMB4_" + hwSource;
    EffTrigPhi1dMB4[iHw]  = new TH1F(hTitle.c_str(),hName.c_str(),100,0,100);
    hTitle                = "TrigPhiEffHHHL_" + hwSource;
    hName                 = "Trigger Efficiency - Phi (" + hwSource + " HH/HL)";
    EffTrigPhiHHHL[iHw]   = new TH2F(hTitle.c_str(),hName.c_str(),14, 1.,15.,26,0.,26.);
    hTitle                = "TrigPhiEff1dHHHL_" + hwSource;
    EffTrigPhiHHHL1d[iHw] = new TH1F(hTitle.c_str(),hName.c_str(),100,0,100);
  }
  // ***************
  // NOW USING TM INFO (Before <2016  DDU)
  // ***************
     sprintf(thehtit,"CorrFractionPhi");
     sprintf(thehistoname,"Fraction of Correlated Triggers - Phi");
     TH2F * CorrFracPhi= new TH2F(thehtit,thehistoname,14, 1.,15.,26,0.,26.);  // Now using TM info
     sprintf(thehtit,"HFractionTheta");
     sprintf(thehistoname,"Fraction of H Triggers - Theta");
     TH2F * HFracTheta= new TH2F(thehtit,thehistoname,14, 1.,15.,26,0.,26.);
  // ************
   
  TH1F * BXMeanPhi[5][4];
  TH1F * LutMean[5][4];
  TH1F * LutRMS[5][4];
  TH1F * LutDirMean[5][4];
  TH1F * LutDirRMS[5][4];
  TH1D * BX1d[5][4][12];
  for(int iw=0;iw<5;iw++) 
    for(int ic=0;ic<4;ic++) 
      for(int is=0;is<12;is++) BX1d[iw][ic][is]=NULL;

  for (int iSt=1;iSt<=4;++iSt) {
    for (int iWh=-2;iWh<=2;++iWh) {
      sprintf(thehtit,"BXMeanPhi%d_%d",iWh,iSt);
      sprintf(thehistoname,"BX Average Wheel %d Station %d",iWh,iSt);
      BXMeanPhi[iWh+2][iSt-1] = new TH1F(thehtit,thehistoname,14, 1.,15.);
      sprintf(thehtit,"LutPosMean%d_%d",iWh,iSt);
      sprintf(thehistoname,"Position LUTs Residuam Mean Wheel %d Station %d (cm)",iWh,iSt);
      LutMean[iWh+2][iSt-1] = new TH1F(thehtit,thehistoname,14, 1.,15.);
      sprintf(thehtit,"LutPosRMS%d_%d",iWh,iSt);
      sprintf(thehistoname,"Position LUTs Residual RMS Wheel  %d Station %d",iWh,iSt);
      LutRMS[iWh+2][iSt-1] = new TH1F(thehtit,thehistoname,14, 1.,15.);
      sprintf(thehtit,"LutDirMean%d_%d",iWh,iSt);
      sprintf(thehistoname,"Direction LUTs Residuam Mean Wheel %d Station %d (cm)",iWh,iSt);
      LutDirMean[iWh+2][iSt-1] = new TH1F(thehtit,thehistoname,14, 1.,15.);
      sprintf(thehtit,"LutDirRMS%d_%d",iWh,iSt);
      sprintf(thehistoname,"Direction LUTs Residual RMS Wheel  %d Station %d",iWh,iSt);
      LutDirRMS[iWh+2][iSt-1] = new TH1F(thehtit,thehistoname,14, 1.,15.);
    }
  }

  // Put Labels
  for(int isec=1;isec<15;isec++)
    {
      char  EWlabel[40];
      sprintf(EWlabel,"Sect-%d",isec);
      //for (int i=0;i<2;++i) {
      for (int i=0;i<iHwMax;++i) {
      EffTrigPhi[i]->GetXaxis()->SetBinLabel(isec,EWlabel);
      EffTrigPhiHHHL[i]->GetXaxis()->SetBinLabel(isec,EWlabel);
      }
      //if(ProcessDDUTrigger)
      {// Now using TM
        CorrFracPhi->GetXaxis()->SetBinLabel(isec,EWlabel); 
        HFracTheta->GetXaxis()->SetBinLabel(isec,EWlabel);
      }
      BXMeanPhi[4][0]->GetXaxis()->SetBinLabel(isec,EWlabel);
      LutMean[4][0]->GetXaxis()->SetBinLabel(isec,EWlabel);
      LutRMS[4][0]->GetXaxis()->SetBinLabel(isec,EWlabel);
      LutDirMean[4][0]->GetXaxis()->SetBinLabel(isec,EWlabel);
      LutDirRMS[4][0]->GetXaxis()->SetBinLabel(isec,EWlabel);
    }

  for(int iw=0;iw<5;iw++)
    for(int imb=1;imb<5;imb++)
      {
    char  EWlabel[40];
    sprintf(EWlabel,"Wh%d MB%d",iw-2,imb);
      //for (int i=0;i<2;++i) {
      for (int i=0;i<iHwMax;++i) {
        EffTrigPhi[i]->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel);
        EffTrigPhiHHHL[i]->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel);
      }
      if(ProcessDDUTrigger)
      {
         CorrFracPhi->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel); 
         HFracTheta->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel);
      }
      else
      { // Now using TM
         CorrFracPhi->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel);
         HFracTheta->GetYaxis()->SetBinLabel((4-iw)*5+imb+1,EWlabel);
      }
      }

  if(ProcessDDUTrigger)
  {
   for (int wh=-2;wh<=2;wh++) {
    stringstream wheel; wheel << wh;
    string recoFolder       = myMainFolder + "04-LocalTrigger-DDU/Wheel" + wheel.str();
    string histoNameCorr    = recoFolder + "/DDU_CorrFractionPhi_W" + wheel.str();
    string histoNameH       = recoFolder + "/DDU_HFractionTheta_W" + wheel.str();

    TH2F *HCorr  = (TH2F*) myFile -> Get(histoNameCorr.c_str());
    if(!HCorr)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameCorr.c_str() << " not found" << endl;
    else {
      for(int stat=1;stat<=4;++stat) {
       for(int sec=1;sec<=12;sec++) {
         CorrFracPhi->SetBinContent(sec,(2-wh)*5+stat+1,HCorr->GetBinContent(sec,stat)*100);
       }
      }
    }
    TH2F * HH  = (TH2F*) myFile -> Get(histoNameH.c_str());
    if(!HH)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameH.c_str() << " not found" << endl;
    else {
      for(int stat=1;stat<=4;++stat) {
    for(int sec=1;sec<=12;sec++) {
     HFracTheta->SetBinContent(sec,(2-wh)*5+stat+1,HH->GetBinContent(sec,stat)*100);
    }
      }
    }
   }
  }
  else // using TM 
  {
   for (int wh=-2;wh<=2;wh++) {
    stringstream wheel; wheel << wh;
    string recoFolder       = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str();
    string histoNameCorr    = recoFolder + "/TM_CorrFractionPhi_W" + wheel.str();
    string histoNameH       = recoFolder + "/TM_HFractionTheta_W" + wheel.str();

    TH2F *HCorr  = (TH2F*) myFile -> Get(histoNameCorr.c_str());
    if(!HCorr)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameCorr.c_str() << " not found" << endl;
    else {
      for(int stat=1;stat<=4;++stat) {
    for(int sec=1;sec<=12;sec++) {
     CorrFracPhi->SetBinContent(sec,(2-wh)*5+stat+1,HCorr->GetBinContent(sec,stat)*100);
    }
      }
    }
    TH2F * HH  = (TH2F*) myFile -> Get(histoNameH.c_str());
    if(!HH)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameH.c_str() << " not found" << endl;
    else {
      for(int stat=1;stat<=4;++stat) {
        for(int sec=1;sec<=12;sec++) {
            HFracTheta->SetBinContent(sec,(2-wh)*5+stat+1,HH->GetBinContent(sec,stat)*100);
        }
      }
    }
   }
  }


  chIt = myMuonGeom->chambers().begin();
  chEnd = myMuonGeom->chambers().end();
  for (; chIt != chEnd; ++chIt) {
    DTChamberId ch = (*chIt)->id();
    stringstream wheel; wheel << ch.wheel();
    stringstream station; station << ch.station();
    stringstream sector; sector << ch.sector();
    int wh  = ch.wheel();
    int sec= ch.sector();
    int stat = ch.station();
    string histoNameLUTs = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/Segment/TM_PhiResidual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
    TH1F * HResidualLUTs  = (TH1F*) myFile -> Get(histoNameLUTs.c_str());
    if(!HResidualLUTs)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameLUTs.c_str() << " not found" << endl;
    else {
      float mean = HResidualLUTs->GetMean();
      mean = fabs(mean)<14. ? mean : 14.*(mean)/fabs(mean);
      float meanErr = fabs(mean)<14. ? HResidualLUTs->GetMeanError() : 0;
      float rms  = HResidualLUTs->GetRMS();
      rms = rms<1.95 ? rms : 1.95 ;
      float rmsErr  = rms<1.95 ? HResidualLUTs->GetRMSError() : 0;
      LutMean[wh+2][stat-1]->SetBinContent(sec,mean);
      LutMean[wh+2][stat-1]->SetBinError(sec,meanErr);
      LutRMS[wh+2][stat-1]->SetBinContent(sec,rms);
      LutRMS[wh+2][stat-1]->SetBinError(sec,rmsErr);
    }
    string histoNameDirLUTs = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/Segment/TM_PhibResidual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
    TH1F * HResidualDirLUTs  = (TH1F*) myFile -> Get(histoNameDirLUTs.c_str());
    if(!HResidualDirLUTs)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameDirLUTs.c_str() << " not found" << endl;
    else {
      float mean = HResidualDirLUTs->GetMean();
      mean = fabs(mean)<1.95 ? mean : 1.95*(mean)/fabs(mean);
      float meanErr = fabs(mean)<2. ? HResidualDirLUTs->GetMeanError() : 0;
      float rms  = HResidualDirLUTs->GetRMS();
      rms = rms<3.95 ? rms : 3.95 ;
      float rmsErr  = rms<3.95 ? HResidualDirLUTs->GetRMSError() : 0;
      LutDirMean[wh+2][stat-1]->SetBinContent(sec,mean);
      LutDirMean[wh+2][stat-1]->SetBinError(sec,meanErr);
      LutDirRMS[wh+2][stat-1]->SetBinContent(sec,rms);
      LutDirRMS[wh+2][stat-1]->SetBinError(sec,rmsErr);
    }
    if (sec>12) continue;
    //string histoNameBXvsQual = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhi/TM_BXvsQual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
    //Changed to move to TM IN (06/10/2016)
    string histoNameBXvsQual = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhiIn/TM_BXvsQual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
    TH2F * HBXvsQual  = (TH2F*) myFile -> Get(histoNameBXvsQual.c_str());
    if(!HBXvsQual)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameBXvsQual.c_str() << " not found" << endl;
    else {
      stringstream hBXId;
      hBXId << "BX1d" << wh << sec << stat;
      BX1d[wh+2][stat-1][sec-1] = HBXvsQual->ProjectionY(hBXId.str().c_str(),1,7);
      BXMeanPhi[wh+2][stat-1]->SetBinContent(sec,BX1d[wh+2][stat-1][sec-1]->GetMean());
      BXMeanPhi[wh+2][stat-1]->SetBinError(sec,BX1d[wh+2][stat-1][sec-1]->GetMeanError());
    }
  }

  for (int wh=-2;wh<=2;wh++) {
    for (int iHw=0;iHw<iHwMax;++iHw) {
      //string hwSource = iHw ? "DDU" : "TM";
      string hwSource;
      if(iHw==0)  hwSource ="TM";
      if(iHw==1)  hwSource ="DDU";
      //string dqmTrigBase = iHw ? "04-LocalTrigger-DDU/" : "03-LocalTrigger-TM/";
      string dqmTrigBase ; 
      if(iHw==0) dqmTrigBase = "03-LocalTrigger-TM/";
      if(iHw==1) dqmTrigBase = "04-LocalTrigger-DDU/";
      stringstream wheel; wheel << wh;
      string recoFolder       = myMainFolder + dqmTrigBase + "Wheel" + wheel.str();
      string histoNameEff     = recoFolder + "/" + hwSource + "_TrigEffPhi_W" + wheel.str();
      string histoNameEffHHHL = recoFolder + "/" + hwSource + "_TrigEffHHHLPhi_W" + wheel.str();

      TH2F * HEff  = (TH2F*) myFile -> Get(histoNameEff.c_str());
      if(!HEff)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameEff.c_str() << " not found" << endl;
      else {
      for(int stat=1;stat<=4;++stat) {
        for(int sec=1;sec<=12;sec++) {
         EffTrigPhi[iHw]->SetBinContent(sec,(2-wh)*5+stat+1,HEff->GetBinContent(sec,stat)*100);
      if (stat<4) {
        EffTrigPhi1d[iHw]->Fill(HEff->GetBinContent(sec,stat)*100);
      } else {
          EffTrigPhi1dMB4[iHw]->Fill(HEff->GetBinContent(sec,stat)*100);
      }
        }
      }
      }
      TH2F * HEffHHHL  = (TH2F*) myFile -> Get(histoNameEffHHHL.c_str());
      if(!HEffHHHL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameEffHHHL.c_str() << " not found" << endl;
      else {
      for(int stat=1;stat<=4;++stat) {
        for(int sec=1;sec<=12;sec++) {
         EffTrigPhiHHHL[iHw]->SetBinContent(sec,(2-wh)*5+stat+1,HEffHHHL->GetBinContent(sec,stat)*100);
         EffTrigPhiHHHL1d[iHw]->Fill(HEffHHHL->GetBinContent(sec,stat)*100);
        }
      }
      }
    }
  }


  // DRAW   HISTOS...................................
  // color palette to have green when efficiency=1 
  gStyle->SetPalette(20,col);

  TCanvas * ct5[28];
  ct5[0]=new TCanvas("ct5_1","Phi Trigg Eff TM",500,0,900,950);
  ct5[1]=new TCanvas("ct5_2","Phi Trigg HHHL Eff TM",500,0,900,950);
  if(ProcessDDUTrigger)
  {
     ct5[2]=new TCanvas("ct5_3","Phi Trigg Eff DDU",500,0,900,950);
     ct5[3]=new TCanvas("ct5_4","Phi Trigg HHHL Eff DDU",500,0,900,950);
     ct5[4]=new TCanvas("ct5_5","Phi Trigg Corr Fraction",500,0,900,950);
     ct5[5]=new TCanvas("ct5_6","Theta Trigg H Fraction",500,0,900,950);
  }
  else
  { // Ussing TM
     ct5[4]=new TCanvas("ct5_5","Phi Trigg Corr Fraction",500,0,900,950);
     ct5[5]=new TCanvas("ct5_6","Theta Trigg H Fraction",500,0,900,950);
  }
  ct5[6]=new TCanvas("ct5_7","BX Average (AnyQual)",500,0,900,950);
  ct5[7]=new TCanvas("ct5_8","LUTs Pos Mean",500,0,900,950);
  ct5[8]=new TCanvas("ct5_9","LUTs Pos RMS",500,0,900,950);
  ct5[9]=new TCanvas("ct5_10","Phi Trigg Eff 1d TM",500,0,900,950);
  if(ProcessDDUTrigger)
       ct5[10]=new TCanvas("ct5_11","Phi Trigg Eff 1d DDU",500,0,900,950);

  for (int wh=-2;wh<=2;++wh) {
    stringstream canvasNameBX, canvasNameQualTM, canvasNameQualDDU, titleBX, titleQualTM, titleQualDDU;
    canvasNameBX << "BX Distribution (TM data) for Wheel " << wh;
    canvasNameQualTM << "Quality Distribution (TM data) for Wheel " << wh;
    canvasNameQualDDU << "Quality Distribution (DDU data) for Wheel " << wh;
    titleBX << "ct5_" << wh+14;
    titleQualTM << "ct_5" << wh+19;
    titleQualDDU << "ct_5" << wh+24;
    ct5[wh+13] = new TCanvas(titleBX.str().c_str(),canvasNameBX.str().c_str(),500,0,900,950);
    ct5[wh+18] = new TCanvas(titleQualTM.str().c_str(),canvasNameQualTM.str().c_str(),500,0,900,950);
    if(ProcessDDUTrigger)
        ct5[wh+23] = new TCanvas(titleQualDDU.str().c_str(),canvasNameQualDDU.str().c_str(),500,0,900,950);
  }
  ct5[26]=new TCanvas("ct5_27","LUTs Dir Mean",500,0,900,950);
  ct5[27]=new TCanvas("ct5_28","LUTs Dir RMS",500,0,900,950);
  

  TPaveLabel* wheelName[5];
  TPaveLabel* canvasTitle[5];
  memset(wheelName,0,5*sizeof(TPaveLabel*));
  memset(canvasTitle,0,3*sizeof(TPaveLabel*));

  TH2F * hTrig[6];
  hTrig[0]=EffTrigPhi[0];
  hTrig[1]=EffTrigPhiHHHL[0];
  if(ProcessDDUTrigger)
  {
    hTrig[2]=EffTrigPhi[1];
    hTrig[3]=EffTrigPhiHHHL[1];
    hTrig[4]=CorrFracPhi; 
    hTrig[5]=HFracTheta;
  }
  else
  {// Now Calculated from TM 
    hTrig[4]=CorrFracPhi; 
    hTrig[5]=HFracTheta;
  }
  int ihMax=5;
  //int ihMax=1;
  //if(ProcessDDUTrigger)ihMax=5;  // histograms 2 to 5 corresponds to DDU
  for(int ih=0;ih<=ihMax;ih++) {
    if(ProcessDDUTrigger || ih<2 || ih>3)   // histograms 2,3 corresponds to DDU
    {
      ct5[ih]->cd()->SetTopMargin(0.10);
      ct5[ih]->cd()->SetRightMargin(0.12);
      ct5[ih]->cd()->SetFillColor(0);
      ct5[ih]->cd()->SetBorderMode(0);
      ct5[ih]->cd()->SetLineColor(4);
      ct5[ih]->cd();
      hTrig[ih]->GetXaxis()->SetLabelSize(0.03);
      hTrig[ih]->GetYaxis()->SetLabelSize(0.03);
      hTrig[ih]->SetMaximum(100);
      hTrig[ih]->SetMinimum(50);
      if(ih==4) hTrig[ih]->SetMinimum(35);
      if(ih==5) hTrig[ih]->SetMinimum(20);
      if(ih==5) hTrig[ih]->SetMinimum(20);
      hTrig[ih]->SetTitleSize(0.025);
      hTrig[ih]->Draw("colz");
    }// End IF
  }// End loop

  for (int ih=6;ih<=8;++ih) {
    ct5[ih]->SetFillColor(0);
    ct5[ih]->SetTopMargin(0.45);
    ct5[ih]->SetBottomMargin(0.25);
    ct5[ih]->SetRightMargin(0.12);
    ct5[ih]->Divide(1,5,0.,0);
    for (int iw=1;iw<6;++iw) {
      ct5[ih]->cd(iw)->SetFillColor(0);
      ct5[ih]->cd(iw)->SetGrid();
      ct5[ih]->cd(iw)->SetBorderMode(0);
    }
  }

  for (int ih=26;ih<=27;++ih) {
    ct5[ih]->SetFillColor(0);
    ct5[ih]->SetTopMargin(0.45);
    ct5[ih]->SetBottomMargin(0.25);
    ct5[ih]->SetRightMargin(0.12);
    ct5[ih]->Divide(1,5,0.,0);
    for (int iw=1;iw<6;++iw) {
      ct5[ih]->cd(iw)->SetFillColor(0);
      ct5[ih]->cd(iw)->SetGrid();
      ct5[ih]->cd(iw)->SetBorderMode(0);
    }
  }

  for (int iHw=0;iHw<iHwMax;++iHw) {
    TH1F *TrigEff1d[2];
    TrigEff1d[0]=EffTrigPhi1d[iHw];
    TrigEff1d[1]=EffTrigPhiHHHL1d[iHw];
    ct5[9+iHw]->SetFillColor(0);
    ct5[9+iHw]->Divide(1,2);
    for (int ip=0;ip<2;++ip) {
      ct5[9+iHw]->cd(ip+1)->SetFillColor(0);
      ct5[9+iHw]->cd(ip+1)->SetGrid();
      ct5[9+iHw]->cd(ip+1)->SetBorderMode(0);
      ct5[9+iHw]->cd(ip+1);
      TrigEff1d[ip]->GetXaxis()->SetTitle("efficiency");
      TrigEff1d[ip]->SetTitleSize(0.025);
      TrigEff1d[ip]->GetXaxis()->SetTitleSize(0.035);
      TrigEff1d[ip]->SetStats(1);
      float mean = TrigEff1d[ip]->GetBinCenter(TrigEff1d[ip]->GetMaximumBin());
      float min = mean>7.  ? mean-7 : 0;
      float max = mean<93. ? mean+7 : 100;
      try { TrigEff1d[ip]->Fit("gaus","Q","",min,max); } 
      catch(const cms::Exception&) {
	edm::LogError("DTDPGSummary") << "[DTDPGCreateWheelSummary]:  Error fitting Mean Trigger Efficiency";
      }
      TrigEff1d[ip]->Draw();
      if (!ip) {
	EffTrigPhi1dMB4[iHw]->SetLineColor(4);
	//EffTrigPhi1dMB4[iHw]->SetFillColor(38);
	EffTrigPhi1dMB4[iHw]->Draw("same");
      }
    }
  }

  TH1F * hTrig1d[5][5][4];
  memcpy(hTrig1d[0],BXMeanPhi,20*sizeof(TH1F*));
  memcpy(hTrig1d[1],LutMean,20*sizeof(TH1F*));
  memcpy(hTrig1d[2],LutRMS,20*sizeof(TH1F*));
  memcpy(hTrig1d[3],LutDirMean,20*sizeof(TH1F*));
  memcpy(hTrig1d[4],LutDirRMS,20*sizeof(TH1F*));
  //int rmin[5]     = {-2,-15,0,-2,0};
  //int rmax[5]     = {2,15,2,2,4};
  int rmin[5]     = {-4,-15,0,-2,0};
  int rmax[5]     = {3,15,2,2,4};
  char cTitle[5][240] = {"BX Average (AnyQual)", "Trigger Pos Residual Mean (Deg)","Trigger Pos residual RMS","Trigger Dir Residual Mean (Deg)","Trigger Dir residual RMS"};
  for(int ih=0;ih<5;ih++) {
    int ican = (ih<=2) ? ih+6 : ih+23;
    for (int iw=0;iw<5;++iw) {
      stringstream ptitle; ptitle << "Wheel " << (iw-2);
      wheelName[iw] = new TPaveLabel(1.005,0.35,1.14,0.65,ptitle.str().c_str(),"NDC");
      wheelName[iw]->SetBorderSize(0);
      wheelName[iw]->SetFillColor(0);
      ct5[ican]->cd(iw+1);
      ct5[ican]->cd(iw+1)->SetGrid();
      for (int ist=0;ist<4;++ist) {
	if (ih>2 && ist==2) continue;
	string opt = ist ? "same" : "";
	hTrig1d[ih][iw][ist]->SetLineWidth(2);
	hTrig1d[ih][iw][ist]->SetLineColor(ist+1);
	hTrig1d[ih][iw][ist]->GetYaxis()->SetRangeUser(rmin[ih],rmax[ih]);
	hTrig1d[ih][iw][ist]->GetYaxis()->SetLabelSize(0.12);
	hTrig1d[ih][iw][ist]->GetYaxis()->SetNdivisions(504);
	hTrig1d[ih][iw][ist]->GetXaxis()->SetTitle("sector");
	hTrig1d[ih][iw][ist]->GetXaxis()->SetTitleSize(0.12);
	hTrig1d[ih][iw][ist]->GetXaxis()->SetLabelSize(0.12);
	hTrig1d[ih][iw][ist]->GetXaxis()->CenterLabels();
	hTrig1d[ih][iw][ist]->GetXaxis()->SetNdivisions(14);
	hTrig1d[ih][iw][ist]->Draw(opt.c_str());
	if (!ist) wheelName[iw]->Draw();
      }
    }
    canvasTitle[ih] = new TPaveLabel(0.1,0.93,0.85,0.99,cTitle[ih]);
    canvasTitle[ih]->SetFillColor(0);
    canvasTitle[ih]->SetTextColor(4);
    ct5[ican]->cd();
    canvasTitle[ih]->Draw();
  }
  for (int wh=-2;wh<=2;++wh) {

    TCanvas* canvasBX   = ct5[wh+13];
    TCanvas* canvasQual[2];
    canvasQual[0] = ct5[wh+18];
    if(ProcessDDUTrigger)
       canvasQual[1] = ct5[wh+23];
    canvasBX->SetFillColor(0);
    canvasBX->Divide(4,3);
    int iimax=1;
    if(ProcessDDUTrigger)iimax=2;
    for (int i=0;i<iimax;++i) {
      canvasQual[i]->SetFillColor(0);
      canvasQual[i]->Divide(4,3);
    }
    for (int sec=1;sec<=12;++sec) {
      for (int stat=1;stat<=4;++stat) {
      stringstream wheel; wheel << wh;
      stringstream station ; station << stat;
      stringstream sector; sector << sec;
      if (BX1d[wh+2][stat-1][sec-1]) {
        canvasBX->cd(sec);
        canvasBX->cd(sec)->SetGridy();
        BX1d[wh+2][stat-1][sec-1]->SetFillColor(0);
        BX1d[wh+2][stat-1][sec-1]->SetLineWidth(2);
        BX1d[wh+2][stat-1][sec-1]->SetLineColor(stat);
        BX1d[wh+2][stat-1][sec-1]->GetXaxis()->SetTitle("BX");
        BX1d[wh+2][stat-1][sec-1]->GetXaxis()->SetTitleSize(0.08);
        BX1d[wh+2][stat-1][sec-1]->GetXaxis()->SetLabelSize(0.07);
        BX1d[wh+2][stat-1][sec-1]->GetYaxis()->SetLabelSize(0.08);
        if (BX1d[wh+2][stat-1][sec-1]->Integral()) {
         BX1d[wh+2][stat-1][sec-1]->Scale(1./BX1d[wh+2][stat-1][sec-1]->Integral());
        }
        if (stat==1) {
          string title = "#splitline{BX Distribution (TM)}{for Wh " + wheel.str() + " Sec " + sector.str() + "}";
          BX1d[wh+2][stat-1][sec-1]->SetTitle(title.c_str());
          BX1d[wh+2][stat-1][sec-1]->GetYaxis()->SetRangeUser(0,1.1);
          BX1d[wh+2][stat-1][sec-1]->Draw();
        }
        else {
          BX1d[wh+2][stat-1][sec-1]->Draw("same");
        }
      }
      for (int iHw=0;iHw<iHwMax;++iHw) {
        //string hwSource = iHw ? "DDU" : "TM";
        string hwSource;
        if(iHw==0) hwSource = "TM";
        if(iHw==1) hwSource = "DDU";
        string dqmTrigBase = iHw ? "04-LocalTrigger-DDU/" : "03-LocalTrigger-TM/";
    //Changed to move to TM IN (06/10/2016)
        //string histoNameBestQual = myMainFolder + dqmTrigBase + "Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhi/" + hwSource + "_BestQual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
        string histoNameBestQual = myMainFolder + dqmTrigBase + "Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhiIn/" + hwSource + "_BestQual_W" + wheel.str() + "_Sec" +sector.str() + "_St" + station.str();
        TH1F * HBestQual  = (TH1F*) myFile -> Get(histoNameBestQual.c_str());
        if(!HBestQual)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameBestQual.c_str() << " not found" << endl;
        else {
          canvasQual[iHw]->cd(sec);
          canvasQual[iHw]->cd(sec)->SetGridy();
          HBestQual->SetFillColor(0);
          HBestQual->SetLineWidth(2);
          HBestQual->SetLineColor(stat);
          HBestQual->GetXaxis()->SetTitle("quality");
          HBestQual->GetXaxis()->SetTitleSize(0.08);
          HBestQual->GetXaxis()->SetLabelSize(0.10);
          HBestQual->GetYaxis()->SetLabelSize(0.10);
          if (HBestQual->Integral()) {
            HBestQual->Scale(1./HBestQual->Integral());
          }
          if (stat==1) {
            string title = " #splitline{Quality of best " + hwSource + " Trigger}{for Wh " + wheel.str() + " Sec " + sector.str() +"}";
            HBestQual->SetTitle(title.c_str());
            HBestQual->GetYaxis()->SetRangeUser(0,1.1);
            HBestQual->Draw();
          }
          else {
            HBestQual->Draw("same");
          }
        }
      }
      }
    }
  }
  gStyle->SetOptTitle(1); // These need the title!
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleBorderSize(0);
  createGifFile("TriggerEffTM",ct5[0]);
  createGifFile("TriggerEffHHHLTM",ct5[1]);
  if(ProcessDDUTrigger)
  {
    createGifFile("TriggerEffDDU",ct5[2]);
    createGifFile("TriggerEffHHHLDDU",ct5[3]);
    createGifFile("TriggerCorrFraction",ct5[4]); 
    createGifFile("TriggerHFraction",ct5[5],true);
  }
  else
  {// Now Calculated from TM 
    createGifFile("TriggerCorrFraction",ct5[4]); 
    createGifFile("TriggerHFraction",ct5[5]);
  }
    
  gStyle->SetOptFit(1);
  gStyle->SetStatX(0.35);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.90);
  gStyle->SetStatH(0.40);
  createGifFile("MeanTriggerEffTM",ct5[9],true);
  if(ProcessDDUTrigger)createGifFile("MeanTriggerEffDDU",ct5[10],true);
  gStyle->SetOptFit(0);
  gStyle->SetTitleW(0.99);
  for (int wh=-2;wh<=2;++wh) {
       createGifFile("TriggerBX",ct5[wh+13],wh);
       createGifFile("TriggerQualityTM",ct5[wh+18],wh);
       if(ProcessDDUTrigger) createGifFile("TriggerQualityDDU",ct5[wh+23],wh);
  }
  gStyle->SetOptTitle(0);
  createGifFile("TriggerBXAverage",ct5[6]);
  createGifFile("TriggerLutPosMean",ct5[7]);
  createGifFile("TriggerLutPosRMS",ct5[8],true);
  createGifFile("TriggerLutDirMean",ct5[26]);
  createGifFile("TriggerLutDirRMS",ct5[27],true);


  TCanvas * ct57b=new TCanvas("ct57b","BX Average (AnyQual) (zoom)",500,0,900,950);
  ct57b->SetFillColor(0);
  ct57b->SetTopMargin(0.45);
  ct57b->SetBottomMargin(0.25);
  ct57b->SetRightMargin(0.12);
  ct57b->Divide(1,5,0.,0);
  for (int iw=1;iw<6;++iw) {
    ct57b->cd(iw)->SetFillColor(0);
    ct57b->cd(iw)->SetGrid();
    ct57b->cd(iw)->SetBorderMode(0);
  }

  for (int iw=0;iw<5;++iw) {
    stringstream ptitle; ptitle << "Wheel " << (iw-2);
    wheelName[iw] = new TPaveLabel(1.005,0.35,1.14,0.65,ptitle.str().c_str(),"NDC");
    wheelName[iw]->SetBorderSize(0);
    wheelName[iw]->SetFillColor(0);
    ct57b->cd(iw+1);
    ct57b->cd(iw+1)->SetGrid();
    for (int ist=0;ist<4;++ist) {
    string opt = ist ? "same" : "";
    hTrig1d[0][iw][ist]->SetLineWidth(2);
    hTrig1d[0][iw][ist]->SetLineColor(ist+1);
    hTrig1d[0][iw][ist]->GetYaxis()->SetRangeUser(-0.50,0.50);
    hTrig1d[0][iw][ist]->GetYaxis()->SetLabelSize(0.12);
    hTrig1d[0][iw][ist]->GetYaxis()->SetNdivisions(504);
    hTrig1d[0][iw][ist]->GetXaxis()->SetTitle("sector");
    hTrig1d[0][iw][ist]->GetXaxis()->SetTitleSize(0.12);
    hTrig1d[0][iw][ist]->GetXaxis()->SetLabelSize(0.12);
    hTrig1d[0][iw][ist]->GetXaxis()->CenterLabels();
    hTrig1d[0][iw][ist]->GetXaxis()->SetNdivisions(14);
    hTrig1d[0][iw][ist]->Draw(opt.c_str());
    if (!ist) wheelName[iw]->Draw();
    }
  }
  TPaveLabel * canvasTitle57b = new TPaveLabel(0.1,0.93,0.85,0.99,"BX Average (AnyQual) - (Zoom)");
  canvasTitle57b->SetFillColor(0);
  canvasTitle57b->SetTextColor(4);
  ct57b->cd();
  canvasTitle57b->Draw();
  createGifFile("TriggerBXAverageZoom",ct57b);




  for (int i=0;i<5;++i) delete wheelName[i];
  for (int i=0;i<3;++i) delete canvasTitle[i];

  gStyle->SetTitleBorderSize(2);



  // >>>>>>>>>>>>>>>>   MUON TRIGGERS      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary Muon Triggers from L1/L1GMT DQM 
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  TCanvas * ct66=new TCanvas("ct66","Muon Triggers",201,81,650,550);

  gStyle->SetOptTitle();
  ct66->SetTopMargin(0.10);
  //ct66->SetLeftMargin(0.05);
  //ct66->SetRightMargin(0.05);
  //ct66->SetBottomMargin(0.20);
  //ct66->SetFillColor(0);


  gStyle->SetOptStat(11);
  gStyle->SetStatX(0.90);
  gStyle->SetStatW(0.22);
  //gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.30);

   string histoNameMuTrig = myL1TFolder + "L1TGMT/Regional_trigger";
   
   //TH1F * MuTrig  = (TH1F*) myFile -> Get(histoNameMuTrig.c_str());
   TH1F * MuTrig  = NULL; 
   MuTrig  = (TH1F*) myFile -> Get(histoNameMuTrig.c_str());
   if(MuTrig==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameMuTrig.c_str() << " not found" << endl;
   else
	{
          MuTrig->GetXaxis()->SetLabelSize(0.04);
          MuTrig->GetYaxis()->SetLabelSize(0.04);
          MuTrig->GetXaxis()->SetTitleOffset(1.50);
          MuTrig->SetFillColor(5);
          MuTrig->Draw();

          createGifFile("MuonTriggers",ct66);  // If the histo doesn't exist and this is create the program crash 
	}

  gStyle->SetOptTitle(0);


  // >>>>>>>>>>>>>>>>   EVENT LENGTH        <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary histos with the event length of each DDU
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  TCanvas * ct6=new TCanvas("ct6","DDU EventLength",500,0,900,950);
  ct6->SetTopMargin(0.36);
  ct6->SetLeftMargin(0.05);
  ct6->SetRightMargin(0.05);
  ct6->SetBottomMargin(0.20);
  ct6->SetFillColor(0);
  //ct6->Divide(2,5,0.,0.); // For 10 DDUs
  ct6->Divide(1,5,0.,0.); // For 5 DDUs


  sprintf(titlename,"Event Length");
  title = new TPaveLabel(0.1,0.945,0.9,0.98,titlename);
  title->SetFillColor(0);
  title->SetTextColor(4);
  title->Draw();

  gStyle->SetOptStat(111111);
  gStyle->SetStatX(0.88);
  gStyle->SetStatW(0.32);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.40);

  for(int iw=0;iw<5;iw++)
    {
      stringstream dduID; dduID << wheelToDDU(iw-2) ; //  
      string histoNameEL = myMainFolder + "00-DataIntegrity/FED" + dduID.str() +  "/FED" + dduID.str()  + "_EventLenght";
      TH2F * EvLength  = (TH2F*) myFile -> Get(histoNameEL.c_str());
      
      if(EvLength==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameEL.c_str() << " not found" << endl;
      else
	{
          //ct6->cd(iw*2+1);                  // for 10 DDUS
          //ct6->cd(iw*2+1)->SetLogy(1);      // for 10 DDUS
          //ct6->cd(iw*2+1)->SetFillColor(0); // for 10 DDUS
          ct6->cd(iw+1);                  // for 5 DDUS
          ct6->cd(iw+1)->SetLogy(1);      // for 5 DDUS
          ct6->cd(iw+1)->SetFillColor(0); // for 5 DDUS
          EvLength->GetXaxis()->SetLabelSize(0.08);
          EvLength->GetYaxis()->SetLabelSize(0.08);
          EvLength->SetLineColor(1);

          EvLength->Draw();

          //sprintf(titlename,"Wheel %d  S1-S6",iw-2);  // for 10 DDUs
          sprintf(titlename,"Wheel %d  ",iw-2);  // for 5 DDUs
          TPaveLabel* WheelName = new TPaveLabel(0.05,0.84,0.40,0.98,titlename,"NDC");
          WheelName->SetFillColor(0);
          WheelName->SetTextColor(4);
          WheelName->Draw();

	}

      // Uncomment next for 10 DDUs
    /*
      stringstream dduID2; dduID2 << wheelToDDU(iw-2)+5 ; //  
      string histoNameEL2 = myMainFolder + "00-DataIntegrity/FED" + dduID2.str() +  "/FED" + dduID2.str()  + "_EventLenght";
      TH2F * EvLength2  = (TH2F*) myFile -> Get(histoNameEL2.c_str());
      
      if(EvLength2==NULL)LogVerbatim("DTDPGSummary") << "[DTDPGCreateWheelSummary]: Histo = " << histoNameEL2.c_str() << " not found" << endl;
      else
	{
          ct6->cd(iw*2+2);
          ct6->cd(iw*2+2)->SetLogy(1);
          ct6->cd(iw*2+2)->SetFillColor(0);
          EvLength2->GetXaxis()->SetLabelSize(0.08);
          EvLength2->GetYaxis()->SetLabelSize(0.08);

          EvLength2->Draw();

          sprintf(titlename,"Wheel %d  S7-S12",iw-2);
          TPaveLabel* WheelName = new TPaveLabel(0.05,0.84,0.40,0.98,titlename,"NDC");
          WheelName->SetFillColor(0);
          WheelName->SetTextColor(4);
          WheelName->Draw();
	}
    */

    }
           
  createGifFile("EventLength",ct6);


  // >>>>>>>>>>>>>>>>   DATA INTEGRITY      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
  // Summary histos with the data integrity of each DDU
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  gStyle->SetOptStat(10);

  TCanvas * ct7[2];
  ct7[0]=new TCanvas("ct7_a","DDU Data Integrity",500,0,1150,950);
  ct7[1]=new TCanvas("ct7_b","DDU Data Integrity",500,0,1150,950);
  //for(int i=0;i<2;i++) // Loop for 10 DDUs, Uncomment this and createGifFile DataIntegrity_b" if 10DDU
  for(int i=0;i<1;i++) // Not loop need for 5 DDU, mantained for simplicity and in case we need to comeback to 10DDU 
  {
   ct7[i]->SetTopMargin(0.36);
   ct7[i]->SetLeftMargin(0.05);
   ct7[i]->SetRightMargin(0.05);
   ct7[i]->SetBottomMargin(0.20);
   ct7[i]->Divide(3,6,0.0001,0.0001);
  
   gStyle->SetOptTitle(1);
   gStyle->SetPalette(1);
  
   gStyle->SetTitleH(0.1);
   gStyle->SetTitleW(0.97);
   gStyle->SetTitleH(0.1);
   gStyle->SetStatX(0.99);
   gStyle->SetStatW(0.29);
   gStyle->SetStatY(0.88);
   gStyle->SetStatH(0.30);
  
  
   for(int iw=0;iw<5;iw++)
     {
       stringstream dduID; dduID << wheelToDDU(iw-2)+i*5;  
       string histoNameInt = myMainFolder + "00-DataIntegrity/FED" + dduID.str() + "/FED" + dduID.str() + "_TTSValues";
       TH2F * DDUInfo    = (TH2F*) myFile -> Get(histoNameInt.c_str());
       histoNameInt = myMainFolder + "00-DataIntegrity/FED" + dduID.str() + "/FED" + dduID.str() + "_ROSStatus";
       TH2F * RosStatus  = (TH2F*) myFile -> Get(histoNameInt.c_str());
       histoNameInt = myMainFolder + "00-DataIntegrity"+ "/FED" + dduID.str() + "_ROSSummary"; 
       TH2F * RosSummary = (TH2F*) myFile -> Get(histoNameInt.c_str());
  
       if(DDUInfo)
       {
         ct7[i]->cd(iw*3+1);
         DDUInfo->GetXaxis()->SetBinLabel(1,"disc.");
         DDUInfo->GetXaxis()->SetBinLabel(2,"overfl");
         DDUInfo->GetXaxis()->SetBinLabel(3,"out sync");
         DDUInfo->GetXaxis()->SetBinLabel(4,"busy");
         DDUInfo->GetXaxis()->SetBinLabel(5,"ready");
         DDUInfo->GetXaxis()->SetBinLabel(6,"error");
         DDUInfo->GetXaxis()->SetBinLabel(7,"disc.");
  
         DDUInfo->GetYaxis()->SetLabelSize(0.08);
         DDUInfo->GetXaxis()->SetLabelSize(0.08);
         DDUInfo->Draw();
  
         sprintf(titlename,"Wheel %d",iw-2);
         TPaveLabel* WheelName = new TPaveLabel(0.10,0.70,0.40,0.88,titlename,"NDC");
         WheelName->Draw();
  
       }
       if(RosStatus)
       {
         ct7[i]->cd(iw*3+2);
         ct7[i]->cd(iw*3+2)->SetBottomMargin(0.15);
         ct7[i]->cd(iw*3+2)->SetGrid();
  
         RosStatus->SetStats(0);
         RosStatus->GetXaxis()->SetBinLabel(1,"#splitline{ch.}{enable}");
         RosStatus->GetXaxis()->SetBinLabel(2,"#splitline{time}{out}");
         RosStatus->GetXaxis()->SetBinLabel(3,"#splitline{trailer}{lost}");
         RosStatus->GetXaxis()->SetBinLabel(4,"#splitline{fiber}{lost}");
         RosStatus->GetXaxis()->SetBinLabel(5,"#splitline{prop.}{err.}");
         RosStatus->GetXaxis()->SetBinLabel(6,"#splitline{patt.}{err.}");
         RosStatus->GetXaxis()->SetBinLabel(7,"#splitline{sign.}{lost}");
         RosStatus->GetXaxis()->SetBinLabel(8,"#splitline{ROS }{err.}");
         RosStatus->GetXaxis()->SetBinLabel(9,"#splitline{ROS }{in ev.}");
         RosStatus->GetXaxis()->SetBinLabel(10,"#splitline{Miss }{ev.}");
         RosStatus->GetXaxis()->SetBinLabel(11,"#splitline{Ev ID}{Miss}");
         RosStatus->GetXaxis()->SetBinLabel(12,"#splitline{BX}{Miss}");
  
         //RosStatus->LabelsOption("h","X");
         RosStatus->GetYaxis()->SetLabelSize(0.08);
         RosStatus->GetXaxis()->SetLabelSize(0.07);
         RosStatus->Draw("colz");
       }
       if(RosSummary)
       {
         ct7[i]->cd(iw*3+3);
	 ct7[i]->cd(iw*3+3)->SetGrid();

  
         RosSummary->SetStats(0);
         /*
         RosSummary->GetXaxis()->SetBinLabel(1,"#splitline{Link}{TimeOut}");
         RosSummary->GetXaxis()->SetBinLabel(2,"#splitline{Ev.Id.}{Mis.}");
         RosSummary->GetXaxis()->SetBinLabel(3,"#splitline{FIFO}{#splitline{almost}{full}}");
         RosSummary->GetXaxis()->SetBinLabel(4,"#splitline{FIFO}{full}");
         RosSummary->GetXaxis()->SetBinLabel(5,"#splitline{Ceros}{Timeout}");
         RosSummary->GetXaxis()->SetBinLabel(6,"#splitline{Max.}{wds}");
         RosSummary->GetXaxis()->SetBinLabel(7,"#splitline{L1A}{FF}");
         //RosSummary->GetXaxis()->SetBinLabel(8,"#splitline{PC}{#splitline{from}{TDC}}");
         RosSummary->GetXaxis()->SetBinLabel(8,"#splitline{TDC}{#splitline{parity}{err.}}");
  
         RosSummary->GetXaxis()->SetBinLabel(9,"#splitline{BX ID}{Mis.}");
         RosSummary->GetXaxis()->SetBinLabel(10,"TXP");
  
         //RosSummary->GetXaxis()->SetBinLabel(11,"#splitline{TDC}{Fatal}");
         //RosSummary->GetXaxis()->SetBinLabel(12,"#splitline{TDC}{#splitline{FIFO}{Ov.}}");
         //RosSummary->GetXaxis()->SetBinLabel(13,"#splitline{L1}{#splitline{Buffer}{Ov.}}");
  
         RosSummary->GetXaxis()->SetBinLabel(11,"#splitline{L1A}{#splitline{almost}{full}}");
         RosSummary->GetXaxis()->SetBinLabel(12,"#splitline{Ch}{blocked}");
         RosSummary->GetXaxis()->SetBinLabel(13,"#splitline{Ev.}{#splitline{Id.}{Mis.}}");
  
         RosSummary->GetXaxis()->SetBinLabel(14,"#splitline{CEROS}{blocked}");
         RosSummary->GetXaxis()->SetBinLabel(15,"#splitline{TDC}{Fatal}");
         RosSummary->GetXaxis()->SetBinLabel(16,"#splitline{TDC}{#splitline{FIFO}{Ov.}}");
         RosSummary->GetXaxis()->SetBinLabel(17,"#splitline{L1}{#splitline{Buffer}{Ov.}}");
         RosSummary->GetXaxis()->SetBinLabel(18,"#splitline{L1A}{#splitline{FIFO}{Ov.}}");
         RosSummary->GetXaxis()->SetBinLabel(19,"#splitline{TDC}{#splitline{hit}{err.}}");
         RosSummary->GetXaxis()->SetBinLabel(20,"#splitline{TDC}{#splitline{hit}{rej.}}");
         //RosSummary->LabelsOption("h","X");
  
         RosSummary->GetYaxis()->SetLabelSize(0.08);
         RosSummary->GetXaxis()->SetLabelSize(0.04);
         */
         RosSummary->GetXaxis()->SetBinLabel(1, "A");
         RosSummary->GetXaxis()->SetBinLabel(2, "B");
         RosSummary->GetXaxis()->SetBinLabel(3, "C");
         RosSummary->GetXaxis()->SetBinLabel(4, "D");
         RosSummary->GetXaxis()->SetBinLabel(5, "E");
         RosSummary->GetXaxis()->SetBinLabel(6, "F");
         RosSummary->GetXaxis()->SetBinLabel(7, "G");
         RosSummary->GetXaxis()->SetBinLabel(8, "H");
         RosSummary->GetXaxis()->SetBinLabel(9, "I");
         RosSummary->GetXaxis()->SetBinLabel(10,"J");
         RosSummary->GetXaxis()->SetBinLabel(11,"K");
         RosSummary->GetXaxis()->SetBinLabel(12,"L");
         RosSummary->GetXaxis()->SetBinLabel(13,"M");
         RosSummary->GetXaxis()->SetBinLabel(14,"N");
         RosSummary->GetXaxis()->SetBinLabel(15,"O");
         RosSummary->GetXaxis()->SetBinLabel(16,"P");
         RosSummary->GetXaxis()->SetBinLabel(17,"Q");
         RosSummary->GetXaxis()->SetBinLabel(18,"R");
         RosSummary->GetXaxis()->SetBinLabel(19,"S");
         RosSummary->GetXaxis()->SetBinLabel(20,"T");
  
         RosSummary->GetYaxis()->SetLabelSize(0.08);
         RosSummary->GetXaxis()->SetLabelSize(0.11);
         RosSummary->Draw("colz");
       }
     }
  
     ct7[i]->cd(18);
     TText L1;
     L1.SetTextSize(0.09);
     L1.SetTextColor(2);
     L1.DrawText(0.02,0.91,"A = "); 
     L1.DrawText(0.02,0.81,"B = "); 
     L1.DrawText(0.02,0.71,"C = "); 
     L1.DrawText(0.02,0.61,"D = "); 
     L1.DrawText(0.02,0.51,"E = "); 
     L1.DrawText(0.02,0.41,"F = "); 
     L1.DrawText(0.02,0.31,"G = "); 
     L1.DrawText(0.02,0.21,"H = "); 
     L1.DrawText(0.02,0.11,"I = "); 
     L1.DrawText(0.02,0.01,"J = "); 
                                 
     L1.DrawText(0.50,0.91,"K = "); 
     L1.DrawText(0.50,0.81,"L = "); 
     L1.DrawText(0.50,0.71,"M = "); 
     L1.DrawText(0.50,0.61,"N = "); 
     L1.DrawText(0.50,0.51,"O = "); 
     L1.DrawText(0.50,0.41,"P = "); 
     L1.DrawText(0.50,0.31,"Q = "); 
     L1.DrawText(0.50,0.21,"R = "); 
     L1.DrawText(0.50,0.11,"S = "); 
     L1.DrawText(0.50,0.01,"T = "); 
       
     TText L2;
     L2.SetTextSize(0.09);
     L2.SetTextColor(4);
     L2.DrawText(0.07,0.91," Link TimeOut");
     L2.DrawText(0.07,0.81," Ev.Id. Missmatch");
     L2.DrawText(0.07,0.71," FIFO almost full");
     L2.DrawText(0.07,0.61," FIFO  full");
     L2.DrawText(0.07,0.51," Ceros TimeOut");
     L2.DrawText(0.07,0.41," Max.  wds");
     L2.DrawText(0.07,0.31," L1A  FF");
     L2.DrawText(0.07,0.21," TDC parity error");
     L2.DrawText(0.07,0.11," BX ID Missmatch");
     L2.DrawText(0.07,0.01," TXP");
  
     L2.DrawText(0.55,0.91," L1A almost full");
     L2.DrawText(0.55,0.81," Channel blocked");
     L2.DrawText(0.55,0.71," Event Id. Missmatch");
     L2.DrawText(0.55,0.61," CEROS  blocked");
     L2.DrawText(0.55,0.51," TDC Fatal");
     L2.DrawText(0.55,0.41," TDC FIFO Overflow");
     L2.DrawText(0.55,0.31," L1  Buffer Overflow");
     L2.DrawText(0.55,0.21," L1A FIFO Overflow");
     L2.DrawText(0.55,0.11," TDC hit error");  
     L2.DrawText(0.55,0.01," TDC hit rej.");
   
   if(i==0)createGifFile("DataIntegrity_a",ct7[i]);
//   if(i==1)createGifFile("DataIntegrity_b",ct7[i]);
  } 


// ==========================================================================================================
 char Whname[5][20]={"Wm2","Wm1","W0","W1","W2"};

 //gStyle->SetStatX(0.4);
 //gStyle->SetStatW(0.29);
 gStyle->SetStatX(0.6);
 gStyle->SetStatW(0.49);
 gStyle->SetStatY(0.88);
 gStyle->SetStatH(0.38);


 TH1F * first_to_paint;
 TH1F * first_to_paint_MB[4][4];

 if(ProcessDDUTrigger)
 {
 //stringstream folder0; folder0<<"DQMData/Run 108478/DT/Run summary/04-LocalTrigger-DDU/";
   stringstream folder0; folder0 << myMainFolder << "04-LocalTrigger-DDU/"; 

 TH1F *QualBestDDUSec[5][12];
 TH1F *QualBestDDUMB[5][4];
 TH1F *QualBestDDUMBTop[5][4];
 TH1F *QualBestDDUMBBottom[5][4];
 TH1F *QualBestDDUMBVertical[5][4];
 for (int iw=0; iw<5;iw++)
  for (int ins=1; ins<13;ins++)
   for(int ic=1;ic<5;ic++){
    //Changed to move to TM IN (06/10/2016)
    // stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhi/" ;
     stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhiIn/" ;
     stringstream hname; hname << "DDU_BestQual_W" << iw-2 << "_Sec" << ins << "_St" << ic;

     char hnamec[240];
     sprintf(hnamec,"%s%s%s",folder0.str().c_str(),folder1.str().c_str(),hname.str().c_str());
     TH1F * theHisto= (TH1F*)myFile->Get(hnamec);

     if(ins==1)
     {
        stringstream hname2; hname2 << "DDU_BestQual_" << Whname[iw] << "_MB" << ic ;
        QualBestDDUMB[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUMB[iw][ic-1]->SetTitle(hname2.str().c_str());

        stringstream hname3; hname3 << "DDU_BestQual_" << Whname[iw] << "_MB" << ic << "_Vertical";
        QualBestDDUMBVertical[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUMBVertical[iw][ic-1]->SetTitle(hname3.str().c_str());

     }
     if(ins==2)
     {
        stringstream hname2; hname2 << "DDU_BestQual_" << Whname[iw] << "_MB" << ic << "_Top"  ;
        QualBestDDUMBTop[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUMBTop[iw][ic-1]->SetTitle(hname2.str().c_str());


     }
     if(ins==8)
     {
        stringstream hname2; hname2 << "DDU_BestQual_" << Whname[iw] << "_MB" << ic << "_Bottom"  ;
        QualBestDDUMBBottom[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUMBBottom[iw][ic-1]->SetTitle(hname2.str().c_str());

     }


     if(ins>2)QualBestDDUMB[iw][ic-1]->Add(theHisto);   // 
     if(ins>2 && ins<7)QualBestDDUMBTop[iw][ic-1]->Add(theHisto);
     if(ins>8)QualBestDDUMBBottom[iw][ic-1]->Add(theHisto);
     if(ins==7)QualBestDDUMBVertical[iw][ic-1]->Add(theHisto);

     if(ic==1)
     {
        stringstream hname2; hname2 << "DDU_BestQual_" << Whname[iw] << "_S" << ins;
        QualBestDDUSec[iw][ins-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUSec[iw][ins-1]->SetTitle(hname2.str().c_str());
     }
     if(ic>1)QualBestDDUSec[iw][ins-1]->Add(theHisto);
   }

 for(int ic=0;ic<4;ic++)
 {
   first_to_paint_MB[ic][0]=QualBestDDUMB[0][ic];
   first_to_paint_MB[ic][1]=QualBestDDUMBTop[0][ic];
   first_to_paint_MB[ic][2]=QualBestDDUMBBottom[0][ic];
   first_to_paint_MB[ic][3]=QualBestDDUMBVertical[0][ic];
   float nbmax[4]={0,0,0,0};
   nbmax[0]=QualBestDDUMB[0][ic]->GetMaximum();
   nbmax[1]=QualBestDDUMBTop[0][ic]->GetMaximum();
   nbmax[2]=QualBestDDUMBBottom[0][ic]->GetMaximum();
   nbmax[3]=QualBestDDUMBVertical[0][ic]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax0[4]={0,0,0,0};
      nbmax0[0]=QualBestDDUMB[iw][ic]->GetMaximum();
      nbmax0[1]=QualBestDDUMBTop[iw][ic]->GetMaximum();
      nbmax0[2]=QualBestDDUMBBottom[iw][ic]->GetMaximum();
      nbmax0[3]=QualBestDDUMBVertical[iw][ic]->GetMaximum();

      if(nbmax0[0]>nbmax[0])
      {
         nbmax[0]=nbmax0[0]; first_to_paint_MB[ic][0]=QualBestDDUMB[iw][ic];
      }
      if(nbmax0[1]>nbmax[1])
      {
         nbmax[1]=nbmax0[1]; first_to_paint_MB[ic][1]=QualBestDDUMBTop[iw][ic];
      }
      if(nbmax0[2]>nbmax[2])
      {
         nbmax[2]=nbmax0[2]; first_to_paint_MB[ic][2]=QualBestDDUMBBottom[iw][ic];
      }
      if(nbmax0[3]>nbmax[3])
      {
         nbmax[3]=nbmax0[3]; first_to_paint_MB[ic][3]=QualBestDDUMBVertical[iw][ic];
      }
   }
 }

 for(int ic=1;ic<5;ic++)
 {
   stringstream hname21; hname21 << "DDU_BestQual_MB" << ic ;
   stringstream hname22; hname22 << "DDU_BestQual_MB" << ic << "_Top"  ;
   stringstream hname23; hname23 << "DDU_BestQual_MB" << ic << "_Bottom"  ;
   stringstream hname24; hname24 << "DDU_BestQual_MB" << ic << "_Vertical"  ;
   first_to_paint_MB[ic-1][0]->SetTitle(hname21.str().c_str());
   first_to_paint_MB[ic-1][1]->SetTitle(hname22.str().c_str());
   first_to_paint_MB[ic-1][2]->SetTitle(hname23.str().c_str());
   first_to_paint_MB[ic-1][3]->SetTitle(hname24.str().c_str());

   for(int ityp=0;ityp<4;ityp++)
   {
     first_to_paint_MB[ic-1][ityp]->GetXaxis()->SetLabelSize(0.07);
     first_to_paint_MB[ic-1][ityp]->GetYaxis()->SetLabelSize(0.05);
     first_to_paint_MB[ic-1][ityp]->SetNdivisions(505);
   }
 }


 TCanvas *QualBestDDUAll = new TCanvas("QualBestDDUAll", "",201,81,999,950);
 QualBestDDUAll->Divide(4,4) ;
 for(int ip=1;ip<17;ip++)
 {
   QualBestDDUAll->cd(ip)->SetFillColor(0) ;
   QualBestDDUAll->cd(ip)->SetFrameFillColor(0) ;
 }
 for(int ic=1;ic<5;ic++)
 {
    first_to_paint_MB[ic-1][0]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetLabelSize(0.12);

    QualBestDDUAll->cd((ic-1)*4+1)->SetFillColor(0); first_to_paint_MB[ic-1][0]->Draw();
    QualBestDDUAll->cd((ic-1)*4+2)->SetFillColor(0); first_to_paint_MB[ic-1][1]->Draw();
    QualBestDDUAll->cd((ic-1)*4+3)->SetFillColor(0); first_to_paint_MB[ic-1][2]->Draw();
    QualBestDDUAll->cd((ic-1)*4+4)->SetFillColor(0); first_to_paint_MB[ic-1][3]->Draw();

    for(int iw=0;iw<5;iw++)
    {
      QualBestDDUMB[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUMBTop[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUMBBottom[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUMBVertical[iw][ic-1]->SetLineColor(iw+1);
      if(iw==4)
       {
         QualBestDDUMB[iw][ic-1]->SetLineColor(6);
         QualBestDDUMBTop[iw][ic-1]->SetLineColor(6);
         QualBestDDUMBBottom[iw][ic-1]->SetLineColor(6);
         QualBestDDUMBVertical[iw][ic-1]->SetLineColor(6);
       }
      QualBestDDUAll->cd((ic-1)*4+1);  QualBestDDUMB[iw][ic-1]->Draw("same");
      QualBestDDUAll->cd((ic-1)*4+2);  QualBestDDUMBTop[iw][ic-1]->Draw("same");
      QualBestDDUAll->cd((ic-1)*4+3);  QualBestDDUMBBottom[iw][ic-1]->Draw("same");
      QualBestDDUAll->cd((ic-1)*4+4);  QualBestDDUMBVertical[iw][ic-1]->Draw("same");
    }
 }

 createGifFile("QualBestAllDDU",QualBestDDUAll);
 delete QualBestDDUAll;


 TCanvas *QualBestDDUWh[5];
 QualBestDDUWh[0] = new TCanvas("QualBestDDUWh-2", "",201,81,999,950);
 QualBestDDUWh[1] = new TCanvas("QualBestDDUWh-1", "",201,81,999,950);
 QualBestDDUWh[2] = new TCanvas("QualBestDDUWh0", "",201,81,999,950);
 QualBestDDUWh[3] = new TCanvas("QualBestDDUWh+1", "",201,81,999,950);
 QualBestDDUWh[4] = new TCanvas("QualBestDDUWh+2", "",201,81,999,950);
 for(int iw=0;iw<5;iw++)
 {
   QualBestDDUWh[iw]->Divide(3,4) ;
   for(int isec=1;isec<13;isec++)
   {
      QualBestDDUWh[iw]->cd(isec);
      QualBestDDUWh[iw]->cd(isec)->SetFillColor(0);
      QualBestDDUWh[iw]->cd(isec)->SetFrameFillColor(0);
      QualBestDDUSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.07);
      QualBestDDUSec[iw][isec-1]->GetYaxis()->SetLabelSize(0.05);
      QualBestDDUSec[iw][isec-1]->SetNdivisions(505);
      QualBestDDUSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.12);
      QualBestDDUSec[iw][isec-1]->Draw();   // 
   }
 }

 for(int iw=0;iw<5;iw++)
 {
   createGifFile("QualBestSecDDU",QualBestDDUWh[iw],(iw-2));
   delete QualBestDDUWh[iw];
 }


 TCanvas *QualBestDDUWhAll;
 QualBestDDUWhAll = new TCanvas("QualBestDDUWhAll", "",201,81,999,950);
 QualBestDDUWhAll->Divide(3,4) ;
 for(int ins=1;ins<13;ins++)
 {
   first_to_paint=QualBestDDUSec[0][ins-1];
   float nbmax1=QualBestDDUSec[0][ins-1]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax2=QualBestDDUSec[iw][ins-1]->GetMaximum();
      if(nbmax2>nbmax1)
      {
         nbmax1=nbmax2; first_to_paint=QualBestDDUSec[iw][ins-1];
      }
    }

    QualBestDDUWhAll->cd(ins) ;
    QualBestDDUWhAll->cd(ins)->SetFillColor(0) ;
    QualBestDDUWhAll->cd(ins)->SetFrameFillColor(0) ;
    stringstream hname; hname << "DDU_BestQual_S" << ins;
    first_to_paint->SetTitle(hname.str().c_str());
    first_to_paint->GetXaxis()->SetLabelSize(0.12);
    first_to_paint->Draw();
    for(int iw=0;iw<5;iw++)
    {

      QualBestDDUSec[iw][ins-1]->SetLineColor(iw+1);
      if(iw==4) QualBestDDUSec[iw][ins-1]->SetLineColor(6);
      QualBestDDUSec[iw][ins-1]->Draw("same");
    }
 }

 createGifFile("QualBestSecAllDDU",QualBestDDUWhAll);
 delete QualBestDDUWhAll;


 TH1F *QualBestDDUThetaSec[5][12];
 TH1F *QualBestDDUThetaMB[5][4];
 TH1F *QualBestDDUThetaMBTop[5][4];
 TH1F *QualBestDDUThetaMBBottom[5][4];
 TH1F *QualBestDDUThetaMBVertical[5][4];
 char QualLabelTheta[7][3]={"  ","L","  ","H","  ","  ","  "};
 for (int iw=0; iw<5;iw++)
  for (int ins=1; ins<13;ins++)
   for(int ic=1;ic<4;ic++){// MB4 has not Theta SL 

     stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerTheta/" ;
     stringstream hname; hname << "DDU_ThetaBestQual_W" << iw-2 << "_Sec" << ins << "_St" << ic;

     char hnamec[240];
     sprintf(hnamec,"%s%s%s",folder0.str().c_str(),folder1.str().c_str(),hname.str().c_str());
     TH1F * theHisto= (TH1F*)myFile->Get(hnamec);

     if(ins==1)
     {
        stringstream hname2; hname2 << "DDU_ThetaBestQual_" << Whname[iw] << "_MB" << ic ;

        QualBestDDUThetaMB[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUThetaMB[iw][ic-1]->SetTitle(hname2.str().c_str());

        stringstream hname3; hname3 << "DDU_ThetaBestQual_" << Whname[iw] << "_MB" << ic << "_Vertical" ;

        QualBestDDUThetaMBVertical[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUThetaMBVertical[iw][ic-1]->SetTitle(hname3.str().c_str());

     }
     if(ins==2)
     {
        stringstream hname2; hname2 << "DDU_ThetaBestQual_" << Whname[iw] << "_MB" << ic << "_Top"  ;

        QualBestDDUThetaMBTop[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUThetaMBTop[iw][ic-1]->SetTitle(hname2.str().c_str());


     }
     if(ins==8)
     {
        stringstream hname2; hname2 << "DDU_ThetaBestQual_" << Whname[iw] << "_MB" << ic << "_Bottom"  ;
        QualBestDDUThetaMBBottom[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUThetaMBBottom[iw][ic-1]->SetTitle(hname2.str().c_str());

     }

     if(ins>2)QualBestDDUThetaMB[iw][ic-1]->Add(theHisto);   // 
     if(ins>2 && ins<7)QualBestDDUThetaMBTop[iw][ic-1]->Add(theHisto);
     if(ins>8)QualBestDDUThetaMBBottom[iw][ic-1]->Add(theHisto);
     if(ins==7)QualBestDDUThetaMBVertical[iw][ic-1]->Add(theHisto);

     if(ic==1) 
     {
        stringstream hname2; hname2 << "DDU_ThetaBestQual_" << Whname[iw] << "_S" << ins;
        QualBestDDUThetaSec[iw][ins-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestDDUThetaSec[iw][ins-1]->SetTitle(hname2.str().c_str());
     }
        if(ic>1)QualBestDDUThetaSec[iw][ins-1]->Add(theHisto);
   }


 for(int ic=0;ic<3;ic++)
 {
   first_to_paint_MB[ic][0]=QualBestDDUThetaMB[0][ic];
   first_to_paint_MB[ic][1]=QualBestDDUThetaMBTop[0][ic];
   first_to_paint_MB[ic][2]=QualBestDDUThetaMBBottom[0][ic];
   first_to_paint_MB[ic][3]=QualBestDDUThetaMBVertical[0][ic];
   float nbmax[4]={0,0,0,0};
   nbmax[0]=QualBestDDUThetaMB[0][ic]->GetMaximum();
   nbmax[1]=QualBestDDUThetaMBTop[0][ic]->GetMaximum();
   nbmax[2]=QualBestDDUThetaMBBottom[0][ic]->GetMaximum();
   nbmax[3]=QualBestDDUThetaMBVertical[0][ic]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax0[4]={0,0,0,0};
      nbmax0[0]=QualBestDDUThetaMB[iw][ic]->GetMaximum();
      nbmax0[1]=QualBestDDUThetaMBTop[iw][ic]->GetMaximum();
      nbmax0[2]=QualBestDDUThetaMBBottom[iw][ic]->GetMaximum();
      nbmax0[3]=QualBestDDUThetaMBVertical[iw][ic]->GetMaximum();

      if(nbmax0[0]>nbmax[0])
      {
         nbmax[0]=nbmax0[0]; first_to_paint_MB[ic][0]=QualBestDDUThetaMB[iw][ic];
      }
      if(nbmax0[1]>nbmax[1])
      {
         nbmax[1]=nbmax0[1]; first_to_paint_MB[ic][1]=QualBestDDUThetaMBTop[iw][ic];
      }
      if(nbmax0[2]>nbmax[2])
      {
         nbmax[2]=nbmax0[2]; first_to_paint_MB[ic][2]=QualBestDDUThetaMBBottom[iw][ic];
      }
      if(nbmax0[3]>nbmax[3])
      {
         nbmax[3]=nbmax0[3]; first_to_paint_MB[ic][3]=QualBestDDUThetaMBVertical[iw][ic];
      }
   }
 }
 for(int ic=1;ic<4;ic++)
 {
   stringstream hname21; hname21 << "DDU_ThetaBestQual_MB" << ic ;
   stringstream hname22; hname22 << "DDU_ThetaBestQual_MB" << ic << "_Top"  ;
   stringstream hname23; hname23 << "DDU_ThetaBestQual_MB" << ic << "_Bottom"  ;
   stringstream hname24; hname24 << "DDU_ThetaBestQual_MB" << ic << "_Vertical"  ;
   first_to_paint_MB[ic-1][0]->SetTitle(hname21.str().c_str());
   first_to_paint_MB[ic-1][1]->SetTitle(hname22.str().c_str());
   first_to_paint_MB[ic-1][2]->SetTitle(hname23.str().c_str());
   first_to_paint_MB[ic-1][3]->SetTitle(hname24.str().c_str());

   for(int ityp=0;ityp<4;ityp++)
   {
     first_to_paint_MB[ic-1][ityp]->GetXaxis()->SetLabelSize(0.07);
     first_to_paint_MB[ic-1][ityp]->GetYaxis()->SetLabelSize(0.05);
     first_to_paint_MB[ic-1][ityp]->SetNdivisions(505);
   }
 }


 TCanvas *ThetaQualBestDDUAll = new TCanvas("QualBestDDUAll", "",201,81,999,950);
 ThetaQualBestDDUAll->Divide(4,4) ;
 for(int ip=1;ip<17;ip++)
 {
    ThetaQualBestDDUAll->cd(ip)->SetFillColor(0) ;
    ThetaQualBestDDUAll->cd(ip)->SetFrameFillColor(0) ;
 }
 for(int ic=1;ic<4;ic++)
 {
   for(int ib=1;ib<8;ib++)
   {

    first_to_paint_MB[ic-1][0]->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);

   } 

    first_to_paint_MB[ic-1][0]->GetXaxis()->SetLabelSize(0.15);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetLabelSize(0.15);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetLabelSize(0.15);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetLabelSize(0.15);

    ThetaQualBestDDUAll->cd((ic-1)*4+1)->SetFillColor(0); first_to_paint_MB[ic-1][0]->Draw();
    ThetaQualBestDDUAll->cd((ic-1)*4+2)->SetFillColor(0); first_to_paint_MB[ic-1][1]->Draw();
    ThetaQualBestDDUAll->cd((ic-1)*4+3)->SetFillColor(0); first_to_paint_MB[ic-1][2]->Draw();
    ThetaQualBestDDUAll->cd((ic-1)*4+4)->SetFillColor(0); first_to_paint_MB[ic-1][3]->Draw();

    for(int iw=0;iw<5;iw++)
    {
      QualBestDDUThetaMB[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUThetaMBTop[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUThetaMBBottom[iw][ic-1]->SetLineColor(iw+1);
      QualBestDDUThetaMBVertical[iw][ic-1]->SetLineColor(iw+1);
      if(iw==4)
       {
         QualBestDDUThetaMB[iw][ic-1]->SetLineColor(6);
         QualBestDDUThetaMBTop[iw][ic-1]->SetLineColor(6);
         QualBestDDUThetaMBBottom[iw][ic-1]->SetLineColor(6);
         QualBestDDUThetaMBVertical[iw][ic-1]->SetLineColor(6);
       }
      ThetaQualBestDDUAll->cd((ic-1)*4+1);  QualBestDDUThetaMB[iw][ic-1]->Draw("same");
      ThetaQualBestDDUAll->cd((ic-1)*4+2);  QualBestDDUThetaMBTop[iw][ic-1]->Draw("same");
      ThetaQualBestDDUAll->cd((ic-1)*4+3);  QualBestDDUThetaMBBottom[iw][ic-1]->Draw("same");
      ThetaQualBestDDUAll->cd((ic-1)*4+4);  QualBestDDUThetaMBVertical[iw][ic-1]->Draw("same");
    }
 }
 
 createGifFile("QualBestThetaAllDDU",ThetaQualBestDDUAll);
 delete ThetaQualBestDDUAll;

 TCanvas *QualBestDDUThetaWh[5];
 QualBestDDUThetaWh[0] = new TCanvas("QualBestDDUWh-2", "",201,81,999,950);
 QualBestDDUThetaWh[1] = new TCanvas("QualBestDDUWh-1", "",201,81,999,950);
 QualBestDDUThetaWh[2] = new TCanvas("QualBestDDUWh0", "",201,81,999,950);
 QualBestDDUThetaWh[3] = new TCanvas("QualBestDDUWh+1", "",201,81,999,950);
 QualBestDDUThetaWh[4] = new TCanvas("QualBestDDUWh+2", "",201,81,999,950);
 for(int iw=0;iw<5;iw++)
 {
   QualBestDDUThetaWh[iw]->Divide(3,4) ;
   for(int isec=1;isec<13;isec++)
   {
      QualBestDDUThetaWh[iw]->cd(isec);
      QualBestDDUThetaWh[iw]->cd(isec)->SetFillColor(0);
      QualBestDDUThetaWh[iw]->cd(isec)->SetFrameFillColor(0);
      QualBestDDUThetaSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.07);
      QualBestDDUThetaSec[iw][isec-1]->GetYaxis()->SetLabelSize(0.05);
      QualBestDDUThetaSec[iw][isec-1]->SetNdivisions(505);
      for(int ib=1;ib<8;ib++)
       QualBestDDUThetaSec[iw][isec-1]->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);
      QualBestDDUThetaSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.15);
      QualBestDDUThetaSec[iw][isec-1]->Draw();   // 
   }
 }
 for(int iw=0;iw<5;iw++)
 {
   createGifFile("QualBestThetaSecDDU",QualBestDDUThetaWh[iw],iw-2,false);
   delete QualBestDDUThetaWh[iw];
 }

 TCanvas *QualBestDDUThetaWhAll;
 QualBestDDUThetaWhAll = new TCanvas("QualBestDDUThetaWhAll", "",201,81,999,950);
 QualBestDDUThetaWhAll->Divide(3,4) ;
 for(int ins=1;ins<13;ins++)
 {
   first_to_paint=QualBestDDUThetaSec[0][ins-1];
   float nbmax1=QualBestDDUThetaSec[0][ins-1]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax2=QualBestDDUThetaSec[iw][ins-1]->GetMaximum();
      if(nbmax2>nbmax1)
      {
         nbmax1=nbmax2; first_to_paint=QualBestDDUThetaSec[iw][ins-1];
      }
    }

    QualBestDDUThetaWhAll->cd(ins) ;
    QualBestDDUThetaWhAll->cd(ins)->SetFillColor(0) ;
    QualBestDDUThetaWhAll->cd(ins)->SetFrameFillColor(0) ;
    stringstream hname; hname << "DDU_BestThetaQual_S" << ins;
    first_to_paint->SetTitle(hname.str().c_str());
    for(int ib=1;ib<8;ib++)
      first_to_paint->GetXaxis()->SetBinLabel(ib,QualLabelTheta[ib-1]);
    first_to_paint->GetXaxis()->SetLabelSize(0.15);
    first_to_paint->Draw();
    for(int iw=0;iw<5;iw++)
    {

      QualBestDDUThetaSec[iw][ins-1]->SetLineColor(iw+1);
      if(iw==4) QualBestDDUThetaSec[iw][ins-1]->SetLineColor(6);
      QualBestDDUThetaSec[iw][ins-1]->Draw("same");
    }
  }


 createGifFile("QualBestThetaSecAllDDU",QualBestDDUThetaWhAll);
 delete QualBestDDUThetaWhAll;
 } // END if(ProcessDDUTrigger)


 stringstream folder0b; folder0b<< myMainFolder << "03-LocalTrigger-TM/";

 TH1F *QualBestTMSec[5][12];
 TH1F *QualBestTMMB[5][4];
 TH1F *QualBestTMMBTop[5][4];
 TH1F *QualBestTMMBBottom[5][4];
 TH1F *QualBestTMMBVertical[5][4];

 for (int iw=0; iw<5;iw++)
  for (int ins=1; ins<13;ins++)
   for(int ic=1;ic<5;ic++){

    //Changed to move to TM IN (06/10/2016)
    //stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhi/" ;
    stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhiIn/" ;
    stringstream hname; hname << "TM_BestQual_W" << iw-2 << "_Sec" << ins << "_St" << ic;

    char hnamec[240];
    sprintf(hnamec,"%s%s%s",folder0b.str().c_str(),folder1.str().c_str(),hname.str().c_str());
    TH1F * theHisto= NULL; 
    theHisto= (TH1F*)myFile->Get(hnamec);

    if(theHisto !=NULL) // Protecting in case the histo doesn't exist 
    {

     if(ins==1)
     {
        stringstream hname2; hname2 << "TM_BestQual_" << Whname[iw] << "_MB" << ic ;
        QualBestTMMB[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestTMMB[iw][ic-1]->SetTitle(hname2.str().c_str());

        stringstream hname3; hname3 << "TM_BestQual_" << Whname[iw] << "_MB" << ic << "_Vertical";
        QualBestTMMBVertical[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestTMMBVertical[iw][ic-1]->SetTitle(hname3.str().c_str());

     }
     if(ins==2)
     {
        stringstream hname2; hname2 << "TM_BestQual_" << Whname[iw] << "_MB" << ic << "_Top"  ;
        QualBestTMMBTop[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestTMMBTop[iw][ic-1]->SetTitle(hname2.str().c_str());
     }
     if(ins==8)
     {
        stringstream hname2; hname2 << "TM_BestQual_" << Whname[iw] << "_MB" << ic << "_Bottom"  ;
        QualBestTMMBBottom[iw][ic-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestTMMBBottom[iw][ic-1]->SetTitle(hname2.str().c_str());

     }
     if(ins>2)QualBestTMMB[iw][ic-1]->Add(theHisto);   // 
     if(ins>2 && ins<7)QualBestTMMBTop[iw][ic-1]->Add(theHisto);
     if(ins>8)QualBestTMMBBottom[iw][ic-1]->Add(theHisto);
     if(ins==7)QualBestTMMBVertical[iw][ic-1]->Add(theHisto);

     if(ic==1)
     {
        stringstream hname2; hname2 << "TM_BestQual_" << Whname[iw] << "_S" << ins;
        QualBestTMSec[iw][ins-1]=(TH1F*)theHisto->Clone(hname.str().c_str());
        QualBestTMSec[iw][ins-1]->SetTitle(hname2.str().c_str());
     }
     if(ic>1)QualBestTMSec[iw][ins-1]->Add(theHisto);
    }
   }

 for(int ic=0;ic<4;ic++)
 {
   first_to_paint_MB[ic][0]=QualBestTMMB[0][ic];
   first_to_paint_MB[ic][1]=QualBestTMMBTop[0][ic];
   first_to_paint_MB[ic][2]=QualBestTMMBBottom[0][ic];
   first_to_paint_MB[ic][3]=QualBestTMMBVertical[0][ic];
   float nbmax[4]={0,0,0,0};
   nbmax[0]=QualBestTMMB[0][ic]->GetMaximum();
   nbmax[1]=QualBestTMMBTop[0][ic]->GetMaximum();
   nbmax[2]=QualBestTMMBBottom[0][ic]->GetMaximum();
   nbmax[3]=QualBestTMMBVertical[0][ic]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax0[4]={0,0,0,0};
      nbmax0[0]=QualBestTMMB[iw][ic]->GetMaximum();
      nbmax0[1]=QualBestTMMBTop[iw][ic]->GetMaximum();
      nbmax0[2]=QualBestTMMBBottom[iw][ic]->GetMaximum();
      nbmax0[3]=QualBestTMMBVertical[iw][ic]->GetMaximum();

      if(nbmax0[0]>nbmax[0])
      {
         nbmax[0]=nbmax0[0]; first_to_paint_MB[ic][0]=QualBestTMMB[iw][ic];
      }
      if(nbmax0[1]>nbmax[1])
      {
         nbmax[1]=nbmax0[1]; first_to_paint_MB[ic][1]=QualBestTMMBTop[iw][ic];
      }
      if(nbmax0[2]>nbmax[2])
      {
         nbmax[2]=nbmax0[2]; first_to_paint_MB[ic][2]=QualBestTMMBBottom[iw][ic];
      }
      if(nbmax0[3]>nbmax[3])
      {
         nbmax[3]=nbmax0[3]; first_to_paint_MB[ic][3]=QualBestTMMBVertical[iw][ic];
      }
   }
 }

 for(int ic=1;ic<5;ic++)
 {
   stringstream hname21; hname21 << "TM_BestQual_MB" << ic ;
   stringstream hname22; hname22 << "TM_BestQual_MB" << ic << "_Top"  ;
   stringstream hname23; hname23 << "TM_BestQual_MB" << ic << "_Bottom"  ;
   stringstream hname24; hname24 << "TM_BestQual_MB" << ic << "_Vertical"  ;
   first_to_paint_MB[ic-1][0]->SetTitle(hname21.str().c_str());
   first_to_paint_MB[ic-1][1]->SetTitle(hname22.str().c_str());
   first_to_paint_MB[ic-1][2]->SetTitle(hname23.str().c_str());
   first_to_paint_MB[ic-1][3]->SetTitle(hname24.str().c_str());

   for(int ityp=0;ityp<4;ityp++)
   {
     first_to_paint_MB[ic-1][ityp]->GetXaxis()->SetLabelSize(0.07);
     first_to_paint_MB[ic-1][ityp]->GetYaxis()->SetLabelSize(0.05);
     first_to_paint_MB[ic-1][ityp]->SetNdivisions(505);
   }
 }


 TCanvas *QualBestTMAll = new TCanvas("QualBestTMAll", "",201,81,999,950);
 QualBestTMAll->Divide(4,4) ;
 for(int ip=1;ip<17;ip++)
 {
   QualBestTMAll->cd(ip)->SetFillColor(0) ;
   QualBestTMAll->cd(ip)->SetFrameFillColor(0) ;
 }
 for(int ic=1;ic<5;ic++)
 {
    first_to_paint_MB[ic-1][0]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetLabelSize(0.12);

    QualBestTMAll->cd((ic-1)*4+1)->SetFillColor(0); first_to_paint_MB[ic-1][0]->Draw();
    QualBestTMAll->cd((ic-1)*4+2)->SetFillColor(0); first_to_paint_MB[ic-1][1]->Draw();
    QualBestTMAll->cd((ic-1)*4+3)->SetFillColor(0); first_to_paint_MB[ic-1][2]->Draw();
    QualBestTMAll->cd((ic-1)*4+4)->SetFillColor(0); first_to_paint_MB[ic-1][3]->Draw();

    for(int iw=0;iw<5;iw++)
    {
      QualBestTMMB[iw][ic-1]->SetLineColor(iw+1);
      QualBestTMMBTop[iw][ic-1]->SetLineColor(iw+1);
      QualBestTMMBBottom[iw][ic-1]->SetLineColor(iw+1);
      QualBestTMMBVertical[iw][ic-1]->SetLineColor(iw+1);
      if(iw==4)
       {
         QualBestTMMB[iw][ic-1]->SetLineColor(6);
         QualBestTMMBTop[iw][ic-1]->SetLineColor(6);
         QualBestTMMBBottom[iw][ic-1]->SetLineColor(6);
         QualBestTMMBVertical[iw][ic-1]->SetLineColor(6);
       }
      QualBestTMAll->cd((ic-1)*4+1);  QualBestTMMB[iw][ic-1]->Draw("same");
      QualBestTMAll->cd((ic-1)*4+2);  QualBestTMMBTop[iw][ic-1]->Draw("same");
      QualBestTMAll->cd((ic-1)*4+3);  QualBestTMMBBottom[iw][ic-1]->Draw("same");
      QualBestTMAll->cd((ic-1)*4+4);  QualBestTMMBVertical[iw][ic-1]->Draw("same");
    }
 }


 createGifFile("QualBestAllTM",QualBestTMAll);
 delete QualBestTMAll;

 TCanvas *QualBestTMWh[5];
 QualBestTMWh[0] = new TCanvas("QualBestTMWh-2", "",201,81,999,950);
 QualBestTMWh[1] = new TCanvas("QualBestTMWh-1", "",201,81,999,950);
 QualBestTMWh[2] = new TCanvas("QualBestTMWh0", "",201,81,999,950);
 QualBestTMWh[3] = new TCanvas("QualBestTMWh+1", "",201,81,999,950);
 QualBestTMWh[4] = new TCanvas("QualBestTMWh+2", "",201,81,999,950);
 for(int iw=0;iw<5;iw++)
 {
   QualBestTMWh[iw]->Divide(3,4) ;
   for(int isec=1;isec<13;isec++)
   {
      QualBestTMWh[iw]->cd(isec);
      QualBestTMWh[iw]->cd(isec)->SetFillColor(0);
      QualBestTMWh[iw]->cd(isec)->SetFrameFillColor(0);
      QualBestTMSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.07);
      QualBestTMSec[iw][isec-1]->GetYaxis()->SetLabelSize(0.05);
      QualBestTMSec[iw][isec-1]->SetNdivisions(505);
      QualBestTMSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.12);
      QualBestTMSec[iw][isec-1]->Draw();   // 
   }
 }



 for(int iw=0;iw<5;iw++)
 {
   createGifFile("QualBestSecTM",QualBestTMWh[iw],iw-2,false);
   delete QualBestTMWh[iw];
 }

 TCanvas *QualBestTMWhAll;
 QualBestTMWhAll = new TCanvas("QualBestTMWhAll", "",201,81,999,950);
 QualBestTMWhAll->Divide(3,4) ;
 for(int ins=1;ins<13;ins++)
 {
   first_to_paint=QualBestTMSec[0][ins-1];
   float nbmax1=QualBestTMSec[0][ins-1]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax2=QualBestTMSec[iw][ins-1]->GetMaximum();
      if(nbmax2>nbmax1)
      {
         nbmax1=nbmax2; first_to_paint=QualBestTMSec[iw][ins-1];
      }
    }

    QualBestTMWhAll->cd(ins) ;
    QualBestTMWhAll->cd(ins)->SetFillColor(0) ;
    QualBestTMWhAll->cd(ins)->SetFrameFillColor(0) ;
    stringstream hname; hname << "TM_BestQual_S" << ins;
    first_to_paint->SetTitle(hname.str().c_str());
    first_to_paint->GetXaxis()->SetLabelSize(0.12);
    first_to_paint->Draw();
    for(int iw=0;iw<5;iw++)
    {
      QualBestTMSec[iw][ins-1]->SetLineColor(iw+1);
      if(iw==4) QualBestTMSec[iw][ins-1]->SetLineColor(6);
      QualBestTMSec[iw][ins-1]->Draw("same");
    }
 }

 createGifFile("QualBestSecAllTM",QualBestTMWhAll);
 delete QualBestTMWhAll;


 TH1F *QualBX0TMSec[5][12];
 TH1F *QualBX0TMMB[5][4];
 TH1F *QualBX0TMMBTop[5][4];
 TH1F *QualBX0TMMBBottom[5][4];
 TH1F *QualBX0TMMBVertical[5][4];
 for (int iw=0; iw<5;iw++)
  for(int ic=1;ic<5;ic++){
     stringstream hname21; hname21 << "TM_BX0Qual_" << Whname[iw] << "_MB" << ic ;
     QualBX0TMMB[iw][ic-1]= new TH1F(hname21.str().c_str(),hname21.str().c_str(),7, 0.5,7.5);

     stringstream hname22; hname22 << "TM_BX0Qual_" << Whname[iw] << "_MB" << ic << "_Top"  ;
     QualBX0TMMBTop[iw][ic-1]= new TH1F(hname22.str().c_str(),hname22.str().c_str(),7, 0.5,7.5);

     stringstream hname23; hname23 << "TM_BX0Qual_" << Whname[iw] << "_MB" << ic << "_Bottom"  ;
     QualBX0TMMBBottom[iw][ic-1]= new TH1F(hname23.str().c_str(),hname23.str().c_str(),7, 0.5,7.5);

     stringstream hname24; hname24 << "TM_BX0Qual_" << Whname[iw] << "_MB" << ic << "_Vertical"  ;
     QualBX0TMMBVertical[iw][ic-1]= new TH1F(hname24.str().c_str(),hname24.str().c_str(),7, 0.5,7.5);
 }

 for (int iw=0; iw<5;iw++)
  for(int isec=1;isec<13;isec++){
     stringstream hname21; hname21 << "TM_BX0Qual_" << Whname[iw] << "_S" << isec ;
     QualBX0TMSec[iw][isec-1]= new TH1F(hname21.str().c_str(),hname21.str().c_str(),7,0.5,7.5);
 }


 for (int iw=0; iw<5;iw++)
  for (int ins=1; ins<13;ins++)
   for(int ic=1;ic<5;ic++){

    //Changed to move to TM IN (06/10/2016)
    // stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhi/" ;
     stringstream folder1; folder1 << "Wheel" << iw-2 << "/Sector" << ins << "/Station" << ic << "/LocalTriggerPhiIn/" ;
     stringstream hname; hname << "TM_BXvsQual_W" << iw-2 << "_Sec" << ins << "_St" << ic;

     char hnamec[240];
     sprintf(hnamec,"%s%s%s",folder0b.str().c_str(),folder1.str().c_str(),hname.str().c_str());
     TH2F * theHisto= (TH2F*)myFile->Get(hnamec);

     int bin0=-1; // Depending on the range of the minBXTM the bin for which BX=0 is different 
     int nbin0=0;
     int totBins=theHisto->GetYaxis()->GetNbins();
     for(int ib=1;ib<totBins; ib++)
     {
         float cc=theHisto->GetYaxis()->GetBinCenter(ib);
         if(int(cc)==0) { bin0=ib; nbin0=0;}
     }
     if(bin0>0)
     {
       if (nbin0>1)
         for(int jjkk=0;jjkk<20;jjkk++)cout << "[DTDPGCreateWheelSummary] WARNING. Bad range. BX=0 found " << nbin0 << " times in histogram " << hnamec << endl;

       for(int iqualb=1;iqualb<8;iqualb++)
        {
          float qualentries=theHisto->GetBinContent(iqualb, bin0);
          QualBX0TMMB[iw][ic-1]->Fill(float(iqualb),qualentries);   // 
          if(ins>1 && ins<7)QualBX0TMMBTop[iw][ic-1]->Fill(float(iqualb),qualentries);
          if(ins>7)QualBX0TMMBBottom[iw][ic-1]->Fill(float(iqualb),qualentries);
          if(ins==1 || ins==7)QualBX0TMMBVertical[iw][ic-1]->Fill(float(iqualb),qualentries);
        
          QualBX0TMSec[iw][ins-1]->Fill(float(iqualb),qualentries);   // 
        }
     }
     else
       for(int jjkk=0;jjkk<20;jjkk++)cout << "[DTDPGCreateWheelSummary] ERROR. Bad range. BX=0 doesn't found in histogram " << hnamec << endl; 
    
   }
 for(int ic=0;ic<4;ic++)
 {  
   first_to_paint_MB[ic][0]=QualBX0TMMB[0][ic];
   first_to_paint_MB[ic][1]=QualBX0TMMBTop[0][ic];
   first_to_paint_MB[ic][2]=QualBX0TMMBBottom[0][ic];
   first_to_paint_MB[ic][3]=QualBX0TMMBVertical[0][ic];
   float nbmax[4]={0,0,0,0};
   nbmax[0]=QualBX0TMMB[0][ic]->GetMaximum();
   nbmax[1]=QualBX0TMMBTop[0][ic]->GetMaximum();
   nbmax[2]=QualBX0TMMBBottom[0][ic]->GetMaximum();
   nbmax[3]=QualBX0TMMBVertical[0][ic]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax0[4]={0,0,0,0};
      nbmax0[0]=QualBX0TMMB[iw][ic]->GetMaximum();
      nbmax0[1]=QualBX0TMMBTop[iw][ic]->GetMaximum();
      nbmax0[2]=QualBX0TMMBBottom[iw][ic]->GetMaximum();
      nbmax0[3]=QualBX0TMMBVertical[iw][ic]->GetMaximum();

      if(nbmax0[0]>nbmax[0])
      {
         nbmax[0]=nbmax0[0]; first_to_paint_MB[ic][0]=QualBX0TMMB[iw][ic];
      }
      if(nbmax0[1]>nbmax[1])
      {
         nbmax[1]=nbmax0[1]; first_to_paint_MB[ic][1]=QualBX0TMMBTop[iw][ic];
      }
      if(nbmax0[2]>nbmax[2])
      {
         nbmax[2]=nbmax0[2]; first_to_paint_MB[ic][2]=QualBX0TMMBBottom[iw][ic];
      }
      if(nbmax0[3]>nbmax[3])
      {
         nbmax[3]=nbmax0[3]; first_to_paint_MB[ic][3]=QualBX0TMMBVertical[iw][ic];
      }
   }
 }

 for(int ic=1;ic<5;ic++)
 {
   stringstream hname21; hname21 << "TM_BX0Qual_MB" << ic ;
   stringstream hname22; hname22 << "TM_BX0Qual_MB" << ic << "_Top"  ;
   stringstream hname23; hname23 << "TM_BX0Qual_MB" << ic << "_Bottom"  ;
   stringstream hname24; hname24 << "TM_BX0Qual_MB" << ic << "_Vertical"  ;

   first_to_paint_MB[ic-1][0]->SetTitle(hname21.str().c_str());
   first_to_paint_MB[ic-1][1]->SetTitle(hname22.str().c_str());
   first_to_paint_MB[ic-1][2]->SetTitle(hname23.str().c_str());
   first_to_paint_MB[ic-1][3]->SetTitle(hname24.str().c_str());

   for(int ityp=0;ityp<4;ityp++)
   {
     first_to_paint_MB[ic-1][ityp]->GetXaxis()->SetLabelSize(0.07);
     first_to_paint_MB[ic-1][ityp]->GetYaxis()->SetLabelSize(0.05);
     first_to_paint_MB[ic-1][ityp]->SetNdivisions(505);
   }
 }

 char QualLabel[7][3]={"LI","LO","HI","HO","LL","HL","HH"};
 TCanvas *QualBX0TMAll = new TCanvas("QualBX0TMAll", "",201,81,999,950);
 QualBX0TMAll->Divide(4,4) ;
 for(int ip=1;ip<17;ip++)
 {
   QualBX0TMAll->cd(ip)->SetFillColor(0) ;
   QualBX0TMAll->cd(ip)->SetFrameFillColor(0) ;
 }
 for(int ic=1;ic<5;ic++)
 {
   for(int ib=1;ib<8;ib++)
   {
    first_to_paint_MB[ic-1][0]->GetXaxis()->SetBinLabel(ib,QualLabel[ib-1]);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetBinLabel(ib,QualLabel[ib-1]);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetBinLabel(ib,QualLabel[ib-1]);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetBinLabel(ib,QualLabel[ib-1]);
   }

    first_to_paint_MB[ic-1][0]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][1]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][2]->GetXaxis()->SetLabelSize(0.12);
    first_to_paint_MB[ic-1][3]->GetXaxis()->SetLabelSize(0.12);

    QualBX0TMAll->cd((ic-1)*4+1)->SetFillColor(0); first_to_paint_MB[ic-1][0]->Draw("histo");
    QualBX0TMAll->cd((ic-1)*4+2)->SetFillColor(0); first_to_paint_MB[ic-1][1]->Draw("histo");
    QualBX0TMAll->cd((ic-1)*4+3)->SetFillColor(0); first_to_paint_MB[ic-1][2]->Draw("histo");
    QualBX0TMAll->cd((ic-1)*4+4)->SetFillColor(0); first_to_paint_MB[ic-1][3]->Draw("histo");

    for(int iw=0;iw<5;iw++)
    {
      QualBX0TMMB[iw][ic-1]->SetLineColor(iw+1);
      QualBX0TMMBTop[iw][ic-1]->SetLineColor(iw+1);
      QualBX0TMMBBottom[iw][ic-1]->SetLineColor(iw+1);
      QualBX0TMMBVertical[iw][ic-1]->SetLineColor(iw+1);
      if(iw==4)
       {
         QualBX0TMMB[iw][ic-1]->SetLineColor(6);
         QualBX0TMMBTop[iw][ic-1]->SetLineColor(6);
         QualBX0TMMBBottom[iw][ic-1]->SetLineColor(6);
         QualBX0TMMBVertical[iw][ic-1]->SetLineColor(6);
       }
      QualBX0TMAll->cd((ic-1)*4+1);  QualBX0TMMB[iw][ic-1]->Draw("histosame");
      QualBX0TMAll->cd((ic-1)*4+2);  QualBX0TMMBTop[iw][ic-1]->Draw("histosame");
      QualBX0TMAll->cd((ic-1)*4+3);  QualBX0TMMBBottom[iw][ic-1]->Draw("histosame");
      QualBX0TMAll->cd((ic-1)*4+4);  QualBX0TMMBVertical[iw][ic-1]->Draw("histosame");
    }
 }

 createGifFile("QualBX0AllTM",QualBX0TMAll);
 delete QualBX0TMAll;


 TCanvas *QualBX0TMWh[5];
 QualBX0TMWh[0] = new TCanvas("QualBX0TMWh-2", "",201,81,999,950);
 QualBX0TMWh[1] = new TCanvas("QualBX0TMWh-1", "",201,81,999,950);
 QualBX0TMWh[2] = new TCanvas("QualBX0TMWh0", "",201,81,999,950);
 QualBX0TMWh[3] = new TCanvas("QualBX0TMWh+1", "",201,81,999,950);
 QualBX0TMWh[4] = new TCanvas("QualBX0TMWh+2", "",201,81,999,950);
 for(int iw=0;iw<5;iw++)
 {
   QualBX0TMWh[iw]->Divide(3,4) ;
   for(int isec=1;isec<13;isec++)
   {
      QualBX0TMWh[iw]->cd(isec);
      QualBX0TMWh[iw]->cd(isec)->SetFillColor(0);
      QualBX0TMWh[iw]->cd(isec)->SetFrameFillColor(0);
      QualBX0TMSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.10);
      QualBX0TMSec[iw][isec-1]->GetYaxis()->SetLabelSize(0.06);
      QualBX0TMSec[iw][isec-1]->SetNdivisions(505);
      for(int ib=1;ib<8;ib++)
         QualBX0TMSec[iw][isec-1]->GetXaxis()->SetBinLabel(ib,QualLabel[ib-1]);
      QualBX0TMSec[iw][isec-1]->GetXaxis()->SetLabelSize(0.12);
      QualBX0TMSec[iw][isec-1]->Draw("histo");   // 
   }
 }
 for(int iw=0;iw<5;iw++)
 {
   createGifFile("QualBX0SecTM",QualBX0TMWh[iw],(iw-2),false);
   delete QualBX0TMWh[iw];
 }

 TCanvas *QualBX0TMWhAll;
 QualBX0TMWhAll = new TCanvas("QualBX0TMWhAll", "",201,81,999,950);
 QualBX0TMWhAll->Divide(3,4) ;
 for(int ins=1;ins<13;ins++)
 {
   first_to_paint=QualBX0TMSec[0][ins-1];
   float nbmax1=QualBX0TMSec[0][ins-1]->GetMaximum();
   for(int iw=1;iw<5;iw++)
   {
      float nbmax2=QualBX0TMSec[iw][ins-1]->GetMaximum();
      if(nbmax2>nbmax1)
      {
         nbmax1=nbmax2; first_to_paint=QualBX0TMSec[iw][ins-1];
      }
    }

    QualBX0TMWhAll->cd(ins) ;
    QualBX0TMWhAll->cd(ins)->SetFillColor(0) ;
    QualBX0TMWhAll->cd(ins)->SetFrameFillColor(0) ;
    stringstream hname; hname << "TM_BX0Qual_S" << ins;
    first_to_paint->SetTitle(hname.str().c_str());
    first_to_paint->Draw("histo");
    for(int iw=0;iw<5;iw++)
    {

      QualBX0TMSec[iw][ins-1]->SetLineColor(iw+1);
      if(iw==4) QualBX0TMSec[iw][ins-1]->SetLineColor(6);
      QualBX0TMSec[iw][ins-1]->GetXaxis()->SetLabelSize(0.12);
      QualBX0TMSec[iw][ins-1]->Draw("histosame");
    }
 }
 

 createGifFile("QualBX0SecAllTM",QualBX0TMWhAll);
 delete QualBX0TMWhAll;


  // delete TCanvases
  delete c0;
  delete ct2;
  delete ct3;
  delete ct6;
  delete ct66;
  for (int i=0;i<2;i++) delete ct7[i];

  if(ProcessDDUTrigger)
     for (int i=0;i<4;i++) delete ct5[i];
   

  for (int i=0;i<=4;i++) {
    delete c1[i];
    delete ct4[i];
  }


  if (myFile) { 
    myFile->Close();
    delete myFile;
  }
  edm::LogVerbatim ("DTDPGSummary") << " closing txt dead channel file ";
  DeadChannelList->close();


}
