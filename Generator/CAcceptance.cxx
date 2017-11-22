#ifndef CACCEPTANCE_CXX
#define CACCEPTANCE_CXX

#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include "CAcceptance.h"
#include "../ManualFunctions.h"

using namespace std;

CAcceptance::CAcceptance (boolFunction func) {
	initFlag_ = 0;
	SetInputFile ("../Generated.root");
	SetOutputFile ("../Acceptance.root");
	SetAcceptanceFunction (func);
}


CAcceptance::CAcceptance (boolFunction func, TString inputFileName, TString outputFileName) {
	initFlag_ = 0;
	SetInputFile (inputFileName);
	SetOutputFile (outputFileName);
	SetAcceptanceFunction (func);
}


void CAcceptance::SetInputFile (TString inputFileName) {
	//if (inputFile_ != 0) delete inputFile_;
	inputFile_ = new TFile (inputFileName, "READ");
}


void CAcceptance::SetOutputFile (TString outputFileName) {
	//if (outputFile_ != 0) delete outputFile_;
	outputFile_ = new TFile (outputFileName, "RECREATE");
}


void CAcceptance::SetAcceptanceFunction (boolFunction func) {
	AcceptanceFunction_ = func;
}


void CAcceptance::Init (){
	initFlag_ = 1;
	gRandom -> SetSeed(0);
	event_ = new CEvent;
	inputTree_ = (TTree*) inputFile_ -> Get ("Tree");
	inputTree_ -> SetBranchAddress ("Event", &event_);
	outputTree_ = new TTree ("Tree", "Events with acceptance applied");
	outputTree_ -> Branch ("Event", &event_, 128000, 4);
	outputFile_ -> cd();
}


void CAcceptance::ProcessTree () {
	if (initFlag_ == 0)
		Init ();
	Int_t mh;
	Long64_t nentries = inputTree_ -> GetEntries ();
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        cout << "\rEvent " << jentry + 1 << " from " << nentries;
		inputTree_ -> GetEntry(jentry);
		mh = event_ -> GetMh ();

		for (Int_t itrack = 1; itrack <= mh; itrack++) {
			track_ = event_ -> GetTrack (itrack);
			if (AcceptanceFunction_ (*track_) == 0){
				event_ -> RemoveTrack (itrack);
				mh--;
				itrack--;
			}
		}
		event_ -> SetCent (GetCentralityClass(mh));

		outputTree_ -> Fill ();
	}

	Finish ();
}


void CAcceptance::Finish () {
	outputFile_ -> Write ();
	outputFile_ -> Close ();
}

#endif // CACCEPTANCE_CXX
