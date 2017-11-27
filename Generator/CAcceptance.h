#ifndef CACCEPTANCE_H
#define CACCEPTANCE_H

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "../MyDataTree/CEvent.h"
#include "../config.h"

class CAcceptance {

public:
	CAcceptance (boolFunction func);
	CAcceptance (boolFunction func, TString inputFileName, TString outputFileName);
	void SetInputFile (TString inputFileName = "Generated.root");
	void SetOutputFile (TString outputFileName = "Acceptance.root");
	void SetAcceptanceFunction (boolFunction func);
	void Init ();
	void ProcessTree ();
	void Finish ();

private:
	Int_t nEvents_;
	TFile* outputFile_;
	TFile* inputFile_;
	TTree* inputTree_;
	TTree* outputTree_;
	CEvent* event_;
	CTrack* track_;
	Bool_t initFlag_;

	boolFunction AcceptanceFunction_;
};

#endif // CACCEPTANCE_H
