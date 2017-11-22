#ifndef CTREEBUILDER_H
#define CTREEBUILDER_H

#include <TMath.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "../MyDataTree/CEvent.h"

using namespace std;

class CTreeBuilder {

public:
	CTreeBuilder ();
	CTreeBuilder (TString outputFileName);
	void SetOutputFileName (TString outputFileName = "../Generated.root");
	void AddEvent (CEvent* ev);
	void Init ();
	void Finish ();

private:
	TString outputFileName_;
	TFile* outputFile_;
	TTree* tree_;
	CEvent* event_;
	Bool_t initFlag_;
};

#endif // CTREEBUILDER_H
