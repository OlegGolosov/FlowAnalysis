#ifndef CTREEBUILDER_CXX
#define CTREEBUILDER_CXX

#include <iostream>
#include "CTreeBuilder.h"

using namespace std;

CTreeBuilder::CTreeBuilder () {
	initFlag_ = 0;
	SetOutputFileName ("../Generated.root");
}


CTreeBuilder::CTreeBuilder (TString outputFileName) {
	initFlag_ = 0;
	SetOutputFileName (outputFileName);
}


void CTreeBuilder::SetOutputFileName (TString outputFileName) {
	initFlag_ = 0;
	outputFileName_ = outputFileName;
} 


void CTreeBuilder::Init (){
	initFlag_ = 1;
	outputFile_ = new TFile (outputFileName_, "RECREATE");
	outputFile_ -> cd();
	tree_ = new TTree ("Tree", "Generated Events Tree");	
	event_ = new CEvent;
	tree_ -> Branch ("Event", &event_, 128000, 4);
}


void CTreeBuilder::AddEvent (CEvent* ev) {
	event_ = ev;
	tree_ -> Fill ();
}


void CTreeBuilder::Finish () {
	outputFile_ -> cd ();
	tree_ -> Write ();
	outputFile_ -> Close ();
}

#endif // CTREEBUILDER_CXX