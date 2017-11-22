#ifndef CEVENTGENERATOR_CXX
#define CEVENTGENERATOR_CXX

#include <cassert>
#include <TMath.h>
#include <TRandom3.h>
#include <TTime.h>
#include <iostream>
#include "CEventGenerator.h"

using namespace std;

CEventGenerator::CEventGenerator (Int_t nHarmonics) {
	nHarmonics_ = nHarmonics;
	variables_ ["nHarmonics"] = nHarmonics;
	initFlag_ = 0;
	activeHarmonics_ = new Bool_t [nHarmonics];
	v_ = new Float_t [nHarmonics];
	psi_ = new Float_t [nHarmonics];
	HarmonicFunctions = new floatFunction [nHarmonics];

	for (Int_t i = 0; i < nHarmonics; i++){
		activeHarmonics_ [i] = 0;
		HarmonicFunctions [i] = 0;
	}

	gRandom -> SetSeed(0);
}


CEventGenerator::~CEventGenerator () {
	delete activeHarmonics_;
	delete v_;
	delete psi_;
}


void CEventGenerator::Init () {
	initFlag_ = 1;
}


void CEventGenerator::SetHarmonicFunction (Int_t n, floatFunction func){
	HarmonicFunctions [n - 1] = func;
	activeHarmonics_ [n - 1] = 1;
}

Float_t CEventGenerator::rho (Float_t phi) {
    Float_t f = 0.5;
    for (Int_t i = 0; i < nHarmonics_; i++) {
        f += v_ [i] * cos ((i + 1) * (phi - psi_ [i]));
    }
    if (f < 0) {
        cout << "Invalid flow coefficients: rho < 0\n";
        assert (0);
    }
    f /= PI;
   return f;
}

CEvent *CEventGenerator::GenerateEvent () {

    //INITIATION CHECK
    if (initFlag_ == 0) {
		cout << endl;
		cout << "=================================================" << endl;
		cout << "         Please, initiate the generator!         " << endl;
		cout << "=================================================" << endl;
		return 0;
	}

	Int_t mh, charge;
	Float_t pt, eta, phi, pid;

	event_.Clear ();
	GetPsi (variables_);
	GetEventVariables (variables_);
	event_.SetCent (variables_ ["cent"]);
	for (Int_t i = 0; i < nHarmonics_; i++) {
		if (activeHarmonics_ [i] == 1) {
			psi_ [i] = variables_ [Form ("psi_%i", i + 1)];
		}
		else psi_ [i] = 0.0;
		event_.SetPsi_n (i + 1, psi_ [i]);
	}

	mh = static_cast <int> (variables_ ["mh"]);

	for (Int_t j = 0; j < mh; j++) {
        GetTrackVariables (variables_);
        for (Int_t i = 0; i < nHarmonics_; i++) {
            if (activeHarmonics_ [i] == 1) {
                v_ [i] = HarmonicFunctions [i] (variables_);
            }
        }
		phi = gRandom -> Rndm () * 2 * PI;
		if (gRandom -> Rndm () < rho (phi))	{
            pt = variables_ ["pt"];
            eta = variables_ ["eta"];
            pid = variables_ ["pid"];
            charge = variables_ ["charge"];
            event_.AddTrack (pt, eta, phi, charge, pid);
		}
		else j--;
	}

	return &event_;
}

void CEventGenerator::SetPsiFunction (voidFunction func) {
	GetPsi = func;
}
void CEventGenerator::SetEventFunction (voidFunction func) {
	GetEventVariables = func;
}
void CEventGenerator::SetTrackFunction (voidFunction func) {
	GetTrackVariables = func;
}

#endif // CEVENTGENERATOR_CXX
