#ifndef CEVENT_CXX
#define CEVENT_CXX

#include <iostream>
#include <TString.h>
#include "CEvent.h"

using namespace std;

CEvent::CEvent () : TObject(), mh (0), tracks (new TClonesArray ("CTrack")) {
    SetNrun (0);
	for (Int_t i = 0; i < MAXNHARMONICS; i++) {
		psi [i] = 0.0;
	}
}


CEvent::~CEvent () {
	Clear ();
	delete tracks;
	//delete psi_;
}


void CEvent::SetNrun (Int_t nRun_) {
	nRun = nRun_;
}

void CEvent::SetEvetoFull (Float_t EvetoFull_) {
	EvetoFull = EvetoFull_;
}

void CEvent::SetEveto (Float_t Eveto_ [4]) {
    for (Int_t i = 0; i < 4; i++) {
        Eveto [i] = Eveto_ [i];
    }
}


void CEvent::SetCent (Float_t cent_) {
	cent = cent_;
}


void CEvent::AddTrack (Float_t pt, Float_t eta, Float_t phi, Int_t charge, Int_t pid) {
	new ((*tracks)[mh]) CTrack (pt, eta, phi, charge, pid);
	mh += 1;
}

void CEvent::Clear () {
	//delete psi_;
	tracks -> Clear ();
	nHarmonics = 0;
	mh = 0;
	cent = 0.0;
}


void CEvent::SetPsi_n (Int_t n, Float_t psi_n) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting psi_%i: %i > MAXNHARMONICS", n, n) << endl;
		psi [n - 1] = psi_n;
		if (n > nHarmonics)
			nHarmonics ++;
}

void CEvent::SetXa (Int_t n, Float_t Xa_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Xa_%i: %i > MAXNHARMONICS", n, n) << endl;
		Xa [n - 1] = Xa_;
		if (n > nHarmonics)
			nHarmonics ++;
}

void CEvent::SetXb (Int_t n, Float_t Xc_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Xa_%i: %i > MAXNHARMONICS", n, n) << endl;
		Xa [n - 1] = Xc_;
}

void CEvent::SetXc (Int_t n, Float_t Xc_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Xa_%i: %i > MAXNHARMONICS", n, n) << endl;
		Xa [n - 1] = Xc_;
}

void CEvent::SetYa (Int_t n, Float_t Ya_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Ya_%i: %i > MAXNHARMONICS", n, n) << endl;
		Ya [n - 1] = Ya_;
}

void CEvent::SetYb (Int_t n, Float_t Yc_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Ya_%i: %i > MAXNHARMONICS", n, n) << endl;
		Ya [n - 1] = Yc_;
}

void CEvent::SetYc (Int_t n, Float_t Yc_) {
		if (n > MAXNHARMONICS)
			cout << Form ("Error setting Ya_%i: %i > MAXNHARMONICS", n, n) << endl;
		Ya [n - 1] = Yc_;
}

Int_t CEvent::GetNHarmonics () {
	return nHarmonics;
}

Float_t CEvent::GetPsi_n (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting psi %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return psi [n - 1];
}

Float_t CEvent::GetXa (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Xa %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Xa [n - 1];
}

Float_t CEvent::GetXb (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Xb %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Xb [n - 1];
}

Float_t CEvent::GetXc (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Xc %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Xc [n - 1];
}

Float_t CEvent::GetYa (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Ya %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Ya [n - 1];
}

Float_t CEvent::GetYb (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Yb %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Yb [n - 1];
}

Float_t CEvent::GetYc (Int_t n) {
	if (n > nHarmonics) {
		cout << endl << Form ("Error getting Yc %i: %i > nHarmonics!", n, n) << endl;
		return - 1.0;
	}
	return Yc [n - 1];
}

Int_t CEvent::GetNrun () {
	return nRun;
}

Int_t CEvent::GetMh () {
	return mh;
}

Float_t CEvent::GetCent (){
	return cent;
}

Float_t CEvent::GetEvetoFull () {
	return EvetoFull;
}

Float_t* CEvent::GetEveto () {
	return Eveto;
}


CTrack* CEvent::GetTrack (Int_t n) const {
	if (n < 0 || n > mh) {
		cout << endl << Form ("Error getting track: %i < 0 or %i > mh!", n, n) << endl;
		return 0;
	}
	return ((CTrack*) tracks -> AddrAt(n - 1));
}

void CEvent::RemoveTrack (Int_t i) {
	tracks -> RemoveAt (i - 1);
	tracks -> Compress ();
	mh--;
}

ClassImp(CEvent);

#endif // CEVENT_CXX
