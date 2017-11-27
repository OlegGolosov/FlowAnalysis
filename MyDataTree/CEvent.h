#ifndef CEVENT_H
#define CEVENT_H

#include <TObject.h>
#include <TClonesArray.h>
#include <vector>
#include "CTrack.h"
#include "../config.h"

using namespace std;

class CEvent : public TObject {
	friend class CAcceptance;
public:
	CEvent();
	virtual ~CEvent();
	void SetPsi_n (Int_t n, Float_t psi_n);
	void SetCent (Float_t cent);
	void SetNrun (Int_t nRun);
	void SetEvetoFull (Float_t EvetoFull);
	void SetEveto (Float_t Eveto [4]);
	void SetXa (Int_t n, Float_t Xa_);
	void SetXb (Int_t n, Float_t Xb_);
	void SetXc (Int_t n, Float_t Xc_);
	void SetYa (Int_t n, Float_t Ya_);
	void SetYb (Int_t n, Float_t Yb_);
	void SetYc (Int_t n, Float_t Yc_);
	void AddTrack (Float_t pt, Float_t eta, Float_t phi, Int_t charge, Int_t pid);
	void Clear ();

	Int_t GetNrun ();
	Int_t GetNHarmonics ();
	Float_t GetPsi_n (Int_t n);
	Int_t GetMh ();
	Float_t GetCent ();
	Float_t GetXa (Int_t n);
	Float_t GetXb (Int_t n);
	Float_t GetXc (Int_t n);
	Float_t GetYa (Int_t n);
	Float_t GetYb (Int_t n);
	Float_t GetYc (Int_t n);
	Float_t GetEvetoFull ();
	Float_t* GetEveto ();
	CTrack* GetTrack (Int_t n) const;
	void RemoveTrack (Int_t i);

private:
	Int_t nHarmonics;
	Int_t nRun;
	Int_t mh;
	Float_t EvetoFull;
	Float_t Eveto [4];
	Float_t cent;
	Float_t Xa [MAXNHARMONICS], Xb [MAXNHARMONICS], Xc [MAXNHARMONICS];
	Float_t Ya [MAXNHARMONICS], Yb [MAXNHARMONICS], Yc [MAXNHARMONICS];
	Float_t psi [MAXNHARMONICS];
	TClonesArray* tracks; //->


	ClassDef (CEvent, 1);
};


#endif // CEVENT_H
