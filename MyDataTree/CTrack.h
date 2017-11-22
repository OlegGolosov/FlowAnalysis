#ifndef CTRACK_H
#define CTRACK_H

#include <TObject.h>

class CTrack : public TObject {
public:
	CTrack ();
	CTrack (Float_t pt, Float_t eta, Float_t phi, Int_t charge, Int_t pid_);
	virtual ~CTrack ();
	void SetP (Float_t p_);
	void SetPt (Float_t pt_);
	void SetEta (Float_t eta_);
	void SetRap (Float_t rap_);
	void SetPhi (Float_t phi_);
	void SetPid (Int_t pid_);
	void SetCharge (Int_t charge_);
	void SetdEdx_full (Float_t dEdx_full_);
	void SetdEdx_MTPC (Float_t dEdx_MTPC_);
	void Set (Float_t pt, Float_t eta, Float_t phi, Int_t charge, Int_t pid_);
	Float_t GetP ();
	Float_t GetPt ();
	Float_t GetEta ();
	Float_t GetRap ();
	Float_t GetPhi ();
	Int_t GetPid ();
	Int_t GetCharge ();
	Float_t GetdEdx_full ();
	Float_t GetdEdx_MTPC ();

private:
	Float_t p;
	Float_t pt;
	Float_t eta;
	Float_t rap;
	Float_t phi;
	Int_t pid;
	Int_t charge;
	Float_t dEdx_full;
	Float_t dEdx_MTPC;

	ClassDef(CTrack, 1);
};



#endif // CTRACK_H
