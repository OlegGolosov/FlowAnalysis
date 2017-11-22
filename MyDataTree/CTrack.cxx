#ifndef CTRACK_CXX
#define CTRACK_CXX

#include "CTrack.h"
#include <iostream>

using namespace std;

CTrack::CTrack ()
	:	TObject (){
}


CTrack::CTrack (Float_t pt_, Float_t eta_, Float_t phi_, Int_t charge_, Int_t pid_)
	:	TObject (),
		pt (pt_),
		eta (eta_),
		phi (phi_),
		charge (charge_),
		pid (pid_) {
}


CTrack::~CTrack () {
	//delete this;
}


void CTrack::Set (Float_t pt_, Float_t eta_, Float_t phi_, Int_t charge_, Int_t pid_) {
	pt = pt_;
	eta = eta_;
	phi = phi_;
	pid = pid_;
	charge = charge_;
}

void CTrack::SetP (Float_t p_) {
	p = p_;
}

void CTrack::SetPt (Float_t pt_) {
	pt = pt_;
}

void CTrack::SetEta (Float_t eta_) {
	eta = eta_;
}

void CTrack::SetRap (Float_t rap_) {
	rap = rap_;
}

void CTrack::SetPhi (Float_t phi_) {
	phi = phi_;
}

void CTrack::SetPid (Int_t pid_) {
	pid = pid_;
}

void CTrack::SetCharge (Int_t charge_) {
	charge = charge_;
}

void CTrack::SetdEdx_full (Float_t dEdx_full_) {
	dEdx_full = dEdx_full_;
}

void CTrack::SetdEdx_MTPC (Float_t dEdx_MTPC_) {
	dEdx_MTPC = dEdx_MTPC_;
}

Float_t CTrack::GetP () {
	return p;
}

Float_t CTrack::GetPt () {
	return pt;
}

Float_t CTrack::GetEta () {
	return eta;
}

Float_t CTrack::GetRap () {
	return rap;
}

Float_t CTrack::GetPhi () {
	return phi;
}

Int_t CTrack::GetPid () {
	return pid;
}

Int_t CTrack::GetCharge () {
	return charge;
}

Float_t CTrack::GetdEdx_full () {
	return dEdx_full;
}

Float_t CTrack::GetdEdx_MTPC () {
	return dEdx_MTPC;
}

ClassImp(CTrack);

#endif // CTRACK_CXX
