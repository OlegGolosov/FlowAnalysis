#ifndef CTREECONVERTERHADES_CXX
#define CTREECONVERTERHADES_CXX

#include <iostream>
#include <TMath.h>
#include <TH1F.h> // QA
#include <TH2F.h> // QA
#include <TProfile.h> // QA
#include <THStack.h> // QA
#include <TCanvas.h> // QA
#include <TLegend.h> // QA
#include <TStyle.h> // QA
#include <TROOT.h> // QA
#include "CTreeConverterHADES.h"

using namespace std;

CTreeConverterHADES::CTreeConverterHADES () {
	initFlag_ = 0;
}


CTreeConverterHADES::CTreeConverterHADES (TString inputFileName, TString outputFileName, TString inputTreeName) {
	initFlag_ = 0;
	SetInputFileName (inputFileName);
	SetOutputFileName (outputFileName);
	SetInputTreeName (inputTreeName);
	SetdEdxSource ();
	SetdEdxSource ();
}

void CTreeConverterHADES::SetInputFileName (TString inputFileName) {
    inputFileName_ = inputFileName;
}

void CTreeConverterHADES::SetOutputFileName (TString outputFileName) {
    outputFileName_ = outputFileName;
}

void CTreeConverterHADES::SetInputTreeName (TString inputTreeName) {
    inputTreeName_ = inputTreeName;
}

Bool_t CTreeConverterHADES::SetInputFile () {
	//if (inputFile_) delete inputFile_;
	inputFile_ = new TFile (inputFileName_ + ".root", "READ");
	if (!inputFile_) {
		cout << "Input file '" << inputFileName_ << "' not found!!!\n";
		return 0;
	}
	return 1;
}


Bool_t CTreeConverterHADES::SetOutputFile () {
	//if (outputFile_) delete outputFile_;
	outputFile_ = new TFile (outputFileName_ + ".root", "RECREATE");
	if (!outputFile_) {
		cout << "Output file '" << outputFileName_ << "' is not accessible!!!\n";
		return 0;
	}
	return 1;
}


Bool_t CTreeConverterHADES::SetInputTree () {
	//if (inputTree_) delete inputTree_;
	inputTree_ = (TTree*) inputFile_ -> Get (inputTreeName_);
	if (!inputTree_) {
		cout << "Input tree '" << inputTreeName_ << "' not found!!!\n";
		return 0;
	}
	return 1;
}

void CTreeConverterHADES::SetMhRange (Int_t mhMin, Int_t mhMax) {
    mhMin_ = mhMin;
    mhMax_ = mhMax;
}

void CTreeConverterHADES::SetCentRange (Float_t centMin, Float_t centMax) {
    centMin_ = centMin;
    centMax_ = centMax;
}

void CTreeConverterHADES::SetNbinsMh (Int_t nBinsMh) {
	nBinsMh_ = nBinsMh;
}

void CTreeConverterHADES::SetNbinsCent (Int_t nBinsCent) {
	nBinsCent_ = nBinsCent;
}

void CTreeConverterHADES::SetdEdxSource (Int_t TPCid) {
    dEdxSource_ = TPCid;
}

void CTreeConverterHADES::SetCentralityMethod (Int_t centMethod) {
    centMethod_ = centMethod;
}

Bool_t CTreeConverterHADES::Init (){
	if (!SetInputFile ()) return 0;
	if (!SetOutputFile ()) return 0;
	if (!SetInputTree ()) return 0;

   inputTree_->SetBranchAddress("nRpcClust", &nRpcClust, &b_nRpcClust);
   inputTree_->SetBranchAddress("nRpcClustCut", &nRpcClustCut, &b_nRpcClustCut);
   inputTree_->SetBranchAddress("nRpcHits", &nRpcHits, &b_nRpcHits);
   inputTree_->SetBranchAddress("nRpcHitsCut", &nRpcHitsCut, &b_nRpcHitsCut);
   inputTree_->SetBranchAddress("nTofHits", &nTofHits, &b_nTofHits);
   inputTree_->SetBranchAddress("nTofHitsCut", &nTofHitsCut, &b_nTofHitsCut);
   inputTree_->SetBranchAddress("primaryTracks", &primaryTracks, &b_primaryTracks);
   inputTree_->SetBranchAddress("selectedTracks", &selectedTracks, &b_selectedTracks);
   inputTree_->SetBranchAddress("psiEPa", &psiEPa, &b_psiEPa);
   inputTree_->SetBranchAddress("psiEPb", &psiEPb, &b_psiEPb);
   inputTree_->SetBranchAddress("nA", &nA, &b_nA);
   inputTree_->SetBranchAddress("nB", &nB, &b_nB);
   inputTree_->SetBranchAddress("trigInd", trigInd, &b_trigInd);
   inputTree_->SetBranchAddress("runId", &runId, &b_runId);
   inputTree_->SetBranchAddress("nWallHitsTot", &nWallHitsTot, &b_nWallHitsTot);
   inputTree_->SetBranchAddress("runTime", &runTime, &b_runTime);
   inputTree_->SetBranchAddress("cuts", cuts, &b_cuts);
   inputTree_->SetBranchAddress("wallModuleIndex", wallModuleIndex, &b_wallModuleIndex);
   inputTree_->SetBranchAddress("wallHitTime", wallHitTime, &b_wallHitTime);
   inputTree_->SetBranchAddress("wallHitCharge", wallHitCharge, &b_wallHitCharge);
   inputTree_->SetBranchAddress("wallHitDistance", wallHitDistance, &b_wallHitDistance);
   inputTree_->SetBranchAddress("wallHitRing", wallHitRing, &b_wallHitRing);
   inputTree_->SetBranchAddress("wallHitPhi", wallHitPhi, &b_wallHitPhi);
   inputTree_->SetBranchAddress("wallHitX", wallHitX, &b_wallHitX);
   inputTree_->SetBranchAddress("wallHitY", wallHitY, &b_wallHitY);
   inputTree_->SetBranchAddress("wallHitZ", wallHitZ, &b_wallHitZ);
   inputTree_->SetBranchAddress("wallHitEta", wallHitEta, &b_wallHitEta);
   inputTree_->SetBranchAddress("isWallHitOk", isWallHitOk, &b_isWallHitOk);
   inputTree_->SetBranchAddress("wallChargeTot_mod", wallChargeTot_mod, &b_wallChargeTot_mod);
   inputTree_->SetBranchAddress("wallChargeTot", &wallChargeTot, &b_wallChargeTot);
   inputTree_->SetBranchAddress("wallChargeTot_ring", wallChargeTot_ring, &b_wallChargeTot_ring);
   inputTree_->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   inputTree_->SetBranchAddress("nProtons", &nProtons, &b_nProtons);
   inputTree_->SetBranchAddress("nPiP", &nPiP, &b_nPiP);
   inputTree_->SetBranchAddress("nPiM", &nPiM, &b_nPiM);
   inputTree_->SetBranchAddress("vX", &vX, &b_vX);
   inputTree_->SetBranchAddress("vY", &vY, &b_vY);
   inputTree_->SetBranchAddress("vZ", &vZ, &b_vZ);
   inputTree_->SetBranchAddress("vChi2", &vChi2, &b_vChi2);
   inputTree_->SetBranchAddress("time", &time, &b_time);
   inputTree_->SetBranchAddress("pid", pid_, &b_pid);
   inputTree_->SetBranchAddress("p", p_, &b_p);
   inputTree_->SetBranchAddress("phi", phi_, &b_phi);
   inputTree_->SetBranchAddress("theta", theta_, &b_theta);
   inputTree_->SetBranchAddress("pt", pt_, &b_pt);
   inputTree_->SetBranchAddress("rapidity", rapidity, &b_rapidity);
   inputTree_->SetBranchAddress("eta", eta_, &b_eta);
   inputTree_->SetBranchAddress("metaBeta", metaBeta, &b_metaBeta);
   inputTree_->SetBranchAddress("metaMass", metaMass, &b_metaMass);
   inputTree_->SetBranchAddress("charge", charge_, &b_charge);
   inputTree_->SetBranchAddress("mdcdEdx", mdcdEdx, &b_mdcdEdx);
   inputTree_->SetBranchAddress("tofdEdx", tofdEdx, &b_tofdEdx);
   inputTree_->SetBranchAddress("DCAxy", DCAxy, &b_DCAxy);
   inputTree_->SetBranchAddress("DCAz", DCAz, &b_DCAz);
   inputTree_->SetBranchAddress("mdcNhits", mdcNhits, &b_mdcNhits);
   inputTree_->SetBranchAddress("mdcNhitsInner", mdcNhitsInner, &b_mdcNhitsInner);
   inputTree_->SetBranchAddress("mdcNhitsOuter", mdcNhitsOuter, &b_mdcNhitsOuter);
   inputTree_->SetBranchAddress("chi2all", chi2all, &b_chi2all);
   inputTree_->SetBranchAddress("chi2inner", chi2inner, &b_chi2inner);
   inputTree_->SetBranchAddress("chi2outer", chi2outer, &b_chi2outer);
   inputTree_->SetBranchAddress("metaQ", metaQ, &b_metaQ);
   inputTree_->SetBranchAddress("metaMatchRadius", metaMatchRadius, &b_metaMatchRadius);
   inputTree_->SetBranchAddress("pCorr", pCorr, &b_pCorr);
   inputTree_->SetBranchAddress("pt_corr", pt_corr, &b_pt_corr);
   inputTree_->SetBranchAddress("rapidity_corr", rapidity_corr, &b_rapidity_corr);
   inputTree_->SetBranchAddress("metaDx", metaDx, &b_metaDx);
   inputTree_->SetBranchAddress("metaDy", metaDy, &b_metaDy);
   inputTree_->SetBranchAddress("mdcSecId", mdcSecId, &b_mdcSecId);

	return 1;
}


Bool_t CTreeConverterHADES::CheckEventCuts () {
    Int_t timeBin = averages_ [0] -> FindBin (runTime);
    vX -= averages_ [0] -> GetBinContent (timeBin);
    vY -= averages_ [1] -> GetBinContent (timeBin);

    for (Int_t i = 0; i < 8; i++){
    	if (!cuts[i]) return 0;
    }
    if (vZ < -60.0 || vZ > 0.0) return 0;
    if (vChi2 < 0.5 || vChi2 > 40.0) return 0;
    if (TMath::Sqrt (vX * vX + vY * vY) > 3.0) return 0;
    if (trigInd [2] != 1) return 0;
    return 1;
}

Float_t CTreeConverterHADES::GetCentralityClass (Int_t mh) { //Kardan centrality classes
//    10%
//    Double_t xAxis1[8] = {0, 20, 60, 88, 121, 160, 250, 280};
//    static const Int_t nCentClasses = 6;
//    const Float_t centMax = 60.0;
//    Float_t centClassLimits [nCentClasses + 1] = {0, 20, 60, 88, 121, 160, 280};

//    5%
    Double_t xAxis1[12] = {0, 39, 60, 74, 88, 104, 121, 140, 160, 182, 250, 280};
    static const Int_t nCentClasses = 10;
    const Float_t centMax = 50.0;
    Float_t centClassLimits [nCentClasses + 1] = {0, 39, 60, 74, 88, 104, 121, 140, 160, 182, 280};

// Sadovsky
//    static const Int_t nCentClasses = 8;
//    const Float_t centMax = 40.0;
//      Float_t centClassLimits [nCentClasses + 1] = {0, 20, 50, 90, 130, 160, 190, 220, 280};
	Float_t centClassWidth = centMax / (Float_t) nCentClasses;

	for (Int_t i = 0; i < nCentClasses; i++) {
        if (mh > centClassLimits [i] && mh <= centClassLimits [i + 1]) return centClassWidth * (nCentClasses - i - 0.5);
	}
	return -1.0;
}

Float_t CTreeConverterHADES::GetCentralityClass (Float_t Eveto) {
 //    Float_t Ebeam = 8.32e3;
 //    Float_t r = Eveto / Ebeam;
 //    static const Int_t nCentClasses = 6;
	// Float_t centClassLimits [nCentClasses - 1] = {0.169, 0.314, 0.509, 0.66, 0.778};
 //    Float_t centClassWidth = 100.0 / (Float_t) nCentClasses;
	// if (r <= centClassLimits [0]) return centClassWidth * 0.5;
	// for (Int_t i = 1; i < nCentClasses - 1; i++) {
 //        if (r > centClassLimits [i - 1] && r <= centClassLimits [i]) return centClassWidth * (i + 0.5);
	// }
	// if (r > centClassLimits [nCentClasses - 2]) return 100.0 - centClassWidth * 0.5;
	return -1.0;
}

void CTreeConverterHADES::SetSNN (Float_t sNN) {
    sNN_ = sNN;
}

Float_t CTreeConverterHADES::GetRapidity (Float_t pt, Float_t eta, Int_t pid) {
 //    Float_t mProton = 938.2720814;
 //    Float_t m = particleMass [pid];
 //    if (pid == 0) return -999.0;
 //    Float_t yBeam = TMath::Log (sNN_ * 1.0e3 / mProton);
 //    Float_t a = TMath::Sqrt (pt * pt * 1.0e6 * TMath::CosH(eta) * TMath::CosH(eta) + m * m);
 //    Float_t b = pt * 1.0e3 * TMath::SinH (eta);
	// return 0.5 * TMath::Log ((a + b) / (a - b)) - yBeam;
	return -999.;
}

Bool_t CTreeConverterHADES::CheckTrackCuts (Int_t itrack) {
  //Kardan's cuts for tracks
    if (TMath::Abs (DCAz [itrack] - vZ) > 15.0) return 0;
    if (TMath::Abs (DCAxy [itrack]) > 15.0) return 0;
//    if (metaBeta [itrack] > 1.0) return 0;
    return 1;
}


Int_t CTreeConverterHADES::GetTrackPid (Int_t itrack) {
    if (pid_[itrack] == 14) return kProton;
    if (pid_[itrack] == 8) return kPionPlus;
    if (pid_[itrack] == 9) return kPionMinus;
    return kNID;
}

void CTreeConverterHADES::GetAverages () {
    ULong64_t tMin = 1e19, tMax = 0., nBinsTime;

    Long64_t nentries = inputTree_ -> GetEntries ();

    tMin = inputTree_ -> GetMinimum ("runTime") / 1000;
    tMax = inputTree_ -> GetMaximum ("runTime") / 1000 + 1;
    tMin *= 1000;
    tMax *= 1000;
    nBinsTime = (tMax - tMin) / 100;

    averages_.push_back (new TProfile ("pVxTime", "#LTX_{vertex}#GT over time; time; #LTX_{vertex}#GT", nBinsTime, tMin, tMax));
    averages_.push_back (new TProfile ("pVyTime", "#LTY_{vertex}#GT over time; time; #LTY_{vertex}#GT", nBinsTime, tMin, tMax));

    cout << "/nGetting averages:/n";
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		cout << "\rEvent " << jentry + 1 << " from " << nentries;
		inputTree_ -> GetEntry (jentry);
		if (vX > -100. && vY > -100.) {
            averages_ [0] -> Fill (runTime, vX);
            averages_ [1] -> Fill (runTime, vY);
        }
    }

}

Bool_t CTreeConverterHADES::ConvertTree () {
	if (!Init ()) return 0;

	CEvent* event_ = new CEvent;
	outputTree_ = new TTree ("Tree", "Converted tree");
	outputTree_ -> Branch ("Event", &event_, 128000, 4);
	outputFile_ -> cd ();

//    QA list
//    Float_t pMin = 0., pMax = 5.0,ptMin = 0.0, ptMax = 4.0, etaMin = 0.0, etaMax = 2.5, rapMin = -1.0, rapMax = 1.5, tofdEdxMin = -2., tofdEdxMax = 10., mdcdEdxMin = 0.0, mdcdEdxMax = 1.2; // raw
    Float_t pMin = 0., pMax = 5.0, ptMin = 0.0, ptMax = 2.5, etaMin = 0.0, etaMax = 2.1, rapMin = -1.0, rapMax = 1.5, tofdEdxMin = 0., tofdEdxMax = 10., mdcdEdxMin = 0.0, mdcdEdxMax = 15., thetaMin = 0.1, thetaMax = 1.6; // cut

	TFile *histFile = new TFile (outputFileName_ + "_QA.root", "recreate");

    gStyle -> SetOptStat (111111);
    gROOT -> ForceStyle ();


    TH1F *hMDC_hits = new TH1F ("hMDC_hits", "Number of MDC hits; N^{MDC}_{hits}; counts", 12, 14, 25);
    TH1F *hMDC_hits_inner = new TH1F ("hMDC_hits_inner", "Number of MDC hits for inner plates; N^{MDC}_{hits}; counts", 10, 3, 13);
    TH1F *hMDC_hits_outer = new TH1F ("hMDC_hits_outer", "Number of MDC hits for outer plates; N^{MDC}_{hits}; counts", 9, 4, 13);
	TH1F *hMDC_tracks = new TH1F ("hMDC_tracks", "Number of MDC tracks; N_{tracks}^{MDC}; Nevents", 100, 0, 100);
	TH1F *hTOF_RPC_hits = new TH1F ("hTOF_RPC_hits", "Number of TOF + RPC hits; N_{hits}^{TOF} + N_{hits}^{RPC}; Nevents", 280, 0, 280);
	TH1F *hTOF_hits = new TH1F ("hTOF_hits", "Number of TOF hits; N_{hits}^{TOF}; Nevents", 100, 0, 100);
	TH1F *hRPC_hits = new TH1F ("hRPC_hits", "Number of RPC hits; N_{hits}^{RPC}; Nevents", 200, 0, 200);
	TH1F *hFW_hits = new TH1F ("hFW_hits", "Number of FW hits; N_{hits}^{FW}; Nevents", 160, 0, 160);
	TH1F *hFW_hits_cut = new TH1F ("hFW_hits_cut", "Number of FW hits (with cuts); N_{hits}^{FW}; Nevents", 40, 0, 40);
	TH1F *hFW_charge = new TH1F ("hFW_charge", "FW charge; FW_{charge}; Nevents", 1000, 0, 11000);
	TH1F *hFW_charge_cut = new TH1F ("hFW_charge_cut", "FW charge (with cuts); FW_{charge}; Nevents", 500, 0, 3500);
	TH1F *hFW_hit_time = new TH1F ("hFW_hit_time", "FW time; FW_{time}; Nhits", 500, -40., 940.);
	TH1F *hFW_hit_time_cut = new TH1F ("hFW_hit_time_cut", "FW time (with cuts); FW_{time}; Nhits", 500, 20., 35.);
	TH1F *hFW_hit_charge = new TH1F ("hFW_hit_charge", "FW charge; FW_{charge}; Nhits", 500, 0., 3000.);
	TH1F *hFW_hit_charge_cut = new TH1F ("hFW_hit_charge_cut", "FW charge (with cuts); FW_{charge}; Nhits", 210, 70., 140.);
	TH2F *h2FW_mod_hits = new TH2F ("h2FW_mod_hits", "FW hit number for separate channels (before cuts); channel ID; N_{hits}^{FW_mod}", 305, 0, 305, 2, 0, 2);
	TH2F *h2FW_mod_hits_cut = new TH2F ("h2FW_mod_hits_cut", "FW hit number for separate channels (after cuts); channel ID; N_{hits}^{FW_mod}", 305, 0, 305, 2, 0, 2);
	TH2F *h2FW_mod_hit_time = new TH2F ("h2FW_mod_hit_time", "FW hit time for separate channels (before cuts); channel ID; hit time", 305, 0, 305, 100, 0, 900);
	TH2F *h2FW_mod_hit_time_cut = new TH2F ("h2FW_mod_hit_time_cut", "FW hit time for separate channels (after cuts); channel ID; hit time", 305, 0, 305, 100, 15, 35);
	TH2F *h2FW_mod_hit_charge = new TH2F ("h2FW_mod_hit_charge", "FW hit charge for separate channels (before cuts); channel ID; hit charge", 305, 0, 305, 500, 0, 3000);
	TH2F *h2FW_mod_hit_charge_cut = new TH2F ("h2FW_mod_hit_charge_cut", "FW hit charge for separate channels (after cuts); channel ID; hit charge", 305, 0, 305, 500, 80, 130);
	TH2F *h2FW_mod_charge = new TH2F ("h2FW_mod_charge", "FW charge for separate channels (before cuts); channel ID; charge", 305, 0, 305, 500, 0, 3000);
	TH2F *h2FW_mod_charge_cut = new TH2F ("h2FW_mod_charge_cut", "FW charge for separate channels (after cuts); channel ID; charge", 305, 0, 305, 500, 0, 150);
	TH2F *h2FW_ring_hits = new TH2F ("h2FW_ring_hits", "FW hit number for separate rings (before cuts); ring ID; N_{hits}^{FW_ring}",9, 0, 9, 25, 0, 25);
	TH2F *h2FW_ring_hits_cut = new TH2F ("h2FW_ring_hits_cut", "FW hit number for separate rings (after cuts); ring ID; N_{hits}^{FW_ring}",9, 0, 9, 15, 0, 15);
	TH2F *h2FW_ring_hit_time = new TH2F ("h2FW_ring_hit_time", "FW hit time for separate rings (before cuts); ring ID; hit time", 9, 0, 9, 100, -40, 1400);
	TH2F *h2FW_ring_hit_time_cut = new TH2F ("h2FW_ring_hit_time_cut", "FW hit time for separate rings (after cuts); ring ID; hit time",9, 0, 9, 100, 15, 35);
	TH2F *h2FW_ring_hit_charge = new TH2F ("h2FW_ring_hit_charge", "FW hit charge for separate rings (before cuts); ring ID; hit charge",9, 0, 9, 500, 0, 3000);
	TH2F *h2FW_ring_hit_charge_cut = new TH2F ("h2FW_ring_hit_charge_cut", "FW hit charge for separate rings (after cuts); ring ID; hit charge",9, 0, 9, 500, 80, 130);
	TH2F *h2FW_ring_charge = new TH2F ("h2FW_ring_charge", "FW charge for separate rings (before cuts); ring ID; charge",9, 0, 9, 500, 0, 4500);
	TH2F *h2FW_ring_charge_cut = new TH2F ("h2FW_ring_charge_cut", "FW charge for separate rings (after cuts); ring ID; charge",9, 0, 9, 500, 0, 1300);
	TH2F *h2MDC_TOF_RPC_hits = new TH2F ("h2MDC_TOF_RPC_hits", "Correlation between number of MDC tracks and TOF+RPC hits; N_{tracks}^{MDC}; N_{hits}^{TOF} + N_{hits}^{RPC}", 100, 0, 100, 280, 0, 280);
	TH2F *h2MDC_RPC_hits = new TH2F ("h2MDC_RPC_hits", "Correlation between number of MDC tracks and RPC hits; N_{tracks}^{MDC}; N_{hits}^{RPC}", 100, 0, 100, 180, 0, 180);
	TH2F *h2MDC_TOF_hits = new TH2F ("h2MDC_TOF_hits", "Correlation between number of MDC tracks and TOF hits; N_{tracks}^{MDC}; N_{hits}^{TOF}", 100, 0, 100, 84, 0, 84);
	TH2F *h2RPC_TOF_hits = new TH2F ("h2RPC_TOF_hits", "Correlation between number of RPC and TOF hits; N_{hits}^{RPC}; N_{hits}^{TOF}", 180, 0, 180, 84, 0, 84);
	TH2F *h2FW_MDC_hits = new TH2F ("h2FW_MDC_hits", "Correlation between number of FW and MDC tracks; N_{hits}^{FW}; N_{tracks}^{MDC}", 160, 0, 160, 100, 0, 100);
	TH2F *h2FW_TOF_RPC_hits = new TH2F ("h2FW_TOF_RPC_hits", "Correlation between number of FW and TOF+RPC hits; N_{hits}^{FW}; N_{hits}^{TOF} + N_{hits}^{RPC}", 160, 0, 160, 280, 0, 280);
	TH2F *h2FW_RPC_hits = new TH2F ("h2FW_RPC_hits", "Correlation between number of FW and RPC hits; N_{hits}^{FW}; N_{hits}^{RPC}", 160, 0, 160, 180, 0, 180);
	TH2F *h2FW_TOF_hits = new TH2F ("h2FW_TOF_hits", "Correlation between number of FW and TOF hits; N_{hits}^{FW}; N_{hits}^{TOF}", 160, 0, 160, 84, 0, 84);
	TH2F *h2FW_MDC_charge = new TH2F ("h2FW_MDC_charge", "Correlation between FW charge and number of MDC tracks; FW_{charge}; N_{tracks}^{MDC}", 500, 0, 10000, 100, 0, 100);
	TH2F *h2FW_TOF_RPC_charge = new TH2F ("h2FW_TOF_RPC_charge", "Correlation between FW charge and number of TOF+RPC hits; FW_{charge}; N_{hits}^{TOF} + N_{hits}^{RPC}", 500, 0, 10000, 280, 0, 280);
	TH2F *h2FW_RPC_charge = new TH2F ("h2FW_RPC_charge", "Correlation between FW charge and number of RPC hits; FW_{charge}; N_{hits}^{RPC}", 500, 0, 10000, 180, 0, 180);
	TH2F *h2FW_TOF_charge = new TH2F ("h2FW_TOF_charge", "Correlation between FW charge and number of TOF hits; FW_{charge}; N_{hits}^{TOF}", 500, 0, 10000, 84, 0, 84);
    TH1F *hCent = new TH1F ("hCent", "Centrality distribution; centrality, %; nEvents", nBinsCent_, centMin_, centMax_);
	TH2F *h2Cent_MDC_tracks = new TH2F ("h2Cent_MDC_tracks", "Centrality vs number of MDC tracks; centrality, %; N_{tracks}^{MDC}", nBinsCent_, centMin_, centMax_, 100, 0, 100);
	TH2F *h2Cent_TOF_RPC_hits = new TH2F ("h2Cent_TOF_RPC_hits", "Centrality vs number of RPC and TOF hits; Centrality, %; nRpcClustCut + nTofHitsCut", nBinsCent_, centMin_, centMax_, 250, 0, 250);
    TH2F *h2VertexXY = new TH2F ("h2VertexXY", "Vertex X and Y distribution;Vertex X;Vertex Y;nEvents", 100, -20., 20., 100, -20., 20.);
    TH1F *hVertexZ = new TH1F ("hVertexZ", "Vertex Z distribution;Vertex Z;nEvents", 500, -60., 0.);
    TH1F *hVertex_chi2 = new TH1F ("hVertex_chi2", "Vertex #chi^{2} / nTracks;#chi^{2} / nTracks;nEvents", 100, 0., 6.);

    TH1F *hMetaMatchRadius = new TH1F ("hMetaMatchRadius", "Meta Match Radius; MetaMatchRadius; counts", 100, 0., 60.);
    TH1F *hMetaDx = new TH1F ("hMetaDx", "Meta Dx; Meta Dx; nCounts", 100, -50., 50.);
    TH1F *hMetaDy = new TH1F ("hMetaDy", "Meta Dy; Meta Dy; nCounts", 100, -50., 50.);
    TH1F *hMetaQ = new TH1F ("hMetaQ", "Meta quality; Meta Q; counts", 100, 0., 3.);
    TH1F *hChi2_nHits_all = new TH1F ("hChi2all", "Track #chi^{2} / N_{hits} for all segments;#chi^{2} / N_{hits};counts", 500, 0.0, 70.0);
    TH1F *hChi2inner = new TH1F ("hChi2inner", "Track #chi^{2} / N_{inner}^{hits} for inner segments;#chi^{2} / N_{inner}^{hits};counts", 100, 0.0, 20.0);
    TH1F *hChi2outer = new TH1F ("hChi2outer", "Track #chi^{2} / N_{outer}^{hits} for outer segments;#chi^{2} / N_{outer}^{hits};counts", 100, 0.0, 50.0);
    TH1F *hDCAxy = new TH1F ("DCAxy", "DCA_{xy} distribution;DCA_{xy};counts", 100, -20.0, 20.0);
    TH1F *hDCAz_Vz = new TH1F ("hDCAz_Vz", "DCA_{z} - Z_{vertex} distribution;DCA_{z} - Z_{vertex};counts", 100, -20.0, 20.0);
    TH1F *hMass = new TH1F ("hMass", "Mass distribution; M[MeV]; counts", 100, -5000., 5000.);
    TH1F *hBeta = new TH1F ("hBeta", "#beta distribution; #beta; counts", 100, -2., 2.);
    TH2F *h2logPqMass = new TH2F ("h2logPqMass", "q*log(20P) vs mass;q*log(20P); M[MeV]", 500, -log (20. * pMax), log (20. * pMax), 500, 0., 5000.);
    TH2F *h2logPqBeta = new TH2F ("h2logPqBeta", "q*log(20P) vs #beta;q*log(20P); #beta", 500, -log (20. * pMax), log (20. * pMax), 500, 0., 2.);
    TH1F *hP = new TH1F ("hP", "P distribution; P[GeV/c]; counts", 100, pMin, pMax);
    TH1F *hPcorr = new TH1F ("hPcorr", "Pcorr distribution; P^{corr}[GeV/c]; counts", 100, pMin, pMax);
    TH1F *hPt = new TH1F ("hPt", "P_{T} distribution; P_{T}[GeV/c]; counts", 100, ptMin, ptMax);
    TH1F *hPtCorr = new TH1F ("hPtCorr", "P_{T}^{corr} distribution; P_{T}^{corr} [GeV/c]; counts", 100, ptMin, ptMax);
	TH1F *hEta = new TH1F ("hEta", "#eta distribution;#eta; counts", 100, etaMin, etaMax);
	TH1F *hRap = new TH1F ("hRap", "Rapidity distribution;Rapidity; counts", 100, rapMin, rapMax);
	TH1F *hRapCorr = new TH1F ("hRapCorr", "Corrected rapidity distribution;Corrected rapidity; counts", 100, rapMin, rapMax);
	TH1F *hTheta = new TH1F ("hTheta", "#theta distribution;#theta; counts", 100, thetaMin, thetaMax);
	TH1F *hPhi = new TH1F ("hPhi", "#phi distribution; #phi; counts", 100, -PI, PI);
    TH2F *h2EtaPt = new TH2F ("h2EtaPt", "#eta vs P_{T} distribution; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 700, ptMin, ptMax);
    TH2F *h2RapPt = new TH2F ("h2RapPt", "Rapidity vs P_{T} distribution; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 700, ptMin, ptMax);
    TH2F *h2PhiPt = new TH2F ("h2PhiPt", "#phi vs P_{T}distribution;#phi; P_{T}[GeV/c]; counts", 500, -PI, PI, 700, ptMin, ptMax);
	TH2F *h2EtaPhi = new TH2F ("h2EtaPhi", "#eta vs #phi distribution;#eta; #phi; counts", 500, etaMin, etaMax, 500, -PI, PI);
	TH2F *h2RapPhi = new TH2F ("h2RapPhi", "Rapidity vs #phi distribution;Rapidity; counts", 500, rapMin, rapMax, 500, -PI, PI);
	TH1F *hdEdx_TOF = new TH1F ("hdEdx_TOF", "TOF dE/dx distribution", 500, tofdEdxMin, tofdEdxMax);
	TH1F *hdEdx_MDC = new TH1F ("hdEdx_MDC", "MDC dE/dx distribution", 500, mdcdEdxMin, mdcdEdxMax);

	TH2F *h2logPqdEdx_TOF = new TH2F ("h2logPqdEdx_TOF", "Energy loss in TOF vs q*log(20P); q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, tofdEdxMin, tofdEdxMax);
    TH2F *h2logPqdEdx_MDC = new TH2F ("h2logPqdEdx_MDC", "Energy loss in MDC vs q*log(20P); q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, mdcdEdxMin, mdcdEdxMax);

    TH1F *hMassPi = new TH1F ("hMassPi", "Mass distribution for #pi- and #pi+; M[MeV]; counts", 100, -600., 600.);
    TH1F *hBetaPi = new TH1F ("hBetaPi", "#beta distribution for #pi- and #pi+; #beta; counts", 100, -2., 2.);
    TH2F *h2logPqMassPi = new TH2F ("h2logPqMassPi", "q*log(20P) vs mass for #pi- and #pi+;q*log(20P); M[MeV]", 500, -log (20. * pMax), log (20. * pMax), 500, 0., 500.);
    TH2F *h2logPqBetaPi = new TH2F ("h2logPqBetaPi", "q*log(20P) vs #beta for #pi- and #pi+;q*log(20P); #beta", 500, -log (20. * pMax), log (20. * pMax), 500, 0., 2.);

	TH2F *h2logPqdEdxPi_TOF = new TH2F ("h2logPqdEdxPi_TOF", "Energy loss in TOF vs q*log(20P) for #pi- and #pi+; q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, tofdEdxMin, tofdEdxMax);
	TH2F *h2logPqdEdxPi_MDC = new TH2F ("h2logPqdEdxPi_MDC", "Energy loss in MDC vs q*log(20P) for #pi- and #pi+; q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, mdcdEdxMin, mdcdEdxMax);

    TH2F *h2EtaPtPi = new TH2F ("h2EtaPtPi", "#eta vs P_{T} distribution for #pi+ and #pi-; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 500, ptMin, ptMax);
    TH2F *h2RapPtPi = new TH2F ("h2RapPtPi", "Rapidity vs P_{T} distribution for #pi+ and #pi-; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 500, ptMin, ptMax);

    TH1F *hPtPiMin = new TH1F ("hPtPiMin", "P_{T} distribution for #pi-; P_{T}[GeV/c]; counts", 100, ptMin, ptMax);
	TH1F *hEtaPiMin = new TH1F ("hEtaPiMin", "#eta distribution for #pi-;#eta; counts", 100, etaMin, etaMax);
	TH1F *hRapPiMin = new TH1F ("hRapPiMin", "Rapidity distribution for #pi-;Rapidity; counts", 100, rapMin, rapMax);
	TH1F *hThetaPiMin = new TH1F ("hThetaPiMin", "#theta distribution for #pi-;#theta; counts", 100, thetaMin, thetaMax);
	TH1F *hPhiPiMin = new TH1F ("hPhiPiMin", "#phi distribution for #pi-; #phi; counts", 100, -PI, PI);
	TH2F *h2PhiPtPiMin = new TH2F ("h2PhiPtPiMin", "#phi vs P_{T}distribution for #pi-;#phi; P_{T}[GeV/c]; counts", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPiMin = new TH2F ("h2EtaPtPiMin", "#eta vs P_{T} distribution for #pi-; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPiMin = new TH2F ("h2EtaPhiPiMin", "#eta vs #phi distribution fPrimaryParticles_fTmeanCharge [] [2]for #pi-;#eta; #phi; counts", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPiMin = new TH2F ("h2RapPtPiMin", "Rapidity vs P_{T} distribution for #pi-; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPiMin = new TH2F ("h2RapPhiPiMin", "Rapidity vs #phi distribution for #pi-;Rapidity; #phi; counts", 500, rapMin, rapMax, 500, -PI, PI);

    TH1F *hPtPiPlus = new TH1F ("hPtPiPlus", "P_{T} distribution for #pi+; P_{T}[GeV/c]; counts", 100, ptMin, ptMax);
	TH1F *hEtaPiPlus = new TH1F ("hEtaPiPlus", "#eta distribution for #pi+;#eta; counts", 100, etaMin, etaMax);
	TH1F *hRapPiPlus = new TH1F ("hRapPiPlus", "Rapidity distribution for #pi+;Rapidity; counts", 100, rapMin, rapMax);
	TH1F *hThetaPiPlus = new TH1F ("hThetaPiPlus", "#theta distribution for #pi+;#theta; counts", 100, thetaMin, thetaMax);
	TH1F *hPhiPiPlus = new TH1F ("hPhiPiPlus", "#phi distribution for #pi+; #phi; counts", 100, -PI, PI);
	TH2F *h2PhiPtPiPlus = new TH2F ("h2PhiPtPiPlus", "#phi vs P_{T}distribution for #pi+;#phi; P_{T}[GeV/c]; counts", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPiPlus = new TH2F ("h2EtaPtPiPlus", "#eta vs P_{T} distribution for #pi+; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPiPlus = new TH2F ("h2EtaPhiPiPlus", "#eta vs #phi distribution for #pi+;#eta; #phi; counts", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPiPlus = new TH2F ("h2RapPtPiPlus", "Rapidity vs P_{T} distribution for #pi+; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPiPlus = new TH2F ("h2RapPhiPiPlus", "Rapidity vs #phi distribution for #pi+;Rapidity; #phi; counts", 500, rapMin, rapMax, 500, -PI, PI);

    TH1F *hMassProton = new TH1F ("hMassProton", "Mass distribution for protons; M[MeV]; counts", 100, 0.0, 1400.);
    TH1F *hBetaProton = new TH1F ("hBetaProton", "#beta distribution for protons; #beta; counts", 100, 0., 1.2);
	TH2F *h2logPqdEdxProton_TOF = new TH2F ("h2logPqdEdxProton_TOF", "Energy loss in TOF vs q*log(20P) for protons; q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, tofdEdxMin, tofdEdxMax);
	TH2F *h2logPqdEdxProton_MDC = new TH2F ("h2logPqdEdxProton_MDC", "Energy loss in MDC vs q*log(20P) for protons; q * log (20P); dE/dx [MIP]", 500, -log (20. * pMax), log (20. * pMax), 500, mdcdEdxMin, mdcdEdxMax);
    TH2F *h2logPqMassProton = new TH2F ("h2logPqMassProton", "q*log(20P) vs mass for protons;q*log(20P); M[MeV]", 500, -log (20. * pMax), log (20. * pMax), 500, 500., 1500.);
    TH2F *h2logPqBetaProton = new TH2F ("h2logPqBetaProton", "q*log(20P) vs #beta for protons;q*log(20P); #beta", 500, -log (20. * pMax), log (20. * pMax), 500, 0., 2.);

	TH1F *hPtProton = new TH1F ("hPtProton", "P_{T} distribution for protons; P_{T}[GeV/c]; counts", 100, ptMin, ptMax);
	TH1F *hEtaProton = new TH1F ("hEtaProton", "#eta distribution for protons;#eta; counts", 100, etaMin, etaMax);
	TH1F *hRapProton = new TH1F ("hRapProton", "Rapidity distribution for protons;Rapidity; counts", 100, rapMin, rapMax);
	TH1F *hThetaProton = new TH1F ("hThetaProton", "#theta distribution for protons;#theta; counts", 100, thetaMin, thetaMax);
	TH1F *hPhiProton = new TH1F ("hPhiProton", "#phi distribution for protons; #phi; counts", 100, -PI, PI);
	TH2F *h2PhiPtProton = new TH2F ("h2PhiPtProton", "#phi vs P_{T}distribution for protons;#phi; P_{T}[GeV/c]; counts", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtProton = new TH2F ("h2EtaPtProton", "#eta vs P_{T} distribution for protons; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiProton = new TH2F ("h2EtaPhiProton", "#eta vs #phi distribution for protons;#eta; #phi; counts", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtProton = new TH2F ("h2RapPtProton", "Rapidity vs P_{T} distribution for protons; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiProton = new TH2F ("h2RapPhiProton", "Rapidity vs #phi distribution for protons;Rapidity; #phi; counts", 500, rapMin, rapMax, 500, -PI, PI);
	TH2F *h2RapEtaProton = new TH2F ("h2RapEtaProton", "Rapidity vs #eta distribution for protons;Rapidity; #eta; counts", 500, rapMin, rapMax, 500, etaMin, etaMax);

//	TH1F *hPtAntiProton = new TH1F ("hPtAntiProton", "P_{T} distribution for antiprotons; P_{T}[GeV/c]; counts", 100, ptMin, ptMax);
//	TH1F *hEtaAntiProton = new TH1F ("hEtaAntiProton", "#eta distribution for antiprotons;#eta; counts", 100, etaMin, etaMax);
//	TH1F *hRapAntiProton = new TH1F ("hRapAntiProton", "Rapidity distribution for antiprotons;Rapidity; counts", 100, rapMin, rapMax);
//	TH1F *hThetaAntiProton = new TH1F ("hThetaAntiProton", "#theta distribution for antiprotons;#theta; counts", 100, thetaMin, thetaMax);
//	TH1F *hPhiAntiProton = new TH1F ("hPhiAntiProton", "#phi distribution for antiprotons; #phi; counts", 100, -PI, PI);
//	TH2F *h2PhiPtAntiProton = new TH2F ("h2PhiPtAntiProton", "#phi vs P_{T}distribution for antiprotons;#phi; P_{T}[GeV/c]; counts", 500, -PI, PI, 500, ptMin, ptMax);
//	TH2F *h2EtaPtAntiProton = new TH2F ("h2EtaPtAntiProton", "#eta vs P_{T} distribution for antiprotons; #eta;P_{T}[GeV/c]; counts", 500, etaMin, etaMax, 500, ptMin, ptMax);
//	TH2F *h2EtaPhiAntiProton = new TH2F ("h2EtaPhiAntiProton", "#eta vs #phi distribution for antiprotons;#eta; #phi; counts", 500, etaMin, etaMax, 500, -PI, PI);
//    TH2F *h2RapPtAntiProton = new TH2F ("h2RapPtAntiProton", "Rapidity vs P_{T} distribution for antiprotons; Rapidity;P_{T}[GeV/c]; counts", 500, rapMin, rapMax, 500, ptMin, ptMax);
//	TH2F *h2RapPhiAntiProton = new TH2F ("h2RapPhiAntiProton", "Rapidity vs #phi distribution for antiprotons;Rapidity; #phi; counts", 500, rapMin, rapMax, 500, -PI, PI);

	TH2F *h2PidCharge = new TH2F ("h2PidCharge", "Particle ID vs charge; pID; charge", kNPartTypes, -0.5, kNPartTypes - 0.5, 3, -1.5, 1.5);

	TProfile *pNtracksTime = new TProfile ("pNtracksTime", "Mean number of MDC tracks during one spill;N_{tracks}^{MDC};runTime", 100, 919800, 919900);
	TProfile *pRPChitsTime = new TProfile ("pRPChitsTime", "Mean number of RPC hits during one spill;N_{hits}^{RPC};runTime", 100, 919800, 919900);
	TProfile *pTOFhitsTime = new TProfile ("pTOFhitsTime", "Mean number of TOF hits during one spill;N_{hits}^{TOF};runTime", 100, 919800, 919900);
	TProfile *pTOFRPChitsTime = new TProfile ("pTOFRPChitsTime", "Mean number of TOF+RPC hits during one spill;N_{hits}^{TOF+RPC};runTime", 100, 919800, 919900);

	Int_t trackIndex, nRun, pid, charge, FWhits_cut, FWhits_mod [304], FWhits_ring [9], FWhits_mod_cut [304], FWhits_ring_cut [9];
	Float_t cent, pt, eta, rap, phi, p, p_corr, theta, FWcharge_cut, wallChargeTot_mod_cut [304], wallChargeTot_ring_cut [9];

    cout << "Input Tree: " << inputFileName_ << endl;
	Long64_t nentries = inputTree_ -> GetEntries ();

    GetAverages ();

    cout << "\nConverting tree:\n";
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		cout << "\rEvent " << jentry + 1 << " from " << nentries;
		inputTree_ -> GetEntry (jentry);
		if (!CheckEventCuts ()) continue;

        //centrality classes - ToF+RPC
		//if (centMethod_ == 1) cent = GetCentralityClass (nRpcClustCut+nTofHitsCut);
		//if (centMethod_ == 2) cent = GetCentralityClass (fEveto);
		cent = GetCentralityClass (nRpcHitsCut + nTofHitsCut);
//		cent = GetCentralityClass (nRpcClustCut + nTofHitsCut);

		if (cent < centMin_ || cent > centMax_) continue;

		nRun = TMath::Abs(runId);
		event_ -> SetNrun (nRun);
		event_ -> SetCent (cent);
		event_ -> SetPsiA (psiEPa);
		event_ -> SetPsiB (psiEPb);
		event_ -> SetNa (nA);
		event_ -> SetNb (nB);

        trackIndex = 0;
		for (Int_t itrack = 0; itrack < nTracks; itrack++) {
		    if (!CheckTrackCuts (itrack)) continue;
		    pid = GetTrackPid (itrack);
//		    if (pid == kNID) continue;
		  	trackIndex++;

		    charge = charge_ [itrack];
		    pt = pt_corr [itrack] / 1000.;
		    eta = eta_ [itrack];
		    rap = rapidity_corr [itrack];
		    phi = phi_ [itrack];
		    p = p_ [itrack] / 1000.;
		    p_corr = pCorr [itrack] / 1000.;
		    theta = theta_ [itrack];
		    if (phi > PI) phi -= 2 * PI;

            event_ -> AddTrack (pt, eta, phi, charge, pid);
            event_ -> GetTrack (trackIndex) -> SetP (p_corr);
            event_ -> GetTrack (trackIndex) -> SetRap (rap);

            hMetaMatchRadius -> Fill (metaMatchRadius [itrack]);
            hMetaDx -> Fill (metaDx [itrack]);
            hMetaDy -> Fill (metaDy [itrack]);
            hMetaQ -> Fill (metaQ [itrack]);
            hChi2_nHits_all -> Fill (chi2all [itrack] / (Float_t) mdcNhits [itrack]);
            hChi2inner -> Fill (chi2inner [itrack] / (Float_t) mdcNhitsInner [itrack]);
            hChi2outer -> Fill (chi2outer [itrack] / (Float_t) mdcNhitsOuter [itrack]);
            hDCAxy -> Fill (DCAxy [itrack]);
            hDCAz_Vz -> Fill (DCAz [itrack] - vZ);
            hMDC_hits -> Fill (mdcNhits [itrack]);
            hMDC_hits_inner -> Fill (mdcNhitsInner [itrack]);
            hMDC_hits_outer -> Fill (mdcNhitsOuter [itrack]);
            hMass -> Fill (metaMass [itrack] * charge);
            hBeta -> Fill (metaBeta [itrack] * charge);
            h2logPqMass -> Fill (TMath::Log(20.0 * p) * charge, metaMass [itrack]);
            h2logPqBeta -> Fill (TMath::Log(20.0 * p) * charge, metaBeta [itrack]);

            hP -> Fill (p);
            hPt -> Fill (pt_ [itrack] / 1000.);
            hEta -> Fill (eta_ [itrack]);
            hRap -> Fill (rapidity [itrack]);
            hPcorr -> Fill (p_corr);
            hPtCorr -> Fill (pt_corr [itrack] / 1000.);
            hRapCorr -> Fill (rapidity_corr [itrack]);
		    hTheta -> Fill (theta);
		    hPhi -> Fill (phi);
		    h2EtaPt -> Fill (eta, pt);
		    h2RapPt -> Fill (rap, pt);
		    h2PhiPt -> Fill (phi, pt);
		    h2EtaPhi -> Fill (eta, phi);
		    h2RapPhi -> Fill (rap, phi);

            hdEdx_TOF -> Fill (tofdEdx [itrack]);
            hdEdx_MDC -> Fill (mdcdEdx [itrack]);
            h2logPqdEdx_TOF -> Fill (TMath::Log(20.0 * p) * charge, tofdEdx [itrack]);
            h2logPqdEdx_MDC -> Fill (TMath::Log(20.0 * p) * charge, mdcdEdx [itrack]);
			h2PidCharge -> Fill (pid, charge);

			if (pid == kPionMinus || pid == kPionPlus) {
                hMassPi -> Fill (metaMass [itrack] * charge);
                hBetaPi -> Fill (metaBeta [itrack] * charge);
                h2logPqMassPi -> Fill (TMath::Log(20.0 * p) * charge, metaMass [itrack]);
                h2logPqBetaPi -> Fill (TMath::Log(20.0 * p) * charge, metaBeta [itrack]);
                h2logPqdEdxPi_TOF -> Fill (TMath::Log(20.0 * p) * charge, mdcdEdx [itrack]);
                h2logPqdEdxPi_MDC -> Fill (TMath::Log(20.0 * p) * charge, tofdEdx [itrack]);
                h2EtaPtPi -> Fill (eta, pt);
                h2RapPtPi -> Fill (rap, pt);
			}

			if (pid == kPionMinus) {
                hPtPiMin -> Fill (pt);
                hEtaPiMin -> Fill (eta);
                hRapPiMin -> Fill (rap);
                hThetaPiMin -> Fill (theta);
                hPhiPiMin -> Fill (phi);
                h2EtaPtPiMin -> Fill (eta, pt);
                h2RapPtPiMin -> Fill (rap, pt);
                h2PhiPtPiMin -> Fill (phi, pt);
                h2EtaPhiPiMin -> Fill (eta, phi);
                h2RapPhiPiMin -> Fill (rap, phi);
			}

			if (pid == kPionPlus) {
                hPtPiPlus -> Fill (pt);
                hEtaPiPlus -> Fill (eta);
                hRapPiPlus -> Fill (rap);
                hThetaPiPlus -> Fill (theta);
                hPhiPiPlus -> Fill (phi);
                h2EtaPtPiPlus -> Fill (eta, pt);
                h2RapPtPiPlus -> Fill (rap, pt);
                h2PhiPtPiPlus -> Fill (phi, pt);
                h2EtaPhiPiPlus -> Fill (eta, phi);
                h2RapPhiPiPlus -> Fill (rap, phi);
			}

			if (pid == kProton || pid == kAntiProton) {
                hMassProton -> Fill (metaMass [itrack] * charge);
                hBetaProton -> Fill (metaBeta [itrack] * charge);
                h2logPqdEdxProton_TOF -> Fill (TMath::Log(20.0 * p) * charge, tofdEdx [itrack]);
                h2logPqdEdxProton_MDC -> Fill (TMath::Log(20.0 * p) * charge, mdcdEdx [itrack]);
                h2logPqMassProton -> Fill (TMath::Log(20.0 * p) * charge, metaMass [itrack]);
                h2logPqBetaProton -> Fill (TMath::Log(20.0 * p) * charge, metaBeta [itrack]);
			}

			if (pid == kProton) {
                hPtProton -> Fill (pt);
                hEtaProton -> Fill (eta);
                hRapProton -> Fill (rap);
                hThetaProton -> Fill (theta);
                hPhiProton -> Fill (phi);
                h2EtaPtProton -> Fill (eta, pt);
                h2RapPtProton -> Fill (rap, pt);
                h2PhiPtProton -> Fill (phi, pt);
                h2EtaPhiProton -> Fill (eta, phi);
                h2RapPhiProton -> Fill (rap, phi);
                h2RapEtaProton -> Fill (rap, eta);
			}

//			if (pid == kAntiProton) {
//                hPtAntiProton -> Fill (pt);
//                hEtaAntiProton -> Fill (eta);
//                hRapAntiProton -> Fill (rap);
//                hThetaAntiProton -> Fill (theta);
//                hPhiAntiProton -> Fill (phi);
//                h2EtaPtAntiProton -> Fill (eta, pt);
//                h2RapPtAntiProton -> Fill (rap, pt);
//                h2PhiPtAntiProton -> Fill (phi, pt);
//                h2EtaPhiAntiProton -> Fill (eta, phi);
//                h2RapPhiAntiProton -> Fill (rap, phi);
//			}
        }
        //END track loop;

        FWhits_cut = 0;
        FWcharge_cut = 0.;
        for (Int_t i = 0; i < 304; i++) {
            FWhits_mod [i] = 0;
            FWhits_mod_cut [i] = 0;
            wallChargeTot_mod_cut [i] = 0.;
        }
        for (Int_t i = 0; i < 9; i++) {
            FWhits_ring [i] = 0;
            FWhits_ring_cut [i] = 0.;
            wallChargeTot_ring_cut [i] = 0.;
        }

        for (Int_t j = 0; j < nWallHitsTot; j++) {
            FWhits_mod [wallModuleIndex [j]] ++;
            FWhits_ring [wallHitRing [j]] ++;
            h2FW_mod_hit_time -> Fill (wallModuleIndex [j], wallHitTime [j]);
            h2FW_ring_hit_time -> Fill (wallHitRing [j], wallHitTime [j]);
            hFW_hit_time -> Fill (wallHitTime [j]);
            h2FW_mod_hit_charge -> Fill (wallModuleIndex [j], wallHitCharge [j]);
            h2FW_ring_hit_charge -> Fill (wallHitRing [j], wallHitCharge [j]);
            hFW_hit_charge -> Fill (wallHitCharge [j]);


            if (!isWallHitOk [j]) continue;
//            if (wallHitCharge [j] > 120.) continue;
            if (wallModuleIndex [j] <= 4 && wallHitCharge [j] < 80.) continue;
            if (wallModuleIndex [j] > 4 && wallModuleIndex [j] <= 6 && wallHitCharge [j] < 84.) continue;
            if (wallModuleIndex [j] > 6 && wallHitCharge [j] < 88.) continue;
            if (wallHitTime [j] < 20. || wallHitTime [j] > 35.) continue;
            if (wallHitPhi [j] > PI) wallHitPhi [j] -= 2 * PI;

            FWhits_cut++;
            FWhits_mod_cut [wallModuleIndex [j]] ++;
            FWhits_ring_cut [wallHitRing [j]] ++;
            FWcharge_cut += wallHitCharge [j];
            wallChargeTot_mod_cut [wallModuleIndex [j]] += wallHitCharge [j];
            wallChargeTot_ring_cut [wallHitRing [j]] += wallHitCharge [j];

            hFW_hit_time_cut -> Fill (wallHitTime [j]);
            hFW_hit_charge_cut -> Fill (wallHitCharge [j]);
            h2FW_mod_hit_time_cut -> Fill (wallModuleIndex [j], wallHitTime [j]);
            h2FW_mod_hit_charge_cut -> Fill (wallModuleIndex [j], wallHitCharge [j]);
            h2FW_ring_hit_time_cut -> Fill (wallHitRing [j], wallHitTime [j]);
            h2FW_ring_hit_charge_cut -> Fill (wallHitRing [j], wallHitCharge [j]);

            event_ -> AddTrack (wallHitTime [j], wallHitDistance [j], wallHitPhi [j], wallHitCharge [j], kFW);
            event_ -> GetTrack (trackIndex + FWhits_cut) -> SetRap (wallHitRing [j]);
        }
        //END FW hits loop

        for (Int_t i = 0; i < 304; i++) {
            h2FW_mod_hits -> Fill (i, FWhits_mod [i]);
            h2FW_mod_hits_cut -> Fill (i, FWhits_mod_cut [i]);
            h2FW_mod_charge -> Fill (i, wallChargeTot_mod [i]);
            h2FW_mod_charge_cut -> Fill (i, wallChargeTot_mod_cut [i]);
        }
        for (Int_t i = 0; i < 9; i++) {
            h2FW_ring_hits -> Fill (i, FWhits_ring [i]);
            h2FW_ring_hits_cut -> Fill (i, FWhits_ring_cut [i]);
            h2FW_ring_charge -> Fill (i, wallChargeTot_ring [i]);
            h2FW_ring_charge_cut -> Fill (i, wallChargeTot_ring_cut [i]);
        }

        hTOF_RPC_hits -> Fill (nRpcHitsCut + nTofHitsCut);
		hMDC_tracks -> Fill (trackIndex);
        hTOF_hits -> Fill (nTofHitsCut);
        hRPC_hits -> Fill (nRpcHitsCut);
        h2RPC_TOF_hits -> Fill (nRpcHitsCut, nTofHitsCut);
        hFW_hits -> Fill (nWallHitsTot);
        hFW_charge -> Fill (wallChargeTot);
        hFW_hits_cut -> Fill (FWhits_cut);
        hFW_charge_cut -> Fill (FWcharge_cut);

        h2FW_MDC_hits -> Fill (nWallHitsTot, trackIndex);
        h2MDC_TOF_RPC_hits -> Fill (trackIndex, nRpcHitsCut + nTofHitsCut);
        h2MDC_RPC_hits -> Fill (trackIndex, nRpcHitsCut);
        h2MDC_TOF_hits -> Fill (trackIndex, nTofHitsCut);
        h2FW_TOF_RPC_hits -> Fill (nWallHitsTot, nRpcHitsCut + nTofHitsCut);
        h2FW_RPC_hits -> Fill (nWallHitsTot, nRpcHitsCut);
        h2FW_TOF_hits -> Fill (nWallHitsTot, nTofHitsCut);
        h2FW_MDC_charge -> Fill (wallChargeTot, trackIndex);
        h2FW_TOF_RPC_charge -> Fill (wallChargeTot, nRpcHitsCut + nTofHitsCut);
        h2FW_RPC_charge -> Fill (wallChargeTot, nRpcHitsCut);
        h2FW_TOF_charge -> Fill (wallChargeTot, nTofHitsCut);
		h2Cent_TOF_RPC_hits -> Fill (cent, nRpcHitsCut + nTofHitsCut);
		h2Cent_MDC_tracks -> Fill (cent, trackIndex);
		hCent -> Fill (cent);
		h2VertexXY -> Fill (vX, vY);
		hVertexZ -> Fill (vZ);
		hVertex_chi2 -> Fill (vChi2 / (Float_t) nTracks);

        pNtracksTime -> Fill (runTime, trackIndex);
        pRPChitsTime -> Fill (runTime, nRpcHitsCut);
        pTOFhitsTime -> Fill (runTime, nTofHitsCut);
        pTOFRPChitsTime -> Fill (runTime, nRpcHitsCut + nTofHitsCut);

		outputTree_ -> Fill ();
		event_ -> Clear ();
	}
	histFile -> Write ();
	histFile -> Close ();
	Finish ();
	return 1;
}


void CTreeConverterHADES::Finish () {
	outputFile_ -> Write ();
	outputFile_ -> Close ();
}

#endif // CTREECONVERTERHADES_CXX
