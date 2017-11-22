#ifndef CTREECONVERTER_CXX
#define CTREECONVERTER_CXX

#include <iostream>
#include <TMath.h>
#include <TH1F.h> // QA
#include <TH2F.h> // QA
#include <TProfile.h> // QA
#include <THStack.h> // QA
#include <TCanvas.h> // QA
#include <TLegend.h> // QA
#include <TStyle.h> // QA
#include "CTreeConverter.h"

using namespace std;

CTreeConverter::CTreeConverter () {
	initFlag_ = 0;
}


CTreeConverter::CTreeConverter (TString inputFileName, TString outputFileName, TString inputTreeName) {
	initFlag_ = 0;
	SetInputFileName (inputFileName);
	SetOutputFileName (outputFileName);
	SetInputTreeName (inputTreeName);
	SetdEdxSource ();
	SetdEdxSource ();
}

void CTreeConverter::SetInputFileName (TString inputFileName) {
    inputFileName_ = inputFileName;
}

void CTreeConverter::SetOutputFileName (TString outputFileName) {
    outputFileName_ = outputFileName;
}

void CTreeConverter::SetInputTreeName (TString inputTreeName) {
    inputTreeName_ = inputTreeName;
}

Bool_t CTreeConverter::SetInputFile () {
	//if (inputFile_) delete inputFile_;
	inputFile_ = new TFile (inputFileName_ + ".root", "READ");
	if (!inputFile_) {
		cout << "Input file '" << inputFileName_ << "' not found!!!\n";
		return 0;
	}
	return 1;
}


Bool_t CTreeConverter::SetOutputFile () {
	//if (outputFile_) delete outputFile_;
	outputFile_ = new TFile (outputFileName_ + ".root", "RECREATE");
	if (!outputFile_) {
		cout << "Output file '" << outputFileName_ << "' not accessible!!!\n";
		return 0;
	}
	return 1;
}


Bool_t CTreeConverter::SetInputTree () {
	//if (inputTree_) delete inputTree_;
	inputTree_ = (TTree*) inputFile_ -> Get (inputTreeName_);
	if (!inputTree_) {
		cout << "Input tree '" << inputTreeName_ << "' not found!!!\n";
		return 0;
	}
	return 1;
}

void CTreeConverter::SetMhRange (Int_t mhMin, Int_t mhMax) {
    mhMin_ = mhMin;
    mhMax_ = mhMax;
}

void CTreeConverter::SetCentRange (Float_t centMin, Float_t centMax) {
    centMin_ = centMin;
    centMax_ = centMax;
}

void CTreeConverter::SetNbinsMh (Int_t nBinsMh) {
	nBinsMh_ = nBinsMh;
}

void CTreeConverter::SetNbinsCent (Int_t nBinsCent) {
	nBinsCent_ = nBinsCent;
}

void CTreeConverter::SetdEdxSource (Int_t TPCid) {
    dEdxSource_ = TPCid;
}

void CTreeConverter::SetCentralityMethod (Int_t centMethod) {
    centMethod_ = centMethod;
}

Bool_t CTreeConverter::Init (){
	if (!SetInputFile ()) return 0;
	if (!SetOutputFile ()) return 0;
	if (!SetInputTree ()) return 0;

    inputTree_ -> SetBranchAddress("fVertex_fPchi2", &fVertex_fPchi2, &b_fVertex_fPchi2);
    inputTree_ -> SetBranchAddress("fVeto_fAdcHadron", fVeto_fAdcHadron, &b_Veto_fAdcHadron);
   inputTree_ -> SetBranchAddress("fNRun", &fNRun, &b_fNRun);
   inputTree_ -> SetBranchAddress("fNEvent", &fNEvent, &b_fNEvent);
   inputTree_ -> SetBranchAddress("fTriggerMask", &fTriggerMask, &b_fTriggerMask);
   inputTree_ -> SetBranchAddress("fDate", &fDate, &b_fDate);
   inputTree_ -> SetBranchAddress("fTime", &fTime, &b_fTime);
   inputTree_ -> SetBranchAddress("fEveto", &fEveto, &b_fEveto);
   inputTree_ -> SetBranchAddress("fVertexX", &fVertexX, &b_fVertexX);
   inputTree_ -> SetBranchAddress("fVertexY", &fVertexY, &b_fVertexY);
   inputTree_ -> SetBranchAddress("fVertexZ", &fVertexZ, &b_fVertexZ);
   inputTree_ -> SetBranchAddress("fWfaNbeam", &fWfaNbeam, &b_fWfaNbeam);
   inputTree_ -> SetBranchAddress("fWfaNinter", &fWfaNinter, &b_fWfaNinter);
   inputTree_ -> SetBranchAddress("fWfaBeamTime", &fWfaBeamTime, &b_fWfaBeamTime);
   inputTree_ -> SetBranchAddress("fWfaInterTime", &fWfaInterTime, &b_fWfaInterTime);
   inputTree_ -> SetBranchAddress("fNPrimaryParticles", &fNPrimaryParticles, &b_fNPrimaryParticles);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fIdDet", fPrimaryParticles_fIdDet, &b_fPrimaryParticles_fIdDet);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fCharge", fPrimaryParticles_fCharge, &b_fPrimaryParticles_fCharge);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fNPoint", fPrimaryParticles_fNPoint, &b_fPrimaryParticles_fNPoint);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fNFitPoint", fPrimaryParticles_fNFitPoint, &b_fPrimaryParticles_fNFitPoint);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fNDedxPoint", fPrimaryParticles_fNDedxPoint, &b_fPrimaryParticles_fNDedxPoint);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fNMaxPoint", fPrimaryParticles_fNMaxPoint, &b_fPrimaryParticles_fNMaxPoint);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTmeanCharge", fPrimaryParticles_fTmeanCharge, &b_fPrimaryParticles_fTmeanCharge);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fPz", fPrimaryParticles_fPz, &b_fPrimaryParticles_fPz);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fEta", fPrimaryParticles_fEta, &b_fPrimaryParticles_fEta);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fPhi", fPrimaryParticles_fPhi, &b_fPrimaryParticles_fPhi);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fPt", fPrimaryParticles_fPt, &b_fPrimaryParticles_fPt);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fSigPx", fPrimaryParticles_fSigPx, &b_fPrimaryParticles_fSigPx);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fSigPy", fPrimaryParticles_fSigPy, &b_fPrimaryParticles_fSigPy);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fBx", fPrimaryParticles_fBx, &b_fPrimaryParticles_fBx);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fBy", fPrimaryParticles_fBy, &b_fPrimaryParticles_fBy);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fPchi2", fPrimaryParticles_fPchi2, &b_fPrimaryParticles_fPchi2);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fXFirst", fPrimaryParticles_fXFirst, &b_fPrimaryParticles_fXFirst);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fYFirst", fPrimaryParticles_fYFirst, &b_fPrimaryParticles_fYFirst);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fZFirst", fPrimaryParticles_fZFirst, &b_fPrimaryParticles_fZFirst);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fXLast", fPrimaryParticles_fXLast, &b_fPrimaryParticles_fXLast);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fYLast", fPrimaryParticles_fYLast, &b_fPrimaryParticles_fYLast);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fZLast", fPrimaryParticles_fZLast, &b_fPrimaryParticles_fZLast);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fLabel ", fPrimaryParticles_fLabel , &b_fPrimaryParticles_fLabel );
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofIflag", fPrimaryParticles_fTofIflag, &b_fPrimaryParticles_fTofIflag);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofId", fPrimaryParticles_fTofId, &b_fPrimaryParticles_fTofId);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofX", fPrimaryParticles_fTofX, &b_fPrimaryParticles_fTofX);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofY", fPrimaryParticles_fTofY, &b_fPrimaryParticles_fTofY);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofPathl", fPrimaryParticles_fTofPathl, &b_fPrimaryParticles_fTofPathl);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofCharge", fPrimaryParticles_fTofCharge, &b_fPrimaryParticles_fTofCharge);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofMass2", fPrimaryParticles_fTofMass2, &b_fPrimaryParticles_fTofMass2);
   inputTree_ -> SetBranchAddress("fPrimaryParticles_fTofSigMass2", fPrimaryParticles_fTofSigMass2, &b_fPrimaryParticles_fTofSigMass2);
   inputTree_ -> SetBranchAddress("fRing_fADChadron", fRing_fADChadron, &b_fRing_fADChadron);

	return 1;
}


Bool_t CTreeConverter::CheckEventCuts () {
    if (TMath::Abs(fVertexX) > 0.5) return 0;
    if (TMath::Abs(fVertexY) > 0.5) return 0;
    if (fVertexZ < - 581.6 || fVertexZ > - 580.6) return 0;
    if (fVertex_fPchi2 <= 0.5 || fVertex_fPchi2 > 1.0) return 0;
    if (fNPrimaryParticles <= 10) return 0;
    return 1;
}

Float_t CTreeConverter::GetCentralityClass (Int_t mh) {
    static const Int_t nCentClasses = 20;
	Float_t centClassLimits [nCentClasses - 1] = {364.2, 307.3, 259.9, 218.5, 183.1, 152.1, 125.0,
        102.1, 81.9, 64.7, 50.2, 38.1, 28.3, 20.4, 14.5, 9.9, 6.6, 4.4, 2.6};
//    static const Int_t nCentClasses = 100;
//	Float_t centClassLimits [nCentClasses - 1] = {
//	    424.02, 406.301, 391.289, 377.754, 364.957, 352.406, 340.348, 329.027, 317.953, 307.371,
//        297.035, 287.438, 277.84, 268.734, 259.629, 250.646, 242.156, 234.035, 225.914, 218.285,
//        210.779, 203.52, 196.383, 189.615, 182.848, 176.326, 170.051, 163.898, 157.746, 151.963,
//        146.303, 140.766, 135.413, 130.122, 125.2, 120.278, 115.418, 110.865, 106.497, 102.129,
//        97.7607, 93.6387, 89.5781, 85.7021, 82.0107, 78.2578, 74.6279, 71.2441, 67.9219, 64.7227,
//        61.6465, 58.5703, 55.6787, 52.9717, 50.2646, 47.6807, 45.1582, 42.728, 40.3901, 38.1138,
//        36.022, 33.9917, 32.0845, 30.1465, 28.3315, 26.5781, 24.9478, 23.3789, 21.887, 20.4719,
//        19.1184, 17.8264, 16.6575, 15.5193, 14.4734, 13.4429, 12.5046, 11.5664, 10.7051, 9.90527,
//        9.18237, 8.47485, 7.82886, 7.21362, 6.65991, 6.12927, 5.65247, 5.17566, 4.76807, 4.38354,
//        3.99133, 3.66064, 3.31458, 2.96466, 2.61859, 2.28021, 1.92645, 1.48041, 1.04205};

    Float_t centClassWidth = 100.0 / (Float_t) nCentClasses;
	if (mh >= centClassLimits [0]) return centClassWidth * 0.5;
	for (Int_t i = 1; i < nCentClasses - 1; i++) {
        if (mh < centClassLimits [i - 1] && mh >= centClassLimits [i]) return centClassWidth * (i + 0.5);
	}
	if (mh < centClassLimits [nCentClasses - 2]) return 100.0 - centClassWidth * 0.5;
	return -1.0;
}

Float_t CTreeConverter::GetCentralityClass (Float_t Eveto) {
    Float_t Ebeam = 8.32e3;
    Float_t r = Eveto / Ebeam;
    static const Int_t nCentClasses = 6;
	Float_t centClassLimits [nCentClasses - 1] = {0.169, 0.314, 0.509, 0.66, 0.778};
    Float_t centClassWidth = 100.0 / (Float_t) nCentClasses;
	if (r <= centClassLimits [0]) return centClassWidth * 0.5;
	for (Int_t i = 1; i < nCentClasses - 1; i++) {
        if (r > centClassLimits [i - 1] && r <= centClassLimits [i]) return centClassWidth * (i + 0.5);
	}
	if (r > centClassLimits [nCentClasses - 2]) return 100.0 - centClassWidth * 0.5;
	return -1.0;
}

void CTreeConverter::SetSNN (Float_t sNN) {
    sNN_ = sNN;
}

Float_t CTreeConverter::GetRapidity (Float_t pt, Float_t eta, Int_t pid) {
    Float_t mProton = 938.2720814;
    Float_t m = particleMass [pid];
    if (pid == 0) return -999.0;
    Float_t yBeam = TMath::Log (sNN_ * 1.0e3 / mProton);
    Float_t a = TMath::Sqrt (pt * pt * 1.0e6 * TMath::CosH(eta) * TMath::CosH(eta) + m * m);
    Float_t b = pt * 1.0e3 * TMath::SinH (eta);
	return 0.5 * TMath::Log ((a + b) / (a - b)) - yBeam;
}

Bool_t CTreeConverter::CheckTrackCuts (Int_t itrack) {
    Float_t fitFraction = 0.0;
    Int_t nVTPC1 = fPrimaryParticles_fNMaxPoint [itrack] [0];
	Int_t nVTPC2 = fPrimaryParticles_fNMaxPoint [itrack] [1];
	Int_t nMTPC = fPrimaryParticles_fNMaxPoint [itrack] [2];
	Float_t Bx = fPrimaryParticles_fBx [itrack];
	Float_t By = fPrimaryParticles_fBy [itrack];


    if (fPrimaryParticles_fNMaxPoint [itrack] [3] > 0)
		fitFraction = (Float_t) fPrimaryParticles_fNFitPoint[itrack][3] / fPrimaryParticles_fNMaxPoint[itrack][3];
    if (fitFraction <= 0.55 || fitFraction > 1.0) return 0;
    if (fPrimaryParticles_fPchi2 [itrack]  < 0.0 || fPrimaryParticles_fPchi2 [itrack] > 10.0) return 0;
//    if (nVTPC1 + nVTPC2 < 20 || nMTPC < 30) return 0;
//    if (Bx * Bx / 4.0 + By * By > 1.0) return 0;
    if (nVTPC1 + nVTPC2 < 20 && nMTPC < 30) return 0;
    if (TMath::Abs(Bx) > 2.0 || TMath::Abs(Bx) > 2.0 ) return 0;
    if (fPrimaryParticles_fPt [itrack] > 2.5) return 0;
    if (fPrimaryParticles_fEta [itrack] < 1.4 || fPrimaryParticles_fEta [itrack] > 5.0) return 0;
    if (fPrimaryParticles_fTmeanCharge [itrack] [dEdxSource_] / 1000.0 < 0.01) return 0; // mean energy loss
    return 1;
}

Int_t CTreeConverter::GetTrackPid (Int_t itrack) {
    Float_t p = TMath::Sqrt (fPrimaryParticles_fPz [itrack] * fPrimaryParticles_fPz [itrack] +
                             fPrimaryParticles_fPt [itrack] * fPrimaryParticles_fPt [itrack]);
    Float_t dEdx = fPrimaryParticles_fTmeanCharge [itrack] [dEdxSource_] / 1000.0;
    Int_t charge = fPrimaryParticles_fCharge [itrack];

    if (dEdx > 0.931 + 0.138 * TMath::Log(p) && dEdx < 1.232 + 0.174 * TMath::Log(p)) {
            if (charge > 0) return kPionPlus;
            else return kPionMinus;
    }

//    if (dEdx < 0.931 + 0.138 * TMath::Log(p) && dEdx > 0.5 && p > 3.0) {
    if (dEdx < 0.931 + 0.138 * TMath::Log(p) && dEdx > 0.5) {
        if (charge > 0) return kProton;
        if (charge < 0) return kAntiProton;
    }

    if (dEdx < 1.7 && dEdx > 1.5 && TMath::Log(p) > 0.2 && TMath::Log(p) < 2.2) {
        if (charge < 0) return kElectron;
        if (charge > 0) return kPositron;
    }
    return kNID;
}

Bool_t CTreeConverter::ConvertTree () {
	if (!Init ()) return 0;
	CEvent* event_ = new CEvent;
	outputTree_ = new TTree ("Tree", "Converted tree");
	outputTree_ -> Branch ("Event", &event_, 128000, 4);
	outputFile_ -> cd();

    // QA list
//    Float_t ptMin = 0.0, ptMax = 5.0, etaMin = 0.0, etaMax = 14.0, rapMin = -2.5, rapMax = 5.0; // raw
    Float_t ptMin = 0.0, ptMax = 2.5, etaMin = 1.0, etaMax = 5.5, rapMin = -2.5, rapMax = 3.0; // cut

	TFile *histFile = new TFile (outputFileName_ + "_QA.root", "recreate");

    TH1F *hMh = new TH1F ("hMh", "Multiplicity distribution; mult; nEvents", nBinsMh_, mhMin_, mhMax_);
    TH1F *hCent = new TH1F ("hCent", "Centrality distribution; cent; nEvents", nBinsCent_, centMin_, centMax_);
    TH1F *hEvetoFull = new TH1F ("hEvetoFull", "E_{VETO}^{full} distribution; E_{VETO}^{full}; nEvents", 100, 0.0, 10000.0);

    TH1F *hEveto [4];
    for (Int_t i = 0; i < 4; i++) {
        hEveto [i] = new TH1F (Form ("hEveto%i", i + 1), Form ("E_{VETO}^{%i} distribution; E_{VETO}^{%i}; nEvents", i + 1, i + 1), 100, 0.0, 1000.0);
    }

    TH2F *h2CentMh = new TH2F ("h2CentMh", "Centrality vs multiplicity distribution", nBinsCent_, centMin_, centMax_, nBinsMh_, mhMin_, mhMax_);
    TH2F *h2CentEveto = new TH2F ("h2CentEveto", "Centrality vs E_{VETO} distribution", nBinsCent_, centMin_, centMax_, 100, 0.0, 10000.0);
    TH1F *hPt = new TH1F ("hPt", "P_{T} distribution; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEta = new TH1F ("hEta", "#eta distribution;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRap = new TH1F ("hRap", "Rapidity distribution;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hTheta = new TH1F ("hTheta", "#theta distribution;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhi = new TH1F ("hPhi", "#phi distribution; #phi; nTracks", 100, -PI, PI);
    TH2F *h2EtaPt = new TH2F ("h2EtaPt", "#eta vs P_{T} distribution; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
    TH2F *h2RapPt = new TH2F ("h2RapPt", "Rapidity vs P_{T} distribution; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
    TH2F *h2PhiPt = new TH2F ("h2PhiPt", "#phi vs P_{T}distribution;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPhi = new TH2F ("h2EtaPhi", "#eta vs #phi distribution;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
	TH2F *h2RapPhi = new TH2F ("h2RapPhi", "Rapidity vs #phi distribution;Rapidity; nTracks", 500, rapMin, rapMax, 500, -PI, PI);
	TH1F *hdEdx = new TH1F ("hdEdx", "dEdx distribution", 500, 0.0, 10.0);
    TProfile *hZvertMh = new TProfile ("hZvertMh", "Z_{vertex} vs multiplicity distribution;mult;Z_{vertex}", 100, 10, 510);
    TProfile *hMhZvert = new TProfile ("hMhZvert", "Multiplicity vs Z_{vertex} distribution;Z_{vertex};mult", 10, -581.6, -580.6);

	TH2F *h2logPdEdx = new TH2F ("h2logPdEdx", "Energy loss vs momentum logarithm; log (P); dE/dx [MIP]", 200, -2.0, 5.0, 500, 0.0, 4.0);
    TH2F *h2PqdEdx = new TH2F ("h2PqdEdx", "Energy loss vs momentum multiplied by charge; P*q [GeV/c]; dE/dx [MIP]", 500, -100.0, 100.0, 500, 0.0, 4.0);

	TH1F *hPtPos = new TH1F ("hPtPos", "P_{T} distribution for positively charged particles; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaPos = new TH1F ("hEtaPos", "#eta distribution for positively charged particles;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapPos = new TH1F ("hRapPos", "Rapidity distribution for positively charged particles;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaPos = new TH1F ("hThetaPos", "#theta distribution for positively charged particles;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiPos = new TH1F ("hPhiPos", "#phi distribution for positively charged particles; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtPos = new TH2F ("h2PhiPtPos", "#phi vs P_{T}distribution for positively charged particles;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPos = new TH2F ("h2EtaPtPos", "#eta vs P_{T} distribution for positively charged particles; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPos = new TH2F ("h2EtaPhiPos", "#eta vs #phi distribution for positively charged particles;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPos = new TH2F ("h2RapPtPos", "Rapidity vs P_{T} distribution for positively charged particles; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPos = new TH2F ("h2RapPhiPos", "Rapidity vs #phi distribution for positively charged particles;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

	TH1F *hPtNeg = new TH1F ("hPtNeg", "P_{T} distribution for negatively charged particles; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaNeg = new TH1F ("hEtaNeg", "#eta distribution for negatively charged particles;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapNeg = new TH1F ("hRapNeg", "Rapidity distribution for negatively charged particles;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaNeg = new TH1F ("hThetaNeg", "#theta distribution for negatively charged particles;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiNeg = new TH1F ("hPhiNeg", "#phi distribution for negatively charged particles; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtNeg = new TH2F ("h2PhiPtNeg", "#phi vs P_{T}distribution for negatively charged particles;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtNeg = new TH2F ("h2EtaPtNeg", "#eta vs P_{T} distribution for negatively charged particles; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiNeg = new TH2F ("h2EtaPhiNeg", "#eta vs #phi distribution for negatively charged particles;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtNeg = new TH2F ("h2RapPtNeg", "Rapidity vs P_{T} distribution for negatively charged particles; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiNeg = new TH2F ("h2RapPhiNeg", "Rapidity vs #phi distribution for negatively charged particles;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

	TH2F *h2logPdEdxPi = new TH2F ("h2logPdEdxPi", "Energy loss vs momentum logarithm for #pi- and #pi+; log (P); dE/dx [MIP]", 200, -2.0, 5.0, 500, 0.0, 4.0);
    TH2F *h2PqdEdxPi = new TH2F ("h2PqdEdxPi", "Energy loss vs momentum multiplied by charge for #pi- and #pi+; P*q [GeV/c]; dE/dx [MIP]", 500, -100.0, 100.0, 500, 0.0, 4.0);
    TH2F *h2EtaPtPi = new TH2F ("h2EtaPtPi", "#eta vs P_{T} distribution for #pi+ and #pi-; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
    TH2F *h2RapPtPi = new TH2F ("h2RapPtPi", "Rapidity vs P_{T} distribution for #pi+ and #pi-; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);

    TH1F *hPtPiMin = new TH1F ("hPtPiMin", "P_{T} distribution for #pi-; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaPiMin = new TH1F ("hEtaPiMin", "#eta distribution for #pi-;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapPiMin = new TH1F ("hRapPiMin", "Rapidity distribution for #pi-;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaPiMin = new TH1F ("hThetaPiMin", "#theta distribution for #pi-;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiPiMin = new TH1F ("hPhiPiMin", "#phi distribution for #pi-; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtPiMin = new TH2F ("h2PhiPtPiMin", "#phi vs P_{T}distribution for #pi-;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPiMin = new TH2F ("h2EtaPtPiMin", "#eta vs P_{T} distribution for #pi-; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPiMin = new TH2F ("h2EtaPhiPiMin", "#eta vs #phi distribution for #pi-;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPiMin = new TH2F ("h2RapPtPiMin", "Rapidity vs P_{T} distribution for #pi-; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPiMin = new TH2F ("h2RapPhiPiMin", "Rapidity vs #phi distribution for #pi-;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

    TH1F *hPtPiPlus = new TH1F ("hPtPiPlus", "P_{T} distribution for #pi+; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaPiPlus = new TH1F ("hEtaPiPlus", "#eta distribution for #pi+;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapPiPlus = new TH1F ("hRapPiPlus", "Rapidity distribution for #pi+;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaPiPlus = new TH1F ("hThetaPiPlus", "#theta distribution for #pi+;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiPiPlus = new TH1F ("hPhiPiPlus", "#phi distribution for #pi+; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtPiPlus = new TH2F ("h2PhiPtPiPlus", "#phi vs P_{T}distribution for #pi+;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPiPlus = new TH2F ("h2EtaPtPiPlus", "#eta vs P_{T} distribution for #pi+; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPiPlus = new TH2F ("h2EtaPhiPiPlus", "#eta vs #phi distribution for #pi+;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPiPlus = new TH2F ("h2RapPtPiPlus", "Rapidity vs P_{T} distribution for #pi+; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPiPlus = new TH2F ("h2RapPhiPiPlus", "Rapidity vs #phi distribution for #pi+;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

    TH2F *h2logPdEdxProton = new TH2F ("h2logPdEdxProton", "Energy loss vs momentum logarithm for protons and antiprotons; log (P); dE/dx [MIP]", 200, -2.0, 5.0, 500, 0.0, 4.0);
    TH2F *h2PqdEdxProton = new TH2F ("h2PqdEdxProton", "Energy loss vs momentum multiplied by charge for protons and antiprotons; P*q [GeV/c]; dE/dx [MIP]", 500, -100.0, 100.0, 500, 0.0, 4.0);
	TH1F *hPtProton = new TH1F ("hPtProton", "P_{T} distribution for protons; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaProton = new TH1F ("hEtaProton", "#eta distribution for protons;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapProton = new TH1F ("hRapProton", "Rapidity distribution for protons;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaProton = new TH1F ("hThetaProton", "#theta distribution for protons;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiProton = new TH1F ("hPhiProton", "#phi distribution for protons; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtProton = new TH2F ("h2PhiPtProton", "#phi vs P_{T}distribution for protons;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtProton = new TH2F ("h2EtaPtProton", "#eta vs P_{T} distribution for protons; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiProton = new TH2F ("h2EtaPhiProton", "#eta vs #phi distribution for protons;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtProton = new TH2F ("h2RapPtProton", "Rapidity vs P_{T} distribution for protons; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiProton = new TH2F ("h2RapPhiProton", "Rapidity vs #phi distribution for protons;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);
	TH2F *h2RapEtaProton = new TH2F ("h2RapEtaProton", "Rapidity vs #eta distribution for protons;Rapidity; #eta; nTracks", 500, rapMin, rapMax, 500, etaMin, etaMax);


	TH1F *hPtAntiProton = new TH1F ("hPtAntiProton", "P_{T} distribution for antiprotons; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaAntiProton = new TH1F ("hEtaAntiProton", "#eta distribution for antiprotons;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapAntiProton = new TH1F ("hRapAntiProton", "Rapidity distribution for antiprotons;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaAntiProton = new TH1F ("hThetaAntiProton", "#theta distribution for antiprotons;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiAntiProton = new TH1F ("hPhiAntiProton", "#phi distribution for antiprotons; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtAntiProton = new TH2F ("h2PhiPtAntiProton", "#phi vs P_{T}distribution for antiprotons;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtAntiProton = new TH2F ("h2EtaPtAntiProton", "#eta vs P_{T} distribution for antiprotons; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiAntiProton = new TH2F ("h2EtaPhiAntiProton", "#eta vs #phi distribution for antiprotons;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtAntiProton = new TH2F ("h2RapPtAntiProton", "Rapidity vs P_{T} distribution for antiprotons; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiAntiProton = new TH2F ("h2RapPhiAntiProton", "Rapidity vs #phi distribution for antiprotons;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

    TH2F *h2logPdEdxElectron = new TH2F ("h2logPdEdxElectron", "Energy loss vs momentum logarithm for electrons and positrons; log (P); dE/dx [MIP]", 200, -2.0, 5.0, 500, 0.0, 4.0);
    TH2F *h2PqdEdxElectron = new TH2F ("h2PqdEdxElectron", "Energy loss vs momentum multiplied by charge for electrons and positrons; P*q [GeV/c]; dE/dx [MIP]", 500, -100.0, 100.0, 500, 0.0, 4.0);
    TH1F *hPtElectron = new TH1F ("hPtElectron", "P_{T} distribution for electrons; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaElectron = new TH1F ("hEtaElectron", "#eta distribution for electrons;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapElectron = new TH1F ("hRapElectron", "Rapidity distribution for electrons;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaElectron = new TH1F ("hThetaElectron", "#theta distribution for electrons;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiElectron = new TH1F ("hPhiElectron", "#phi distribution for electrons; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtElectron = new TH2F ("h2PhiPtElectron", "#phi vs P_{T}distribution for electrons;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtElectron = new TH2F ("h2EtaPtElectron", "#eta vs P_{T} distribution for electrons; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiElectron = new TH2F ("h2EtaPhiElectron", "#eta vs #phi distribution for electrons;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtElectron = new TH2F ("h2RapPtElectron", "Rapidity vs P_{T} distribution for electrons; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiElectron = new TH2F ("h2RapPhiElectron", "Rapidity vs #phi distribution for electrons;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

    TH1F *hPtPositron = new TH1F ("hPtPositron", "P_{T} distribution for positrons; P_{T}[GeV/c]; nTracks", 100, ptMin, ptMax);
	TH1F *hEtaPositron = new TH1F ("hEtaPositron", "#eta distribution for positrons;#eta; nTracks", 100, etaMin, etaMax);
	TH1F *hRapPositron = new TH1F ("hRapPositron", "Rapidity distribution for positrons;Rapidity; nTracks", 100, rapMin, rapMax);
	TH1F *hThetaPositron = new TH1F ("hThetaPositron", "#theta distribution for positrons;#theta; nTracks", 100, 0.0, 0.5);
	TH1F *hPhiPositron = new TH1F ("hPhiPositron", "#phi distribution for positrons; #phi; nTracks", 100, -PI, PI);
	TH2F *h2PhiPtPositron = new TH2F ("h2PhiPtPositron", "#phi vs P_{T}distribution for positrons;#phi; P_{T}[GeV/c]; nTracks", 500, -PI, PI, 500, ptMin, ptMax);
	TH2F *h2EtaPtPositron = new TH2F ("h2EtaPtPositron", "#eta vs P_{T} distribution for positrons; #eta;P_{T}[GeV/c]; nTracks", 500, etaMin, etaMax, 500, ptMin, ptMax);
	TH2F *h2EtaPhiPositron = new TH2F ("h2EtaPhiPositron", "#eta vs #phi distribution for positrons;#eta; #phi; nTracks", 500, etaMin, etaMax, 500, -PI, PI);
    TH2F *h2RapPtPositron = new TH2F ("h2RapPtPositron", "Rapidity vs P_{T} distribution for positrons; Rapidity;P_{T}[GeV/c]; nTracks", 500, rapMin, rapMax, 500, ptMin, ptMax);
	TH2F *h2RapPhiPositron = new TH2F ("h2RapPhiPositron", "Rapidity vs #phi distribution for positrons;Rapidity; #phi; nTracks", 500, rapMin, rapMax, 500, -PI, PI);

	TH2F *h2MTPCtot = new TH2F ("h2MTPCtot", "dE/dx: MTPC vs total; MTPC; total; dE/dx [MIP]", 500, 0.0, 4.0, 500, 0.0, 4.0);
	TH2F *h2PidCharge = new TH2F ("h2PidCharge", "Particle ID vs charge; pID; charge", kNPartTypes, -0.5, kNPartTypes - 0.5, 3, -1.5, 1.5);

	// END QA list

	Int_t trackIndex, nRun, mh, mh_cut, pid, charge;
	Float_t cent, pt, pz, p, eta, rap, theta, phi, dEdx, m;
//	Float_t phiVeto [4] = {0.25 * PI, 0.75 * PI, 1.25 * PI, 1.75 * PI}; // 1 2
//                                                                      // 4 3
//	Float_t phiVeto [4] = {0.25 * PI, 0.75 * PI, 1.75 * PI, 1.25 * PI}; // 1 2
//                                                                      // 3 4
	Float_t phiVeto [4] = {0.25 * PI, 1.75 * PI, 0.75 * PI, 1.25 * PI}; // 1 4
                                                                        // 2 3
//	Float_t phiVeto [4] = {1.75 * PI, 1.25 * PI, 0.75 * PI, 0.25 * PI}; // 1 3
                                                                     // // 2 4
	Float_t ptVeto = 1.0, etaVeto = 1.0, rapVeto = 2.23, chargeVeto = 0.0;

    cout << "Converting Tree: " << inputFileName_ << endl;
	Long64_t nentries = inputTree_ -> GetEntries ();
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		cout << "\rEvent " << jentry + 1 << " from " << nentries;
		inputTree_ -> GetEntry (jentry);
		if (!CheckEventCuts ()) continue;
		mh = fNPrimaryParticles;
		mh_cut = 0;

		for (Int_t itrack = 0; itrack < mh; itrack++) {
		    if (!CheckTrackCuts (itrack)) continue;
            mh_cut++;
		}

		if (mh_cut < 10) continue;

		if (centMethod_ == 1) cent = GetCentralityClass (mh_cut);
		if (centMethod_ == 2) cent = GetCentralityClass (fEveto);
		if (cent < centMin_ || cent > centMax_) continue;
		trackIndex = 0;
		nRun = fNRun;
		event_ -> SetNrun (nRun);
		event_ -> SetCent (cent);
		event_ -> SetEvetoFull (fEveto);
		event_ -> SetEveto (fVeto_fAdcHadron);
		hCent -> Fill (cent);
		hMh -> Fill (mh_cut);
		hZvertMh -> Fill (mh_cut, fVertexZ);
		hMhZvert -> Fill (fVertexZ, mh_cut);
		h2CentMh -> Fill (cent, mh_cut);
		hEvetoFull -> Fill (fEveto);
		for (Int_t i = 0; i < 4; i++) hEveto [i] -> Fill (fVeto_fAdcHadron [i]);
		h2CentEveto -> Fill (cent, fEveto);


		// QA fill
		for (Int_t itrack = 0; itrack < mh; itrack++) {
		    if (!CheckTrackCuts (itrack)) continue;
		    trackIndex++;
		    pt = fPrimaryParticles_fPt [itrack];
		    pz = fPrimaryParticles_fPz [itrack];
		    p = TMath::Sqrt (pt * pt + pz * pz);
		    eta = fPrimaryParticles_fEta [itrack];
		    theta = 2 * TMath::ATan (TMath::Exp (-eta));
		    phi = fPrimaryParticles_fPhi [itrack];
		    if (phi > PI) phi -= 2 * PI;
		    dEdx = fPrimaryParticles_fTmeanCharge [itrack] [dEdxSource_] / 1000.0;
		    charge = fPrimaryParticles_fCharge [itrack];
		    pid = GetTrackPid (itrack);
            rap = GetRapidity (pt, eta, pid);

		    hdEdx -> Fill (dEdx); // QA fill
		    h2logPdEdx -> Fill (TMath::Log(p), dEdx); // QA fill
			h2PqdEdx -> Fill (p * charge, dEdx); // QA fill
			h2PidCharge -> Fill (pid, charge); // QA fill
		    h2MTPCtot -> Fill (fPrimaryParticles_fTmeanCharge [itrack] [2] / 1000.0, fPrimaryParticles_fTmeanCharge [itrack] [3] / 1000.0);


            event_ -> AddTrack (pt, eta, phi, charge, pid);
            event_ -> GetTrack (trackIndex) -> SetdEdx_MTPC (fPrimaryParticles_fTmeanCharge [itrack] [2] / 1000.0);
            event_ -> GetTrack (trackIndex) -> SetdEdx_full (fPrimaryParticles_fTmeanCharge [itrack] [3] / 1000.0);
            event_ -> GetTrack (trackIndex) -> SetP (p);
            event_ -> GetTrack (trackIndex) -> SetRap (rap);

		    hPt -> Fill (pt);
		    hEta -> Fill (eta);
		    hTheta -> Fill (theta);
		    hPhi -> Fill (phi);
		    h2EtaPt -> Fill (eta, pt);
		    h2RapPt -> Fill (rap, pt);
		    h2PhiPt -> Fill (phi, pt);
		    h2EtaPhi -> Fill (eta, phi);
		    h2RapPhi -> Fill (rap, phi);

			if (charge > 0.0) {
                hPtPos -> Fill (pt);
                hEtaPos -> Fill (eta);
                hRapPos -> Fill (rap);
                hThetaPos -> Fill (theta);
                hPhiPos -> Fill (phi);
                h2EtaPtPos -> Fill (eta, pt);
                h2RapPtPos -> Fill (rap, pt);
                h2PhiPtPos -> Fill (phi, pt);
                h2EtaPhiPos -> Fill (eta, phi);
                h2RapPhiPos -> Fill (rap, phi);
			}

			if (charge < 0.0) {
                hPtNeg -> Fill (pt);
                hEtaNeg -> Fill (eta);
                hRapNeg -> Fill (rap);
                hThetaNeg -> Fill (theta);
                hPhiNeg -> Fill (phi);
                h2EtaPtNeg -> Fill (eta, pt);
                h2RapPtNeg -> Fill (rap, pt);
                h2PhiPtNeg -> Fill (phi, pt);
                h2EtaPhiNeg -> Fill (eta, phi);
                h2RapPhiNeg -> Fill (rap, phi);
			}

			if (pid == kPionMinus || pid == kPionPlus) {
                h2logPdEdxPi -> Fill (TMath::Log(p), dEdx);
                h2PqdEdxPi -> Fill (p * charge, dEdx);
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

			if (pid == kElectron || pid == kPositron) {
                h2logPdEdxElectron -> Fill (TMath::Log(p), dEdx);
                h2PqdEdxElectron -> Fill (p * charge, dEdx);
			}

			if (pid == kElectron) {
                hPtElectron -> Fill (pt);
                hEtaElectron -> Fill (eta);
                hRapElectron -> Fill (rap);
                hThetaElectron -> Fill (theta);
                hPhiElectron -> Fill (phi);
                h2EtaPtElectron -> Fill (eta, pt);
                h2RapPtElectron -> Fill (rap, pt);
                h2PhiPtElectron -> Fill (phi, pt);
                h2EtaPhiElectron -> Fill (eta, phi);
                h2RapPhiElectron -> Fill (rap, phi);
			}

			if (pid == kPositron) {
                hPtPositron -> Fill (pt);
                hEtaPositron -> Fill (eta);
                hRapPositron -> Fill (rap);
                hThetaPositron -> Fill (theta);
                hPhiPositron -> Fill (phi);
                h2EtaPtPositron -> Fill (eta, pt);
                h2RapPtPositron -> Fill (rap, pt);
                h2PhiPtPositron -> Fill (phi, pt);
                h2EtaPhiPositron -> Fill (eta, phi);
                h2RapPhiPositron -> Fill (rap, phi);
			}

			if (pid == kProton || pid == kAntiProton) {
                h2logPdEdxProton -> Fill (TMath::Log(p), dEdx);
                h2PqdEdxProton -> Fill (p * charge, dEdx);
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

			if (pid == kAntiProton) {
                hPtAntiProton -> Fill (pt);
                hEtaAntiProton -> Fill (eta);
                hRapAntiProton -> Fill (rap);
                hThetaAntiProton -> Fill (theta);
                hPhiAntiProton -> Fill (phi);
                h2EtaPtAntiProton -> Fill (eta, pt);
                h2RapPtAntiProton -> Fill (rap, pt);
                h2PhiPtAntiProton -> Fill (phi, pt);
                h2EtaPhiAntiProton -> Fill (eta, phi);
                h2RapPhiAntiProton -> Fill (rap, phi);
			}
		}
        // END QA FILL

        event_ -> AddTrack (ptVeto, etaVeto, phiVeto [0], chargeVeto, kVeto1);
        event_ -> GetTrack (trackIndex + 1) -> SetRap (rapVeto);
        h2PidCharge -> Fill (kVeto1, chargeVeto);
        event_ -> AddTrack (ptVeto, etaVeto, phiVeto [1], chargeVeto, kVeto2);
        event_ -> GetTrack (trackIndex + 2) -> SetRap (rapVeto);
        h2PidCharge -> Fill (kVeto2, chargeVeto);
        event_ -> AddTrack (ptVeto, etaVeto, phiVeto [2], chargeVeto, kVeto3);
        event_ -> GetTrack (trackIndex + 3) -> SetRap (rapVeto);
        h2PidCharge -> Fill (kVeto3, chargeVeto);
        event_ -> AddTrack (ptVeto, etaVeto, phiVeto [3], chargeVeto, kVeto4);
        event_ -> GetTrack (trackIndex + 4) -> SetRap (rapVeto);
        h2PidCharge -> Fill (kVeto4, chargeVeto);

		outputTree_ -> Fill ();
		event_ -> Clear ();
	}
/*
	TCanvas *c1 = new TCanvas ("c1", "c1", 800, 600);
    gStyle -> SetLegendBorderSize (0);

	TLegend *leg = new TLegend (0.7, 0.2, 0.85, 0.5);
	hPtNeg -> SetLineColor (kRed);
	hEtaNeg -> SetLineColor (kRed);
	hPhiNeg -> SetLineColor (kRed);
	leg -> SetFillColor (0);
    leg -> SetTextSize(0.04);
	leg -> AddEntry (hPtPos, "positive particles", "l");
	leg -> AddEntry (hPtNeg, "negative particles", "l");
    THStack *hsPt = new THStack ("hsPt", "P_{T} distribution for positive and negative particles;P_{T};nTracks");
    THStack *hsEta = new THStack ("hsEta", "#eta distribution for positive and negative particles;#eta;nTracks");
    THStack *hsPhi = new THStack ("hsPhi", "#phi distribution for positive and negative particles;#phi;nTracks");

    c1 -> SetLogy ();
    hsPt -> Add (hPtPos);
    hsPt -> Add (hPtNeg);
    hsPt -> Draw ("nostack");
    leg -> Draw ();
    c1 -> Write ("hPtPosNeg");

    c1 -> SetLogy (0);
    hsEta -> Add (hEtaPos);
    hsEta -> Add (hEtaNeg);
    hsEta -> Draw ("nostack");
    leg -> Draw ();
    c1 -> Write ("hEtaPosNeg");

    hsPhi -> Add (hPhiPos);
    hsPhi -> Add (hPhiNeg);
    hsPhi -> Draw ("nostack");
    leg -> Draw ();
    c1 -> Write ("hPhiPosNeg");

    TH2F *histList [4];
	TCanvas *c2 = new TCanvas ("c2", "c2", 800, 600);
    c2 -> Divide (2, 2);

    histList [0] = h2EtaPtPos;
    histList [1] = h2EtaPtPiPlus;
    histList [2] = h2EtaPtProton;
    histList [3] = h2EtaPtPositron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        gPad -> SetLogz ();
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("EtaPtPos");

    histList [0] = h2EtaPtNeg;
    histList [1] = h2EtaPtPiMin;
    histList [2] = h2EtaPtAntiProton;
    histList [3] = h2EtaPtElectron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("EtaPtNeg");

    histList [0] = h2PhiPtPos;
    histList [1] = h2PhiPtPiPlus;
    histList [2] = h2PhiPtProton;
    histList [3] = h2PhiPtPositron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("PhiPtPos");

    histList [0] = h2PhiPtNeg;
    histList [1] = h2PhiPtPiMin;
    histList [2] = h2PhiPtAntiProton;
    histList [3] = h2PhiPtElectron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("PhiPtNeg");

    histList [0] = h2EtaPhiPos;
    histList [1] = h2EtaPhiPiPlus;
    histList [2] = h2EtaPhiProton;
    histList [3] = h2EtaPhiPositron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("EtaPhiPos");

    histList [0] = h2EtaPhiNeg;
    histList [1] = h2EtaPhiPiMin;
    histList [2] = h2EtaPhiAntiProton;
    histList [3] = h2EtaPhiElectron;

    for (Int_t j = 0; j < 4; j++) {
        c2 -> cd (j + 1);
        histList [j] -> Draw ("colz");
    }
    c2 -> Write ("EtaPhiNeg");
*/
	histFile -> Write ();
	Finish ();
	return 1;
}


void CTreeConverter::Finish () {
	outputFile_ -> Write ();
	outputFile_ -> Close ();
}

#endif // CTREECONVERTER_CXX
