#ifndef CFLOWRECONSTRUCTOR_CXX
#define CFLOWRECONSTRUCTOR_CXX

#include <TMath.h>
#include <Riostream.h>
#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <climits>
#include <vector>
#include <float.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <THStack.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <cassert>

#include "../QnCorrections/QnCorrectionsLog.h"
#include "../QnCorrections/QnCorrectionsEventClassVariablesSet.h"
#include "../QnCorrections/QnCorrectionsCutAbove.h"
#include "../QnCorrections/QnCorrectionsCutBelow.h"
#include "../QnCorrections/QnCorrectionsCutOutside.h"
#include "../QnCorrections/QnCorrectionsCutsSet.h"
#include "../QnCorrections/QnCorrectionsCutSetBit.h"
#include "../QnCorrections/QnCorrectionsCutValue.h"
#include "../QnCorrections/QnCorrectionsCutWithin.h"
#include "../QnCorrections/QnCorrectionsHistogram.h"
#include "../QnCorrections/QnCorrectionsHistogramChannelized.h"
#include "../QnCorrections/QnCorrectionsProfile.h"
#include "../QnCorrections/QnCorrectionsProfile3DCorrelations.h"
#include "../QnCorrections/QnCorrectionsProfileChannelized.h"
#include "../QnCorrections/QnCorrectionsProfileChannelizedIngress.h"
#include "../QnCorrections/QnCorrectionsProfileComponents.h"
#include "../QnCorrections/QnCorrectionsProfileCorrelationComponents.h"
#include "../QnCorrections/QnCorrectionsProfileCorrelationComponentsHarmonics.h"
#include "../QnCorrections/QnCorrectionsDataVector.h"
#include "../QnCorrections/QnCorrectionsDataVectorChannelized.h"
#include "../QnCorrections/QnCorrectionsQnVector.h"
#include "../QnCorrections/QnCorrectionsQnVectorBuild.h"
#include "../QnCorrections/QnCorrectionsCorrectionsSetOnInputData.h"
#include "../QnCorrections/QnCorrectionsCorrectionsSetOnQvector.h"
#include "../QnCorrections/QnCorrectionsCorrectionOnInputData.h"
#include "../QnCorrections/QnCorrectionsCorrectionOnQvector.h"
#include "../QnCorrections/QnCorrectionsDetector.h"
#include "../QnCorrections/QnCorrectionsDetectorConfigurationsSet.h"
#include "../QnCorrections/QnCorrectionsDetectorConfigurationChannels.h"
#include "../QnCorrections/QnCorrectionsDetectorConfigurationTracks.h"
#include "../QnCorrections/QnCorrectionsInputGainEqualization.h"
#include "../QnCorrections/QnCorrectionsQnVectorRecentering.h"
#include "../QnCorrections/QnCorrectionsQnVectorAlignment.h"
#include "../QnCorrections/QnCorrectionsQnVectorTwistAndRescale.h"

#include "../MyDataTree/CEvent.h"
#include "../MyDataTree/CTrack.h"
#include "../Accessory/accessory.h"
#include "CFlowReconstructor.h"

using namespace std;


CFlowReconstructor::CFlowReconstructor () {
	nHarmonics = 0;
	SetNbinsMh (10);
	SetNbinsCent (10);
	SetNbinsPt (10);
	SetNbinsEta (10);
	SetNonUniformInputFileName ("../Acceptance");
	SetComment ("");
	SetNbinsEtaRefl (0);
	SetNrunRange (0, 0);
	harmonicFunctionSet = 0;
	uniformSet = 0;
	useAutoHistRanges_ = 0;
	useZeroSubevents_ = 0;
	SetVariable ();
	resChargeSet = 0;
	mhRangeForFlowSet_ = 0;
	centRangeForFlowSet_ = 0;
	propagateResolutionSign_ = 0;
	samplingMethod_ = kNoSampling;
	SetResolutionMethod (kThreeSubevents);
	SetNsteps (3);
}

void CFlowReconstructor::SetHistFileName (TString histFileName) {
    histFileName_ = histFileName;
}

void CFlowReconstructor::SetNsteps (Int_t nSteps) {
    nSteps_ = nSteps;
}

void CFlowReconstructor::SetComment (TString comment) {
    comment_ = comment;
}

bool CFlowReconstructor::AddHarmonic (Int_t n) {
	Int_t* temp;
	temp = harmonicsMap;
	harmonicsMap = new Int_t [nHarmonics + 1];
	for (int i = 0; i < nHarmonics; i++) {
        if (temp [i] == n) {
            cout << "Harmonic was already added!";
            return 0;
        }
		harmonicsMap [i] = temp [i];
	}
	harmonicsMap [nHarmonics] = n;
	nHarmonics++;
	//delete [] temp; // HELP!!!
	etaLim_.push_back (new Float_t [6]);
	ptLim_.push_back (new Float_t [6]);
	resSign_.push_back (new Int_t [3]);
	SetResolutionSigns(n, 1, 1, 1);
    resParticles.push_back (new vector <int> [3]);
	etaAveragingRange_.push_back (new Float_t [2]);
	ptAveragingRange_.push_back (new Float_t [2]);
	SetReferenceOption (n, "");
	return 1;
}

void CFlowReconstructor::SetVariable (TString var) {
    if (var == "y") varName_ = "#it{y}";
    else varName_ = "#eta";
}

void CFlowReconstructor::UseZeroSubevents (Bool_t useZeroSubevents) {
    useZeroSubevents_ = useZeroSubevents;
}

void CFlowReconstructor::SetSamplingMethod (Int_t samplingMethod) {
    samplingMethod_ = samplingMethod;
}

void CFlowReconstructor::PropagateResolutionSign (Bool_t propagateResolutionSign) {
    propagateResolutionSign_ = propagateResolutionSign;
}

void CFlowReconstructor::UseAutoHistRanges (Bool_t useAutoHistRanges) {
    useAutoHistRanges_ = useAutoHistRanges;
}

void CFlowReconstructor::ExcludeRun (Int_t nRun) {
    excludedRuns_.push_back (nRun);
}

void CFlowReconstructor::SetNrunRange (Int_t nRunMin, Int_t nRunMax) {
    nRunMin_ = nRunMin;
    nRunMax_ = nRunMax;
    nRuns_ = nRunMax - nRunMin + 1;
}

void CFlowReconstructor::SetMhRange (Int_t mhMin, Int_t mhMax) {
    mhMin_ = mhMin;
    mhMax_ = mhMax;
}

void CFlowReconstructor::SetMhRangeForFlow (Int_t mhLow, Int_t mhHigh) {
    mhRangeForFlowSet_ = 1;
    mhLow_ = mhLow;
    mhHigh_ = mhHigh;
}

void CFlowReconstructor::SetCentRange (Float_t centMin, Float_t centMax) {
    centMin_ = centMin;
    centMax_ = centMax;
}

void CFlowReconstructor::SetCentRangeForFlow (Float_t centLow, Float_t centHigh) {
    centRangeForFlowSet_ = 1;
    centLow_ = centLow;
    centHigh_ = centHigh;
}

void CFlowReconstructor::SetPtRange (Float_t ptMin, Float_t ptMax) {
    ptMin_ = ptMin;
    ptMax_ = ptMax;
}

void CFlowReconstructor::SetEtaRange (Float_t etaMin, Float_t etaMax) {
    etaMin_ = etaMin;
    etaMax_ = etaMax;
}

void CFlowReconstructor::SetResolutionCharge (Int_t charge) {
    resCharge = charge;
    resChargeSet = 1;
}

void CFlowReconstructor::AddFlowParticle (Int_t pid) {
    flowParticles.push_back (pid);
}

bool CFlowReconstructor::AddResolutionParticle (Int_t harmonic, Int_t subevent, Int_t pid) {
    Int_t n = -1;
    for (Int_t i = 0; i < nHarmonics; i++) {
        if (harmonic == harmonicsMap [i]) n = i;
    }
    if (n == -1) {
        cout << "Harmonic number " << harmonic << " was not introduced!\n";
        return 0;
    }
    resParticles [n][subevent - 1].push_back (pid);
    return 1;
}


bool CFlowReconstructor::SetResolutionSigns (Int_t harmonic, Int_t signA, Int_t signB, Int_t signC) {
    if (TMath::Abs (signA) != 1 || TMath::Abs (signB) != 1 || TMath::Abs (signC) != 1) {
        cout << "One of sign's absolute value is not equal to 1.0 for harmonic " << harmonic << endl;
        return 0;
    }
    Int_t n = -1;
    for (Int_t i = 0; i < nHarmonics; i++) {
        if (harmonic == harmonicsMap [i]) n = i;
    }
    if (n == -1) {
        cout << "Harmonic number " << harmonic << " was not introduced!\n";
        return 0;
    }
    resSign_ [n][0] = signA;
    resSign_ [n][1] = signB;
    resSign_ [n][2] = signC;
    return 1;
}


void CFlowReconstructor::SetNbinsBS (Int_t nBinsBS) {
    nBinsBS_ = nBinsBS;
}

void CFlowReconstructor::BuildSampleTree (TTree *inputTree) {
    cout << "Building sample tree!\n";
    TRandom3 r (0);
	Long64_t n, nEvents;
	Int_t W [1000];
	nEvents = inputTree -> GetEntries ();
	vector <Int_t*> events;

	TFile *sampleFile = new TFile (histFileName_ + "_sample.root", "RECREATE");
	TH2F *h2nEventSampleWeight = new TH2F ("h2nEventSampleWeight", "Event weights in samples;Nevent;Sample", nEvents, 0, nEvents, nBinsBS_, 0, nBinsBS_);
    TTree *sampleTree = new TTree ("samples", "samples");
    sampleTree -> Branch ("W", &W, "W[1000]/I");

	for (Long64_t i = 0; i < nEvents; i++) { // initialize with zeroes
        events.push_back (new Int_t [1000]);
        for (Long64_t j = 0; j < 1000; j++) {
            events [i][j] = 0;
        }
	}

    for (Int_t i = 0; i < nBinsBS_; i++) { // fill matrix
        for (Long64_t j = 0; j < nEvents; j++) {
            n = r.Rndm () * nEvents;
            events [n][i] ++;
        }
	}

	for (Long64_t i = 0; i < nEvents; i++) { // fill tree
        for (Int_t j = 0; j < 1000; j++) {
            W [j] = events [i][j];
        }
        for (Int_t j = 0; j < nBinsBS_; j++) {
            h2nEventSampleWeight -> Fill (i, j, W [j]);
        }
        sampleTree -> Fill ();
	}

	sampleFile -> Write ();
	sampleFile -> Close ();
	cout << "Finished building sample tree!\n";
}

void CFlowReconstructor::AnalyzeTree () {
//	Int_t mh, n, sign [nHarmonics][nBinsEta_], *flagPtEta, *flagPt, *flagEta, currFlagPtEta, currFlagPt, currFlagEta;
//	Int_t pid, charge;
//	Float_t cent, pt, eta, phi, x, y, X, Y;
//	TH1F *hMh, *hCent, *hPt, *hEta, *hPhi, *hMhPt, *hMhEta;
//	TH2F *h2PtEta, *h2MhPtEta, *h2PidCharge;
//	TH1F *hXnPt [nHarmonics], *hXnEta [nHarmonics], *hYnPt [nHarmonics], *hYnEta [nHarmonics];
//	TProfile *pxnXnPt [nHarmonics], *pynYnPt [nHarmonics], *pxnXnEta [nHarmonics], *pynYnEta [nHarmonics];
//	TH2F *h2XnPtEta [nHarmonics], *h2YnPtEta [nHarmonics];
//	TProfile2D *p2xnXnPtEta [nHarmonics], *p2ynYnPtEta [nHarmonics];
//
//	event = new CEvent;
//	inputFile = new TFile (nonUniformInputFileName + ".root", "READ");
//	inputTree = (TTree*) inputFile -> Get ("Tree");
//	inputTree -> SetBranchAddress ("Event", &event);
//	histFile = new TFile (histFileName_ + "_distr.root", "RECREATE");
//	GetVariableRanges ();
//
//	hMh = new TH1F ("hMh", "Multiplicity distribution; mh; nEvents", mhMax_, mhMin_, mhMax_);
//	hCent = new TH1F ("hCent", "Centrality distribution; cent; nEvents", 13, 0.0, 0.65); // temporarily to be configured manually!!!
//	hPt = new TH1F ("hPt", "P_{T} distribution; P_{T}[GeV]; nTracks", 100, ptMin_, ptMax_);
//	hEta = new TH1F ("hEta", varName_ + " distribution; " + varName_ + "; nTracks", 100, etaMin_, etaMax_);
//	h2PtEta = new TH2F ("h2PtEta", "P_{T} vs " + varName_ + " distribution; P_{T}[GeV]; " + varName_ + "; nTracks", 100, ptMin_, ptMax_, 100, etaMin_, etaMax_);
//	hPhi = new TH1F ("hPhi", "#phi distribution; #phi; nTracks", 100, 0.0, 2 * PI);
//	hMhPt = new TH1F ("hMhPt", "Multiplicity distribution versus P_{T}; P_{T}[GeV]; mh", nBinsPt_, ptMin_, ptMax_);
//	hMhEta = new TH1F ("hMhEta", "Multiplicity distribution versus " + varName_ + "; " + varName_ + "; mh", nBinsEta_, etaMin_, etaMax_);
//	h2MhPtEta = new TH2F ("h2MhPtEta", "Multiplicity distribution versus P_{T} and " + varName_ + "; P_{T}[GeV]; " + varName_ + "; mh", nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
//	h2PidCharge = new TH2F ("h2PidCharge", "Particle ID vs charge; pID; charge", kNPartTypes, -0.5, kNPartTypes - 0.5, 3, -1.5, 1.5);
//	for (int i = 0; i < nHarmonics; i++) {
//		n = harmonicsMap [i];
//		pxnXnPt [i] = new TProfile (Form ("px%iX%iPt", n, n), Form ("x_{%i}X_{%i} versus P_{T}; P_{T}[GeV]", n, n), nBinsPt_, ptMin_, ptMax_);
//		pynYnPt [i] = new TProfile (Form ("py%iY%iPt", n, n), Form ("y_{%i}Y_{%i} versus P_{T}; P_{T}[GeV]", n, n), nBinsPt_, ptMin_, ptMax_);
//		pxnXnEta [i] = new TProfile (Form ("px%iX%iEta", n, n), Form ("x_{%i}X_{%i} versus ", n, n) + varName_ + "; " + varName_, nBinsEta_, etaMin_, etaMax_);
//		pynYnEta [i] = new TProfile (Form ("py%iY%iEta", n, n), Form ("y_{%i}Y_{%i} versus ", n, n) + varName_ + "; " + varName_, nBinsEta_, etaMin_, etaMax_);
//		p2xnXnPtEta [i] = new TProfile2D (Form ("p2x%iX%iPtEta", n, n), Form ("x_{%i}X_{%i} versus P_{T} and ", n, n) + varName_ + "; P_{T}[GeV]; " + varName_, nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
//		p2ynYnPtEta [i] = new TProfile2D (Form ("p2y%iY%iPtEta", n, n), Form ("y_{%i}Y_{%i} versus P_{T} and ", n, n) + varName_ + "; P_{T}[GeV]; " + varName_, nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
//		h2XnPtEta [i] = new TH2F (Form ("h2X%iPtEta", n), Form ("X_{%i} versus P_{T} and ", n) + varName_ + "; P_{T}[GeV]; " + varName_, nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
//		h2YnPtEta [i] = new TH2F (Form ("h2Y%iPtEta", n), Form ("Y_{%i} versus P_{T} and ", n) + varName_ + "; P_{T}[GeV]; " + varName_, nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
//		hXnPt [i] = new TH1F (Form ("hX%iPt", n), Form ("X_{%i} versus P_{T}; P_{T}[GeV]; X_{%i}", n, n), nBinsPt_, ptMin_, ptMax_);
//		hYnPt [i] = new TH1F (Form ("hY%iPt", n), Form ("Y_{%i} versus P_{T}; P_{T}[GeV]; Y_{%i}", n, n), nBinsPt_, ptMin_, ptMax_);
//		hXnEta [i] = new TH1F (Form ("hX%iEta", n), Form ("X_{%i} versus ", n) + varName_ + "; " + varName_ + Form ("; X_{%i}", n), nBinsEta_, etaMin_, etaMax_);
//		hYnEta [i] = new TH1F (Form ("hY%iEta", n), Form ("Y_{%i} versus ", n) + varName_ + "; " + varName_ + Form ("; Y_{%i}", n), nBinsEta_, etaMin_, etaMax_);
//	}
//
//	flagPtEta = new Int_t [mhMax_];
//	flagPt = new Int_t [mhMax_];
//	flagEta = new Int_t [mhMax_];
//
//	Float_t h = (etaMax_ - etaMin_) / nBinsEta_;
//	eta = etaMin_ + 0.5 * h;
//	for (Int_t i = 0; i < nHarmonics; i++) {
//		for (Int_t k = 0; k < nBinsEta_; k++) {
//			sign [i][k] = fabs (eta) / eta;
//			eta += h;
//		}
//	}
//
//    cout << endl;
//
//	Long64_t nEvents = inputTree -> GetEntries ();
//	for (Long64_t jentry = 0; jentry < nEvents; jentry++) {
//		inputTree -> GetEntry (jentry);
//		cout << "\rEvent " << jentry + 1 << " from " << nEvents;
//		mh = event -> GetMh ();
//		cent = event -> GetCent ();
//		hMh -> Fill (mh);
//		hCent -> Fill (cent);
//
//		for (Int_t itrack = 1; itrack <= mh; itrack++) {
//			track = event -> GetTrack (itrack);
//			pt = track -> GetPt ();
//			eta = track -> GetEta ();
//			phi = track -> GetPhi ();
//			charge = track -> GetCharge ();
//			pid = track -> GetPid ();
//			flagPtEta [itrack - 1] = h2XnPtEta [0] -> FindBin (pt, eta);
//			flagPt [itrack - 1] = hXnPt [0] -> FindBin (pt);
//			flagEta [itrack - 1] = hXnEta [0] -> FindBin (eta);
//
//            hPt -> Fill (pt);
//            hEta -> Fill (eta);
//			h2PtEta -> Fill (pt, eta);
//            hPhi -> Fill (phi);
//			hMhPt -> Fill (pt);
//			hMhEta -> Fill (eta);
//			h2MhPtEta -> Fill (pt, eta);
//            h2PidCharge -> Fill (pid, charge);
//
//			for (Int_t i = 0; i < nHarmonics; i++) {
//				n = harmonicsMap [i];
//				x = TMath::Cos (n * phi);
//				y = TMath::Sin (n * phi);
//				h2XnPtEta [i] -> Fill (pt, eta, x);
//				h2YnPtEta [i] -> Fill (pt, eta, y);
//				hXnPt [i] -> Fill (pt, x);
//				hYnPt [i] -> Fill (pt, y);
//				hXnEta [i] -> Fill (eta, x);
//				hYnEta [i] -> Fill (eta, y);
//			}
//		}
//
//		for (Int_t i = 0; i < nHarmonics; i++) {
//			n = harmonicsMap [i];
//
//			for (Int_t j = 1; j <= nBinsPt_; j++) {
//				for (Int_t k = 1; k <= nBinsEta_; k++) {
//					X = h2XnPtEta [i] -> GetBinContent (j, k);
//					Y = h2YnPtEta [i] -> GetBinContent (j, k);
//					currFlagPtEta = h2XnPtEta [i] -> GetBin (j, k);
//					for (Int_t itrack = 1; itrack <= mh; itrack++) {
//						track = event -> GetTrack (itrack);
//						pt = track -> GetPt ();
//						eta = track -> GetEta ();
//						phi = track -> GetPhi ();
//						x = TMath::Cos (n * phi);
//						y = TMath::Sin (n * phi);
//						if (flagPtEta [itrack - 1] == currFlagPtEta) {
//							p2xnXnPtEta [i] -> Fill (pt, eta, sign [i][k] * x * (X - x));
//							p2ynYnPtEta [i] -> Fill (pt, eta, sign [i][k] * y * (Y - y));
//						}
//						else {
//							p2xnXnPtEta [i] -> Fill (pt, eta, sign [i][k] * x * X);
//							p2ynYnPtEta [i] -> Fill (pt, eta, sign [i][k] * y * Y);
//						}
//					}
//				}
//			}
//
//			for (Int_t j = 1; j <= nBinsPt_; j++) {
//				X = hXnPt [i] -> GetBinContent (j);
//				Y = hYnPt [i] -> GetBinContent (j);
//				currFlagPt = hXnPt [i] -> GetBin (j);
//				for (Int_t itrack = 1; itrack <= mh; itrack++) {
//					track = event -> GetTrack (itrack);
//					pt = track -> GetPt ();
//					eta = track -> GetEta ();
//					phi = track -> GetPhi ();
//					x = TMath::Cos (n * phi);
//					y = TMath::Sin (n * phi);
//					if (flagPt [itrack - 1] == currFlagPt) {
//						pxnXnPt [i] -> Fill (pt, x * (X - x));
//						pynYnPt [i] -> Fill (pt, y * (Y - y));
//					}
//					else {
//						pxnXnPt [i] -> Fill (pt, x * X);
//						pynYnPt [i] -> Fill (pt, y * Y);
//					}
//				}
//			}
//
//			for (Int_t k = 1; k <= nBinsEta_; k++) {
//				X = hXnEta [i] -> GetBinContent (k);
//				Y = hYnEta [i] -> GetBinContent (k);
//				currFlagEta = hXnEta [i] -> GetBin (k);
//				for (Int_t itrack = 1; itrack <= mh; itrack++) {
//					track = event -> GetTrack (itrack);
//					pt = track -> GetPt ();
//					eta = track -> GetEta ();
//					phi = track -> GetPhi ();
//					x = TMath::Cos (n * phi);
//					y = TMath::Sin (n * phi);
//					if (flagEta [itrack - 1] == currFlagEta) {
//						pxnXnEta [i] -> Fill (eta, sign [i][k] * x * (X - x));
//						pynYnEta [i] -> Fill (eta, sign [i][k] * y * (Y - y));
//					}
//					else {
//						pxnXnEta [i] -> Fill (eta, sign [i][k] * x * X);
//						pynYnEta [i] -> Fill (eta, sign [i][k] * y * Y);
//					}
//				}
//			}
//			h2XnPtEta [i] -> Reset ("ICESM");
//			h2YnPtEta [i] -> Reset ("ICESM");
//			hXnPt [i] -> Reset ("ICESM");
//			hYnPt [i] -> Reset ("ICESM");
//			hXnEta [i] -> Reset ("ICESM");
//			hYnEta [i] -> Reset ("ICESM");
//		}
//	}
//
//	for (Int_t i = 0; i < nHarmonics; i++) {
//		delete h2XnPtEta [i];
//		delete h2YnPtEta [i];
//		delete hXnPt [i];
//		delete hYnPt [i];
//		delete hXnEta [i];
//		delete hYnEta [i];
//	}
//	histFile -> Write ();
//	histFile -> Close ();
//	delete [] flagPtEta;
//	delete [] flagPt;
//	delete [] flagEta;
}

bool CFlowReconstructor::SetEtaSubeventsLimits (Int_t harmonic, Float_t lim1, Float_t lim2, Float_t lim3, Float_t lim4, Float_t lim5, Float_t lim6) {
	Int_t harmPointer = -1;
	for (int i = 0; i < nHarmonics; i++) {
		if (harmonicsMap [i] == harmonic) harmPointer = i;
	}
	if (harmPointer < 0) {
		cout << harmonic << "'th harmonic was not introduced!!!\n";
		return 0;
	}
	etaLim_ [harmPointer] [0] = lim1;
	etaLim_ [harmPointer] [1] = lim2;
	etaLim_ [harmPointer] [2] = lim3;
	etaLim_ [harmPointer] [3] = lim4;
	etaLim_ [harmPointer] [4] = lim5;
	etaLim_ [harmPointer] [5] = lim6;
	return 1;
}

bool CFlowReconstructor::SetEtaAveragingRange (Int_t harmonic, Float_t etaMin, Float_t etaMax) {
	Int_t harmPointer = -1;
	for (int i = 0; i < nHarmonics; i++) {
		if (harmonicsMap [i] == harmonic) harmPointer = i;
	}
	if (harmPointer < 0) {
		cout << harmonic << "'th harmonic was not introduced!!!\n";
		return 0;
	}
	etaAveragingRange_ [harmPointer] [0] = etaMin;
	etaAveragingRange_ [harmPointer] [1] = etaMax;
	return 1;
}

bool CFlowReconstructor::SetPtSubeventsLimits (Int_t harmonic, Float_t lim1, Float_t lim2, Float_t lim3, Float_t lim4, Float_t lim5, Float_t lim6) {
	Int_t harmPointer = -1;
	for (int i = 0; i < nHarmonics; i++) {
		if (harmonicsMap [i] == harmonic) harmPointer = i;
	}
	if (harmPointer < 0) {
		cout << harmonic << "'th harmonic was not introduced!!!\n";
		return 0;
	}
	ptLim_ [harmPointer] [0] = lim1;
	ptLim_ [harmPointer] [1] = lim2;
	ptLim_ [harmPointer] [2] = lim3;
	ptLim_ [harmPointer] [3] = lim4;
	ptLim_ [harmPointer] [4] = lim5;
	ptLim_ [harmPointer] [5] = lim6;
	return 1;
}

bool CFlowReconstructor::SetPtAveragingRange (Int_t harmonic, Float_t etaMin, Float_t etaMax) {
	Int_t harmPointer = -1;
	for (int i = 0; i < nHarmonics; i++) {
		if (harmonicsMap [i] == harmonic) harmPointer = i;
	}
	if (harmPointer < 0) {
		cout << harmonic << "'th harmonic was not introduced!!!\n";
		return 0;
	}
	ptAveragingRange_ [harmPointer] [0] = etaMin;
	ptAveragingRange_ [harmPointer] [1] = etaMax;
	return 1;
}

void CFlowReconstructor::SetNbinsMh (Int_t nBinsMh) {
	nBinsMh_ = nBinsMh;
}

void CFlowReconstructor::SetNbinsCent (Int_t nBinsCent) {
	nBinsCent_ = nBinsCent;
}

void CFlowReconstructor::SetNbinsPt (Int_t nBinsPt) {
	nBinsPt_ = nBinsPt;
}

void CFlowReconstructor::SetNbinsEta (Int_t nBinsEta) {
	nBinsEta_ = nBinsEta;
}

void CFlowReconstructor::SetNbinsEtaRefl (Int_t nBinsEtaRefl) {
    nBinsEtaRefl_ = nBinsEtaRefl;
}

void CFlowReconstructor::SetHarmonicFunction (Int_t n, floatFunction func){
	harmonicFunctions [n - 1] = func;
	harmonicFunctionSet = 1;
}

void CFlowReconstructor::GetVariableRanges (TTree *inputTree) {
	mhMin_ = INT_MAX;
	mhMax_ = 0;
	ptMin_ = FLT_MAX;
	ptMax_ = 0.0;
	etaMin_ = FLT_MAX;
	etaMax_ = -FLT_MAX;
	centMin_ = 1.0;
	centMax_ = 0.0;
	Float_t cent;
	Int_t mh;
	Float_t pt, eta;
	CEvent *event;
	CTrack *track;
	Long64_t nEvents = inputTree -> GetEntries ();
	for (Long64_t jentry = 0; jentry < nEvents; jentry++) {
		inputTree -> GetEntry (jentry);
		mh = event -> GetMh ();
		cent = event -> GetCent ();
		if (mh > mhMax_)
			mhMax_ = mh;
		if (mh < mhMin_)
			mhMin_ = mh;
		if (cent > centMax_)
			centMax_ = cent;
		if (cent < centMin_)
			centMin_ = cent;
		for (Int_t itrack = 1; itrack <= mh; itrack++) {
			track = event -> GetTrack (itrack);
			pt = track -> GetPt ();
			eta = track -> GetEta ();
			if (pt > ptMax_)
				ptMax_ = pt;
			if (pt < ptMin_)
				ptMin_ = pt;
			if (eta > etaMax_)
				etaMax_ = eta;
			if (eta < etaMin_)
				etaMin_ = eta;
		}
	}
	cout << "ptMin_ = " << ptMin_ << "\tetaMin_ = " << etaMin_ << "\tptMax_ = " << ptMax_ << "\tetaMax_ = " << etaMax_ << endl;
}

void CFlowReconstructor::HistShift (TH1* h, Float_t shiftX, Float_t shiftY) {
	Float_t xMin = h -> GetXaxis() -> GetXmin (), xMax = h -> GetXaxis() -> GetXmax ();
	Float_t deltaX = shiftX * h -> GetBinWidth (1);
	TF1 *f = new TF1 ("f", Form ("1 + %f * x", shiftY), xMin, xMax);
	h -> Multiply (f, 1);
	h -> GetXaxis () -> SetLimits (xMin + deltaX, xMax + deltaX);
}

void CFlowReconstructor::ShiftArray (Float_t *arr, Int_t n, Float_t shift) {
	Float_t delta = shift * (arr [n - 1] - arr [0]);
	for (Int_t i = 0; i < n; i++) arr [i] += delta;
}

void CFlowReconstructor::SetUniformInputFileName (TString name) {
    if (name == "") return;
    uniformInputFileName = name;
    uniformSet = 1;
}

void CFlowReconstructor::SetNonUniformInputFileName (TString name) {
	nonUniformInputFileName = name;
}

void CFlowReconstructor::TH2toTH1withSampling (TH2 *h2In, TH1 *hOut, TDirectory *dir) {
    Int_t nBinsDistr, j, k, maxBin, counter, nBinsFit = 20, nEmptyBins = 200;
    Float_t R, Rerr, xMin, xMax, xLow, xHigh, valuesFraction = 0.00;
    TDirectory *newDir;
    TH1F *hTemp;
    Double_t par [3];
    vector <Float_t> values;
    Int_t vsize;
    Int_t nBinsX = h2In -> GetNbinsX ();
    Int_t nBinsBS = h2In -> GetNbinsY ();
    if (dir != 0) newDir = dir -> mkdir (h2In -> GetName ());
    for (Int_t i = 1; i <= nBinsX; i++) {
        values.erase(values.begin(), values.end());

        for (j = 1; j <= nBinsBS; j++) {
            float d = h2In -> GetBinContent (i, j);
            if (h2In -> GetBinContent (i, j) > 0.0 || h2In -> GetBinContent (i, j) <= 0.0)
                values.push_back (h2In -> GetBinContent (i, j));
        }
        sort (values.begin(), values.end());
        vsize = values.size ();
        if (vsize == 0) {
            hOut -> SetBinContent (i, 1.0);
            hOut -> SetBinError (i, 1.0);
            continue;
        }
        xMin = values [0];
        xMax = values [vsize - 1];
        if (TMath::Abs (xMin) == TMath::Abs (xMax)) continue;

        xLow = values [(vsize - 1) * valuesFraction];
        xHigh = values [(vsize - 1) * (1.0 - valuesFraction)];
        nBinsDistr = (xMax - xMin) / (xHigh - xLow) * nBinsFit;

        hTemp = new TH1F (Form ("Bin_%i", i), Form ("Bin %i distribution", i), nBinsDistr, xMin, xMax);
        for (j = 1; j <= nBinsBS; j++) {
            hTemp -> Fill (h2In -> GetBinContent (i, j));
        }

//        hTemp -> Fit ("gaus", "", "", xLow, xHigh);
//        fit = hTemp -> GetFunction ("gaus");
//        fit -> GetParameters (&par [0]);
//        hOut -> SetBinContent (i, par[1]);
//        hOut -> SetBinError (i, par[2]);

        if (dir != 0) { // test
            newDir -> cd ();
            hTemp -> Write ();
        }

        counter = 0;
        maxBin = hTemp -> GetMaximumBin ();
        for (j = 1; counter < nEmptyBins && maxBin - j > 0; j++) {
            if (hTemp -> GetBinContent (maxBin - j) == 0) counter++;
            else counter = 0;
        }

        counter = 0;
        for (k = 1; counter < nEmptyBins && maxBin + k <= nBinsDistr; k++) {
            if (hTemp -> GetBinContent (maxBin + k) == 0) counter++;
            else counter = 0;
        }

        hTemp -> GetXaxis () -> SetRange (maxBin - j, maxBin + k);
        hOut -> SetBinContent (i, hTemp -> GetMean ());
        hOut -> SetBinError (i, hTemp -> GetRMS ());

        if (dir != 0) {
            newDir -> cd ();
            hTemp -> Write ();
        }
        delete hTemp;
    }
}


void CFlowReconstructor::GetCorrelations () {
    for (Int_t i = 0; i < nSteps_; i++) {
        GetCorrelationsLoop (i);
    }
}


void CFlowReconstructor::SetupQnCorrectionsManager (TFile *qnInputFile, TFile *qnPtInputFile, TFile *qnEtaInputFile, QnCorrectionsManager *QnMan, QnCorrectionsManager *QnManPt, QnCorrectionsManager *QnManEta) {
	const Int_t nEventClassesDimensions = 2;
	Int_t n;
	QnMan -> SetCalibrationHistogramsList (qnInputFile);
	QnManPt -> SetCalibrationHistogramsList (qnPtInputFile);
	QnManEta -> SetCalibrationHistogramsList (qnEtaInputFile);
	QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet (nEventClassesDimensions);
    CorrEventClasses -> Add (new QnCorrectionsEventClassVariable (kNrun, VarNames [kNrun], nRuns_, nRunMin_ - 0.5, nRunMax_ + 0.5));
	CorrEventClasses -> Add (new QnCorrectionsEventClassVariable (kCent, VarNames [kCent], nBinsCent_, centMin_, centMax_));

	QnCorrectionsDetector *myDetectorOne;

	vector <QnCorrectionsDetector**> uDetPt, uDetEta;
	vector <QnCorrectionsDetectorConfigurationTracks**> uDetPtConf, uDetEtaConf;
	vector <QnCorrectionsQnVectorTwistAndRescale**> twScalePt, twScaleEta;

	vector <QnCorrectionsDetector*> myDetectorOneA, myDetectorOneB, myDetectorOneC;
	QnCorrectionsDetectorConfigurationTracks *myDetectorOneConf;
	vector <QnCorrectionsDetectorConfigurationTracks*> myDetectorOneAConf, myDetectorOneBConf, myDetectorOneCConf;
	QnCorrectionsQnVectorTwistAndRescale *twScale1;
	vector <QnCorrectionsQnVectorTwistAndRescale*> twScale1A, twScale1B, twScale1C;

	myDetectorOne = new QnCorrectionsDetector (DetectorNames [kDetector1], kDetector1);
	myDetectorOneConf = new QnCorrectionsDetectorConfigurationTracks (
          "D1", CorrEventClasses, nHarmonics, harmonicsMap);
	myDetectorOneConf -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
	myDetectorOneConf -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
	twScale1 = new QnCorrectionsQnVectorTwistAndRescale();
	twScale1 -> SetApplyTwist (kTRUE);
	twScale1 -> SetApplyRescale (kTRUE);
	twScale1 -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
	myDetectorOneConf -> AddCorrectionOnQnVector (twScale1);
	myDetectorOne -> AddDetectorConfiguration (myDetectorOneConf);
	QnMan -> AddDetector (myDetectorOne);

	for (Int_t i = 0; i < nHarmonics; i++) {
		n = harmonicsMap [i];

		myDetectorOneA.push_back (new QnCorrectionsDetector (DetectorNames [kDetector1A] + Form ("_%i", n), kNDetectors * i + kDetector1A));
		myDetectorOneAConf.push_back (new QnCorrectionsDetectorConfigurationTracks (
			  Form ("D1A_%i", n), CorrEventClasses, nHarmonics, harmonicsMap));
		myDetectorOneAConf [i] -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
		myDetectorOneAConf [i] -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
		twScale1A.push_back (new QnCorrectionsQnVectorTwistAndRescale());
		twScale1A [i] -> SetApplyTwist (kTRUE);
		twScale1A [i] -> SetApplyRescale (kTRUE);
		twScale1A [i] -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
		myDetectorOneAConf [i] -> AddCorrectionOnQnVector (twScale1A [i]);
		myDetectorOneA [i] -> AddDetectorConfiguration (myDetectorOneAConf [i]);
		QnMan -> AddDetector (myDetectorOneA [i]);

		myDetectorOneB.push_back (new QnCorrectionsDetector (DetectorNames [kDetector1B] + Form ("_%i", n), kNDetectors * i + kDetector1B));
		myDetectorOneBConf.push_back (new QnCorrectionsDetectorConfigurationTracks (
			  Form ("D1B_%i", n), CorrEventClasses, nHarmonics, harmonicsMap));
		myDetectorOneBConf [i] -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
		myDetectorOneBConf [i] -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
		twScale1B.push_back (new QnCorrectionsQnVectorTwistAndRescale());
		twScale1B [i] -> SetApplyTwist (kTRUE);
		twScale1B [i] -> SetApplyRescale (kTRUE);
		twScale1B [i] -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
		myDetectorOneBConf [i] -> AddCorrectionOnQnVector (twScale1B [i]);
		myDetectorOneB [i] -> AddDetectorConfiguration (myDetectorOneBConf [i]);
		QnMan -> AddDetector (myDetectorOneB [i]);

		myDetectorOneC.push_back (new QnCorrectionsDetector (DetectorNames [kDetector1C] + Form ("_%i", n), kNDetectors * i + kDetector1C));
		myDetectorOneCConf.push_back (new QnCorrectionsDetectorConfigurationTracks (
			  Form ("D1C_%i", n), CorrEventClasses, nHarmonics, harmonicsMap));
		myDetectorOneCConf [i] -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
		myDetectorOneCConf [i] -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
		twScale1C.push_back (new QnCorrectionsQnVectorTwistAndRescale());
		twScale1C [i] -> SetApplyTwist (kTRUE);
		twScale1C [i] -> SetApplyRescale (kTRUE);
		twScale1C [i] -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
		myDetectorOneCConf [i] -> AddCorrectionOnQnVector (twScale1C [i]);
		myDetectorOneC [i] -> AddDetectorConfiguration (myDetectorOneCConf [i]);
		QnMan -> AddDetector (myDetectorOneC [i]);

        uDetPt.push_back (new QnCorrectionsDetector* [nBinsPt_]);
        uDetPtConf.push_back (new QnCorrectionsDetectorConfigurationTracks* [nBinsPt_]);
        twScalePt.push_back (new QnCorrectionsQnVectorTwistAndRescale* [nBinsPt_]);
        uDetEta.push_back (new QnCorrectionsDetector* [nBinsEta_]);
        uDetEtaConf.push_back (new QnCorrectionsDetectorConfigurationTracks* [nBinsEta_]);
        twScaleEta.push_back (new QnCorrectionsQnVectorTwistAndRescale* [nBinsEta_]);

        for (Int_t j = 0; j < nBinsPt_; j++) {
            uDetPt [i][j] = new QnCorrectionsDetector (Form ("DetPt%i", nBinsPt_ * i + j), nBinsPt_ * i + j);
            uDetPtConf [i][j] = new QnCorrectionsDetectorConfigurationTracks (
                  Form ("uDetPt%i", j + 1), CorrEventClasses, nHarmonics, harmonicsMap);
            uDetPtConf [i][j] -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
            uDetPtConf [i][j] -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
            twScalePt [i][j] = new QnCorrectionsQnVectorTwistAndRescale ();
            twScalePt [i][j] -> SetApplyTwist (kTRUE);
            twScalePt [i][j] -> SetApplyRescale (kTRUE);
            twScalePt [i][j] -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
            uDetPtConf [i][j] -> AddCorrectionOnQnVector (twScalePt [i][j]);
            uDetPt [i][j] -> AddDetectorConfiguration (uDetPtConf [i][j]);
            QnManPt -> AddDetector (uDetPt [i][j]);
        }

        for (Int_t j = 0; j < nBinsEta_; j++) {
            uDetEta [i][j] = new QnCorrectionsDetector (Form ("DetEta%i", nBinsEta_ * i + j), nBinsEta_ * i + j);
            uDetEtaConf [i][j] = new QnCorrectionsDetectorConfigurationTracks (
                  Form ("uDetEta%i", j + 1), CorrEventClasses, nHarmonics, harmonicsMap);
            uDetEtaConf [i][j] -> SetQVectorNormalizationMethod (QnCorrectionsQnVector::QVNORM_QoverM);
            uDetEtaConf [i][j] -> AddCorrectionOnQnVector (new QnCorrectionsQnVectorRecentering ());
            twScaleEta [i][j] = new QnCorrectionsQnVectorTwistAndRescale();
            twScaleEta [i][j] -> SetApplyTwist (kTRUE);
            twScaleEta [i][j] -> SetApplyRescale (kTRUE);
            twScaleEta [i][j] -> SetTwistAndRescaleMethod (QnCorrectionsQnVectorTwistAndRescale::TWRESCALE_doubleHarmonic);
            uDetEtaConf [i][j] -> AddCorrectionOnQnVector (twScaleEta [i][j]);
            uDetEta [i][j] -> AddDetectorConfiguration (uDetEtaConf [i][j]);
            QnManEta -> AddDetector (uDetEta [i][j]);
        }
	}

	cout << "\n================ CONFIGURED ================\n\n";

	/* order the appropriate output */
	QnMan -> SetShouldFillQAHistograms (kTRUE);
	QnMan -> SetShouldFillNveQAHistograms (kTRUE);
	QnMan -> SetShouldFillOutputHistograms (kTRUE);
	QnManPt -> SetShouldFillQAHistograms (kTRUE);
	QnManPt -> SetShouldFillNveQAHistograms (kTRUE);
	QnManPt -> SetShouldFillOutputHistograms (kTRUE);
	QnManEta -> SetShouldFillQAHistograms (kTRUE);
	QnManEta -> SetShouldFillNveQAHistograms (kTRUE);
	QnManEta -> SetShouldFillOutputHistograms (kTRUE);

	cout << "\n================ FINISH SETUP ================\n\n";

	QnMan -> InitializeQnCorrectionsFramework();
	QnManPt -> InitializeQnCorrectionsFramework();
	QnManEta -> InitializeQnCorrectionsFramework();

	cout << "\n================ INITIALIZED ================\n\n";

	QnMan -> SetCurrentProcessListName ("Qn");
	QnManPt -> SetCurrentProcessListName ("QnPt");
	QnManEta -> SetCurrentProcessListName ("QnEta");
}

void CFlowReconstructor::FinalizeQnCorrectionsManager (TFile *qnInputFile, TFile *qnOutputFile, QnCorrectionsManager *QnMan) {
	if (qnInputFile) qnInputFile -> Close ();
	QnMan -> FinalizeQnCorrectionsFramework ();
	qnOutputFile -> cd ();
	QnMan -> GetOutputHistogramsList () -> Write (QnMan -> GetOutputHistogramsList () -> GetName (), TObject::kSingleKey);
	QnMan -> GetQAHistogramsList () -> Write (QnMan -> GetQAHistogramsList () -> GetName (), TObject::kSingleKey);
	QnMan -> GetNveQAHistogramsList () -> Write (QnMan -> GetNveQAHistogramsList () -> GetName (), TObject::kSingleKey);
	qnOutputFile -> Close ();
	delete QnMan;
}


void CFlowReconstructor::GetCorrelationsLoop (Int_t step) {
    TString option [4] = {"recreate", "update", "update", "update"};
	TFile *inputFile = new TFile (nonUniformInputFileName + ".root", "READ");
	if (uniformSet) inputFile = new TFile (uniformInputFileName + ".root", "READ");
	TTree *inputTree = (TTree*) inputFile -> Get ("Tree");
	CEvent *event = new CEvent;
	CTrack* track;
	inputTree -> SetBranchAddress ("Event", &event);
	if (useAutoHistRanges_ == 1) GetVariableRanges (inputTree);
	TFile *histFile = new TFile (histFileName_ + "_corr.root", option [step]);
	if (samplingMethod_ == kBootStrapping) BuildSampleTree (inputTree);

	TFile *qnInputFile = new TFile (histFileName_ + Form ("_%i.root", step), "read");
	TFile *qnPtInputFile = new TFile (histFileName_ + Form ("Pt_%i.root", step), "read");
	TFile *qnEtaInputFile = new TFile (histFileName_ + Form ("Eta_%i.root", step), "read");
	if (uniformSet) {
        qnInputFile = new TFile (histFileName_ + Form ("_%i.root", 0), "read");
        qnPtInputFile = new TFile (histFileName_ + Form ("_%i.root", 0), "read");
        qnEtaInputFile = new TFile (histFileName_ + Form ("_%i.root", 0), "read");
    }
	TFile *qnOutputFile = new TFile (histFileName_ + Form ("_%i.root", step + 1), "recreate");
	TFile *qnPtOutputFile = new TFile (histFileName_ + Form ("Pt_%i.root", step + 1), "recreate");
	TFile *qnEtaOutputFile = new TFile (histFileName_ + Form ("Eta_%i.root", step + 1), "recreate");

	QnCorrectionsManager *QnMan = new QnCorrectionsManager ();
	QnCorrectionsManager *QnManPt = new QnCorrectionsManager ();
	QnCorrectionsManager *QnManEta = new QnCorrectionsManager ();

	SetupQnCorrectionsManager (qnInputFile, qnPtInputFile, qnEtaInputFile, QnMan, QnManPt, QnManEta);
	TDirectory *histDir = histFile -> mkdir (dirName [step]);
	TDirectory *servHistDir = histDir -> mkdir ("Source Histograms");

	Long64_t nEvents = inputTree -> GetEntries ();
    Int_t centBin;
    Bool_t skipFlag;
	Int_t *mha = new Int_t [nHarmonics], *mhb = new Int_t [nHarmonics], *mhc = new Int_t [nHarmonics];
	Int_t **subeventFlag;
	Int_t n, nRun, mh, charge, pid, nFlowParts = flowParticles.size (), bsIndex, sMax;
	Float_t cent, pt, eta, phi, m, x, y, sign = 1.0, subeventIndex, weight, *Eveto, summEveto;
//	Float_t ptAvg [nHarmonics][nBinsCent_], ptAvgA [nHarmonics][nBinsCent_], ptAvgB [nHarmonics][nBinsCent_], ptAvgC [nHarmonics][nBinsCent_];
	TRandom3 r (0);

	TFile *sampleFile;
	TTree *sampleTree;
    Int_t W [1000]; // maximum number of samples = 1000
    if (samplingMethod_ == kBootStrapping) {
        sampleFile = new TFile (histFileName_ + "_sample.root", "READ");
        sampleTree = (TTree*) sampleFile -> Get ("samples");
        sampleTree -> SetBranchAddress ("W", &W);
        sMax = nBinsBS_;
    }
    else sMax = 1;

	Float_t *X = new Float_t [nHarmonics], *Y = new Float_t [nHarmonics];
	Float_t *Xa = new Float_t [nHarmonics], *Ya = new Float_t [nHarmonics];
	Float_t *Xb = new Float_t [nHarmonics], *Yb = new Float_t [nHarmonics];
	Float_t *Xc = new Float_t [nHarmonics], *Yc = new Float_t [nHarmonics];
	Float_t *Q = new Float_t [nHarmonics], *Qa = new Float_t [nHarmonics], *Qb = new Float_t [nHarmonics], *Qc = new Float_t [nHarmonics];
	Float_t *psiEP = new Float_t [nHarmonics], *psiEPa = new Float_t [nHarmonics], *psiEPb = new Float_t [nHarmonics], *psiEPc = new Float_t [nHarmonics];
	Float_t *XRP = new Float_t [nHarmonics], *YRP = new Float_t [nHarmonics], *psiRP = new Float_t [nHarmonics];


    // test
        TFile *testFile = new TFile (histFileName_ + Form ("_Q_%i.root", step), "RECREATE");
        TTree *testTree = new TTree ("treeQ", "Test tree");
        testTree -> Branch ("nRun", &nRun);
        testTree -> Branch ("Xa", &Xa, "Xa[2]/F");
        testTree -> Branch ("Xb", &Xb, "Xb[2]/F");
        testTree -> Branch ("Xc", &Xc, "Xc[2]/F");
        testTree -> Branch ("Ya", &Ya, "Ya[2]/F");
        testTree -> Branch ("Yb", &Yb, "Yb[2]/F");
        testTree -> Branch ("Yc", &Yc, "Yc[2]/F");
        testTree -> Branch ("mha", &mha, "mha[2]/I");
        testTree -> Branch ("mhb", &mhb, "mhb[2]/I");
        testTree -> Branch ("mhc", &mhc, "mhc[2]/I");
        testTree -> Branch ("Nsub", &bsIndex, 32000, 4);
        testTree -> Branch ("cent", &cent, 32000, 4);
    // test


    servHistDir -> cd ();
    vector <TH2F*> h2nEventSampleWeight;
	vector <TH1F*> hMh, hMha, hMhb, hMhc;
	vector <TProfile*> pPtCent, pPtCentA, pPtCentB, pPtCentC;

	vector <TH1F*> hNeventsBS, hNtracksBS, hNtracksBSa, hNtracksBSb, hNtracksBSc;
	vector <TH2F*> h2mhCent, h2mhaCent, h2mhbCent, h2mhcCent;
	vector <TH2F*> h2mhMult, h2mhaMult, h2mhbMult, h2mhcMult;
    vector <TH2F*> h2PtEta, h2PtEtaA, h2PtEtaB, h2PtEtaC;
	vector <TProfile*> pPtEta, pPtEtaA, pPtEtaB, pPtEtaC;

    vector <TH2F*> hXaXRPCent_SP, hYaYRPCent_SP, hXbXRPCent_SP, hYbYRPCent_SP, hXcXRPCent_SP, hYcYRPCent_SP;
    vector <TH2F*> hXaXRPCent_EP, hYaYRPCent_EP, hXbXRPCent_EP, hYbYRPCent_EP, hXcXRPCent_EP, hYcYRPCent_EP;

    vector <TH2F*> hXaXRPMult_SP, hYaYRPMult_SP, hXbXRPMult_SP, hYbYRPMult_SP, hXcXRPMult_SP, hYcYRPMult_SP;
    vector <TH2F*> hXaXRPMult_EP, hYaYRPMult_EP, hXbXRPMult_EP, hYbYRPMult_EP, hXcXRPMult_EP, hYcYRPMult_EP;

	vector <TProfile*> pxXCent_RP;
	vector <TProfile*> pxXaCent_SP, pxXbCent_SP, pxXcCent_SP, pyYaCent_SP, pyYbCent_SP, pyYcCent_SP;
	vector <TProfile*> pyXaCent_SP, pyXbCent_SP, pyXcCent_SP, pxYaCent_SP, pxYbCent_SP, pxYcCent_SP;
	vector <TProfile*> pxXaCent_EP, pxXbCent_EP, pxXcCent_EP, pyYaCent_EP, pyYbCent_EP, pyYcCent_EP;
	vector <TProfile*> pyXaCent_EP, pyXbCent_EP, pyXcCent_EP, pxYaCent_EP, pxYbCent_EP, pxYcCent_EP;
	vector <TProfile2D*> p2xXaCent_SP, p2xXbCent_SP, p2xXcCent_SP, p2yYaCent_SP, p2yYbCent_SP, p2yYcCent_SP;
	vector <TProfile2D*> p2yXaCent_SP, p2yXbCent_SP, p2yXcCent_SP, p2xYaCent_SP, p2xYbCent_SP, p2xYcCent_SP;
	vector <TProfile2D*> p2xXaCent_EP, p2xXbCent_EP, p2xXcCent_EP, p2yYaCent_EP, p2yYbCent_EP, p2yYcCent_EP;
	vector <TProfile2D*> p2yXaCent_EP, p2yXbCent_EP, p2yXcCent_EP, p2xYaCent_EP, p2xYbCent_EP, p2xYcCent_EP;

	vector <TProfile*> pxXMult_RP;
	vector <TProfile*> pxXaMult_SP, pxXbMult_SP, pxXcMult_SP, pyYaMult_SP, pyYbMult_SP, pyYcMult_SP;
	vector <TProfile*> pyXaMult_SP, pyXbMult_SP, pyXcMult_SP, pxYaMult_SP, pxYbMult_SP, pxYcMult_SP;
	vector <TProfile*> pxXaMult_EP, pxXbMult_EP, pxXcMult_EP, pyYaMult_EP, pyYbMult_EP, pyYcMult_EP;
	vector <TProfile*> pyXaMult_EP, pyXbMult_EP, pyXcMult_EP, pxYaMult_EP, pxYbMult_EP, pxYcMult_EP;
	vector <TProfile2D*> p2xXaMult_SP, p2xXbMult_SP, p2xXcMult_SP, p2yYaMult_SP, p2yYbMult_SP, p2yYcMult_SP;
	vector <TProfile2D*> p2yXaMult_SP, p2yXbMult_SP, p2yXcMult_SP, p2xYaMult_SP, p2xYbMult_SP, p2xYcMult_SP;
	vector <TProfile2D*> p2xXaMult_EP, p2xXbMult_EP, p2xXcMult_EP, p2yYaMult_EP, p2yYbMult_EP, p2yYcMult_EP;
	vector <TProfile2D*> p2yXaMult_EP, p2yXbMult_EP, p2yXcMult_EP, p2xYaMult_EP, p2xYbMult_EP, p2xYcMult_EP;

	vector <TProfile*> pXaXbCent_SP, pXaXcCent_SP, pXbXcCent_SP, pYaYbCent_SP, pYaYcCent_SP, pYbYcCent_SP;
	vector <TProfile*> pXaYbCent_SP, pXaYcCent_SP, pXbYcCent_SP, pYaXbCent_SP, pYaXcCent_SP, pYbXcCent_SP;
	vector <TProfile*> pXaXbCent_EP, pXaXcCent_EP, pXbXcCent_EP, pYaYbCent_EP, pYaYcCent_EP, pYbYcCent_EP;
	vector <TProfile*> pXaYbCent_EP, pXaYcCent_EP, pXbYcCent_EP, pYaXbCent_EP, pYaXcCent_EP, pYbXcCent_EP;

	vector <TProfile*> pXaXbMult_SP, pXaXcMult_SP, pXbXcMult_SP, pYaYbMult_SP, pYaYcMult_SP, pYbYcMult_SP;
	vector <TProfile*> pXaYbMult_SP, pXaYcMult_SP, pXbYcMult_SP, pYaXbMult_SP, pYaXcMult_SP, pYbXcMult_SP;
	vector <TProfile*> pXaXbMult_EP, pXaXcMult_EP, pXbXcMult_EP, pYaYbMult_EP, pYaYcMult_EP, pYbYcMult_EP;
	vector <TProfile*> pXaYbMult_EP, pXaYcMult_EP, pXbYcMult_EP, pYaXbMult_EP, pYaXcMult_EP, pYbXcMult_EP;

	vector <TProfile2D*> p2XaXbCent_SP, p2XaXcCent_SP, p2XbXcCent_SP, p2YaYbCent_SP, p2YaYcCent_SP, p2YbYcCent_SP;
	vector <TProfile2D*> p2XaYbCent_SP, p2XaYcCent_SP, p2XbYcCent_SP, p2YaXbCent_SP, p2YaXcCent_SP, p2YbXcCent_SP;
	vector <TProfile2D*> p2XaXbCent_EP, p2XaXcCent_EP, p2XbXcCent_EP, p2YaYbCent_EP, p2YaYcCent_EP, p2YbYcCent_EP;
	vector <TProfile2D*> p2XaYbCent_EP, p2XaYcCent_EP, p2XbYcCent_EP, p2YaXbCent_EP, p2YaXcCent_EP, p2YbXcCent_EP;

	vector <TProfile2D*> p2XaXbMult_SP, p2XaXcMult_SP, p2XbXcMult_SP, p2YaYbMult_SP, p2YaYcMult_SP, p2YbYcMult_SP;
	vector <TProfile2D*> p2XaYbMult_SP, p2XaYcMult_SP, p2XbYcMult_SP, p2YaXbMult_SP, p2YaXcMult_SP, p2YbXcMult_SP;
	vector <TProfile2D*> p2XaXbMult_EP, p2XaXcMult_EP, p2XbXcMult_EP, p2YaYbMult_EP, p2YaYcMult_EP, p2YbYcMult_EP;
	vector <TProfile2D*> p2XaYbMult_EP, p2XaYcMult_EP, p2XbYcMult_EP, p2YaXbMult_EP, p2YaXcMult_EP, p2YbXcMult_EP;

    vector <TProfile2D*> p2xXaPtCent_SP, p2yYaPtCent_SP, p2yXaPtCent_SP, p2xYaPtCent_SP;
    vector <TProfile2D*> p2xXbPtCent_SP, p2yYbPtCent_SP, p2yXbPtCent_SP, p2xYbPtCent_SP;
    vector <TProfile2D*> p2xXcPtCent_SP, p2yYcPtCent_SP, p2yXcPtCent_SP, p2xYcPtCent_SP;
    vector <TProfile2D*> p2xXaPtCent_EP, p2yYaPtCent_EP, p2yXaPtCent_EP, p2xYaPtCent_EP;
    vector <TProfile2D*> p2xXbPtCent_EP, p2yYbPtCent_EP, p2yXbPtCent_EP, p2xYbPtCent_EP;
    vector <TProfile2D*> p2xXcPtCent_EP, p2yYcPtCent_EP, p2yXcPtCent_EP, p2xYcPtCent_EP;

    vector <TProfile3D*> p3xXaPtCent_SP, p3yYaPtCent_SP, p3yXaPtCent_SP, p3xYaPtCent_SP;
    vector <TProfile3D*> p3xXbPtCent_SP, p3yYbPtCent_SP, p3yXbPtCent_SP, p3xYbPtCent_SP;
    vector <TProfile3D*> p3xXcPtCent_SP, p3yYcPtCent_SP, p3yXcPtCent_SP, p3xYcPtCent_SP;
    vector <TProfile3D*> p3xXaPtCent_EP, p3yYaPtCent_EP, p3yXaPtCent_EP, p3xYaPtCent_EP;
    vector <TProfile3D*> p3xXbPtCent_EP, p3yYbPtCent_EP, p3yXbPtCent_EP, p3xYbPtCent_EP;
    vector <TProfile3D*> p3xXcPtCent_EP, p3yYcPtCent_EP, p3yXcPtCent_EP, p3xYcPtCent_EP;

    vector <TProfile2D*> p2xXaPtMult_SP, p2yYaPtMult_SP, p2yXaPtMult_SP, p2xYaPtMult_SP;
    vector <TProfile2D*> p2xXbPtMult_SP, p2yYbPtMult_SP, p2yXbPtMult_SP, p2xYbPtMult_SP;
    vector <TProfile2D*> p2xXcPtMult_SP, p2yYcPtMult_SP, p2yXcPtMult_SP, p2xYcPtMult_SP;
    vector <TProfile2D*> p2xXaPtMult_EP, p2yYaPtMult_EP, p2yXaPtMult_EP, p2xYaPtMult_EP;
    vector <TProfile2D*> p2xXbPtMult_EP, p2yYbPtMult_EP, p2yXbPtMult_EP, p2xYbPtMult_EP;
    vector <TProfile2D*> p2xXcPtMult_EP, p2yYcPtMult_EP, p2yXcPtMult_EP, p2xYcPtMult_EP;

    vector <TProfile3D*> p3xXaPtMult_SP, p3yYaPtMult_SP, p3yXaPtMult_SP, p3xYaPtMult_SP;
    vector <TProfile3D*> p3xXbPtMult_SP, p3yYbPtMult_SP, p3yXbPtMult_SP, p3xYbPtMult_SP;
    vector <TProfile3D*> p3xXcPtMult_SP, p3yYcPtMult_SP, p3yXcPtMult_SP, p3xYcPtMult_SP;
    vector <TProfile3D*> p3xXaPtMult_EP, p3yYaPtMult_EP, p3yXaPtMult_EP, p3xYaPtMult_EP;
    vector <TProfile3D*> p3xXbPtMult_EP, p3yYbPtMult_EP, p3yXbPtMult_EP, p3xYbPtMult_EP;
    vector <TProfile3D*> p3xXcPtMult_EP, p3yYcPtMult_EP, p3yXcPtMult_EP, p3xYcPtMult_EP;

    vector <TProfile2D*> p2xXaEtaCent_SP, p2yYaEtaCent_SP, p2yXaEtaCent_SP, p2xYaEtaCent_SP;
    vector <TProfile2D*> p2xXbEtaCent_SP, p2yYbEtaCent_SP, p2yXbEtaCent_SP, p2xYbEtaCent_SP;
    vector <TProfile2D*> p2xXcEtaCent_SP, p2yYcEtaCent_SP, p2yXcEtaCent_SP, p2xYcEtaCent_SP;
    vector <TProfile2D*> p2xXaEtaCent_EP, p2yYaEtaCent_EP, p2yXaEtaCent_EP, p2xYaEtaCent_EP;
    vector <TProfile2D*> p2xXbEtaCent_EP, p2yYbEtaCent_EP, p2yXbEtaCent_EP, p2xYbEtaCent_EP;
    vector <TProfile2D*> p2xXcEtaCent_EP, p2yYcEtaCent_EP, p2yXcEtaCent_EP, p2xYcEtaCent_EP;

    vector <TProfile3D*> p3xXaEtaCent_SP, p3yYaEtaCent_SP, p3yXaEtaCent_SP, p3xYaEtaCent_SP;
    vector <TProfile3D*> p3xXbEtaCent_SP, p3yYbEtaCent_SP, p3yXbEtaCent_SP, p3xYbEtaCent_SP;
    vector <TProfile3D*> p3xXcEtaCent_SP, p3yYcEtaCent_SP, p3yXcEtaCent_SP, p3xYcEtaCent_SP;
    vector <TProfile3D*> p3xXaEtaCent_EP, p3yYaEtaCent_EP, p3yXaEtaCent_EP, p3xYaEtaCent_EP;
    vector <TProfile3D*> p3xXbEtaCent_EP, p3yYbEtaCent_EP, p3yXbEtaCent_EP, p3xYbEtaCent_EP;
    vector <TProfile3D*> p3xXcEtaCent_EP, p3yYcEtaCent_EP, p3yXcEtaCent_EP, p3xYcEtaCent_EP;

    vector <TProfile2D*> p2xXaEtaMult_SP, p2yYaEtaMult_SP, p2yXaEtaMult_SP, p2xYaEtaMult_SP;
    vector <TProfile2D*> p2xXbEtaMult_SP, p2yYbEtaMult_SP, p2yXbEtaMult_SP, p2xYbEtaMult_SP;
    vector <TProfile2D*> p2xXcEtaMult_SP, p2yYcEtaMult_SP, p2yXcEtaMult_SP, p2xYcEtaMult_SP;
    vector <TProfile2D*> p2xXaEtaMult_EP, p2yYaEtaMult_EP, p2yXaEtaMult_EP, p2xYaEtaMult_EP;
    vector <TProfile2D*> p2xXbEtaMult_EP, p2yYbEtaMult_EP, p2yXbEtaMult_EP, p2xYbEtaMult_EP;
    vector <TProfile2D*> p2xXcEtaMult_EP, p2yYcEtaMult_EP, p2yXcEtaMult_EP, p2xYcEtaMult_EP;

    vector <TProfile3D*> p3xXaEtaMult_SP, p3yYaEtaMult_SP, p3yXaEtaMult_SP, p3xYaEtaMult_SP;
    vector <TProfile3D*> p3xXbEtaMult_SP, p3yYbEtaMult_SP, p3yXbEtaMult_SP, p3xYbEtaMult_SP;
    vector <TProfile3D*> p3xXcEtaMult_SP, p3yYcEtaMult_SP, p3yXcEtaMult_SP, p3xYcEtaMult_SP;
    vector <TProfile3D*> p3xXaEtaMult_EP, p3yYaEtaMult_EP, p3yXaEtaMult_EP, p3xYaEtaMult_EP;
    vector <TProfile3D*> p3xXbEtaMult_EP, p3yYbEtaMult_EP, p3yXbEtaMult_EP, p3xYbEtaMult_EP;
    vector <TProfile3D*> p3xXcEtaMult_EP, p3yYcEtaMult_EP, p3yXcEtaMult_EP, p3xYcEtaMult_EP;

	subeventFlag = new Int_t* [nHarmonics];

	for (Int_t i = 0; i < nHarmonics; i++) { // create histograms
		n = harmonicsMap [i];
		subeventFlag [i] = new Int_t [mhMax_];

		h2nEventSampleWeight.push_back (new TH2F (Form ("h2nEventSampleWeight_%i", n), Form ("Event weights in samples (n = %i);Nevent;Sample", n), nEvents, 0, nEvents, nBinsBS_, 0, nBinsBS_));
		pPtCent.push_back (new TProfile (Form ("pPtCent%i", n), Form ("#LTP_{T}#GT versus centrality, harmonic %i; #LTP_{T}#GT; centrality", n), nBinsCent_, centMin_, centMax_));
		pPtCentA.push_back (new TProfile (Form ("pPtCentA%i", n), Form ("#LTP_{T}#GT versus centrality, subevent a_{%i}; #LTP_{T}#GT; centrality", n), nBinsCent_, centMin_, centMax_));
		pPtCentB.push_back (new TProfile (Form ("pPtCentB%i", n), Form ("#LTP_{T}#GT versus centrality, subevent b_{%i}; #LTP_{T}#GT; centrality", n), nBinsCent_, centMin_, centMax_));
		pPtCentC.push_back (new TProfile (Form ("pPtCentC%i", n), Form ("#LTP_{T}#GT versus centrality, subevent c_{%i}; #LTP_{T}#GT; centrality", n), nBinsCent_, centMin_, centMax_));
		hMh.push_back (new TH1F (Form ("hMh%i", n), Form ("Multiplicity distribution, harmonic %i; mh; nEvents", n), 100, 0, mhMax_));
		hMha.push_back (new TH1F (Form ("hMha%i", n), Form ("Multiplicity distribution for subevent 'a_{%i}'; mh; nEvents", n), 100, 0, mhMax_));
		hMhb.push_back (new TH1F (Form ("hMhb%i", n), Form ("Multiplicity distribution for subevent 'b_{%i}'; mh; nEvents", n), 100, 0, mhMax_));
		hMhc.push_back (new TH1F (Form ("hMhc%i", n), Form ("Multiplicity distribution for subevent 'c_{%i}'; mh; nEvents", n), 100, 0, mhMax_));

        hNeventsBS.push_back (new TH1F (Form ("hNeventsBS%i", n), "Subsample distribution of Nevents;subsample;Nevents", nBinsBS_, 0, nBinsBS_));
        hNtracksBS.push_back (new TH1F (Form ("hNtracksBS%i", n), "Subsample distribution of Ntracks;subsample;Ntracks", nBinsBS_, 0, nBinsBS_));
        hNtracksBSa.push_back (new TH1F (Form ("hNtracksBSa%i", n), "Subsample distribution of Ntracks (subevent A);subsample;Ntracks", nBinsBS_, 0, nBinsBS_));
        hNtracksBSb.push_back (new TH1F (Form ("hNtracksBSb%i", n), "Subsample distribution of Ntracks (subevent B);subsample;Ntracks", nBinsBS_, 0, nBinsBS_));
        hNtracksBSc.push_back (new TH1F (Form ("hNtracksBSc%i", n), "Subsample distribution of Ntracks (subevent C);subsample;Ntracks", nBinsBS_, 0, nBinsBS_));

        h2mhCent.push_back (new TH2F (Form ("h2mhCent%i", n), "Ntracks vs centrality distribution;Ntracks;centrality;Nevents", 100, 0, mhMax_, (centMax_ - centMin_) * 100, centMin_, centMax_));
        h2mhaCent.push_back (new TH2F (Form ("h2mhaCent%i", n), "Ntracks vs centrality distribution (subevent A);Ntracks;centrality;Nevents", 500, 0, mhMax_, (centMax_ - centMin_) * 100, centMin_, centMax_));
        h2mhbCent.push_back (new TH2F (Form ("h2mhbCent%i", n), "Ntracks vs centrality distribution (subevent B);Ntracks;centrality;Nevents", 500, 0, mhMax_, (centMax_ - centMin_) * 100, centMin_, centMax_));
        h2mhcCent.push_back (new TH2F (Form ("h2mhcCent%i", n), "Ntracks vs centrality distribution (subevent C);Ntracks;centrality;Nevents", 500, 0, mhMax_, (centMax_ - centMin_) * 100, centMin_, centMax_));

        h2mhMult.push_back (new TH2F (Form ("h2mhMult%i", n), "Ntracks vs multiplicity distribution;Ntracks;multiplicity;Nevents", 100, 0, mhMax_, 100, mhMin_, mhMax_));
        h2mhaMult.push_back (new TH2F (Form ("h2mhaMult%i", n), "Ntracks vs multiplicity distribution (subevent A);Ntracks;multiplicity;Nevents", 500, 0, mhMax_, 100, mhMin_, mhMax_));
        h2mhbMult.push_back (new TH2F (Form ("h2mhbMult%i", n), "Ntracks vs multiplicity distribution (subevent B);Ntracks;multiplicity;Nevents", 500, 0, mhMax_, 100, mhMin_, mhMax_));
        h2mhcMult.push_back (new TH2F (Form ("h2mhcMult%i", n), "Ntracks vs multiplicity distribution (subevent C);Ntracks;multiplicity;Nevents", 500, 0, mhMax_, 100, mhMin_, mhMax_));

        h2PtEta.push_back (new TH2F (Form ("h2PtEta_%i", n), "P_{T} and " + varName_ + " distribution;P_{T}; " + varName_ + "; nTracks", 100, ptMin_, ptMax_, 100, etaMin_, etaMax_));
        h2PtEtaA.push_back (new TH2F (Form ("h2PtEtaA_%i", n), "P_{T} and " + varName_ + " distribution (subevent A);P_{T}; " + varName_ + "; nTracks", 100, ptMin_, ptMax_, 100, etaMin_, etaMax_));
        h2PtEtaB.push_back (new TH2F (Form ("h2PtEtaB_%i", n), "P_{T} and " + varName_ + " distribution (subevent B);P_{T}; " + varName_ + "; nTracks", 100, ptMin_, ptMax_, 100, etaMin_, etaMax_));
        h2PtEtaC.push_back (new TH2F (Form ("h2PtEtaC_%i", n), "P_{T} and " + varName_ + " distribution (subevent C);P_{T}; " + varName_ + "; nTracks", 100, ptMin_, ptMax_, 100, etaMin_, etaMax_));

        pPtEta.push_back (new TProfile (Form ("pPtEta_%i", n), "#LTP_{T}#GT versus " + varName_ + ";" + varName_ + ";#LTP_{T}#GT", 50, etaMin_, etaMax_));
        pPtEtaA.push_back (new TProfile (Form ("pPtEtaA_%i", n), "#LTP_{T}#GT versus " + varName_ + " (subevent A);" + varName_ + ";#LTP_{T}#GT", 50, etaMin_, etaMax_));
        pPtEtaB.push_back (new TProfile (Form ("pPtEtaB_%i", n), "#LTP_{T}#GT versus " + varName_ + " (subevent B);" + varName_ + ";#LTP_{T}#GT", 50, etaMin_, etaMax_));
        pPtEtaC.push_back (new TProfile (Form ("pPtEtaC_%i", n), "#LTP_{T}#GT versus " + varName_ + " (subevent C);" + varName_ + ";#LTP_{T}#GT", 50, etaMin_, etaMax_));
        if (uniformSet) {
            hXaXRPCent_SP.push_back (new TH2F (Form ("hX%iaX%iRPCent_SP", n, n), Form ("X_{%i}^{a}X_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hXbXRPCent_SP.push_back (new TH2F (Form ("hX%ibX%iRPCent_SP", n, n), Form ("X_{%i}^{b}X_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hXcXRPCent_SP.push_back (new TH2F (Form ("hX%icX%iRPCent_SP", n, n), Form ("X_{%i}^{c}X_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYaYRPCent_SP.push_back (new TH2F (Form ("hY%iaY%iRPCent_SP", n, n), Form ("Y_{%i}^{a}Y_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYbYRPCent_SP.push_back (new TH2F (Form ("hY%ibY%iRPCent_SP", n, n), Form ("Y_{%i}^{b}Y_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYcYRPCent_SP.push_back (new TH2F (Form ("hY%icY%iRPCent_SP", n, n), Form ("Y_{%i}^{c}Y_{%i}^{RP} (SP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hXaXRPCent_EP.push_back (new TH2F (Form ("hX%iaX%iRPCent_EP", n, n), Form ("X_{%i}^{a}X_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hXbXRPCent_EP.push_back (new TH2F (Form ("hX%ibX%iRPCent_EP", n, n), Form ("X_{%i}^{b}X_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hXcXRPCent_EP.push_back (new TH2F (Form ("hX%icX%iRPCent_EP", n, n), Form ("X_{%i}^{c}X_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYaYRPCent_EP.push_back (new TH2F (Form ("hY%iaY%iRPCent_EP", n, n), Form ("Y_{%i}^{a}Y_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYbYRPCent_EP.push_back (new TH2F (Form ("hY%ibY%iRPCent_EP", n, n), Form ("Y_{%i}^{b}Y_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));
            hYcYRPCent_EP.push_back (new TH2F (Form ("hY%icY%iRPCent_EP", n, n), Form ("Y_{%i}^{c}Y_{%i}^{RP} (EP);cent", n, n), nBinsCent_, centMin_, centMax_, 100, -1.0, 1.0));

            hXaXRPMult_SP.push_back (new TH2F (Form ("hX%iaX%iRPMult_SP", n, n), Form ("X_{%i}^{a}X_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hXbXRPMult_SP.push_back (new TH2F (Form ("hX%ibX%iRPMult_SP", n, n), Form ("X_{%i}^{b}X_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hXcXRPMult_SP.push_back (new TH2F (Form ("hX%icX%iRPMult_SP", n, n), Form ("X_{%i}^{c}X_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYaYRPMult_SP.push_back (new TH2F (Form ("hY%iaY%iRPMult_SP", n, n), Form ("Y_{%i}^{a}Y_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYbYRPMult_SP.push_back (new TH2F (Form ("hY%ibY%iRPMult_SP", n, n), Form ("Y_{%i}^{b}Y_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYcYRPMult_SP.push_back (new TH2F (Form ("hY%icY%iRPMult_SP", n, n), Form ("Y_{%i}^{c}Y_{%i}^{RP} (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hXaXRPMult_EP.push_back (new TH2F (Form ("hX%iaX%iRPMult_EP", n, n), Form ("X_{%i}^{a}X_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hXbXRPMult_EP.push_back (new TH2F (Form ("hX%ibX%iRPMult_EP", n, n), Form ("X_{%i}^{b}X_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hXcXRPMult_EP.push_back (new TH2F (Form ("hX%icX%iRPMult_EP", n, n), Form ("X_{%i}^{c}X_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYaYRPMult_EP.push_back (new TH2F (Form ("hY%iaY%iRPMult_EP", n, n), Form ("Y_{%i}^{a}Y_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYbYRPMult_EP.push_back (new TH2F (Form ("hY%ibY%iRPMult_EP", n, n), Form ("Y_{%i}^{b}Y_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
            hYcYRPMult_EP.push_back (new TH2F (Form ("hY%icY%iRPMult_EP", n, n), Form ("Y_{%i}^{c}Y_{%i}^{RP} (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_, 100, -1.0, 1.0));
        }

		pxXaCent_SP.push_back (new TProfile (Form ("px%iX%iaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxXbCent_SP.push_back (new TProfile (Form ("px%iX%ibCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxXcCent_SP.push_back (new TProfile (Form ("px%iX%icCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYaCent_SP.push_back (new TProfile (Form ("py%iY%iaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYbCent_SP.push_back (new TProfile (Form ("py%iY%ibCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYcCent_SP.push_back (new TProfile (Form ("py%iY%icCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXaCent_SP.push_back (new TProfile (Form ("py%iX%iaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXbCent_SP.push_back (new TProfile (Form ("py%iX%ibCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXcCent_SP.push_back (new TProfile (Form ("py%iX%icCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYaCent_SP.push_back (new TProfile (Form ("px%iY%iaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYbCent_SP.push_back (new TProfile (Form ("px%iY%ibCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYcCent_SP.push_back (new TProfile (Form ("px%iY%icCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));

		pxXaCent_EP.push_back (new TProfile (Form ("px%iX%iaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxXbCent_EP.push_back (new TProfile (Form ("px%iX%ibCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxXcCent_EP.push_back (new TProfile (Form ("px%iX%icCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYaCent_EP.push_back (new TProfile (Form ("py%iY%iaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYbCent_EP.push_back (new TProfile (Form ("py%iY%ibCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyYcCent_EP.push_back (new TProfile (Form ("py%iY%icCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXaCent_EP.push_back (new TProfile (Form ("py%iX%iaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXbCent_EP.push_back (new TProfile (Form ("py%iX%ibCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pyXcCent_EP.push_back (new TProfile (Form ("py%iX%icCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYaCent_EP.push_back (new TProfile (Form ("px%iY%iaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYbCent_EP.push_back (new TProfile (Form ("px%iY%ibCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
		pxYcCent_EP.push_back (new TProfile (Form ("px%iY%icCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));

		p2xXaCent_SP.push_back (new TProfile2D (Form ("p2x%iX%iaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xXbCent_SP.push_back (new TProfile2D (Form ("p2x%iX%ibCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xXcCent_SP.push_back (new TProfile2D (Form ("p2x%iX%icCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYaCent_SP.push_back (new TProfile2D (Form ("p2y%iY%iaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYbCent_SP.push_back (new TProfile2D (Form ("p2y%iY%ibCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYcCent_SP.push_back (new TProfile2D (Form ("p2y%iY%icCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXaCent_SP.push_back (new TProfile2D (Form ("p2y%iX%iaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXbCent_SP.push_back (new TProfile2D (Form ("p2y%iX%ibCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXcCent_SP.push_back (new TProfile2D (Form ("p2y%iX%icCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYaCent_SP.push_back (new TProfile2D (Form ("p2x%iY%iaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYbCent_SP.push_back (new TProfile2D (Form ("p2x%iY%ibCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYcCent_SP.push_back (new TProfile2D (Form ("p2x%iY%icCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

		p2xXaCent_EP.push_back (new TProfile2D (Form ("p2x%iX%iaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xXbCent_EP.push_back (new TProfile2D (Form ("p2x%iX%ibCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xXcCent_EP.push_back (new TProfile2D (Form ("p2x%iX%icCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYaCent_EP.push_back (new TProfile2D (Form ("p2y%iY%iaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYbCent_EP.push_back (new TProfile2D (Form ("p2y%iY%ibCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yYcCent_EP.push_back (new TProfile2D (Form ("p2y%iY%icCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXaCent_EP.push_back (new TProfile2D (Form ("p2y%iX%iaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXbCent_EP.push_back (new TProfile2D (Form ("p2y%iX%ibCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2yXcCent_EP.push_back (new TProfile2D (Form ("p2y%iX%icCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYaCent_EP.push_back (new TProfile2D (Form ("p2x%iY%iaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYbCent_EP.push_back (new TProfile2D (Form ("p2x%iY%ibCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2xYcCent_EP.push_back (new TProfile2D (Form ("p2x%iY%icCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

		pxXaMult_SP.push_back (new TProfile (Form ("px%iX%iaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxXbMult_SP.push_back (new TProfile (Form ("px%iX%ibMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxXcMult_SP.push_back (new TProfile (Form ("px%iX%icMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYaMult_SP.push_back (new TProfile (Form ("py%iY%iaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYbMult_SP.push_back (new TProfile (Form ("py%iY%ibMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYcMult_SP.push_back (new TProfile (Form ("py%iY%icMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXaMult_SP.push_back (new TProfile (Form ("py%iX%iaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXbMult_SP.push_back (new TProfile (Form ("py%iX%ibMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXcMult_SP.push_back (new TProfile (Form ("py%iX%icMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYaMult_SP.push_back (new TProfile (Form ("px%iY%iaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYbMult_SP.push_back (new TProfile (Form ("px%iY%ibMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYcMult_SP.push_back (new TProfile (Form ("px%iY%icMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));

		pxXaMult_EP.push_back (new TProfile (Form ("px%iX%iaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxXbMult_EP.push_back (new TProfile (Form ("px%iX%ibMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxXcMult_EP.push_back (new TProfile (Form ("px%iX%icMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYaMult_EP.push_back (new TProfile (Form ("py%iY%iaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYbMult_EP.push_back (new TProfile (Form ("py%iY%ibMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyYcMult_EP.push_back (new TProfile (Form ("py%iY%icMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXaMult_EP.push_back (new TProfile (Form ("py%iX%iaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXbMult_EP.push_back (new TProfile (Form ("py%iX%ibMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pyXcMult_EP.push_back (new TProfile (Form ("py%iX%icMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYaMult_EP.push_back (new TProfile (Form ("px%iY%iaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYbMult_EP.push_back (new TProfile (Form ("px%iY%ibMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));
		pxYcMult_EP.push_back (new TProfile (Form ("px%iY%icMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_));

		p2xXaMult_SP.push_back (new TProfile2D (Form ("p2x%iX%iaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xXbMult_SP.push_back (new TProfile2D (Form ("p2x%iX%ibMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xXcMult_SP.push_back (new TProfile2D (Form ("p2x%iX%icMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYaMult_SP.push_back (new TProfile2D (Form ("p2y%iY%iaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYbMult_SP.push_back (new TProfile2D (Form ("p2y%iY%ibMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYcMult_SP.push_back (new TProfile2D (Form ("p2y%iY%icMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXaMult_SP.push_back (new TProfile2D (Form ("p2y%iX%iaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXbMult_SP.push_back (new TProfile2D (Form ("p2y%iX%ibMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXcMult_SP.push_back (new TProfile2D (Form ("p2y%iX%icMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYaMult_SP.push_back (new TProfile2D (Form ("p2x%iY%iaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYbMult_SP.push_back (new TProfile2D (Form ("p2x%iY%ibMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYcMult_SP.push_back (new TProfile2D (Form ("p2x%iY%icMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

		p2xXaMult_EP.push_back (new TProfile2D (Form ("p2x%iX%iaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xXbMult_EP.push_back (new TProfile2D (Form ("p2x%iX%ibMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xXcMult_EP.push_back (new TProfile2D (Form ("p2x%iX%icMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYaMult_EP.push_back (new TProfile2D (Form ("p2y%iY%iaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYbMult_EP.push_back (new TProfile2D (Form ("p2y%iY%ibMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yYcMult_EP.push_back (new TProfile2D (Form ("p2y%iY%icMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXaMult_EP.push_back (new TProfile2D (Form ("p2y%iX%iaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXbMult_EP.push_back (new TProfile2D (Form ("p2y%iX%ibMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2yXcMult_EP.push_back (new TProfile2D (Form ("p2y%iX%icMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYaMult_EP.push_back (new TProfile2D (Form ("p2x%iY%iaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYbMult_EP.push_back (new TProfile2D (Form ("p2x%iY%ibMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2xYcMult_EP.push_back (new TProfile2D (Form ("p2x%iY%icMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        pXaXbCent_SP.push_back (new TProfile (Form ("pX%iaX%ibCent_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaXcCent_SP.push_back (new TProfile (Form ("pX%iaX%icCent_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXbXcCent_SP.push_back (new TProfile (Form ("pX%ibX%icCent_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaYbCent_SP.push_back (new TProfile (Form ("pX%iaY%ibCent_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaYcCent_SP.push_back (new TProfile (Form ("pX%iaY%icCent_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXbYcCent_SP.push_back (new TProfile (Form ("pX%ibY%icCent_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaYbCent_SP.push_back (new TProfile (Form ("pY%iaY%ibCent_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaYcCent_SP.push_back (new TProfile (Form ("pY%iaY%icCent_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYbYcCent_SP.push_back (new TProfile (Form ("pY%ibY%icCent_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaXbCent_SP.push_back (new TProfile (Form ("pY%iaX%ibCent_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaXcCent_SP.push_back (new TProfile (Form ("pY%iaX%icCent_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYbXcCent_SP.push_back (new TProfile (Form ("pY%ibX%icCent_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));

        pXaXbCent_EP.push_back (new TProfile (Form ("pX%iaX%ibCent_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaXcCent_EP.push_back (new TProfile (Form ("pX%iaX%icCent_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXbXcCent_EP.push_back (new TProfile (Form ("pX%ibX%icCent_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaYbCent_EP.push_back (new TProfile (Form ("pX%iaY%ibCent_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXaYcCent_EP.push_back (new TProfile (Form ("pX%iaY%icCent_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pXbYcCent_EP.push_back (new TProfile (Form ("pX%ibY%icCent_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaYbCent_EP.push_back (new TProfile (Form ("pY%iaY%ibCent_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaYcCent_EP.push_back (new TProfile (Form ("pY%iaY%icCent_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYbYcCent_EP.push_back (new TProfile (Form ("pY%ibY%icCent_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaXbCent_EP.push_back (new TProfile (Form ("pY%iaX%ibCent_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYaXcCent_EP.push_back (new TProfile (Form ("pY%iaX%icCent_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
		pYbXcCent_EP.push_back (new TProfile (Form ("pY%ibX%icCent_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));

        pXaXbMult_SP.push_back (new TProfile (Form ("pX%iaX%ibMult_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaXcMult_SP.push_back (new TProfile (Form ("pX%iaX%icMult_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXbXcMult_SP.push_back (new TProfile (Form ("pX%ibX%icMult_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaYbMult_SP.push_back (new TProfile (Form ("pX%iaY%ibMult_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaYcMult_SP.push_back (new TProfile (Form ("pX%iaY%icMult_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXbYcMult_SP.push_back (new TProfile (Form ("pX%ibY%icMult_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaYbMult_SP.push_back (new TProfile (Form ("pY%iaY%ibMult_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaYcMult_SP.push_back (new TProfile (Form ("pY%iaY%icMult_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYbYcMult_SP.push_back (new TProfile (Form ("pY%ibY%icMult_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaXbMult_SP.push_back (new TProfile (Form ("pY%iaX%ibMult_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaXcMult_SP.push_back (new TProfile (Form ("pY%iaX%icMult_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYbXcMult_SP.push_back (new TProfile (Form ("pY%ibX%icMult_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));

        pXaXbMult_EP.push_back (new TProfile (Form ("pX%iaX%ibMult_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaXcMult_EP.push_back (new TProfile (Form ("pX%iaX%icMult_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXbXcMult_EP.push_back (new TProfile (Form ("pX%ibX%icMult_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaYbMult_EP.push_back (new TProfile (Form ("pX%iaY%ibMult_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXaYcMult_EP.push_back (new TProfile (Form ("pX%iaY%icMult_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pXbYcMult_EP.push_back (new TProfile (Form ("pX%ibY%icMult_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaYbMult_EP.push_back (new TProfile (Form ("pY%iaY%ibMult_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaYcMult_EP.push_back (new TProfile (Form ("pY%iaY%icMult_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYbYcMult_EP.push_back (new TProfile (Form ("pY%ibY%icMult_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaXbMult_EP.push_back (new TProfile (Form ("pY%iaX%ibMult_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYaXcMult_EP.push_back (new TProfile (Form ("pY%iaX%icMult_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
		pYbXcMult_EP.push_back (new TProfile (Form ("pY%ibX%icMult_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));

        p2XaXbCent_SP.push_back (new TProfile2D (Form ("p2X%iaX%ibCent_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaXcCent_SP.push_back (new TProfile2D (Form ("p2X%iaX%icCent_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XbXcCent_SP.push_back (new TProfile2D (Form ("p2X%ibX%icCent_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYbCent_SP.push_back (new TProfile2D (Form ("p2X%iaY%ibCent_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYcCent_SP.push_back (new TProfile2D (Form ("p2X%iaY%icCent_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XbYcCent_SP.push_back (new TProfile2D (Form ("p2X%ibY%icCent_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYbCent_SP.push_back (new TProfile2D (Form ("p2Y%iaY%ibCent_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYcCent_SP.push_back (new TProfile2D (Form ("p2Y%iaY%icCent_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YbYcCent_SP.push_back (new TProfile2D (Form ("p2Y%ibY%icCent_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXbCent_SP.push_back (new TProfile2D (Form ("p2Y%iaX%ibCent_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXcCent_SP.push_back (new TProfile2D (Form ("p2Y%iaX%icCent_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YbXcCent_SP.push_back (new TProfile2D (Form ("p2Y%ibX%icCent_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p2XaXbCent_EP.push_back (new TProfile2D (Form ("p2X%iaX%ibCent_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaXcCent_EP.push_back (new TProfile2D (Form ("p2X%iaX%icCent_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XbXcCent_EP.push_back (new TProfile2D (Form ("p2X%ibX%icCent_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYbCent_EP.push_back (new TProfile2D (Form ("p2X%iaY%ibCent_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYcCent_EP.push_back (new TProfile2D (Form ("p2X%iaY%icCent_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2XbYcCent_EP.push_back (new TProfile2D (Form ("p2X%ibY%icCent_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYbCent_EP.push_back (new TProfile2D (Form ("p2Y%iaY%ibCent_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYcCent_EP.push_back (new TProfile2D (Form ("p2Y%iaY%icCent_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YbYcCent_EP.push_back (new TProfile2D (Form ("p2Y%ibY%icCent_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXbCent_EP.push_back (new TProfile2D (Form ("p2Y%iaX%ibCent_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXcCent_EP.push_back (new TProfile2D (Form ("p2Y%iaX%icCent_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p2YbXcCent_EP.push_back (new TProfile2D (Form ("p2Y%ibX%icCent_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p2XaXbMult_SP.push_back (new TProfile2D (Form ("p2X%iaX%ibMult_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaXcMult_SP.push_back (new TProfile2D (Form ("p2X%iaX%icMult_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XbXcMult_SP.push_back (new TProfile2D (Form ("p2X%ibX%icMult_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYbMult_SP.push_back (new TProfile2D (Form ("p2X%iaY%ibMult_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYcMult_SP.push_back (new TProfile2D (Form ("p2X%iaY%icMult_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XbYcMult_SP.push_back (new TProfile2D (Form ("p2X%ibY%icMult_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYbMult_SP.push_back (new TProfile2D (Form ("p2Y%iaY%ibMult_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYcMult_SP.push_back (new TProfile2D (Form ("p2Y%iaY%icMult_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YbYcMult_SP.push_back (new TProfile2D (Form ("p2Y%ibY%icMult_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXbMult_SP.push_back (new TProfile2D (Form ("p2Y%iaX%ibMult_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXcMult_SP.push_back (new TProfile2D (Form ("p2Y%iaX%icMult_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YbXcMult_SP.push_back (new TProfile2D (Form ("p2Y%ibX%icMult_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        p2XaXbMult_EP.push_back (new TProfile2D (Form ("p2X%iaX%ibMult_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaXcMult_EP.push_back (new TProfile2D (Form ("p2X%iaX%icMult_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XbXcMult_EP.push_back (new TProfile2D (Form ("p2X%ibX%icMult_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYbMult_EP.push_back (new TProfile2D (Form ("p2X%iaY%ibMult_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XaYcMult_EP.push_back (new TProfile2D (Form ("p2X%iaY%icMult_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2XbYcMult_EP.push_back (new TProfile2D (Form ("p2X%ibY%icMult_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYbMult_EP.push_back (new TProfile2D (Form ("p2Y%iaY%ibMult_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaYcMult_EP.push_back (new TProfile2D (Form ("p2Y%iaY%icMult_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YbYcMult_EP.push_back (new TProfile2D (Form ("p2Y%ibY%icMult_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXbMult_EP.push_back (new TProfile2D (Form ("p2Y%iaX%ibMult_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YaXcMult_EP.push_back (new TProfile2D (Form ("p2Y%iaX%icMult_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p2YbXcMult_EP.push_back (new TProfile2D (Form ("p2Y%ibX%icMult_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        p2xXaPtCent_SP.push_back (new TProfile2D (Form ("p2x%iX%iaPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXaPtCent_SP.push_back (new TProfile2D (Form ("p2y%iX%iaPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xXbPtCent_SP.push_back (new TProfile2D (Form ("p2x%iX%ibPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXbPtCent_SP.push_back (new TProfile2D (Form ("p2y%iX%ibPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xXcPtCent_SP.push_back (new TProfile2D (Form ("p2x%iX%icPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXcPtCent_SP.push_back (new TProfile2D (Form ("p2y%iX%icPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYaPtCent_SP.push_back (new TProfile2D (Form ("p2y%iY%iaPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYaPtCent_SP.push_back (new TProfile2D (Form ("p2x%iY%iaPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYbPtCent_SP.push_back (new TProfile2D (Form ("p2y%iY%ibPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYbPtCent_SP.push_back (new TProfile2D (Form ("p2x%iY%ibPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYcPtCent_SP.push_back (new TProfile2D (Form ("p2y%iY%icPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYcPtCent_SP.push_back (new TProfile2D (Form ("p2x%iY%icPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));

        p2xXaPtCent_EP.push_back (new TProfile2D (Form ("p2x%iX%iaPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXaPtCent_EP.push_back (new TProfile2D (Form ("p2y%iX%iaPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xXbPtCent_EP.push_back (new TProfile2D (Form ("p2x%iX%ibPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXbPtCent_EP.push_back (new TProfile2D (Form ("p2y%iX%ibPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xXcPtCent_EP.push_back (new TProfile2D (Form ("p2x%iX%icPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yXcPtCent_EP.push_back (new TProfile2D (Form ("p2y%iX%icPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYaPtCent_EP.push_back (new TProfile2D (Form ("p2y%iY%iaPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYaPtCent_EP.push_back (new TProfile2D (Form ("p2x%iY%iaPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYbPtCent_EP.push_back (new TProfile2D (Form ("p2y%iY%ibPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYbPtCent_EP.push_back (new TProfile2D (Form ("p2x%iY%ibPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2yYcPtCent_EP.push_back (new TProfile2D (Form ("p2y%iY%icPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
		p2xYcPtCent_EP.push_back (new TProfile2D (Form ("p2x%iY%icPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));

        p3xXaPtCent_SP.push_back (new TProfile3D (Form ("p3x%iX%iaPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaPtCent_SP.push_back (new TProfile3D (Form ("p3y%iX%iaPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbPtCent_SP.push_back (new TProfile3D (Form ("p3x%iX%ibPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbPtCent_SP.push_back (new TProfile3D (Form ("p3y%iX%ibPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcPtCent_SP.push_back (new TProfile3D (Form ("p3x%iX%icPtCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcPtCent_SP.push_back (new TProfile3D (Form ("p3y%iX%icPtCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaPtCent_SP.push_back (new TProfile3D (Form ("p3y%iY%iaPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaPtCent_SP.push_back (new TProfile3D (Form ("p3x%iY%iaPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbPtCent_SP.push_back (new TProfile3D (Form ("p3y%iY%ibPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbPtCent_SP.push_back (new TProfile3D (Form ("p3x%iY%ibPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcPtCent_SP.push_back (new TProfile3D (Form ("p3y%iY%icPtCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcPtCent_SP.push_back (new TProfile3D (Form ("p3x%iY%icPtCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p3xXaPtCent_EP.push_back (new TProfile3D (Form ("p3x%iX%iaPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaPtCent_EP.push_back (new TProfile3D (Form ("p3y%iX%iaPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbPtCent_EP.push_back (new TProfile3D (Form ("p3x%iX%ibPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbPtCent_EP.push_back (new TProfile3D (Form ("p3y%iX%ibPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcPtCent_EP.push_back (new TProfile3D (Form ("p3x%iX%icPtCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcPtCent_EP.push_back (new TProfile3D (Form ("p3y%iX%icPtCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaPtCent_EP.push_back (new TProfile3D (Form ("p3y%iY%iaPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaPtCent_EP.push_back (new TProfile3D (Form ("p3x%iY%iaPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbPtCent_EP.push_back (new TProfile3D (Form ("p3y%iY%ibPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbPtCent_EP.push_back (new TProfile3D (Form ("p3x%iY%ibPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcPtCent_EP.push_back (new TProfile3D (Form ("p3y%iY%icPtCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcPtCent_EP.push_back (new TProfile3D (Form ("p3x%iY%icPtCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p2xXaPtMult_SP.push_back (new TProfile2D (Form ("p2x%iX%iaPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXaPtMult_SP.push_back (new TProfile2D (Form ("p2y%iX%iaPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXbPtMult_SP.push_back (new TProfile2D (Form ("p2x%iX%ibPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXbPtMult_SP.push_back (new TProfile2D (Form ("p2y%iX%ibPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXcPtMult_SP.push_back (new TProfile2D (Form ("p2x%iX%icPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXcPtMult_SP.push_back (new TProfile2D (Form ("p2y%iX%icPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYaPtMult_SP.push_back (new TProfile2D (Form ("p2y%iY%iaPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYaPtMult_SP.push_back (new TProfile2D (Form ("p2x%iY%iaPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYbPtMult_SP.push_back (new TProfile2D (Form ("p2y%iY%ibPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYbPtMult_SP.push_back (new TProfile2D (Form ("p2x%iY%ibPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYcPtMult_SP.push_back (new TProfile2D (Form ("p2y%iY%icPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYcPtMult_SP.push_back (new TProfile2D (Form ("p2x%iY%icPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));

        p2xXaPtMult_EP.push_back (new TProfile2D (Form ("p2x%iX%iaPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXaPtMult_EP.push_back (new TProfile2D (Form ("p2y%iX%iaPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXbPtMult_EP.push_back (new TProfile2D (Form ("p2x%iX%ibPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXbPtMult_EP.push_back (new TProfile2D (Form ("p2y%iX%ibPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXcPtMult_EP.push_back (new TProfile2D (Form ("p2x%iX%icPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXcPtMult_EP.push_back (new TProfile2D (Form ("p2y%iX%icPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYaPtMult_EP.push_back (new TProfile2D (Form ("p2y%iY%iaPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYaPtMult_EP.push_back (new TProfile2D (Form ("p2x%iY%iaPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYbPtMult_EP.push_back (new TProfile2D (Form ("p2y%iY%ibPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYbPtMult_EP.push_back (new TProfile2D (Form ("p2x%iY%ibPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYcPtMult_EP.push_back (new TProfile2D (Form ("p2y%iY%icPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYcPtMult_EP.push_back (new TProfile2D (Form ("p2x%iY%icPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));

        p3xXaPtMult_SP.push_back (new TProfile3D (Form ("p3x%iX%iaPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaPtMult_SP.push_back (new TProfile3D (Form ("p3y%iX%iaPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbPtMult_SP.push_back (new TProfile3D (Form ("p3x%iX%ibPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbPtMult_SP.push_back (new TProfile3D (Form ("p3y%iX%ibPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcPtMult_SP.push_back (new TProfile3D (Form ("p3x%iX%icPtMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcPtMult_SP.push_back (new TProfile3D (Form ("p3y%iX%icPtMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaPtMult_SP.push_back (new TProfile3D (Form ("p3y%iY%iaPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaPtMult_SP.push_back (new TProfile3D (Form ("p3x%iY%iaPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbPtMult_SP.push_back (new TProfile3D (Form ("p3y%iY%ibPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbPtMult_SP.push_back (new TProfile3D (Form ("p3x%iY%ibPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcPtMult_SP.push_back (new TProfile3D (Form ("p3y%iY%icPtMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcPtMult_SP.push_back (new TProfile3D (Form ("p3x%iY%icPtMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        p3xXaPtMult_EP.push_back (new TProfile3D (Form ("p3x%iX%iaPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaPtMult_EP.push_back (new TProfile3D (Form ("p3y%iX%iaPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbPtMult_EP.push_back (new TProfile3D (Form ("p3x%iX%ibPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbPtMult_EP.push_back (new TProfile3D (Form ("p3y%iX%ibPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcPtMult_EP.push_back (new TProfile3D (Form ("p3x%iX%icPtMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcPtMult_EP.push_back (new TProfile3D (Form ("p3y%iX%icPtMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaPtMult_EP.push_back (new TProfile3D (Form ("p3y%iY%iaPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaPtMult_EP.push_back (new TProfile3D (Form ("p3x%iY%iaPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbPtMult_EP.push_back (new TProfile3D (Form ("p3y%iY%ibPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbPtMult_EP.push_back (new TProfile3D (Form ("p3x%iY%ibPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcPtMult_EP.push_back (new TProfile3D (Form ("p3y%iY%icPtMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcPtMult_EP.push_back (new TProfile3D (Form ("p3x%iY%icPtMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        p2xXaEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iX%iaEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXaEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iX%iaEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xXbEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iX%ibEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXbEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iX%ibEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xXcEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iX%icEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXcEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iX%icEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYaEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iY%iaEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYaEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iY%iaEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYbEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iY%ibEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYbEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iY%ibEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYcEtaCent_SP.push_back (new TProfile2D (Form ("p2y%iY%icEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYcEtaCent_SP.push_back (new TProfile2D (Form ("p2x%iY%icEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));

        p2xXaEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iX%iaEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXaEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iX%iaEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xXbEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iX%ibEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXbEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iX%ibEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xXcEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iX%icEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yXcEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iX%icEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYaEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iY%iaEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYaEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iY%iaEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYbEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iY%ibEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYbEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iY%ibEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2yYcEtaCent_EP.push_back (new TProfile2D (Form ("p2y%iY%icEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
		p2xYcEtaCent_EP.push_back (new TProfile2D (Form ("p2x%iY%icEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));

        p3xXaEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iX%iaEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iX%iaEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iX%ibEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iX%ibEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iX%icEtaCent_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iX%icEtaCent_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iY%iaEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iY%iaEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iY%ibEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iY%ibEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcEtaCent_SP.push_back (new TProfile3D (Form ("p3y%iY%icEtaCent_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcEtaCent_SP.push_back (new TProfile3D (Form ("p3x%iY%icEtaCent_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p3xXaEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iX%iaEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iX%iaEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iX%ibEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iX%ibEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iX%icEtaCent_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iX%icEtaCent_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iY%iaEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iY%iaEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iY%ibEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iY%ibEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcEtaCent_EP.push_back (new TProfile3D (Form ("p3y%iY%icEtaCent_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcEtaCent_EP.push_back (new TProfile3D (Form ("p3x%iY%icEtaCent_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

        p2xXaEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iX%iaEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXaEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iX%iaEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXbEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iX%ibEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXbEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iX%ibEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXcEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iX%icEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXcEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iX%icEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYaEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iY%iaEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYaEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iY%iaEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYbEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iY%ibEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYbEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iY%ibEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYcEtaMult_SP.push_back (new TProfile2D (Form ("p2y%iY%icEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYcEtaMult_SP.push_back (new TProfile2D (Form ("p2x%iY%icEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));

        p2xXaEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iX%iaEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXaEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iX%iaEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXbEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iX%ibEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXbEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iX%ibEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xXcEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iX%icEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yXcEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iX%icEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYaEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iY%iaEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYaEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iY%iaEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYbEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iY%ibEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYbEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iY%ibEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2yYcEtaMult_EP.push_back (new TProfile2D (Form ("p2y%iY%icEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
		p2xYcEtaMult_EP.push_back (new TProfile2D (Form ("p2x%iY%icEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));

        p3xXaEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iX%iaEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iX%iaEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iX%ibEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iX%ibEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iX%icEtaMult_SP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iX%icEtaMult_SP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iY%iaEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iY%iaEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iY%ibEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iY%ibEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcEtaMult_SP.push_back (new TProfile3D (Form ("p3y%iY%icEtaMult_SP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcEtaMult_SP.push_back (new TProfile3D (Form ("p3x%iY%icEtaMult_SP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (SP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

        p3xXaEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iX%iaEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXaEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iX%iaEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXbEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iX%ibEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXbEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iX%ibEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xXcEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iX%icEtaMult_EP", n, n), Form ("#LTx_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yXcEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iX%icEtaMult_EP", n, n), Form ("#LTy_{%i}X_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYaEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iY%iaEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYaEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iY%iaEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{a}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYbEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iY%ibEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYbEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iY%ibEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{b}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3yYcEtaMult_EP.push_back (new TProfile3D (Form ("p3y%iY%icEtaMult_EP", n, n), Form ("#LTy_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
		p3xYcEtaMult_EP.push_back (new TProfile3D (Form ("p3x%iY%icEtaMult_EP", n, n), Form ("#LTx_{%i}Y_{%i}^{c}#GT (EP);" + varName_ + " ;mh;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
	}

	int total = 0;

	for (Long64_t jentry = 0; jentry < nEvents; jentry++) { // loop over events
        cout << "\rEvent " << jentry + 1 << " from " << nEvents;
		QnMan -> ClearEvent ();
		inputTree -> GetEntry (jentry);
		if (samplingMethod_ == kBootStrapping) sampleTree -> GetEntry (jentry);
        else {
            bsIndex = gRandom -> Rndm () * nBinsBS_; // subsampling
            W [bsIndex] = 1; // subsampling
        }
		nRun = event -> GetNrun ();
        skipFlag = 0;

		for (UInt_t i = 0; i < excludedRuns_.size (); i++) {
            if (nRun == excludedRuns_ [i]) skipFlag = 1;
		}
		if (skipFlag == 1) continue;
		mh = event -> GetMh ();
		cent = event -> GetCent ();

		Eveto = event -> GetEveto ();
		summEveto = Eveto [0] + Eveto [1] + Eveto [2] + Eveto [3];

		centBin = p3xXaEtaCent_SP [0] -> GetYaxis () -> FindBin (cent);
        //centrality classes
        if (cent < centMin_ || cent > centMax_) {
            continue;
        }
		if (mh < mhMin_ || mh > mhMax_) {
            continue;
        }
		QnMan -> GetDataContainer () [kNrun] = nRun;
		//QnMan -> GetDataContainer () [kMult] = mh;
		QnMan -> GetDataContainer () [kCent] = cent;

		for (Int_t i = 0; i < nHarmonics; i++) { // zero subevent multiplicities
			mha [i] = 0;
			mhb [i] = 0;
			mhc [i] = 0;
		}

		for (Int_t itrack = 1; itrack <= mh; itrack++) { // loop over tracks
			track = event -> GetTrack (itrack);
            pt = track -> GetPt ();
			if (varName_ == "#it{y}") eta = track -> GetRap ();
			else eta = track -> GetEta ();
			phi = track -> GetPhi ();
			charge = track -> GetCharge ();
			pid = track -> GetPid ();

            if (pid == kVeto1) weight = Eveto [0] / summEveto;
            else if (pid == kVeto2) weight = Eveto [1] / summEveto;
            else if (pid == kVeto3) weight = Eveto [2] / summEveto;
            else if (pid == kVeto4) weight = Eveto [3] / summEveto;
            else if (pid == kFW) weight = charge;
            else weight = 1.0;


            subeventIndex = r.Rndm (); // random subevent

			QnMan -> AddDataVector (kDetector1, phi);

			for (Int_t i = 0; i < nHarmonics; i++) {
                n = harmonicsMap [i];
                subeventFlag [i][itrack - 1] = 0;
                h2PtEta [i] -> Fill (pt, eta);
                pPtEta [i] -> Fill (eta, pt);
                pPtCent [i] -> Fill (cent, pt, eta);
                if (resChargeSet && charge * resCharge < 0) continue; // resolution from differently charged particles

                for (UInt_t j = 0; j < resParticles [i][0].size(); j++) {
                    if (resMethod_ == kRandomSubevent && subeventIndex >= 0.5) break;
                    if (pid == resParticles [i][0][j]) {
                        if (eta > etaLim_ [i][0] && eta < etaLim_ [i][1] && pt > ptLim_ [i][0] && pt < ptLim_ [i][1]) {
                            if (n == 1 && pid == kProton) phi += TMath::Pi (); // patch
                            QnMan -> AddDataVector (kNDetectors * i + kDetector1A, phi, weight);
                            subeventFlag [i][itrack - 1] = 1;
                            mha [i] ++;
                            h2PtEtaA [i] -> Fill (pt, eta);
                            pPtEtaA [i] -> Fill (eta, pt);
                            pPtCentA [i] -> Fill (cent, pt, eta);
                            break;
                        }
                    }
				}

				for (UInt_t j = 0; j < resParticles [i][1].size(); j++) {
                    if (resMethod_ == kRandomSubevent && subeventIndex < 0.5) break;
                    if (pid == resParticles [i][1][j]) {
                        if (eta > etaLim_ [i][2] && eta < etaLim_ [i][3] && pt > ptLim_ [i][2] && pt < ptLim_ [i][3]) {
                            if (n == 1 && pid == kProton) phi += TMath::Pi (); // patch
                            QnMan -> AddDataVector (kNDetectors * i + kDetector1B, phi, weight);
                            subeventFlag [i][itrack - 1] = 2;
                            mhb [i] ++;
                            h2PtEtaB [i] -> Fill (pt, eta);
                            pPtEtaB [i] -> Fill (eta, pt);
                            pPtCentB [i] -> Fill (cent, pt, eta);
                            break;
                        }
                    }
				}

				for (UInt_t j = 0; j < resParticles [i][2].size(); j++) {
                    if (resMethod_ == kRandomSubevent) break;
                    if (pid == resParticles [i][2][j]) {
                        if (eta > etaLim_ [i][4] && eta < etaLim_ [i][5] && pt > ptLim_ [i][4] && pt < ptLim_ [i][5]) {
                            if (n == 1 && pid == kProton) phi += TMath::Pi (); // patch
                            QnMan -> AddDataVector (kNDetectors * i + kDetector1C, phi, weight);
                            subeventFlag [i][itrack - 1] = 3;
                            mhc [i] ++;
                            h2PtEtaC [i] -> Fill (pt, eta);
                            pPtEtaC [i] -> Fill (eta, pt);
                            pPtCentC [i] -> Fill (cent, pt, eta);
                            break;
                        }
                    }
				}
			}
		}

		QnMan -> ProcessEvent ();

		for (Int_t i = 0; i < nHarmonics; i++) {
			n = harmonicsMap [i];
//            for (Int_t j = 0; j < nBinsCent_; j++) { // momentum conservation
//                ptAvg [i][j] = pPtCent [i] -> GetBinContent (j + 1);
//                ptAvgA [i][j] = pPtCentA [i] -> GetBinContent (j + 1);
//                ptAvgB [i][j] = pPtCentB [i] -> GetBinContent (j + 1);
//                ptAvgC [i][j] = pPtCentC [i] -> GetBinContent (j + 1);
//            }

			if (QnMan -> GetDetectorQnVector ("D1")) {
                X [i] = QnMan -> GetDetectorQnVector ("D1", "latest") -> Qx (n);
                Y [i] = QnMan -> GetDetectorQnVector ("D1", "latest") -> Qy (n);
			}
			else {
				total++;
				printf ("\nD1_%i FAILED!!!\tmh = %i\ttotal = %i\n", n, mh, total);
				if (!useZeroSubevents_) continue;
				X [i] = 0.0;
				Y [i] = 0.0;
			}

			if (QnMan -> GetDetectorQnVector (Form ("D1A_%i", n))) {
                Xa [i] = QnMan -> GetDetectorQnVector (Form ("D1A_%i", n), "latest") -> Qx (n);
                Ya [i] = QnMan -> GetDetectorQnVector (Form ("D1A_%i", n), "latest") -> Qy (n);
			}
			else {
				total++;
				printf ("\nD1A_%i FAILED!!!\tmha = %i\ttotal = %i\n", n, mha [i], total);
				if (!useZeroSubevents_) continue;
                Xa [i] = 0.0;
                Ya [i] = 0.0;
			}

			if (QnMan -> GetDetectorQnVector (Form ("D1B_%i", n))) {
                Xb [i] = QnMan -> GetDetectorQnVector (Form ("D1B_%i", n), "latest") -> Qx (n);
                Yb [i] = QnMan -> GetDetectorQnVector (Form ("D1B_%i", n), "latest") -> Qy (n);
			}
			else {
				total++;
				printf ("\nD1B_%i FAILED!!!\tmhb = %i\ttotal = %i\n", n, mhb [i], total);
				if (!useZeroSubevents_) continue;
                Xb [i] = 0.0;
                Yb [i] = 0.0;
			}

			if (QnMan -> GetDetectorQnVector (Form ("D1C_%i", n))) {
                Xc [i] = QnMan -> GetDetectorQnVector (Form ("D1C_%i", n), "latest") -> Qx (n);
                Yc [i] = QnMan -> GetDetectorQnVector (Form ("D1C_%i", n), "latest") -> Qy (n);
			}
			else {
				total++;
				if (resMethod_ == kThreeSubevents)
                    printf ("\nD1C_%i FAILED!!!\tmhc = %i\ttotal = %i\n", n, mhc [i], total);
				if (!useZeroSubevents_ && resMethod_ == kThreeSubevents) continue;
                Xc [i] = 0.0;
                Yc [i] = 0.0;
			}

			hMh [i] -> Fill (mh);
			hMha [i] -> Fill (mha [i]);
			hMhb [i] -> Fill (mhb [i]);
			hMhc [i] -> Fill (mhc [i]);

            h2mhCent [i] -> Fill (mh, cent);
            h2mhaCent [i] -> Fill (mha [i], cent);
            h2mhbCent [i] -> Fill (mhb [i], cent);
            h2mhcCent [i] -> Fill (mhc [i], cent);

            h2mhMult [i] -> Fill (mh, mh);
            h2mhaMult [i] -> Fill (mha [i], mh);
            h2mhbMult [i] -> Fill (mhb [i], mh);
            h2mhcMult [i] -> Fill (mhc [i], mh);

			//psiEP [i] = TMath::ATan2 (Y [i], X [i]) / n;
			psiEPa [i] = TMath::ATan2 (Ya [i], Xa [i]) / n;
			psiEPb [i] = TMath::ATan2 (Yb [i], Xb [i]) / n;
			psiEPc [i] = TMath::ATan2 (Yc [i], Xc [i]) / n;

//			Q [i] = TMath::Sqrt (X [i] * X [i] + Y [i] * Y [i]); // candidate
//			Qa [i] = TMath::Sqrt (Xa [i] * Xa [i] + Ya [i] * Ya [i]);
//			Qb [i] = TMath::Sqrt (Xb [i] * Xb [i] + Yb [i] * Yb [i]);
//			Qc [i] = TMath::Sqrt (Xc [i] * Xc [i] + Yc [i] * Yc [i]);

			if (uniformSet) {
				psiRP [i] = event -> GetPsi_n (n);
				XRP [i] = TMath::Cos (n * psiRP [i]);
				YRP [i] = TMath::Sin (n * psiRP [i]);
			}

            //no sampling
            pXaXbCent_SP [i] -> Fill (cent, Xa [i] * Xb [i]);
            pXaXcCent_SP [i] -> Fill (cent, Xa [i] * Xc [i]);
            pXbXcCent_SP [i] -> Fill (cent, Xb [i] * Xc [i]);
            pXaYbCent_SP [i] -> Fill (cent, Xa [i] * Yb [i]);
            pXaYcCent_SP [i] -> Fill (cent, Xa [i] * Yc [i]);
            pXbYcCent_SP [i] -> Fill (cent, Xb [i] * Yc [i]);
            pYaYbCent_SP [i] -> Fill (cent, Ya [i] * Yb [i]);
            pYaYcCent_SP [i] -> Fill (cent, Ya [i] * Yc [i]);
            pYbYcCent_SP [i] -> Fill (cent, Yb [i] * Yc [i]);
            pYaXbCent_SP [i] -> Fill (cent, Ya [i] * Xb [i]);
            pYaXcCent_SP [i] -> Fill (cent, Ya [i] * Xc [i]);
            pYbXcCent_SP [i] -> Fill (cent, Yb [i] * Xc [i]);

            pXaXbCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]));
            pXaXcCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]));
            pXbXcCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]));
            pXaYbCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]));
            pXaYcCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]));
            pXbYcCent_EP [i] -> Fill (cent, TMath::Cos (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]));
            pYaYbCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]));
            pYaYcCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]));
            pYbYcCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]));
            pYaXbCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]));
            pYaXcCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]));
            pYbXcCent_EP [i] -> Fill (cent, TMath::Sin (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]));

            pXaXbMult_SP [i] -> Fill (mh, Xa [i] * Xb [i]);
            pXaXcMult_SP [i] -> Fill (mh, Xa [i] * Xc [i]);
            pXbXcMult_SP [i] -> Fill (mh, Xb [i] * Xc [i]);
            pXaYbMult_SP [i] -> Fill (mh, Xa [i] * Yb [i]);
            pXaYcMult_SP [i] -> Fill (mh, Xa [i] * Yc [i]);
            pXbYcMult_SP [i] -> Fill (mh, Xb [i] * Yc [i]);
            pYaYbMult_SP [i] -> Fill (mh, Ya [i] * Yb [i]);
            pYaYcMult_SP [i] -> Fill (mh, Ya [i] * Yc [i]);
            pYbYcMult_SP [i] -> Fill (mh, Yb [i] * Yc [i]);
            pYaXbMult_SP [i] -> Fill (mh, Ya [i] * Xb [i]);
            pYaXcMult_SP [i] -> Fill (mh, Ya [i] * Xc [i]);
            pYbXcMult_SP [i] -> Fill (mh, Yb [i] * Xc [i]);

            pXaXbMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]));
            pXaXcMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]));
            pXbXcMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]));
            pXaYbMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]));
            pXaYcMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]));
            pXbYcMult_EP [i] -> Fill (mh, TMath::Cos (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]));
            pYaYbMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]));
            pYaYcMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]));
            pYbYcMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]));
            pYaXbMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]));
            pYaXcMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]));
            pYbXcMult_EP [i] -> Fill (mh, TMath::Sin (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]));

            //sampling
            for (Int_t s = 0; s < sMax; s++) { // loop over samples
                if (samplingMethod_ == kBootStrapping) bsIndex = s;

                hNeventsBS [i] -> Fill (bsIndex, W [bsIndex]);
                hNtracksBS [i] -> Fill (bsIndex, mh);
                hNtracksBSa [i] -> Fill (bsIndex, mha [i]);
                hNtracksBSb [i] -> Fill (bsIndex, mhb [i]);
                hNtracksBSc [i] -> Fill (bsIndex, mhc [i]);

                h2nEventSampleWeight [i] -> Fill (jentry, bsIndex, W [bsIndex]);
                p2XaXbCent_SP [i] -> Fill (cent, bsIndex, Xa [i] * Xb [i], W [bsIndex]);
                p2XaXcCent_SP [i] -> Fill (cent, bsIndex, Xa [i] * Xc [i], W [bsIndex]);
                p2XbXcCent_SP [i] -> Fill (cent, bsIndex, Xb [i] * Xc [i], W [bsIndex]);
                p2XaYbCent_SP [i] -> Fill (cent, bsIndex, Xa [i] * Yb [i], W [bsIndex]);
                p2XaYcCent_SP [i] -> Fill (cent, bsIndex, Xa [i] * Yc [i], W [bsIndex]);
                p2XbYcCent_SP [i] -> Fill (cent, bsIndex, Xb [i] * Yc [i], W [bsIndex]);
                p2YaYbCent_SP [i] -> Fill (cent, bsIndex, Ya [i] * Yb [i], W [bsIndex]);
                p2YaYcCent_SP [i] -> Fill (cent, bsIndex, Ya [i] * Yc [i], W [bsIndex]);
                p2YbYcCent_SP [i] -> Fill (cent, bsIndex, Yb [i] * Yc [i], W [bsIndex]);
                p2YaXbCent_SP [i] -> Fill (cent, bsIndex, Ya [i] * Xb [i], W [bsIndex]);
                p2YaXcCent_SP [i] -> Fill (cent, bsIndex, Ya [i] * Xc [i], W [bsIndex]);
                p2YbXcCent_SP [i] -> Fill (cent, bsIndex, Yb [i] * Xc [i], W [bsIndex]);

                p2XaXbCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                p2XaXcCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2XbXcCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2XaYbCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                p2XaYcCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2XbYcCent_EP [i] -> Fill (cent, bsIndex, TMath::Cos (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YaYbCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                p2YaYcCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YbYcCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YaXbCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                p2YaXcCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2YbXcCent_EP [i] -> Fill (cent, bsIndex, TMath::Sin (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);

                p2XaXbMult_SP [i] -> Fill (mh, bsIndex, Xa [i] * Xb [i], W [bsIndex]);
                p2XaXcMult_SP [i] -> Fill (mh, bsIndex, Xa [i] * Xc [i], W [bsIndex]);
                p2XbXcMult_SP [i] -> Fill (mh, bsIndex, Xb [i] * Xc [i], W [bsIndex]);
                p2XaYbMult_SP [i] -> Fill (mh, bsIndex, Xa [i] * Yb [i], W [bsIndex]);
                p2XaYcMult_SP [i] -> Fill (mh, bsIndex, Xa [i] * Yc [i], W [bsIndex]);
                p2XbYcMult_SP [i] -> Fill (mh, bsIndex, Xb [i] * Yc [i], W [bsIndex]);
                p2YaYbMult_SP [i] -> Fill (mh, bsIndex, Ya [i] * Yb [i], W [bsIndex]);
                p2YaYcMult_SP [i] -> Fill (mh, bsIndex, Ya [i] * Yc [i], W [bsIndex]);
                p2YbYcMult_SP [i] -> Fill (mh, bsIndex, Yb [i] * Yc [i], W [bsIndex]);
                p2YaXbMult_SP [i] -> Fill (mh, bsIndex, Ya [i] * Xb [i], W [bsIndex]);
                p2YaXcMult_SP [i] -> Fill (mh, bsIndex, Ya [i] * Xc [i], W [bsIndex]);
                p2YbXcMult_SP [i] -> Fill (mh, bsIndex, Yb [i] * Xc [i], W [bsIndex]);

                p2XaXbMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                p2XaXcMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2XbXcMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2XaYbMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                p2XaYcMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2XbYcMult_EP [i] -> Fill (mh, bsIndex, TMath::Cos (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YaYbMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                p2YaYcMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YbYcMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPb [i]) * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                p2YaXbMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                p2YaXcMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPa [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                p2YbXcMult_EP [i] -> Fill (mh, bsIndex, TMath::Sin (n * psiEPb [i]) * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
            }
		}

		for (Int_t itrack = 1; itrack <= mh; itrack++) {
			skipFlag = 1;
            track = event -> GetTrack (itrack);
            pid = track -> GetPid ();
            for (Int_t i = 0; i < nFlowParts; i++) { // check pid of the tracks used for unQn
                if (pid == flowParticles [i]) skipFlag = 0;
            }
            if (skipFlag == 1) continue;
			pt = track -> GetPt ();
            if (varName_ == "#it{y}") eta = track -> GetRap ();
			else eta = track -> GetEta ();
			phi = track -> GetPhi ();


			for (Int_t i = 0; i < nHarmonics; i++) {
				n = harmonicsMap [i];
				x = TMath::Cos (n * phi);
				y = TMath::Sin (n * phi);
                if (varName_ == "#it{y}" && eta < 0.0 && n % 2 == 1) {
                    //eta = TMath::Abs (eta);
                    sign = -1.0;
                }
                else sign = 1.0;

				if (uniformSet) {
                    // Fill Monte-Carlo HERE
				}

//				if (Xa [i] != 0.0) Xa [i] += ptAvgA [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] += ptAvgB [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] += ptAvgC [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] += ptAvgFixedA [centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] += ptAvgFixedB [centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] += ptAvgFixedC [centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] += ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] += ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] += ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] += ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] += ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] += ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

				if (subeventFlag [i][itrack - 1] == 1) {
					Xa [i] -= x / mha [i];
					Ya [i] -= y / mha [i];
					psiEPa [i] = TMath::ATan2 (Ya [i], Xa [i]) / n;
				}

				if (subeventFlag [i][itrack - 1] == 2) {
					Xb [i] -= x / mhb [i];
					Yb [i] -= y / mhb [i];
					psiEPb [i] = TMath::ATan2 (Yb [i], Xb [i]) / n;
				}

				if (subeventFlag [i][itrack - 1] == 3) {
					Xc [i] -= x / mhc [i];
					Yc [i] -= y / mhc [i];
					psiEPc [i] = TMath::ATan2 (Yc [i], Xc [i]) / n;
				}

                // no sampling
				if (pt >= ptAveragingRange_ [i][0] && pt <= ptAveragingRange_ [i][1]
                    && eta >= etaAveragingRange_ [i][0] && eta <= etaAveragingRange_ [i][1]) {

                    pxXaCent_SP [i] -> Fill (cent, sign * x * Xa [i]);
                    pxXbCent_SP [i] -> Fill (cent, sign * x * Xb [i]);
                    pxXcCent_SP [i] -> Fill (cent, sign * x * Xc [i]);
                    pyXaCent_SP [i] -> Fill (cent, sign * y * Xa [i]);
                    pyXbCent_SP [i] -> Fill (cent, sign * y * Xb [i]);
                    pyXcCent_SP [i] -> Fill (cent, sign * y * Xc [i]);
                    pyYaCent_SP [i] -> Fill (cent, sign * y * Ya [i]);
                    pyYbCent_SP [i] -> Fill (cent, sign * y * Yb [i]);
                    pyYcCent_SP [i] -> Fill (cent, sign * y * Yc [i]);
                    pxYaCent_SP [i] -> Fill (cent, sign * x * Ya [i]);
                    pxYbCent_SP [i] -> Fill (cent, sign * x * Yb [i]);
                    pxYcCent_SP [i] -> Fill (cent, sign * x * Yc [i]);

                    pxXaCent_EP [i] -> Fill (cent, sign * x * TMath::Cos (n * psiEPa [i]));
                    pxXbCent_EP [i] -> Fill (cent, sign * x * TMath::Cos (n * psiEPb [i]));
                    pxXcCent_EP [i] -> Fill (cent, sign * x * TMath::Cos (n * psiEPc [i]));
                    pxYaCent_EP [i] -> Fill (cent, sign * x * TMath::Sin (n * psiEPa [i]));
                    pxYbCent_EP [i] -> Fill (cent, sign * x * TMath::Sin (n * psiEPb [i]));
                    pxYcCent_EP [i] -> Fill (cent, sign * x * TMath::Sin (n * psiEPc [i]));
                    pyYaCent_EP [i] -> Fill (cent, sign * y * TMath::Sin (n * psiEPa [i]));
                    pyYbCent_EP [i] -> Fill (cent, sign * y * TMath::Sin (n * psiEPb [i]));
                    pyYcCent_EP [i] -> Fill (cent, sign * y * TMath::Sin (n * psiEPc [i]));
                    pyXaCent_EP [i] -> Fill (cent, sign * y * TMath::Cos (n * psiEPa [i]));
                    pyXbCent_EP [i] -> Fill (cent, sign * y * TMath::Cos (n * psiEPb [i]));
                    pyXcCent_EP [i] -> Fill (cent, sign * y * TMath::Cos (n * psiEPc [i]));

                    pxXaMult_SP [i] -> Fill (mh, sign * x * Xa [i]);
                    pxXbMult_SP [i] -> Fill (mh, sign * x * Xb [i]);
                    pxXcMult_SP [i] -> Fill (mh, sign * x * Xc [i]);
                    pyXaMult_SP [i] -> Fill (mh, sign * y * Xa [i]);
                    pyXbMult_SP [i] -> Fill (mh, sign * y * Xb [i]);
                    pyXcMult_SP [i] -> Fill (mh, sign * y * Xc [i]);
                    pyYaMult_SP [i] -> Fill (mh, sign * y * Ya [i]);
                    pyYbMult_SP [i] -> Fill (mh, sign * y * Yb [i]);
                    pyYcMult_SP [i] -> Fill (mh, sign * y * Yc [i]);
                    pxYaMult_SP [i] -> Fill (mh, sign * x * Ya [i]);
                    pxYbMult_SP [i] -> Fill (mh, sign * x * Yb [i]);
                    pxYcMult_SP [i] -> Fill (mh, sign * x * Yc [i]);

                    pxXaMult_EP [i] -> Fill (mh, sign * x * TMath::Cos (n * psiEPa [i]));
                    pxXbMult_EP [i] -> Fill (mh, sign * x * TMath::Cos (n * psiEPb [i]));
                    pxXcMult_EP [i] -> Fill (mh, sign * x * TMath::Cos (n * psiEPc [i]));
                    pxYaMult_EP [i] -> Fill (mh, sign * x * TMath::Sin (n * psiEPa [i]));
                    pxYbMult_EP [i] -> Fill (mh, sign * x * TMath::Sin (n * psiEPb [i]));
                    pxYcMult_EP [i] -> Fill (mh, sign * x * TMath::Sin (n * psiEPc [i]));
                    pyYaMult_EP [i] -> Fill (mh, sign * y * TMath::Sin (n * psiEPa [i]));
                    pyYbMult_EP [i] -> Fill (mh, sign * y * TMath::Sin (n * psiEPb [i]));
                    pyYcMult_EP [i] -> Fill (mh, sign * y * TMath::Sin (n * psiEPc [i]));
                    pyXaMult_EP [i] -> Fill (mh, sign * y * TMath::Cos (n * psiEPa [i]));
                    pyXbMult_EP [i] -> Fill (mh, sign * y * TMath::Cos (n * psiEPb [i]));
                    pyXcMult_EP [i] -> Fill (mh, sign * y * TMath::Cos (n * psiEPc [i]));
                }

                if (pt >= ptMin_ && pt <= ptMax_) {
                    if (eta >= etaAveragingRange_ [i][0] && eta <= etaAveragingRange_ [i][1]) {
                        p2xXaPtCent_SP [i] -> Fill (pt, cent, sign * x * Xa [i]);
                        p2xXbPtCent_SP [i] -> Fill (pt, cent, sign * x * Xb [i]);
                        p2xXcPtCent_SP [i] -> Fill (pt, cent, sign * x * Xc [i]);
                        p2yXaPtCent_SP [i] -> Fill (pt, cent, sign * y * Xa [i]);
                        p2yXbPtCent_SP [i] -> Fill (pt, cent, sign * y * Xb [i]);
                        p2yXcPtCent_SP [i] -> Fill (pt, cent, sign * y * Xc [i]);
                        p2yYaPtCent_SP [i] -> Fill (pt, cent, sign * y * Ya [i]);
                        p2yYbPtCent_SP [i] -> Fill (pt, cent, sign * y * Yb [i]);
                        p2yYcPtCent_SP [i] -> Fill (pt, cent, sign * y * Yc [i]);
                        p2xYaPtCent_SP [i] -> Fill (pt, cent, sign * x * Ya [i]);
                        p2xYbPtCent_SP [i] -> Fill (pt, cent, sign * x * Yb [i]);
                        p2xYcPtCent_SP [i] -> Fill (pt, cent, sign * x * Yc [i]);

                        p2xXaPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Cos (n * psiEPa [i]));
                        p2xXbPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Cos (n * psiEPb [i]));
                        p2xXcPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Cos (n * psiEPc [i]));
                        p2yXaPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Cos (n * psiEPa [i]));
                        p2yXbPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Cos (n * psiEPb [i]));
                        p2yXcPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Cos (n * psiEPc [i]));
                        p2yYaPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Sin (n * psiEPa [i]));
                        p2yYbPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Sin (n * psiEPb [i]));
                        p2yYcPtCent_EP [i] -> Fill (pt, cent, sign * y *TMath::Sin (n * psiEPc [i]));
                        p2xYaPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Sin (n * psiEPa [i]));
                        p2xYbPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Sin (n * psiEPb [i]));
                        p2xYcPtCent_EP [i] -> Fill (pt, cent, sign * x * TMath::Sin (n * psiEPc [i]));

                        p2xXaPtMult_SP [i] -> Fill (pt, mh, sign * x * Xa [i]);
                        p2xXbPtMult_SP [i] -> Fill (pt, mh, sign * x * Xb [i]);
                        p2xXcPtMult_SP [i] -> Fill (pt, mh, sign * x * Xc [i]);
                        p2yXaPtMult_SP [i] -> Fill (pt, mh, sign * y * Xa [i]);
                        p2yXbPtMult_SP [i] -> Fill (pt, mh, sign * y * Xb [i]);
                        p2yXcPtMult_SP [i] -> Fill (pt, mh, sign * y * Xc [i]);
                        p2yYaPtMult_SP [i] -> Fill (pt, mh, sign * y * Ya [i]);
                        p2yYbPtMult_SP [i] -> Fill (pt, mh, sign * y * Yb [i]);
                        p2yYcPtMult_SP [i] -> Fill (pt, mh, sign * y * Yc [i]);
                        p2xYaPtMult_SP [i] -> Fill (pt, mh, sign * x * Ya [i]);
                        p2xYbPtMult_SP [i] -> Fill (pt, mh, sign * x * Yb [i]);
                        p2xYcPtMult_SP [i] -> Fill (pt, mh, sign * x * Yc [i]);

                        p2xXaPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Cos (n * psiEPa [i]));
                        p2xXbPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Cos (n * psiEPb [i]));
                        p2xXcPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Cos (n * psiEPc [i]));
                        p2yXaPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Cos (n * psiEPa [i]));
                        p2yXbPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Cos (n * psiEPb [i]));
                        p2yXcPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Cos (n * psiEPc [i]));
                        p2yYaPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Sin (n * psiEPa [i]));
                        p2yYbPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Sin (n * psiEPb [i]));
                        p2yYcPtMult_EP [i] -> Fill (pt, mh, sign * y *TMath::Sin (n * psiEPc [i]));
                        p2xYaPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Sin (n * psiEPa [i]));
                        p2xYbPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Sin (n * psiEPb [i]));
                        p2xYcPtMult_EP [i] -> Fill (pt, mh, sign * x * TMath::Sin (n * psiEPc [i]));
                    }
                }

                if (eta >= etaMin_ && eta <= etaMax_) {
                    if (pt >= ptAveragingRange_ [i][0] && pt <= ptAveragingRange_ [i][1]) {
                        p2xXaEtaCent_SP [i] -> Fill (eta, cent, x * Xa [i]);
                        p2xXbEtaCent_SP [i] -> Fill (eta, cent, x * Xb [i]);
                        p2xXcEtaCent_SP [i] -> Fill (eta, cent, x * Xc [i]);
                        p2yXaEtaCent_SP [i] -> Fill (eta, cent, y * Xa [i]);
                        p2yXbEtaCent_SP [i] -> Fill (eta, cent, y * Xb [i]);
                        p2yXcEtaCent_SP [i] -> Fill (eta, cent, y * Xc [i]);
                        p2yYaEtaCent_SP [i] -> Fill (eta, cent, y * Ya [i]);
                        p2yYbEtaCent_SP [i] -> Fill (eta, cent, y * Yb [i]);
                        p2yYcEtaCent_SP [i] -> Fill (eta, cent, y * Yc [i]);
                        p2xYaEtaCent_SP [i] -> Fill (eta, cent, x * Ya [i]);
                        p2xYbEtaCent_SP [i] -> Fill (eta, cent, x * Yb [i]);
                        p2xYcEtaCent_SP [i] -> Fill (eta, cent, x * Yc [i]);

                        p2xXaEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Cos (n * psiEPa [i]));
                        p2xXbEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Cos (n * psiEPb [i]));
                        p2xXcEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Cos (n * psiEPc [i]));
                        p2yXaEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Cos (n * psiEPa [i]));
                        p2yXbEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Cos (n * psiEPb [i]));
                        p2yXcEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Cos (n * psiEPc [i]));
                        p2yYaEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Sin (n * psiEPa [i]));
                        p2yYbEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Sin (n * psiEPb [i]));
                        p2yYcEtaCent_EP [i] -> Fill (eta, cent, y * TMath::Sin (n * psiEPc [i]));
                        p2xYaEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Sin (n * psiEPa [i]));
                        p2xYbEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Sin (n * psiEPb [i]));
                        p2xYcEtaCent_EP [i] -> Fill (eta, cent, x * TMath::Sin (n * psiEPc [i]));

                        p2xXaEtaMult_SP [i] -> Fill (eta, mh, x * Xa [i]);
                        p2xXbEtaMult_SP [i] -> Fill (eta, mh, x * Xb [i]);
                        p2xXcEtaMult_SP [i] -> Fill (eta, mh, x * Xc [i]);
                        p2yXaEtaMult_SP [i] -> Fill (eta, mh, y * Xa [i]);
                        p2yXbEtaMult_SP [i] -> Fill (eta, mh, y * Xb [i]);
                        p2yXcEtaMult_SP [i] -> Fill (eta, mh, y * Xc [i]);
                        p2yYaEtaMult_SP [i] -> Fill (eta, mh, y * Ya [i]);
                        p2yYbEtaMult_SP [i] -> Fill (eta, mh, y * Yb [i]);
                        p2yYcEtaMult_SP [i] -> Fill (eta, mh, y * Yc [i]);
                        p2xYaEtaMult_SP [i] -> Fill (eta, mh, x * Ya [i]);
                        p2xYbEtaMult_SP [i] -> Fill (eta, mh, x * Yb [i]);
                        p2xYcEtaMult_SP [i] -> Fill (eta, mh, x * Yc [i]);

                        p2xXaEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Cos (n * psiEPa [i]));
                        p2xXbEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Cos (n * psiEPb [i]));
                        p2xXcEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Cos (n * psiEPc [i]));
                        p2yXaEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Cos (n * psiEPa [i]));
                        p2yXbEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Cos (n * psiEPb [i]));
                        p2yXcEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Cos (n * psiEPc [i]));
                        p2yYaEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Sin (n * psiEPa [i]));
                        p2yYbEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Sin (n * psiEPb [i]));
                        p2yYcEtaMult_EP [i] -> Fill (eta, mh, y * TMath::Sin (n * psiEPc [i]));
                        p2xYaEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Sin (n * psiEPa [i]));
                        p2xYbEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Sin (n * psiEPb [i]));
                        p2xYcEtaMult_EP [i] -> Fill (eta, mh, x * TMath::Sin (n * psiEPc [i]));
                    }
                }

                // sampling
                for (Int_t s = 0; s < sMax; s++) { // loop over samples
                    if (samplingMethod_ == kBootStrapping) bsIndex = s;

                    if (pt >= ptAveragingRange_ [i][0] && pt <= ptAveragingRange_ [i][1]
                        && eta >= etaAveragingRange_ [i][0] && eta <= etaAveragingRange_ [i][1]) {

                        p2xXaCent_SP [i] -> Fill (cent, bsIndex, sign * x * Xa [i], W [bsIndex]);
                        p2xXbCent_SP [i] -> Fill (cent, bsIndex, sign * x * Xb [i], W [bsIndex]);
                        p2xXcCent_SP [i] -> Fill (cent, bsIndex, sign * x * Xc [i], W [bsIndex]);
                        p2yXaCent_SP [i] -> Fill (cent, bsIndex, sign * y * Xa [i], W [bsIndex]);
                        p2yXbCent_SP [i] -> Fill (cent, bsIndex, sign * y * Xb [i], W [bsIndex]);
                        p2yXcCent_SP [i] -> Fill (cent, bsIndex, sign * y * Xc [i], W [bsIndex]);
                        p2yYaCent_SP [i] -> Fill (cent, bsIndex, sign * y * Ya [i], W [bsIndex]);
                        p2yYbCent_SP [i] -> Fill (cent, bsIndex, sign * y * Yb [i], W [bsIndex]);
                        p2yYcCent_SP [i] -> Fill (cent, bsIndex, sign * y * Yc [i], W [bsIndex]);
                        p2xYaCent_SP [i] -> Fill (cent, bsIndex, sign * x * Ya [i], W [bsIndex]);
                        p2xYbCent_SP [i] -> Fill (cent, bsIndex, sign * x * Yb [i], W [bsIndex]);
                        p2xYcCent_SP [i] -> Fill (cent, bsIndex, sign * x * Yc [i], W [bsIndex]);

                        p2xXaCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                        p2xXbCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                        p2xXcCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                        p2xYaCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                        p2xYbCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                        p2xYcCent_EP [i] -> Fill (cent, bsIndex, sign * x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        p2yYaCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                        p2yYbCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                        p2yYcCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        p2yXaCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                        p2yXbCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                        p2yXcCent_EP [i] -> Fill (cent, bsIndex, sign * y * TMath::Cos (n * psiEPc [i]), W [bsIndex]);

                        p2xXaMult_SP [i] -> Fill (mh, bsIndex, sign * x * Xa [i], W [bsIndex]);
                        p2xXbMult_SP [i] -> Fill (mh, bsIndex, sign * x * Xb [i], W [bsIndex]);
                        p2xXcMult_SP [i] -> Fill (mh, bsIndex, sign * x * Xc [i], W [bsIndex]);
                        p2yXaMult_SP [i] -> Fill (mh, bsIndex, sign * y * Xa [i], W [bsIndex]);
                        p2yXbMult_SP [i] -> Fill (mh, bsIndex, sign * y * Xb [i], W [bsIndex]);
                        p2yXcMult_SP [i] -> Fill (mh, bsIndex, sign * y * Xc [i], W [bsIndex]);
                        p2yYaMult_SP [i] -> Fill (mh, bsIndex, sign * y * Ya [i], W [bsIndex]);
                        p2yYbMult_SP [i] -> Fill (mh, bsIndex, sign * y * Yb [i], W [bsIndex]);
                        p2yYcMult_SP [i] -> Fill (mh, bsIndex, sign * y * Yc [i], W [bsIndex]);
                        p2xYaMult_SP [i] -> Fill (mh, bsIndex, sign * x * Ya [i], W [bsIndex]);
                        p2xYbMult_SP [i] -> Fill (mh, bsIndex, sign * x * Yb [i], W [bsIndex]);
                        p2xYcMult_SP [i] -> Fill (mh, bsIndex, sign * x * Yc [i], W [bsIndex]);

                        p2xXaMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                        p2xXbMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                        p2xXcMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                        p2xYaMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                        p2xYbMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                        p2xYcMult_EP [i] -> Fill (mh, bsIndex, sign * x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        p2yYaMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                        p2yYbMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                        p2yYcMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        p2yXaMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                        p2yXbMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                        p2yXcMult_EP [i] -> Fill (mh, bsIndex, sign * y * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                    }

                    if (pt >= ptMin_ && pt <= ptMax_) {
                        if (eta >= etaAveragingRange_ [i][0] && eta <= etaAveragingRange_ [i][1]) {
                            p3xXaPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Xa [i], W [bsIndex]);
                            p3xXbPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Xb [i], W [bsIndex]);
                            p3xXcPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Xc [i], W [bsIndex]);
                            p3yXaPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Xa [i], W [bsIndex]);
                            p3yXbPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Xb [i], W [bsIndex]);
                            p3yXcPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Xc [i], W [bsIndex]);
                            p3yYaPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Ya [i], W [bsIndex]);
                            p3yYbPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Yb [i], W [bsIndex]);
                            p3yYcPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * y * Yc [i], W [bsIndex]);
                            p3xYaPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Ya [i], W [bsIndex]);
                            p3xYbPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Yb [i], W [bsIndex]);
                            p3xYcPtCent_SP [i] -> Fill (pt, cent, bsIndex, sign * x * Yc [i], W [bsIndex]);

                            p3xXaPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3xXbPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3xXcPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yXaPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3yXbPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3yXcPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yYaPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3yYbPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3yYcPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * y *TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                            p3xYaPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3xYbPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3xYcPtCent_EP [i] -> Fill (pt, cent, bsIndex, sign * x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);

                            p3xXaPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Xa [i], W [bsIndex]);
                            p3xXbPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Xb [i], W [bsIndex]);
                            p3xXcPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Xc [i], W [bsIndex]);
                            p3yXaPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Xa [i], W [bsIndex]);
                            p3yXbPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Xb [i], W [bsIndex]);
                            p3yXcPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Xc [i], W [bsIndex]);
                            p3yYaPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Ya [i], W [bsIndex]);
                            p3yYbPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Yb [i], W [bsIndex]);
                            p3yYcPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * y * Yc [i], W [bsIndex]);
                            p3xYaPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Ya [i], W [bsIndex]);
                            p3xYbPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Yb [i], W [bsIndex]);
                            p3xYcPtMult_SP [i] -> Fill (pt, mh, bsIndex, sign * x * Yc [i], W [bsIndex]);

                            p3xXaPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3xXbPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3xXcPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yXaPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3yXbPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3yXcPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yYaPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3yYbPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3yYcPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * y *TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                            p3xYaPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3xYbPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3xYcPtMult_EP [i] -> Fill (pt, mh, bsIndex, sign * x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        }
                    }

                    if (eta >= etaMin_ && eta <= etaMax_) {
                        if (pt >= ptAveragingRange_ [i][0] && pt <= ptAveragingRange_ [i][1]) {
                            p3xXaEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Xa [i], W [bsIndex]);
                            p3xXbEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Xb [i], W [bsIndex]);
                            p3xXcEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Xc [i], W [bsIndex]);
                            p3yXaEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Xa [i], W [bsIndex]);
                            p3yXbEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Xb [i], W [bsIndex]);
                            p3yXcEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Xc [i], W [bsIndex]);
                            p3yYaEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Ya [i], W [bsIndex]);
                            p3yYbEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Yb [i], W [bsIndex]);
                            p3yYcEtaCent_SP [i] -> Fill (eta, cent, bsIndex, y * Yc [i], W [bsIndex]);
                            p3xYaEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Ya [i], W [bsIndex]);
                            p3xYbEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Yb [i], W [bsIndex]);
                            p3xYcEtaCent_SP [i] -> Fill (eta, cent, bsIndex, x * Yc [i], W [bsIndex]);

                            p3xXaEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3xXbEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3xXcEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yXaEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3yXbEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3yXcEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yYaEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3yYbEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3yYcEtaCent_EP [i] -> Fill (eta, cent, bsIndex, y * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                            p3xYaEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3xYbEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3xYcEtaCent_EP [i] -> Fill (eta, cent, bsIndex, x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);

                            p3xXaEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Xa [i], W [bsIndex]);
                            p3xXbEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Xb [i], W [bsIndex]);
                            p3xXcEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Xc [i], W [bsIndex]);
                            p3yXaEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Xa [i], W [bsIndex]);
                            p3yXbEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Xb [i], W [bsIndex]);
                            p3yXcEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Xc [i], W [bsIndex]);
                            p3yYaEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Ya [i], W [bsIndex]);
                            p3yYbEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Yb [i], W [bsIndex]);
                            p3yYcEtaMult_SP [i] -> Fill (eta, mh, bsIndex, y * Yc [i], W [bsIndex]);
                            p3xYaEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Ya [i], W [bsIndex]);
                            p3xYbEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Yb [i], W [bsIndex]);
                            p3xYcEtaMult_SP [i] -> Fill (eta, mh, bsIndex, x * Yc [i], W [bsIndex]);

                            p3xXaEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3xXbEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3xXcEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yXaEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Cos (n * psiEPa [i]), W [bsIndex]);
                            p3yXbEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Cos (n * psiEPb [i]), W [bsIndex]);
                            p3yXcEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Cos (n * psiEPc [i]), W [bsIndex]);
                            p3yYaEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3yYbEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3yYcEtaMult_EP [i] -> Fill (eta, mh, bsIndex, y * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                            p3xYaEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Sin (n * psiEPa [i]), W [bsIndex]);
                            p3xYbEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Sin (n * psiEPb [i]), W [bsIndex]);
                            p3xYcEtaMult_EP [i] -> Fill (eta, mh, bsIndex, x * TMath::Sin (n * psiEPc [i]), W [bsIndex]);
                        }
                    }
                }

                if (subeventFlag [i][itrack - 1] == 1) {
					Xa [i] += x / mha [i];
					Ya [i] += y / mha [i];
					psiEPa [i] = TMath::ATan2 (Ya [i], Xa [i]) / n;
				}

				if (subeventFlag [i][itrack - 1] == 2) {
					Xb [i] += x / mhb [i];
					Yb [i] += y / mhb [i];
					psiEPb [i] = TMath::ATan2 (Yb [i], Xb [i]) / n;
				}

				if (subeventFlag [i][itrack - 1] == 3) {
					Xc [i] += x / mhc [i];
					Yc [i] += y / mhc [i];
					psiEPc [i] = TMath::ATan2 (Yc [i], Xc [i]) / n;
				}

//				if (Xa [i] != 0.0) Xa [i] -= ptAvgA [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] -= ptAvgB [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] -= ptAvgC [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] -= ptAvgFixedA [centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] -= ptAvgFixedB [centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] -= ptAvgFixedC [centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] -= ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] -= ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] -= ptAvg [i][centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];

//				if (Xa [i] != 0.0) Xa [i] -= ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mha [i];
//				if (Xb [i] != 0.0) Xb [i] -= ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mhb [i];
//				if (Xc [i] != 0.0) Xc [i] -= ptAvgFixed [centBin - 1] / pt2Avg [centBin - 1] * pt / mhc [i];
//
//				pPtCent [i] -> Reset ("ICESM");
//				pPtCentA [i] -> Reset ("ICESM");
//				pPtCentB[i] -> Reset ("ICESM");
//				pPtCentC [i] -> Reset ("ICESM");
			}
		}
        testTree -> Fill (); // test
        if (samplingMethod_ == kSubsampling) W [bsIndex] = 0; // subsampling
	}
    cout << endl;

	inputFile -> Close ();
	servHistDir -> Write ();
	histFile -> Close ();

    testFile -> cd (); // test
    testTree -> Write (); // test
    testFile -> Close (); // test
    if (samplingMethod_ == kBootStrapping) sampleFile -> Close ();

//    for (Int_t i = 0; i < nHarmonics; i++) {
//	    QnMan -> GetDetectorQnVectorList (Form ("D1A_%i", n)) -> Print ("",-1);
//	    QnMan -> GetDetectorQnVectorList (Form ("D1B_%i", n)) -> Print ("",-1);
//	    QnMan -> GetDetectorQnVectorList (Form ("D1C_%i", n)) -> Print ("",-1);
//	}

	FinalizeQnCorrectionsManager (qnInputFile, qnOutputFile, QnMan);
	FinalizeQnCorrectionsManager (qnPtInputFile, qnPtOutputFile, QnManPt);
	FinalizeQnCorrectionsManager (qnEtaInputFile, qnEtaOutputFile, QnManEta);
}


void CFlowReconstructor::CalculateCorrelationsWithSampling (TProfile2D *p2XaXb, TProfile *pXaXb) {
    Int_t nBinsX = p2XaXb -> GetNbinsX ();
    Float_t x;
    for (Int_t j = 1; j <= nBinsX; j++) {
        x = pXaXb -> GetXaxis () -> GetBinCenter (j);
        for (Int_t k = 1; k <= nBinsBS_; k++) {
            pXaXb -> Fill (x, p2XaXb -> GetBinContent (j, k));
        }
    }
}

void CFlowReconstructor::CalculateResolutionWithSampling (TProfile2D *p2XaXb, TProfile2D *p2XaXc, TProfile2D *p2XbXc, TH2F *h2Ra, TH2F *h2Rb, TH2F *h2Rc) {
    Int_t l, m;
    Float_t R;
    Int_t nBinsX = p2XaXb -> GetNbinsX ();
    for (Int_t j = 1; j <= nBinsX; j++) {
        for (Int_t k = 1; k <= nBinsBS_; k++) {
            l = k + 3;
            m = k + 7;
            if (l > nBinsBS_) l -= nBinsBS_;
            if (m > nBinsBS_) m -= nBinsBS_;
            if (resMethod_ == kThreeSubevents) {
                R = 2.0 * p2XaXb -> GetBinContent (j, k) *
                        p2XaXc -> GetBinContent (j, k) /
                                  p2XbXc -> GetBinContent (j, k);
    //                    R = 2.0 * p2XaXb -> GetBinContent (j, k) *
    //                              p2XaXc -> GetBinContent (j, l) /
    //                              p2XbXc -> GetBinContent (j, m);
                        h2Ra -> SetBinContent (j, k, R);
    //                    h2Ra -> SetBinContent (j, k, TMath::Abs(R) / R * TMath::Sqrt (TMath::Abs (R)));

                        R = 2.0 * p2XaXb -> GetBinContent (j, k) *
                                  p2XbXc -> GetBinContent (j, k) /
                                  p2XaXc -> GetBinContent (j, k);
    //                    R = 2.0 * p2XaXb -> GetBinContent (j, k) *
    //                              p2XbXc -> GetBinContent (j, l) /
    //                              p2XaXc -> GetBinContent (j, m);
                        h2Rb -> SetBinContent (j, k, R);
    //                    h2Rb -> SetBinContent (j, k, TMath::Abs(R) / R * TMath::Sqrt (TMath::Abs (R)));
    //
                        R = 2.0 * p2XaXc -> GetBinContent (j, k) *
                                  p2XbXc -> GetBinContent (j, k) /
                                  p2XaXb -> GetBinContent (j, k);
    //                    R = 2.0 * p2XaXc -> GetBinContent (j, k) *
    //                              p2XbXc -> GetBinContent (j, l) /
    //                              p2XaXb -> GetBinContent (j, m);
                        h2Rc -> SetBinContent (j, k, R);
    //                    h2Rc -> SetBinContent (j, k, TMath::Abs(R) / R * TMath::Sqrt (TMath::Abs (R)));
                    }

                    if (resMethod_ == kRandomSubevent) {
                        R = 2.0 * p2XaXb -> GetBinContent (j, k);
                        h2Ra -> SetBinContent (j, k, R);
                        h2Rb -> SetBinContent (j, k, R);
    //                    h2Ra -> SetBinContent (j, k, TMath::Abs(R) / R * TMath::Sqrt (TMath::Abs (R)));
    //                    h2Rb -> SetBinContent (j, k, TMath::Abs(R) / R * TMath::Sqrt (TMath::Abs (R)));
                    }
                }
    }
}

void CFlowReconstructor::CalculateResolutionNoSampling (TProfile *pXaXb, TProfile *pXaXc, TProfile *pXbXc, TH1F *hRa, TH1F *hRb, TH1F *hRc) {
    if (resMethod_ == kThreeSubevents) {
        pXaXb -> Sumw2 ();
        pXaXc -> Sumw2 ();
        pXbXc -> Sumw2 ();
        hRa -> Multiply (pXaXb, pXaXc, 2.0);
        hRa -> Divide (pXbXc);
        Acc::CheckHistogram(hRa);
        hRb -> Multiply (pXaXb, pXbXc, 2.0);
        hRb -> Divide (pXaXc);
        Acc::CheckHistogram(hRb);
        hRc -> Multiply (pXaXc, pXbXc, 2.0);
        hRc -> Divide (pXaXb);
        Acc::CheckHistogram(hRc);
    }
    if (resMethod_ == kRandomSubevent) {
        Int_t nBins = hRa -> GetNbinsX ();
        Float_t R, Rerr;
        for (Int_t i = 1; i <= nBins; i++) {
            R = pXaXb -> GetBinContent (i);
            Rerr = pXaXb -> GetBinError (i);
            hRa -> SetBinContent (i, R);
            hRa -> SetBinError (i, Rerr);
            hRb -> SetBinContent (i, R);
            hRb -> SetBinError (i, Rerr);
        }
        hRa -> Scale (2.0);
        hRb -> Scale (2.0);
        Acc::CheckHistogram(hRa);
        Acc::CheckHistogram(hRb);
    }
}

void CFlowReconstructor::Sqrt (TH1 *hR) {
    Int_t nBinsX = hR -> GetNbinsX ();
    Float_t R, Rerr;
    for (Int_t i = 1; i <= nBinsX; i++) {
        R = TMath::Sqrt (TMath::Abs(hR -> GetBinContent (i)));
        Rerr = TMath::Abs (hR -> GetBinError (i) / hR -> GetBinContent (i) * 0.5 * R);
        if (propagateResolutionSign_) R *= (hR -> GetBinContent (i) / TMath::Abs (hR -> GetBinContent (i)));
        hR -> SetBinContent (i, R);
        hR -> SetBinError (i, Rerr);
    }
}

void CFlowReconstructor::Sqrt (TH2 *h2R) {
    Int_t nBinsX = h2R -> GetNbinsX ();
    Int_t nBinsY = h2R -> GetNbinsY ();
    Float_t R, Rerr;
    for (Int_t i = 1; i <= nBinsX; i++) {
        for (Int_t j = 1; j <= nBinsY; j++) {
            R = TMath::Sqrt (TMath::Abs(h2R -> GetBinContent (i, j)));
            Rerr = TMath::Abs (h2R -> GetBinError (i, j) / h2R -> GetBinContent (i, j) * 0.5 * R);
            if (propagateResolutionSign_) R *= (h2R -> GetBinContent (i, j) / TMath::Abs (h2R -> GetBinContent (i, j)));
            h2R -> SetBinContent (i, j, R);
            h2R -> SetBinError (i, j, Rerr);
        }
    }
}


void CFlowReconstructor::CalculateFlowNoSampling (TProfile *pxX, TH1 *hR, TH1 *hV, Int_t sign) {
    pxX -> Sumw2 ();
    hR -> Sumw2 ();
    hV -> Divide (pxX, hR, 2.0);
    hV -> Scale (sign);
}


void CFlowReconstructor::CalculateFlowNoSampling (TProfile2D *p2xX, TH1 *hR, TH1 *hV, Int_t lowerBin, Int_t higherBin, Int_t sign) {
    Int_t nBinsX = p2xX -> GetNbinsX ();
    Int_t nClasses = higherBin - lowerBin + 1;
    Float_t xX, R, V, sample, xXerr, Rerr, Verr;

    for (Int_t i = 1; i <= nBinsX; i++) {
        V = 0.0;
        Verr = 0.0;
        for (Int_t j = lowerBin; j <= higherBin; j++) { // Event classes
            xX = p2xX -> GetBinContent (i, j);
            xXerr = p2xX -> GetBinError (i, j);
            R = hR -> GetBinContent (j);
            Rerr = hR -> GetBinError (j);
            V += 2.0 * xX / R;
            Verr += V * V * (xXerr * xXerr + Rerr * Rerr);
        }
        V /= (Float_t) nClasses;
        Verr = TMath::Sqrt (Verr);
        hV -> SetBinContent (i, V);
        hV -> SetBinError (i, Verr);
    }
    hV -> Scale (sign);
}

void CFlowReconstructor::CalculateFlowWithSampling (TProfile2D *p2xX, TH2 *h2R, TProfile2D *p2V, Int_t sign) {
    Int_t nBinsX = p2xX -> GetNbinsX ();
    Int_t nSamples = p2xX -> GetNbinsY ();
    Float_t xX, R, V, x, sample;

    for (Int_t i = 1; i <= nBinsX; i++) {
        x = p2xX -> GetXaxis () -> GetBinCenter (i);
        for (Int_t j = 1; j <= nSamples; j++) {
            sample = p2xX -> GetZaxis () -> GetBinCenter (j);
            xX = p2xX -> GetBinContent (i, j);
            R = h2R -> GetBinContent (i, j);
            V = 2.0 * xX / R;
            p2V -> Fill (x, sample, V);
        }
    }
    p2V -> Scale (sign);
}


void CFlowReconstructor::CalculateFlowWithSampling (TProfile3D *p3xX, TH2 *h2R, TProfile2D *p2V, Int_t lowerBin, Int_t higherBin, Int_t sign) {
    Int_t nBinsX = p3xX -> GetNbinsX ();
    Int_t nSamples = p3xX -> GetNbinsZ ();
    Float_t xX, R, V, x, sample;

    for (Int_t i = 1; i <= nBinsX; i++) {
        x = p3xX -> GetXaxis () -> GetBinCenter (i);
        for (Int_t j = lowerBin; j <= higherBin; j++) { // Event classes
            for (Int_t k = 1; k <= nSamples; k++) {
                sample = p3xX -> GetZaxis () -> GetBinCenter (k);
                xX = p3xX -> GetBinContent (i, j, k);
                R = h2R -> GetBinContent (j, k);
                V = 2.0 * xX / R;
                p2V -> Fill (x, sample, V);
            }
        }
    }
    p2V -> Scale (sign);
}


void CFlowReconstructor::ReflectRapidity (TH1F *hVEta, TH1F *hVEtaRefl, Int_t nHarmonic) {
    if (nBinsEtaRefl_ == 0) return;
    Int_t sign = TMath::Power (-1.0, nHarmonic);
    for (Int_t i = 0; i < nBinsEtaRefl_; i++) {
        hVEtaRefl -> SetBinContent (i + 1, sign * hVEta -> GetBinContent (nBinsEta_ - i));
        hVEtaRefl -> SetBinError (i + 1, hVEta -> GetBinError (nBinsEta_ - i));
    }
}

void CFlowReconstructor::PlotKinematics (TFile *corrFile, TDirectory *outputDir, Int_t nHarmonic, Int_t step) {
    outputDir -> cd ();
    TCanvas *c;
    TH1 *hList [4];
    TH2 *h2List [4];

    h2List [0] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2mhCent%i", nHarmonic));
    h2List [1] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2mhaCent%i", nHarmonic));
    h2List [2] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2mhbCent%i", nHarmonic));
    h2List [3] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2mhcCent%i", nHarmonic));

    c = new TCanvas (Form ("cMhSubCent%i", nHarmonic), Form ("cMhSubCent%i", nHarmonic), 800, 600);
    c -> Divide (2, 2);
    for (Int_t j = 0; j < 4; j++) {
        c -> cd (j + 1);
        h2List [j] -> Draw ("colz");
        gStyle -> SetOptStat (0);
        gROOT -> ForceStyle();
        gPad -> SetLogz ();
    }
    c -> Write ();
    delete c;

    h2List [0] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2PtEta_%i", nHarmonic));
    h2List [1] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2PtEtaA_%i", nHarmonic));
    h2List [2] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2PtEtaB_%i", nHarmonic));
    h2List [3] = (TH2F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("h2PtEtaC_%i", nHarmonic));

    c = new TCanvas (Form ("cPtEta2D_%i", nHarmonic), Form ("cPtEta2D_%i", nHarmonic), 800, 600);
    c -> Divide (2, 2);
    for (Int_t j = 0; j < 4; j++) {
        c -> cd (j + 1);
        h2List [j] -> Draw ("colz");
        gPad -> SetLogz ();
    }
    c -> Write ();
    delete c;

    hList [0] = (TH1F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pPtEta_%i", nHarmonic));
    hList [1] = (TH1F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pPtEtaA_%i", nHarmonic));
    hList [2] = (TH1F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pPtEtaB_%i", nHarmonic));
    hList [3] = (TH1F*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pPtEtaC_%i", nHarmonic));

    c = new TCanvas (Form ("PtEta1D_%i", nHarmonic), Form ("PtEta1D_%i", nHarmonic), 800, 600);
    gStyle -> SetOptStat (1);
    c -> Divide (2, 2);
    for (Int_t j = 0; j < 4; j++) {
        c -> cd (j + 1);
        hList [j] -> Draw ();
    }
    c -> Write ();
    delete c;
}

void CFlowReconstructor::SetResolutionMethod (Int_t resMethod) {
    resMethod_ = resMethod;
}

void CFlowReconstructor::PlotResolution (TH1 *hList1 [12], TH1 *hList2 [12], TH1 *hList3 [9], TH1 *hList4 [9], Int_t nHist, TDirectory *dir) {
//    Int_t markerColors [8] = {1, 2, 3, 4, 6, 8, 9, 46};
//    Int_t markerStyles [8] = {24, 25, 26, 27, 28, 30, 32, 5};
    Int_t markerColors [7] = {1, 2, 3, 4, 2, 12, 28};
    Int_t markerStyles1 [14] = {24, 25, 26, 32, 30, 27, 28, 20, 21, 22, 23, 29, 33, 34};
    Int_t markerStyles2 [14] = {20, 21, 22, 23, 29, 33, 34, 24, 25, 26, 32, 30, 27, 28};
    char sub [3] = {'a', 'b', 'c'}, sub1 [3] = {'a', 'a', 'b'}, sub2 [3] = {'b', 'c', 'c'};
    char q [] = {'x', 'y', 'q'}, Q [] = {'X', 'Y', 'Q'};
    dir -> cd ();

    Int_t k, nHarmonic;
    TString method, eventClassVariable, histName = hList1 [0] -> GetName ();
    if (histName.Contains ("SP")) method = "SP";
    if (histName.Contains ("EP")) method = "EP";
    if (histName.Contains ("Cent")) eventClassVariable = "Cent";
    if (histName.Contains ("Mult")) eventClassVariable = "Mult";

    for (Int_t i = 0; i < nHarmonics; i++) {
        nHarmonic = harmonicsMap [i];
        if (histName.Contains (Form ("%i", nHarmonic))) break;
    }

    Float_t shift;
    Int_t iMax = 3;
    if (resMethod_ == kRandomSubevent) iMax = 1;
    TCanvas *c;
    gStyle -> SetOptStat (0);
    gStyle -> SetLegendBorderSize (0);
    TString xAxisTitle = hList2 [0] -> GetXaxis () -> GetTitle ();
    THStack *hs, *hsa [4];


    for (Int_t i = 0; i < iMax; i++) { // plot three * four QQ combinations
        c = new TCanvas (Form ("cQ%i%cQ%i%c", nHarmonic, sub1 [i], nHarmonic, sub2 [i]) + eventClassVariable + "_" + method, Form ("cQ%i%cQ%i%c", nHarmonic, sub1 [i], nHarmonic, sub2 [i]) + eventClassVariable + "_" + method, 800, 600);
        hs = new THStack (Form ("hsQ%i%cQ%i%c", nHarmonic, sub1 [i], nHarmonic, sub2 [i]) + eventClassVariable + "_" + method, Form ("#LTQ_{%i}^{%c}Q_{%i}^{%c}#GT (", nHarmonic, sub1 [i], nHarmonic, sub2 [i]) + method + ");" + xAxisTitle);

        for (Int_t j = 0; j < 4; j++) {
            k = i * 4 + j;
            shift = -0.175 + 0.1 * j; // -0.175 -0.075 0.025 0.125
            hList1 [k] -> SetMarkerStyle (markerStyles1 [j]);
            hList1 [k] -> SetMarkerColor (markerColors [j]);
            hList1 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList1 [k], shift);
            hs -> Add (hList1 [k]);
//            HistShift (hList1 [k], -shift);

            // sampling
            hList2 [k] -> SetMarkerStyle (markerStyles2 [j]);
            hList2 [k] -> SetMarkerColor (markerColors [j]);
            hList2 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList2 [k], shift + 0.05);
            hs -> Add (hList2 [k]);
//            HistShift (hList2 [k], -shift - 0.05);
        }
        hs -> Draw ("nostack p e1X0");
        gPad -> BuildLegend () -> SetFillColor (0);
//        gPad -> BuildLegend (0.5, 0.2, 0.85, 0.5) -> SetFillColor (0);
        c -> Write ();
        delete c;
        delete hs;
        gStyle -> SetOptStat (1);
    }

    // plot resolution correction factors
    iMax = 3;
    if (resMethod_ == kRandomSubevent) iMax = 2;
    for (Int_t i = 0; i < iMax; i++) {
        c = new TCanvas (Form ("cR%i%c", nHarmonic, sub [i]) + eventClassVariable + "_" + method, Form ("cR%i%c", nHarmonic, sub [i]) + eventClassVariable + "_" + method, 800, 600);
        gStyle -> SetLegendBorderSize (0);
        hs = new THStack (Form ("R%i%c", nHarmonic, sub [i]) + eventClassVariable + "_" + method, Form ("R_{%i}^{%c,", nHarmonic, sub [i]) + method + "};" + xAxisTitle + ";" + Form ("R_{%i}^{%c,", nHarmonic, sub [i]) + method + "}");
        for (Int_t j = 0; j < 3; j++) {
            k = i * 3 + j;
            shift = -0.125 + 0.1 * j; // -0.125 -0.025 0.075
            hList3 [k] -> SetMarkerStyle (markerStyles1 [j]);
            hList3 [k] -> SetMarkerColor (markerColors [j]);
            hList3 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList3 [k], shift);
            hs -> Add (hList3 [k]);
//            HistShift (hList3 [k], -shift);

            // sampling
            hList4 [k] -> SetMarkerStyle (markerStyles2 [j]);
            hList4 [k] -> SetMarkerColor (markerColors [j]);
            hList4 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList4 [k], shift + 0.05);
            hs -> Add (hList4 [k]);
//            HistShift (hList4 [k], -shift - 0.05);
        }
        hs -> Draw ("nostack p e1X0");
        gPad -> BuildLegend () -> SetFillColor (0);
    //    gPad -> BuildLegend (0.5, 0.2, 0.85, 0.5) -> SetFillColor (0);
        c -> Write ();
        delete c;
        delete hs;
    }

    for (Int_t i = 0; i < 2; i++) { // plot correlations and resolution correction factors (4 pads)
        c = new TCanvas (Form ("cR%i%c%c", nHarmonic, Q [i], Q [i]) + eventClassVariable + "_" + method, Form ("R%i%c%c", nHarmonic, Q [i], Q [i]) + eventClassVariable + "_" + method, 800, 600);
        c -> Divide (2, 2);
        shift = -0.025; // -0.025
        for (Int_t j = 0; j < 3; j++) {
            c -> cd (j + 1);
            hsa [j] = new THStack (Form ("%c%i%c%c%i%c", Q [i], nHarmonic, sub1 [j], Q [i], nHarmonic, sub2 [j]) + eventClassVariable + "_" + method, Form ("#LT%c_{%i}^{%c}%c_{%i}^{%c}#GT (", Q [i], nHarmonic, sub1 [j], Q [i], nHarmonic, sub2 [j]) + method + ");" + xAxisTitle);
            k = j * 4 + i;
            hList1 [k] -> SetMarkerStyle (markerStyles1 [j]);
            hList1 [k] -> SetMarkerColor (markerColors [j]);
            hList1 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList1 [k], shift);
            hsa [j] -> Add (hList1 [k]);
//            HistShift (hList1 [k], - shift);

            // sampling
            hList2 [k] -> SetMarkerStyle (markerStyles2 [j]);
            hList2 [k] -> SetMarkerColor (markerColors [j]);
            hList2 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList2 [k], shift + 0.05);
            hsa [j] -> Add (hList2 [k]);
//            HistShift (hList2 [k], - shift - 0.05);

            hsa [j] -> Draw ("nostack p e1X0");
            gPad -> BuildLegend () -> SetFillColor (0);
//            gPad -> BuildLegend (0.5, 0.2, 0.85, 0.5) -> SetFillColor (0);
        }
        c -> cd (4);
        hsa [3] = new THStack (Form ("R%i%c", nHarmonic, q [i]) + eventClassVariable + "_" + method, Form ("R_{%i, x}^{", nHarmonic) + method + "};" + xAxisTitle + ";" + Form ("R_{%i, x}^{", nHarmonic) + method + "}");
        for (Int_t j = 0; j < iMax; j++) {
            k = j * 3 + i;
            hList3 [k] -> SetMarkerStyle (markerStyles1 [j]);
            hList3 [k] -> SetMarkerColor (markerColors [j]);
            hList3 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList3 [k], shift);
            hsa [3] -> Add (hList3 [k]);
//            HistShift (hList3 [k], - shift);

            // sampling
            hList4 [k] -> SetMarkerStyle (markerStyles2 [j]);
            hList4 [k] -> SetMarkerColor (markerColors [j]);
            hList4 [k] -> SetLineColor (markerColors [j]);
            HistShift (hList4 [k], shift + 0.05);
            hsa [3] -> Add (hList4 [k]);
//            HistShift (hList4 [k], - shift - 0.05);
        }
        hsa [3] -> Draw ("nostack p e1X0");
        gPad -> BuildLegend () -> SetFillColor (0);
//      gPad -> BuildLegend (0.5, 0.2, 0.85, 0.5) -> SetFillColor (0);
        c -> Write ();
        for (Int_t j = 0; j < 4; j++) delete hsa [j];
        delete c;
    }
    gStyle -> SetOptStat (1);
}

void CFlowReconstructor::CombineSubevents (TH2 *pa, TH2 *pb, TH2 *pc, TH2 *p) {
    p -> Add (pa, pb);
    if (resMethod_ == kThreeSubevents) p -> Add (pc);
}

void CFlowReconstructor::CombineSubevents (TH1 *ha, TH1 *hb, TH1 *hc, TH1 *h) {
    if (resMethod_ == kThreeSubevents) {
        h -> Add (ha, hb, 0.33, 0.33);
        h -> Add (hc, 0.33);
    }

    if (resMethod_ == kRandomSubevent) {
        h -> Add (ha, hb, 0.5, 0.5);
    }
}

void CFlowReconstructor::PlotFlow (TH1 *hList1 [5], TH1 *hList2 [5], TH1 *hList3 [5], TH1 *hList4 [5]) {
    Int_t markerColors [7] = {1, 2, 3, 4, 2, 12, 28};
    Int_t markerStyles1 [14] = {24, 25, 26, 32, 30, 27, 28, 20, 21, 22, 23, 29, 33, 34};
    Int_t markerStyles2 [14] = {20, 21, 22, 23, 29, 33, 34, 24, 25, 26, 32, 30, 27, 28};

    TString histName = hList1 [0] -> GetName (), histTitle = hList1 [0] -> GetTitle ();
    TString xAxisTitle = hList1 [0] -> GetXaxis () -> GetTitle ();
    TString yAxisTitle = hList1 [0] -> GetYaxis () -> GetTitle ();
    if (histTitle.Contains ("x+y")) {
        for (Int_t i = 0; i < histTitle.Length (); i++) {
            if (histTitle [i] == 'x') {
                histTitle.Remove (i, 5);
            }
        }
    }
    histName.Remove (0, 1);
    TCanvas *c = new TCanvas ("c" + histName, histTitle, 800, 600);
    gStyle -> SetLegendBorderSize (0);
    gStyle -> SetOptStat (0);
    THStack *hs = new THStack ("hs" + histName, histTitle + ";" + xAxisTitle + ";" + yAxisTitle);

    Int_t iMax = 5;
    Float_t shift;

    if (hList1 [4] == 0) {
        if (resMethod_ == kThreeSubevents) iMax = 4;
        else iMax = 3;
    }

    for (Int_t i = 0; i < iMax; i++) {
        shift = -0.2 + 0.1 * i; // -0.2 -0.1 0.0 0.1 0.2
        hList1 [i] -> SetMarkerStyle (markerStyles1 [i]);
        hList1 [i] -> SetMarkerColor (markerColors [i]);
        hList1 [i] -> SetLineColor (markerColors [i]);
        HistShift (hList1 [i], shift);
        hs -> Add (hList1 [i]);
//        HistShift (hList1 [i], -shift);

        if (hList3 != 0) { // reflected
            hList3 [i] -> SetMarkerStyle (markerStyles1 [i]);
            hList3 [i] -> SetMarkerColor (markerColors [i]);
            hList3 [i] -> SetLineColor (markerColors [i]);
            HistShift (hList3 [i], shift);
            hs -> Add (hList3 [i]);
    //        HistShift (hList3 [i], -shift);
        }

        // sampling
        hList2 [i] -> SetMarkerStyle (markerStyles2 [i]);
        hList2 [i] -> SetMarkerColor (markerColors [i]);
        hList2 [i] -> SetLineColor (markerColors [i]);
        HistShift (hList2 [i], shift + 0.05);
        hs -> Add (hList2 [i]);
//        HistShift (hList2 [i], -shift - 0.05);

        if (hList4 != 0) { // reflected
            hList4 [i] -> SetMarkerStyle (markerStyles2 [i]);
            hList4 [i] -> SetMarkerColor (markerColors [i]);
            hList4 [i] -> SetLineColor (markerColors [i]);
            HistShift (hList4 [i], shift + 0.05);
            hs -> Add (hList4 [i]);
    //        HistShift (hList4 [i], -shift - 0.05);
        }
    }

    hs -> Draw ("nostack p e1X0");
//    gPad -> BuildLegend () -> SetFillColor (0);
    gPad -> BuildLegend (0.5, 0.2, 0.85, 0.5) -> SetFillColor (0);
    c -> Write ();

    delete c;
    delete hs;
    gStyle -> SetOptStat (1);
}

void CFlowReconstructor::GetFlow () {
    for (Int_t i = 0; i < nSteps_; i++) {
        GetFlowLoop (i);
    }
}
void CFlowReconstructor::WritePreviousResults (TDirectory *dir) { // put here anything you like
    dir -> cd ();
//    // NA49 40 AGeV 0-20%
//    static const Int_t nyN = 9;
//    double y1N [nyN] = {0.090000, 0.290000, 0.490000, 0.690000, 0.890000, 1.090000, 1.290000, 1.490000, 1.690000};
//    double vy1N [nyN] = {-0.002917, -0.001064, -0.005188, 0.000165, -0.006057, -0.010585, -0.016514, -0.021671, -0.021830};
//    double vey1N [nyN] = {0.002751, 0.002762, 0.002888, 0.003165, 0.003426, 0.003532, 0.003687, 0.004000, 0.004559};
//    double pt1N [13] = {0.060000, 0.160000, 0.260000, 0.360000, 0.460000, 0.560000, 0.660000, 0.760000, 0.910000, 1.110000, 1.310000, 1.510000, 1.810000};
//    double vpt1N [13] = {-0.024160, -0.013253, -0.004783, -0.000505, -0.005080, -0.009384, -0.002377, -0.004053, -0.001248, 0.022148, 0.013643, 0.007168, 0.011767};
//    double vept1N [13] = {0.003133, 0.002263, 0.002610, 0.002996, 0.003476, 0.004130, 0.004980, 0.006008, 0.005742, 0.008518, 0.012902, 0.019462, 0.022513};
//    double y2N [nyN] = {0.090000, 0.290000, 0.490000, 0.690000, 0.890000, 1.090000, 1.290000, 1.490000, 1.690000};
//    double vy2N [nyN] = {-0.009240, -0.000071, -0.000733, -0.009338, -0.010038, 0.007246, 0.045612, 0.051558, -0.026392};
//    double vey2N [nyN] = {0.016173, 0.016304, 0.016799, 0.017105, 0.018528, 0.020074, 0.022294, 0.025597, 0.030496};
//    double pt2N [13] = {0.060000, 0.160000, 0.260000, 0.360000, 0.460000, 0.560000, 0.660000, 0.760000, 0.910000, 1.110000, 1.310000, 1.510000, 1.810000};
//    double vpt2N [13] = {0.013361, -0.013973, -0.008927, 0.003440, -0.018276, -0.008425, 0.016685, 0.040448, 0.040749, 0.032528, 0.124015, 0.095045, 0.133742};
//    double vept2N [13] = {0.021908, 0.013895, 0.015526, 0.016826, 0.019069, 0.022166, 0.026685, 0.032319, 0.030887, 0.046341, 0.070777, 0.107107, 0.112914};
//
//    //STAR 7.7 GeV 0-10%
//    static const Int_t nPtS = 4;
//    double y1S [10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
//    double vy1S [10] = {0.000441316,0.00182928,-0.000142038,-0.00089236,0.000221004,-0.00526669,-0.00246784,-0.00131499,0.000838314,-0.00750076};
//    double vey1S [10] = {0.00435667,0.00256068,0.00221022,0.00216492,0.00218349,0.00217873,0.00215614,0.00220043,0.00258115,0.00451173};
//    double pt2S [nPtS] = {0.375,0.735,1.125,1.545};
//    double vpt2S [nPtS] = {0.0177992,0.0314537,0.0367388,0.0391578};
//    double vept2S [nPtS] = {0.000866997,0.00168551,0.00375154,0.00859142};

     // NA49 40 AGeV 20-40%
    static const Int_t nyN = 9;
    double y1N [nyN] = {0.060000, 0.260000, 0.460000, 0.660000, 0.860000, 1.060000, 1.260000, 1.460000, 1.660000};
    double vy1N [nyN] = {-0.001103, -0.001252, -0.004370, -0.008522, -0.011762, -0.020268, -0.029440, -0.035648, -0.043193};
    double vey1N [nyN] = {0.002042, 0.002036, 0.002102, 0.002293, 0.002483, 0.002568, 0.002685, 0.002905, 0.003290};
    double pt1N [13] = {0.050000, 0.150000, 0.250000, 0.350000, 0.450000, 0.550000, 0.650000, 0.750000, 0.900000, 1.100000, 1.300000, 1.500000, 1.800000};
    double vpt1N [13] = {-0.020437, -0.020135, -0.007285, -0.013533, -0.015285, -0.012257, -0.017684, -0.011000, 0.002180, 0.005412, 0.025699, 0.044994, 0.062799};
    double vept1N [13] = {0.002234, 0.001629, 0.001890, 0.002187, 0.002557, 0.003063, 0.003726, 0.004541, 0.004365, 0.006611, 0.010104, 0.015539, 0.018450};
    double y2N [nyN] = {0.060000, 0.260000, 0.460000, 0.660000, 0.860000, 1.060000, 1.260000, 1.460000, 1.660000};
    double vy2N [nyN] = {0.024770, 0.023833, 0.023160, 0.019666, 0.032268, 0.024772, 0.031206, 0.022143, 0.030117};
    double vey2N [nyN] = {0.004201, 0.004208, 0.004344, 0.004440, 0.004779, 0.005149, 0.005708, 0.006510, 0.007703};
    double pt2N [13] = {0.050000, 0.150000, 0.250000, 0.350000, 0.450000, 0.550000, 0.650000, 0.750000, 0.900000, 1.100000, 1.300000, 1.500000, 1.800000};
    double vpt2N [13] = {-0.003180, 0.010014, 0.016026, 0.024243, 0.032966, 0.041051, 0.049668, 0.053051, 0.065806, 0.093092, 0.117910, 0.113181, 0.130355};
    double vept2N [13] = {0.005445, 0.003514, 0.003980, 0.004329, 0.004934, 0.005813, 0.007037, 0.008615, 0.008376, 0.012768, 0.019785, 0.030081, 0.040669};

    //STAR 7.7 GeV 10-40%
    static const Int_t nPtS = 7;
    double y1S [10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
    double vy1S [10] = {0.0177322,0.0148518,0.0097621,0.00484279,0.00296891,-0.000666419,-0.0057781,-0.00903913,-0.0132229,-0.0207987};
    double vey1S [10] = {0.00112423,0.000661336,0.000570038,0.000561605,0.000566271,0.000564669,0.000558678,0.000565402,0.000656118,0.00111938};
    double pt2S [nPtS] = {0.285,0.495,0.675,0.885,1.155,1.545,1.935};
    double vpt2S [nPtS] = {0.0264438,0.0453791,0.0613248,0.0752108,0.0886761,0.110444,0.110866};
    double vept2S [nPtS] = {0.000459202,0.000602041,0.00087003,0.00129944,0.00164998,0.00383577,0.00917259};


    // NA49 40-60%
//    static const Int_t nyN = 9;
//    double y1N [9] = {0.030000, 0.230000, 0.430000, 0.630000, 0.830000, 1.030000, 1.230000, 1.430000, 1.630000};
//    double vy1N [9] = {0.000573, -0.004717, -0.012663, -0.018926, -0.033673, -0.045617, -0.053355, -0.073629, -0.086166};
//    double vey1N [9] = {0.002462, 0.002437, 0.002501, 0.002700, 0.002936, 0.003021, 0.003147, 0.003398, 0.003793};
//    double pt1N [13] = {0.040000, 0.140000, 0.240000, 0.340000, 0.440000, 0.540000, 0.640000, 0.740000, 0.890000, 1.090000, 1.290000, 1.490000, 1.790000};
//    double vpt1N [13] = {-0.021820, -0.030928, -0.026938, -0.032392, -0.041378, -0.036585, -0.029139, -0.021672, -0.022864, 0.004461, -0.005396, 0.031558, 0.003868};
//    double vept1N [13] = {0.002614, 0.001889, 0.002207, 0.002577, 0.003042, 0.003682, 0.004526, 0.005563, 0.005500, 0.008517, 0.013315, 0.021167, 0.011700};
//    double y2N [9] = {0.030000, 0.230000, 0.430000, 0.630000, 0.830000, 1.030000, 1.230000, 1.430000, 1.630000};
//    double vy2N [9] = {0.061216, 0.043069, 0.041216, 0.025346, 0.031566, 0.034724, 0.045622, 0.009114, 0.017752};
//    double vey2N [9] = {0.007401, 0.007365, 0.007576, 0.007737, 0.008329, 0.008876, 0.009790, 0.011048, 0.012956};
//    double pt2N [13] = {0.040000, 0.140000, 0.240000, 0.340000, 0.440000, 0.540000, 0.640000, 0.740000, 0.890000, 1.090000, 1.290000, 1.490000, 1.790000};
//    double vpt2N [13] = {0.013888, 0.009422, 0.044918, 0.017056, 0.058014, 0.067450, 0.053726, 0.057203, 0.090791, 0.096066, 0.186829, 0.095682, 0.042780};
//    double vept2N [13] = {0.009278, 0.006000, 0.006881, 0.007520, 0.008602, 0.010205, 0.012495, 0.015477, 0.015341, 0.023959, 0.037999, 0.058463, 0.029478};
//
//    //STAR  7.7 GeV 40-80%
//    static const Int_t nPtS = 6;
//    double y1S [10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
//    double vy1S [10] = {0.051107,0.0363774,0.0267159,0.0184652,0.00620654,-0.0041068,-0.0167095,-0.0267019,-0.0376915,-0.0526559};
//    double vey1S [10] = {0.00227639,0.00134638,0.0011662,0.00115529,0.0011687,0.00116132,0.00114483,0.00115022,0.00132452,0.0022337};
//    double pt2S [nPtS] = {0.285,0.495,0.675,0.885,1.125,1.545};
//    double vpt2S [nPtS] = {0.0257936,0.0513509,0.0728293,0.0895198,0.110976,0.100925};
//    double vept2S [nPtS] = {0.00177321,0.00235618,0.00345922,0.00532486,0.00703612,0.0174535};


    double y1Nrefl [nyN * 2], vy1Nrefl [nyN * 2], vey1Nrefl [nyN * 2];
    double y2Nrefl [nyN * 2], vy2Nrefl [nyN * 2], vey2Nrefl [nyN * 2];

    for (Int_t j = 0; j < nyN; j++) {
        y1Nrefl [j] = (-1) * y1N [nyN - 1 - j];
        y2Nrefl [j] = (-1) * y2N [nyN - 1 - j];
        vy1Nrefl [j] = (-1) * vy1N [nyN - 1 - j];
        vy2Nrefl [j] = vy2N [nyN - 1 - j];
        vey1Nrefl [j] = vey1N [nyN - 1 - j];
        vey2Nrefl [j] = vey2N [nyN - 1 - j];
    }
    for (Int_t j = 0; j < nyN; j++) {
        y1Nrefl [j + nyN] = y1N [j];
        y2Nrefl [j + nyN] = y2N [j];
        vy1Nrefl [j + nyN] = vy1N [j];
        vy2Nrefl [j + nyN] = vy2N [j];
        vey1Nrefl [j + nyN] = vey1N [j];
        vey2Nrefl [j + nyN] = vey2N [j];
    }

    TGraphErrors *gV1yN = new TGraphErrors (9, y1N, vy1N, 0, vey1N);
    TGraphErrors *gV1yNrefl = new TGraphErrors (18, y1Nrefl, vy1Nrefl, 0, vey1Nrefl);
    TGraphErrors *gV1yS = new TGraphErrors (10, y1S, vy1S, 0, vey1S);
    TGraphErrors *gV1PtN = new TGraphErrors (13, pt1N, vpt1N, 0, vept1N);
    TGraphErrors *gV2yN = new TGraphErrors (9, y2N, vy2N, 0, vey2N);
    TGraphErrors *gV2yNrefl = new TGraphErrors (18, y2Nrefl, vy2Nrefl, 0, vey2Nrefl);
    TGraphErrors *gV2PtN = new TGraphErrors (13, pt2N, vpt2N, 0, vept2N);
    TGraphErrors *gV2PtS = new TGraphErrors (nPtS, pt2S, vpt2S, 0, vept2S);

    TGraphErrors *graphList [8];

    graphList [0] = gV1yN;
    graphList [1] = gV1PtN;
    graphList [2] = gV2yN;
    graphList [3] = gV2PtN;
    graphList [4] = gV1yNrefl;
    graphList [5] = gV2yNrefl;
    graphList [6] = gV1yS;
    graphList [7] = gV2PtS;

    for (Int_t j = 0; j < 6; j++) {
        graphList [j] -> SetMarkerStyle (28);
        graphList [j] -> SetMarkerSize (1.5);
        graphList [j] -> SetMarkerColor (4);
        graphList [j] -> SetLineColor (4);
    }

    for (Int_t j = 6; j < 8; j++) {
        graphList [j] -> SetMarkerStyle (30);
        graphList [j] -> SetMarkerSize (1.5);
        graphList [j] -> SetMarkerColor (1);
        graphList [j] -> SetLineColor (1);
    }

    gV1PtN -> Write ("gV1PtN");
    gV1yN -> Write ("gV1yN");
    gV1yNrefl -> Write ("gV1yNrefl");
    gV1yS -> Write ("gV1yS");
    gV2PtN -> Write ("gV2PtN");
    gV2PtS -> Write ("gV2PtS");
    gV2yN -> Write ("gV2yN");
    gV2yNrefl -> Write ("gV2yNrefl");
}

void CFlowReconstructor::GetFlowLoop (Int_t step) {
    Int_t nProfs = 7; // plot xY and yX flow
//    Int_t nProfs = 4; // do not plot xY and yX flow
    TString option [4] = {"recreate", "update", "update", "update"};
    TFile *corrFile = new TFile (histFileName_ + "_corr.root", "read");
    TFile *flowFile = new TFile (histFileName_ + "_flow.root", option [step]);
    if (step == 0) WritePreviousResults (flowFile); // for comparison with existing data
    TDirectory *stepDir = flowFile -> mkdir (dirName [step]);
    TDirectory *distrDir = stepDir -> mkdir ("Distributions");

    if (mhRangeForFlowSet_ == 0) {
        mhLow_ = mhMin_;
        mhHigh_ = mhMax_;
    }
    if (centRangeForFlowSet_ == 0) {
        centLow_ = centMin_;
        centHigh_ = centMax_;
    }
    if (harmonicFunctionSet) FillReferenceHist ();
	Int_t n, mh;
	Float_t R, Rerr, cent, pt, eta, bsIndex, shift, sign;

// test
    Float_t QQ [3][2][2]; // [AB, AC, BC][x, y][x, y]
	TFile *testFile = new TFile (histFileName_ + "QQ.root", "RECREATE");
    TTree *testTree = new TTree ("testTree", "Test tree");
	testTree -> Branch ("QQ", &QQ, "QQ[3][2][2]/F");
	testTree -> Branch ("Nsub", &nBinsBS_, 32000, 4);
	testTree -> Branch ("Cent", &cent, 32000, 4);
	testTree -> Branch ("Harm", &n, 32000, 4);
// test

    TH1 *hList1 [12], *hList2 [12], *hList3 [12], *hList4 [12];

	vector <TProfile*> pxXCent_RP, pxXMult_RP, pxXPt_RP, pxXEta_RP;

	vector <TProfile*> pxXaCent_SP, pxXbCent_SP, pxXcCent_SP, pyYaCent_SP, pyYbCent_SP, pyYcCent_SP;
	vector <TProfile*> pxYaCent_SP, pxYbCent_SP, pxYcCent_SP, pyXaCent_SP, pyXbCent_SP, pyXcCent_SP;
	vector <TProfile*> pqQaCent_SP, pqQbCent_SP, pqQcCent_SP;
	vector <TProfile*> pxXaCent_EP, pxXbCent_EP, pxXcCent_EP, pyYaCent_EP, pyYbCent_EP, pyYcCent_EP;
	vector <TProfile*> pxYaCent_EP, pxYbCent_EP, pxYcCent_EP, pyXaCent_EP, pyXbCent_EP, pyXcCent_EP;
	vector <TProfile*> pqQaCent_EP, pqQbCent_EP, pqQcCent_EP;

	vector <TProfile*> pxXaMult_SP, pxXbMult_SP, pxXcMult_SP, pyYaMult_SP, pyYbMult_SP, pyYcMult_SP;
	vector <TProfile*> pxYaMult_SP, pxYbMult_SP, pxYcMult_SP, pyXaMult_SP, pyXbMult_SP, pyXcMult_SP;
	vector <TProfile*> pqQaMult_SP, pqQbMult_SP, pqQcMult_SP;
	vector <TProfile*> pxXaMult_EP, pxXbMult_EP, pxXcMult_EP, pyYaMult_EP, pyYbMult_EP, pyYcMult_EP;
	vector <TProfile*> pxYaMult_EP, pxYbMult_EP, pxYcMult_EP, pyXaMult_EP, pyXbMult_EP, pyXcMult_EP;
	vector <TProfile*> pqQaMult_EP, pqQbMult_EP, pqQcMult_EP;

	vector <TProfile2D*> p2xXaCent_SP, p2xXbCent_SP, p2xXcCent_SP, p2yYaCent_SP, p2yYbCent_SP, p2yYcCent_SP;
	vector <TProfile2D*> p2xYaCent_SP, p2xYbCent_SP, p2xYcCent_SP, p2yXaCent_SP, p2yXbCent_SP, p2yXcCent_SP;
	vector <TProfile2D*> p2qQaCent_SP, p2qQbCent_SP, p2qQcCent_SP;
	vector <TProfile2D*> p2xXaCent_EP, p2xXbCent_EP, p2xXcCent_EP, p2yYaCent_EP, p2yYbCent_EP, p2yYcCent_EP;
	vector <TProfile2D*> p2xYaCent_EP, p2xYbCent_EP, p2xYcCent_EP, p2yXaCent_EP, p2yXbCent_EP, p2yXcCent_EP;
	vector <TProfile2D*> p2qQaCent_EP, p2qQbCent_EP, p2qQcCent_EP;

	vector <TProfile2D*> p2xXaMult_SP, p2xXbMult_SP, p2xXcMult_SP, p2yYaMult_SP, p2yYbMult_SP, p2yYcMult_SP;
	vector <TProfile2D*> p2xYaMult_SP, p2xYbMult_SP, p2xYcMult_SP, p2yXaMult_SP, p2yXbMult_SP, p2yXcMult_SP;
	vector <TProfile2D*> p2qQaMult_SP, p2qQbMult_SP, p2qQcMult_SP;
	vector <TProfile2D*> p2xXaMult_EP, p2xXbMult_EP, p2xXcMult_EP, p2yYaMult_EP, p2yYbMult_EP, p2yYcMult_EP;
	vector <TProfile2D*> p2xYaMult_EP, p2xYbMult_EP, p2xYcMult_EP, p2yXaMult_EP, p2yXbMult_EP, p2yXcMult_EP;
	vector <TProfile2D*> p2qQaMult_EP, p2qQbMult_EP, p2qQcMult_EP;

	vector <TProfile2D*> p2xXaPtCent_SP, p2xXbPtCent_SP, p2xXcPtCent_SP, p2yYaPtCent_SP, p2yYbPtCent_SP, p2yYcPtCent_SP;
	vector <TProfile2D*> p2xYaPtCent_SP, p2xYbPtCent_SP, p2xYcPtCent_SP, p2yXaPtCent_SP, p2yXbPtCent_SP, p2yXcPtCent_SP;
	vector <TProfile2D*> p2qQaPtCent_SP, p2qQbPtCent_SP, p2qQcPtCent_SP;
	vector <TProfile2D*> p2xXaPtCent_EP, p2xXbPtCent_EP, p2xXcPtCent_EP, p2yYaPtCent_EP, p2yYbPtCent_EP, p2yYcPtCent_EP;
	vector <TProfile2D*> p2xYaPtCent_EP, p2xYbPtCent_EP, p2xYcPtCent_EP, p2yXaPtCent_EP, p2yXbPtCent_EP, p2yXcPtCent_EP;
	vector <TProfile2D*> p2qQaPtCent_EP, p2qQbPtCent_EP, p2qQcPtCent_EP;

	vector <TProfile2D*> p2xXaPtMult_SP, p2xXbPtMult_SP, p2xXcPtMult_SP, p2yYaPtMult_SP, p2yYbPtMult_SP, p2yYcPtMult_SP;
	vector <TProfile2D*> p2xYaPtMult_SP, p2xYbPtMult_SP, p2xYcPtMult_SP, p2yXaPtMult_SP, p2yXbPtMult_SP, p2yXcPtMult_SP;
	vector <TProfile2D*> p2qQaPtMult_SP, p2qQbPtMult_SP, p2qQcPtMult_SP;
	vector <TProfile2D*> p2xXaPtMult_EP, p2xXbPtMult_EP, p2xXcPtMult_EP, p2yYaPtMult_EP, p2yYbPtMult_EP, p2yYcPtMult_EP;
	vector <TProfile2D*> p2xYaPtMult_EP, p2xYbPtMult_EP, p2xYcPtMult_EP, p2yXaPtMult_EP, p2yXbPtMult_EP, p2yXcPtMult_EP;
	vector <TProfile2D*> p2qQaPtMult_EP, p2qQbPtMult_EP, p2qQcPtMult_EP;

	vector <TProfile2D*> p2xXaEtaCent_SP, p2xXbEtaCent_SP, p2xXcEtaCent_SP, p2yYaEtaCent_SP, p2yYbEtaCent_SP, p2yYcEtaCent_SP;
	vector <TProfile2D*> p2xYaEtaCent_SP, p2xYbEtaCent_SP, p2xYcEtaCent_SP, p2yXaEtaCent_SP, p2yXbEtaCent_SP, p2yXcEtaCent_SP;
	vector <TProfile2D*> p2qQaEtaCent_SP, p2qQbEtaCent_SP, p2qQcEtaCent_SP;
	vector <TProfile2D*> p2xXaEtaCent_EP, p2xXbEtaCent_EP, p2xXcEtaCent_EP, p2yYaEtaCent_EP, p2yYbEtaCent_EP, p2yYcEtaCent_EP;
	vector <TProfile2D*> p2xYaEtaCent_EP, p2xYbEtaCent_EP, p2xYcEtaCent_EP, p2yXaEtaCent_EP, p2yXbEtaCent_EP, p2yXcEtaCent_EP;
	vector <TProfile2D*> p2qQaEtaCent_EP, p2qQbEtaCent_EP, p2qQcEtaCent_EP;

	vector <TProfile2D*> p2xXaEtaMult_SP, p2xXbEtaMult_SP, p2xXcEtaMult_SP, p2yYaEtaMult_SP, p2yYbEtaMult_SP, p2yYcEtaMult_SP;
	vector <TProfile2D*> p2xYaEtaMult_SP, p2xYbEtaMult_SP, p2xYcEtaMult_SP, p2yXaEtaMult_SP, p2yXbEtaMult_SP, p2yXcEtaMult_SP;
	vector <TProfile2D*> p2qQaEtaMult_SP, p2qQbEtaMult_SP, p2qQcEtaMult_SP;
	vector <TProfile2D*> p2xXaEtaMult_EP, p2xXbEtaMult_EP, p2xXcEtaMult_EP, p2yYaEtaMult_EP, p2yYbEtaMult_EP, p2yYcEtaMult_EP;
	vector <TProfile2D*> p2xYaEtaMult_EP, p2xYbEtaMult_EP, p2xYcEtaMult_EP, p2yXaEtaMult_EP, p2yXbEtaMult_EP, p2yXcEtaMult_EP;
	vector <TProfile2D*> p2qQaEtaMult_EP, p2qQbEtaMult_EP, p2qQcEtaMult_EP;

	vector <TProfile3D*> p3xXaPtCent_SP, p3xXbPtCent_SP, p3xXcPtCent_SP, p3yYaPtCent_SP, p3yYbPtCent_SP, p3yYcPtCent_SP;
	vector <TProfile3D*> p3xYaPtCent_SP, p3xYbPtCent_SP, p3xYcPtCent_SP, p3yXaPtCent_SP, p3yXbPtCent_SP, p3yXcPtCent_SP;
	vector <TProfile3D*> p3qQaPtCent_SP, p3qQbPtCent_SP, p3qQcPtCent_SP;
	vector <TProfile3D*> p3xXaPtCent_EP, p3xXbPtCent_EP, p3xXcPtCent_EP, p3yYaPtCent_EP, p3yYbPtCent_EP, p3yYcPtCent_EP;
	vector <TProfile3D*> p3xYaPtCent_EP, p3xYbPtCent_EP, p3xYcPtCent_EP, p3yXaPtCent_EP, p3yXbPtCent_EP, p3yXcPtCent_EP;
	vector <TProfile3D*> p3qQaPtCent_EP, p3qQbPtCent_EP, p3qQcPtCent_EP;

	vector <TProfile3D*> p3xXaPtMult_SP, p3xXbPtMult_SP, p3xXcPtMult_SP, p3yYaPtMult_SP, p3yYbPtMult_SP, p3yYcPtMult_SP;
	vector <TProfile3D*> p3xYaPtMult_SP, p3xYbPtMult_SP, p3xYcPtMult_SP, p3yXaPtMult_SP, p3yXbPtMult_SP, p3yXcPtMult_SP;
	vector <TProfile3D*> p3qQaPtMult_SP, p3qQbPtMult_SP, p3qQcPtMult_SP;
	vector <TProfile3D*> p3xXaPtMult_EP, p3xXbPtMult_EP, p3xXcPtMult_EP, p3yYaPtMult_EP, p3yYbPtMult_EP, p3yYcPtMult_EP;
	vector <TProfile3D*> p3xYaPtMult_EP, p3xYbPtMult_EP, p3xYcPtMult_EP, p3yXaPtMult_EP, p3yXbPtMult_EP, p3yXcPtMult_EP;
	vector <TProfile3D*> p3qQaPtMult_EP, p3qQbPtMult_EP, p3qQcPtMult_EP;

	vector <TProfile3D*> p3xXaEtaCent_SP, p3xXbEtaCent_SP, p3xXcEtaCent_SP, p3yYaEtaCent_SP, p3yYbEtaCent_SP, p3yYcEtaCent_SP;
	vector <TProfile3D*> p3xYaEtaCent_SP, p3xYbEtaCent_SP, p3xYcEtaCent_SP, p3yXaEtaCent_SP, p3yXbEtaCent_SP, p3yXcEtaCent_SP;
	vector <TProfile3D*> p3qQaEtaCent_SP, p3qQbEtaCent_SP, p3qQcEtaCent_SP;
	vector <TProfile3D*> p3xXaEtaCent_EP, p3xXbEtaCent_EP, p3xXcEtaCent_EP, p3yYaEtaCent_EP, p3yYbEtaCent_EP, p3yYcEtaCent_EP;
	vector <TProfile3D*> p3xYaEtaCent_EP, p3xYbEtaCent_EP, p3xYcEtaCent_EP, p3yXaEtaCent_EP, p3yXbEtaCent_EP, p3yXcEtaCent_EP;
	vector <TProfile3D*> p3qQaEtaCent_EP, p3qQbEtaCent_EP, p3qQcEtaCent_EP;

	vector <TProfile3D*> p3xXaEtaMult_SP, p3xXbEtaMult_SP, p3xXcEtaMult_SP, p3yYaEtaMult_SP, p3yYbEtaMult_SP, p3yYcEtaMult_SP;
	vector <TProfile3D*> p3xYaEtaMult_SP, p3xYbEtaMult_SP, p3xYcEtaMult_SP, p3yXaEtaMult_SP, p3yXbEtaMult_SP, p3yXcEtaMult_SP;
	vector <TProfile3D*> p3qQaEtaMult_SP, p3qQbEtaMult_SP, p3qQcEtaMult_SP;
	vector <TProfile3D*> p3xXaEtaMult_EP, p3xXbEtaMult_EP, p3xXcEtaMult_EP, p3yYaEtaMult_EP, p3yYbEtaMult_EP, p3yYcEtaMult_EP;
	vector <TProfile3D*> p3xYaEtaMult_EP, p3xYbEtaMult_EP, p3xYcEtaMult_EP, p3yXaEtaMult_EP, p3yXbEtaMult_EP, p3yXcEtaMult_EP;
	vector <TProfile3D*> p3qQaEtaMult_EP, p3qQbEtaMult_EP, p3qQcEtaMult_EP;

	vector <TProfile2D*> p2XaXbCent_SP, p2YaYbCent_SP, p2XaYbCent_SP, p2YaXbCent_SP, p2QaQbCent_SP;
	vector <TProfile2D*> p2XaXcCent_SP, p2YaYcCent_SP, p2XaYcCent_SP, p2YaXcCent_SP, p2QaQcCent_SP;
	vector <TProfile2D*> p2XbXcCent_SP, p2YbYcCent_SP, p2XbYcCent_SP, p2YbXcCent_SP, p2QbQcCent_SP;
	vector <TProfile2D*> p2XaXbCent_EP, p2YaYbCent_EP, p2XaYbCent_EP, p2YaXbCent_EP, p2QaQbCent_EP;
	vector <TProfile2D*> p2XaXcCent_EP, p2YaYcCent_EP, p2XaYcCent_EP, p2YaXcCent_EP, p2QaQcCent_EP;
	vector <TProfile2D*> p2XbXcCent_EP, p2YbYcCent_EP, p2XbYcCent_EP, p2YbXcCent_EP, p2QbQcCent_EP;

	vector <TProfile2D*> p2XaXbMult_SP, p2YaYbMult_SP, p2XaYbMult_SP, p2YaXbMult_SP, p2QaQbMult_SP;
	vector <TProfile2D*> p2XaXcMult_SP, p2YaYcMult_SP, p2XaYcMult_SP, p2YaXcMult_SP, p2QaQcMult_SP;
	vector <TProfile2D*> p2XbXcMult_SP, p2YbYcMult_SP, p2XbYcMult_SP, p2YbXcMult_SP, p2QbQcMult_SP;
	vector <TProfile2D*> p2XaXbMult_EP, p2YaYbMult_EP, p2XaYbMult_EP, p2YaXbMult_EP, p2QaQbMult_EP;
	vector <TProfile2D*> p2XaXcMult_EP, p2YaYcMult_EP, p2XaYcMult_EP, p2YaXcMult_EP, p2QaQcMult_EP;
	vector <TProfile2D*> p2XbXcMult_EP, p2YbYcMult_EP, p2XbYcMult_EP, p2YbXcMult_EP, p2QbQcMult_EP;

    vector <TH2F*> h2RxaCent_SP, h2RyaCent_SP, h2RaCent_SP, h2RxbCent_SP, h2RybCent_SP, h2RbCent_SP;
    vector <TH2F*> h2RxcCent_SP, h2RycCent_SP, h2RcCent_SP, h2RxaCent_EP, h2RyaCent_EP, h2RaCent_EP;
    vector <TH2F*> h2RxbCent_EP, h2RybCent_EP, h2RbCent_EP, h2RxcCent_EP, h2RycCent_EP, h2RcCent_EP;

    vector <TH2F*> h2RxaMult_SP, h2RyaMult_SP, h2RaMult_SP, h2RxbMult_SP, h2RybMult_SP, h2RbMult_SP;
    vector <TH2F*> h2RxcMult_SP, h2RycMult_SP, h2RcMult_SP, h2RxaMult_EP, h2RyaMult_EP, h2RaMult_EP;
    vector <TH2F*> h2RxbMult_EP, h2RybMult_EP, h2RbMult_EP, h2RxcMult_EP, h2RycMult_EP, h2RcMult_EP;

	vector <TProfile*> pXaXbCent_SP, pYaYbCent_SP, pXaYbCent_SP, pYaXbCent_SP, pQaQbCent_SP;
	vector <TProfile*> pXaXcCent_SP, pYaYcCent_SP, pXaYcCent_SP, pYaXcCent_SP, pQaQcCent_SP;
	vector <TProfile*> pXbXcCent_SP, pYbYcCent_SP, pXbYcCent_SP, pYbXcCent_SP, pQbQcCent_SP;
	vector <TProfile*> pXaXbCent_EP, pYaYbCent_EP, pXaYbCent_EP, pYaXbCent_EP, pQaQbCent_EP;
	vector <TProfile*> pXaXcCent_EP, pYaYcCent_EP, pXaYcCent_EP, pYaXcCent_EP, pQaQcCent_EP;
	vector <TProfile*> pXbXcCent_EP, pYbYcCent_EP, pXbYcCent_EP, pYbXcCent_EP, pQbQcCent_EP;

	vector <TProfile*> pXaXbMult_SP, pYaYbMult_SP, pXaYbMult_SP, pYaXbMult_SP, pQaQbMult_SP;
	vector <TProfile*> pXaXcMult_SP, pYaYcMult_SP, pXaYcMult_SP, pYaXcMult_SP, pQaQcMult_SP;
	vector <TProfile*> pXbXcMult_SP, pYbYcMult_SP, pXbYcMult_SP, pYbXcMult_SP, pQbQcMult_SP;
	vector <TProfile*> pXaXbMult_EP, pYaYbMult_EP, pXaYbMult_EP, pYaXbMult_EP, pQaQbMult_EP;
	vector <TProfile*> pXaXcMult_EP, pYaYcMult_EP, pXaYcMult_EP, pYaXcMult_EP, pQaQcMult_EP;
	vector <TProfile*> pXbXcMult_EP, pYbYcMult_EP, pXbYcMult_EP, pYbXcMult_EP, pQbQcMult_EP;

    vector <TH1F*> hRxaCent_SP, hRyaCent_SP, hRaCent_SP, hRxbCent_SP, hRybCent_SP, hRbCent_SP;
    vector <TH1F*> hRxcCent_SP, hRycCent_SP, hRcCent_SP;
    vector <TH1F*> hRxaCent_EP, hRyaCent_EP, hRaCent_EP, hRxbCent_EP, hRybCent_EP, hRbCent_EP;
    vector <TH1F*> hRxcCent_EP, hRycCent_EP, hRcCent_EP;

    vector <TH1F*> hRxaMult_SP, hRyaMult_SP, hRaMult_SP, hRxbMult_SP, hRybMult_SP, hRbMult_SP;
    vector <TH1F*> hRxcMult_SP, hRycMult_SP, hRcMult_SP;
    vector <TH1F*> hRxaMult_EP, hRyaMult_EP, hRaMult_EP, hRxbMult_EP, hRybMult_EP, hRbMult_EP;
    vector <TH1F*> hRxcMult_EP, hRycMult_EP, hRcMult_EP;

	vector <TProfile*> pXaXbCentBS_SP, pYaYbCentBS_SP, pXaYbCentBS_SP, pYaXbCentBS_SP, pQaQbCentBS_SP;
	vector <TProfile*> pXaXcCentBS_SP, pYaYcCentBS_SP, pXaYcCentBS_SP, pYaXcCentBS_SP, pQaQcCentBS_SP;
	vector <TProfile*> pXbXcCentBS_SP, pYbYcCentBS_SP, pXbYcCentBS_SP, pYbXcCentBS_SP, pQbQcCentBS_SP;
	vector <TProfile*> pXaXbCentBS_EP, pYaYbCentBS_EP, pXaYbCentBS_EP, pYaXbCentBS_EP, pQaQbCentBS_EP;
	vector <TProfile*> pXaXcCentBS_EP, pYaYcCentBS_EP, pXaYcCentBS_EP, pYaXcCentBS_EP, pQaQcCentBS_EP;
	vector <TProfile*> pXbXcCentBS_EP, pYbYcCentBS_EP, pXbYcCentBS_EP, pYbXcCentBS_EP, pQbQcCentBS_EP;

	vector <TProfile*> pXaXbMultBS_SP, pYaYbMultBS_SP, pXaYbMultBS_SP, pYaXbMultBS_SP, pQaQbMultBS_SP;
	vector <TProfile*> pXaXcMultBS_SP, pYaYcMultBS_SP, pXaYcMultBS_SP, pYaXcMultBS_SP, pQaQcMultBS_SP;
	vector <TProfile*> pXbXcMultBS_SP, pYbYcMultBS_SP, pXbYcMultBS_SP, pYbXcMultBS_SP, pQbQcMultBS_SP;
	vector <TProfile*> pXaXbMultBS_EP, pYaYbMultBS_EP, pXaYbMultBS_EP, pYaXbMultBS_EP, pQaQbMultBS_EP;
	vector <TProfile*> pXaXcMultBS_EP, pYaYcMultBS_EP, pXaYcMultBS_EP, pYaXcMultBS_EP, pQaQcMultBS_EP;
	vector <TProfile*> pXbXcMultBS_EP, pYbYcMultBS_EP, pXbYcMultBS_EP, pYbXcMultBS_EP, pQbQcMultBS_EP;

    vector <TH1F*> hRxaCentBS_SP, hRyaCentBS_SP, hRaCentBS_SP, hRxbCentBS_SP, hRybCentBS_SP, hRbCentBS_SP;
    vector <TH1F*> hRxcCentBS_SP, hRycCentBS_SP, hRcCentBS_SP, hRxaCentBS_EP, hRyaCentBS_EP, hRaCentBS_EP;
    vector <TH1F*> hRxbCentBS_EP, hRybCentBS_EP, hRbCentBS_EP, hRxcCentBS_EP, hRycCentBS_EP, hRcCentBS_EP;

    vector <TH1F*> hRxaMultBS_SP, hRyaMultBS_SP, hRaMultBS_SP, hRxbMultBS_SP, hRybMultBS_SP, hRbMultBS_SP;
    vector <TH1F*> hRxcMultBS_SP, hRycMultBS_SP, hRcMultBS_SP;
    vector <TH1F*> hRxaMultBS_EP, hRyaMultBS_EP, hRaMultBS_EP, hRxbMultBS_EP, hRybMultBS_EP, hRbMultBS_EP;
    vector <TH1F*> hRxcMultBS_EP, hRycMultBS_EP, hRcMultBS_EP;

    vector <TProfile2D*> p2VxaCent_SP, p2VyaCent_SP, p2VaCent_SP, p2VxYaCent_SP, p2VyXaCent_SP;
    vector <TProfile2D*> p2VxbCent_SP, p2VybCent_SP, p2VbCent_SP, p2VxYbCent_SP, p2VyXbCent_SP;
    vector <TProfile2D*> p2VxcCent_SP, p2VycCent_SP, p2VcCent_SP, p2VxYcCent_SP, p2VyXcCent_SP;
    vector <TProfile2D*> p2VxCent_SP, p2VyCent_SP, p2VCent_SP, p2VxYCent_SP, p2VyXCent_SP;
    vector <TProfile2D*> p2VxaCent_EP, p2VyaCent_EP, p2VaCent_EP, p2VxYaCent_EP, p2VyXaCent_EP;
    vector <TProfile2D*> p2VxbCent_EP, p2VybCent_EP, p2VbCent_EP, p2VxYbCent_EP, p2VyXbCent_EP;
    vector <TProfile2D*> p2VxcCent_EP, p2VycCent_EP, p2VcCent_EP, p2VxYcCent_EP, p2VyXcCent_EP;
    vector <TProfile2D*> p2VxCent_EP, p2VyCent_EP, p2VCent_EP, p2VxYCent_EP, p2VyXCent_EP;

    vector <TProfile2D*> p2VxaMult_SP, p2VyaMult_SP, p2VaMult_SP, p2VxYaMult_SP, p2VyXaMult_SP;
    vector <TProfile2D*> p2VxbMult_SP, p2VybMult_SP, p2VbMult_SP, p2VxYbMult_SP, p2VyXbMult_SP;
    vector <TProfile2D*> p2VxcMult_SP, p2VycMult_SP, p2VcMult_SP, p2VxYcMult_SP, p2VyXcMult_SP;
    vector <TProfile2D*> p2VxMult_SP, p2VyMult_SP, p2VMult_SP, p2VxYMult_SP, p2VyXMult_SP;
    vector <TProfile2D*> p2VxaMult_EP, p2VyaMult_EP, p2VaMult_EP, p2VxYaMult_EP, p2VyXaMult_EP;
    vector <TProfile2D*> p2VxbMult_EP, p2VybMult_EP, p2VbMult_EP, p2VxYbMult_EP, p2VyXbMult_EP;
    vector <TProfile2D*> p2VxcMult_EP, p2VycMult_EP, p2VcMult_EP, p2VxYcMult_EP, p2VyXcMult_EP;
    vector <TProfile2D*> p2VxMult_EP, p2VyMult_EP, p2VMult_EP, p2VxYMult_EP, p2VyXMult_EP;

    vector <TProfile2D*> p2VxaPtCent_SP, p2VyaPtCent_SP, p2VaPtCent_SP, p2VxYaPtCent_SP, p2VyXaPtCent_SP;
    vector <TProfile2D*> p2VxbPtCent_SP, p2VybPtCent_SP, p2VbPtCent_SP, p2VxYbPtCent_SP, p2VyXbPtCent_SP;
    vector <TProfile2D*> p2VxcPtCent_SP, p2VycPtCent_SP, p2VcPtCent_SP, p2VxYcPtCent_SP, p2VyXcPtCent_SP;
    vector <TProfile2D*> p2VxPtCent_SP, p2VyPtCent_SP, p2VPtCent_SP, p2VxYPtCent_SP, p2VyXPtCent_SP;
    vector <TProfile2D*> p2VxaPtCent_EP, p2VyaPtCent_EP, p2VaPtCent_EP, p2VxYaPtCent_EP, p2VyXaPtCent_EP;
    vector <TProfile2D*> p2VxbPtCent_EP, p2VybPtCent_EP, p2VbPtCent_EP, p2VxYbPtCent_EP, p2VyXbPtCent_EP;
    vector <TProfile2D*> p2VxcPtCent_EP, p2VycPtCent_EP, p2VcPtCent_EP, p2VxYcPtCent_EP, p2VyXcPtCent_EP;
    vector <TProfile2D*> p2VxPtCent_EP, p2VyPtCent_EP, p2VPtCent_EP, p2VxYPtCent_EP, p2VyXPtCent_EP;

    vector <TProfile2D*> p2VxaEtaCent_SP, p2VyaEtaCent_SP, p2VaEtaCent_SP, p2VxYaEtaCent_SP, p2VyXaEtaCent_SP;
    vector <TProfile2D*> p2VxbEtaCent_SP, p2VybEtaCent_SP, p2VbEtaCent_SP, p2VxYbEtaCent_SP, p2VyXbEtaCent_SP;
    vector <TProfile2D*> p2VxcEtaCent_SP, p2VycEtaCent_SP, p2VcEtaCent_SP, p2VxYcEtaCent_SP, p2VyXcEtaCent_SP;
    vector <TProfile2D*> p2VxEtaCent_SP, p2VyEtaCent_SP, p2VEtaCent_SP, p2VxYEtaCent_SP, p2VyXEtaCent_SP;
    vector <TProfile2D*> p2VxaEtaCent_EP, p2VyaEtaCent_EP, p2VaEtaCent_EP, p2VxYaEtaCent_EP, p2VyXaEtaCent_EP;
    vector <TProfile2D*> p2VxbEtaCent_EP, p2VybEtaCent_EP, p2VbEtaCent_EP, p2VxYbEtaCent_EP, p2VyXbEtaCent_EP;
    vector <TProfile2D*> p2VxcEtaCent_EP, p2VycEtaCent_EP, p2VcEtaCent_EP, p2VxYcEtaCent_EP, p2VyXcEtaCent_EP;
    vector <TProfile2D*> p2VxEtaCent_EP, p2VyEtaCent_EP, p2VEtaCent_EP, p2VxYEtaCent_EP, p2VyXEtaCent_EP;

    vector <TProfile2D*> p2VxaPtMult_SP, p2VyaPtMult_SP, p2VaPtMult_SP, p2VxYaPtMult_SP, p2VyXaPtMult_SP;
    vector <TProfile2D*> p2VxbPtMult_SP, p2VybPtMult_SP, p2VbPtMult_SP, p2VxYbPtMult_SP, p2VyXbPtMult_SP;
    vector <TProfile2D*> p2VxcPtMult_SP, p2VycPtMult_SP, p2VcPtMult_SP, p2VxYcPtMult_SP, p2VyXcPtMult_SP;
    vector <TProfile2D*> p2VxPtMult_SP, p2VyPtMult_SP, p2VPtMult_SP, p2VxYPtMult_SP, p2VyXPtMult_SP;
    vector <TProfile2D*> p2VxaPtMult_EP, p2VyaPtMult_EP, p2VaPtMult_EP, p2VxYaPtMult_EP, p2VyXaPtMult_EP;
    vector <TProfile2D*> p2VxbPtMult_EP, p2VybPtMult_EP, p2VbPtMult_EP, p2VxYbPtMult_EP, p2VyXbPtMult_EP;
    vector <TProfile2D*> p2VxcPtMult_EP, p2VycPtMult_EP, p2VcPtMult_EP, p2VxYcPtMult_EP, p2VyXcPtMult_EP;
    vector <TProfile2D*> p2VxPtMult_EP, p2VyPtMult_EP, p2VPtMult_EP, p2VxYPtMult_EP, p2VyXPtMult_EP;

    vector <TProfile2D*> p2VxaEtaMult_SP, p2VyaEtaMult_SP, p2VaEtaMult_SP, p2VxYaEtaMult_SP, p2VyXaEtaMult_SP;
    vector <TProfile2D*> p2VxbEtaMult_SP, p2VybEtaMult_SP, p2VbEtaMult_SP, p2VxYbEtaMult_SP, p2VyXbEtaMult_SP;
    vector <TProfile2D*> p2VxcEtaMult_SP, p2VycEtaMult_SP, p2VcEtaMult_SP, p2VxYcEtaMult_SP, p2VyXcEtaMult_SP;
    vector <TProfile2D*> p2VxEtaMult_SP, p2VyEtaMult_SP, p2VEtaMult_SP, p2VxYEtaMult_SP, p2VyXEtaMult_SP;
    vector <TProfile2D*> p2VxaEtaMult_EP, p2VyaEtaMult_EP, p2VaEtaMult_EP, p2VxYaEtaMult_EP, p2VyXaEtaMult_EP;
    vector <TProfile2D*> p2VxbEtaMult_EP, p2VybEtaMult_EP, p2VbEtaMult_EP, p2VxYbEtaMult_EP, p2VyXbEtaMult_EP;
    vector <TProfile2D*> p2VxcEtaMult_EP, p2VycEtaMult_EP, p2VcEtaMult_EP, p2VxYcEtaMult_EP, p2VyXcEtaMult_EP;
    vector <TProfile2D*> p2VxEtaMult_EP, p2VyEtaMult_EP, p2VEtaMult_EP, p2VxYEtaMult_EP, p2VyXEtaMult_EP;

    vector <TH1F*> hVxaCent_SP, hVyaCent_SP, hVaCent_SP, hVxYaCent_SP, hVyXaCent_SP;
    vector <TH1F*> hVxbCent_SP, hVybCent_SP, hVbCent_SP, hVxYbCent_SP, hVyXbCent_SP;
    vector <TH1F*> hVxcCent_SP, hVycCent_SP, hVcCent_SP, hVxYcCent_SP, hVyXcCent_SP;
    vector <TH1F*> hVxCent_SP, hVyCent_SP, hVCent_SP, hVxYCent_SP, hVyXCent_SP;
    vector <TH1F*> hVxaCent_EP, hVyaCent_EP, hVaCent_EP, hVxYaCent_EP, hVyXaCent_EP;
    vector <TH1F*> hVxbCent_EP, hVybCent_EP, hVbCent_EP, hVxYbCent_EP, hVyXbCent_EP;
    vector <TH1F*> hVxcCent_EP, hVycCent_EP, hVcCent_EP, hVxYcCent_EP, hVyXcCent_EP;
    vector <TH1F*> hVxCent_EP, hVyCent_EP, hVCent_EP, hVxYCent_EP, hVyXCent_EP;

    vector <TH1F*> hVxaMult_SP, hVyaMult_SP, hVaMult_SP, hVxYaMult_SP, hVyXaMult_SP;
    vector <TH1F*> hVxbMult_SP, hVybMult_SP, hVbMult_SP, hVxYbMult_SP, hVyXbMult_SP;
    vector <TH1F*> hVxcMult_SP, hVycMult_SP, hVcMult_SP, hVxYcMult_SP, hVyXcMult_SP;
    vector <TH1F*> hVxMult_SP, hVyMult_SP, hVMult_SP, hVxYMult_SP, hVyXMult_SP;
    vector <TH1F*> hVxaMult_EP, hVyaMult_EP, hVaMult_EP, hVxYaMult_EP, hVyXaMult_EP;
    vector <TH1F*> hVxbMult_EP, hVybMult_EP, hVbMult_EP, hVxYbMult_EP, hVyXbMult_EP;
    vector <TH1F*> hVxcMult_EP, hVycMult_EP, hVcMult_EP, hVxYcMult_EP, hVyXcMult_EP;
    vector <TH1F*> hVxMult_EP, hVyMult_EP, hVMult_EP, hVxYMult_EP, hVyXMult_EP;

    vector <TH1F*> hVxaCentBS_SP, hVyaCentBS_SP, hVaCentBS_SP, hVxYaCentBS_SP, hVyXaCentBS_SP;
    vector <TH1F*> hVxbCentBS_SP, hVybCentBS_SP, hVbCentBS_SP, hVxYbCentBS_SP, hVyXbCentBS_SP;
    vector <TH1F*> hVxcCentBS_SP, hVycCentBS_SP, hVcCentBS_SP, hVxYcCentBS_SP, hVyXcCentBS_SP;
    vector <TH1F*> hVxCentBS_SP, hVyCentBS_SP, hVCentBS_SP, hVxYCentBS_SP, hVyXCentBS_SP;
    vector <TH1F*> hVxaCentBS_EP, hVyaCentBS_EP, hVaCentBS_EP, hVxYaCentBS_EP, hVyXaCentBS_EP;
    vector <TH1F*> hVxbCentBS_EP, hVybCentBS_EP, hVbCentBS_EP, hVxYbCentBS_EP, hVyXbCentBS_EP;
    vector <TH1F*> hVxcCentBS_EP, hVycCentBS_EP, hVcCentBS_EP, hVxYcCentBS_EP, hVyXcCentBS_EP;
    vector <TH1F*> hVxCentBS_EP, hVyCentBS_EP, hVCentBS_EP, hVxYCentBS_EP, hVyXCentBS_EP;

    vector <TH1F*> hVxaMultBS_SP, hVyaMultBS_SP, hVaMultBS_SP, hVxYaMultBS_SP, hVyXaMultBS_SP;
    vector <TH1F*> hVxbMultBS_SP, hVybMultBS_SP, hVbMultBS_SP, hVxYbMultBS_SP, hVyXbMultBS_SP;
    vector <TH1F*> hVxcMultBS_SP, hVycMultBS_SP, hVcMultBS_SP, hVxYcMultBS_SP, hVyXcMultBS_SP;
    vector <TH1F*> hVxMultBS_SP, hVyMultBS_SP, hVMultBS_SP, hVxYMultBS_SP, hVyXMultBS_SP;
    vector <TH1F*> hVxaMultBS_EP, hVyaMultBS_EP, hVaMultBS_EP, hVxYaMultBS_EP, hVyXaMultBS_EP;
    vector <TH1F*> hVxbMultBS_EP, hVybMultBS_EP, hVbMultBS_EP, hVxYbMultBS_EP, hVyXbMultBS_EP;
    vector <TH1F*> hVxcMultBS_EP, hVycMultBS_EP, hVcMultBS_EP, hVxYcMultBS_EP, hVyXcMultBS_EP;
    vector <TH1F*> hVxMultBS_EP, hVyMultBS_EP, hVMultBS_EP, hVxYMultBS_EP, hVyXMultBS_EP;

    vector <TH1F*> hVxaPtCent_SP, hVyaPtCent_SP, hVaPtCent_SP, hVxYaPtCent_SP, hVyXaPtCent_SP;
    vector <TH1F*> hVxbPtCent_SP, hVybPtCent_SP, hVbPtCent_SP, hVxYbPtCent_SP, hVyXbPtCent_SP;
    vector <TH1F*> hVxcPtCent_SP, hVycPtCent_SP, hVcPtCent_SP, hVxYcPtCent_SP, hVyXcPtCent_SP;
    vector <TH1F*> hVxPtCent_SP, hVyPtCent_SP, hVPtCent_SP, hVxYPtCent_SP, hVyXPtCent_SP;
    vector <TH1F*> hVxaPtCent_EP, hVyaPtCent_EP, hVaPtCent_EP, hVxYaPtCent_EP, hVyXaPtCent_EP;
    vector <TH1F*> hVxbPtCent_EP, hVybPtCent_EP, hVbPtCent_EP, hVxYbPtCent_EP, hVyXbPtCent_EP;
    vector <TH1F*> hVxcPtCent_EP, hVycPtCent_EP, hVcPtCent_EP, hVxYcPtCent_EP, hVyXcPtCent_EP;
    vector <TH1F*> hVxPtCent_EP, hVyPtCent_EP, hVPtCent_EP, hVxYPtCent_EP, hVyXPtCent_EP;

    vector <TH1F*> hVxaPtMult_SP, hVyaPtMult_SP, hVaPtMult_SP, hVxYaPtMult_SP, hVyXaPtMult_SP;
    vector <TH1F*> hVxbPtMult_SP, hVybPtMult_SP, hVbPtMult_SP, hVxYbPtMult_SP, hVyXbPtMult_SP;
    vector <TH1F*> hVxcPtMult_SP, hVycPtMult_SP, hVcPtMult_SP, hVxYcPtMult_SP, hVyXcPtMult_SP;
    vector <TH1F*> hVxPtMult_SP, hVyPtMult_SP, hVPtMult_SP, hVxYPtMult_SP, hVyXPtMult_SP;
    vector <TH1F*> hVxaPtMult_EP, hVyaPtMult_EP, hVaPtMult_EP, hVxYaPtMult_EP, hVyXaPtMult_EP;
    vector <TH1F*> hVxbPtMult_EP, hVybPtMult_EP, hVbPtMult_EP, hVxYbPtMult_EP, hVyXbPtMult_EP;
    vector <TH1F*> hVxcPtMult_EP, hVycPtMult_EP, hVcPtMult_EP, hVxYcPtMult_EP, hVyXcPtMult_EP;
    vector <TH1F*> hVxPtMult_EP, hVyPtMult_EP, hVPtMult_EP, hVxYPtMult_EP, hVyXPtMult_EP;

    vector <TH1F*> hVxaPtCentBS_SP, hVyaPtCentBS_SP, hVaPtCentBS_SP, hVxYaPtCentBS_SP, hVyXaPtCentBS_SP;
    vector <TH1F*> hVxbPtCentBS_SP, hVybPtCentBS_SP, hVbPtCentBS_SP, hVxYbPtCentBS_SP, hVyXbPtCentBS_SP;
    vector <TH1F*> hVxcPtCentBS_SP, hVycPtCentBS_SP, hVcPtCentBS_SP, hVxYcPtCentBS_SP, hVyXcPtCentBS_SP;
    vector <TH1F*> hVxPtCentBS_SP, hVyPtCentBS_SP, hVPtCentBS_SP, hVxYPtCentBS_SP, hVyXPtCentBS_SP;
    vector <TH1F*> hVxaPtCentBS_EP, hVyaPtCentBS_EP, hVaPtCentBS_EP, hVxYaPtCentBS_EP, hVyXaPtCentBS_EP;
    vector <TH1F*> hVxbPtCentBS_EP, hVybPtCentBS_EP, hVbPtCentBS_EP, hVxYbPtCentBS_EP, hVyXbPtCentBS_EP;
    vector <TH1F*> hVxcPtCentBS_EP, hVycPtCentBS_EP, hVcPtCentBS_EP, hVxYcPtCentBS_EP, hVyXcPtCentBS_EP;
    vector <TH1F*> hVxPtCentBS_EP, hVyPtCentBS_EP, hVPtCentBS_EP, hVxYPtCentBS_EP, hVyXPtCentBS_EP;

    vector <TH1F*> hVxaPtMultBS_SP, hVyaPtMultBS_SP, hVaPtMultBS_SP, hVxYaPtMultBS_SP, hVyXaPtMultBS_SP;
    vector <TH1F*> hVxbPtMultBS_SP, hVybPtMultBS_SP, hVbPtMultBS_SP, hVxYbPtMultBS_SP, hVyXbPtMultBS_SP;
    vector <TH1F*> hVxcPtMultBS_SP, hVycPtMultBS_SP, hVcPtMultBS_SP, hVxYcPtMultBS_SP, hVyXcPtMultBS_SP;
    vector <TH1F*> hVxPtMultBS_SP, hVyPtMultBS_SP, hVPtMultBS_SP, hVxYPtMultBS_SP, hVyXPtMultBS_SP;
    vector <TH1F*> hVxaPtMultBS_EP, hVyaPtMultBS_EP, hVaPtMultBS_EP, hVxYaPtMultBS_EP, hVyXaPtMultBS_EP;
    vector <TH1F*> hVxbPtMultBS_EP, hVybPtMultBS_EP, hVbPtMultBS_EP, hVxYbPtMultBS_EP, hVyXbPtMultBS_EP;
    vector <TH1F*> hVxcPtMultBS_EP, hVycPtMultBS_EP, hVcPtMultBS_EP, hVxYcPtMultBS_EP, hVyXcPtMultBS_EP;
    vector <TH1F*> hVxPtMultBS_EP, hVyPtMultBS_EP, hVPtMultBS_EP, hVxYPtMultBS_EP, hVyXPtMultBS_EP;

    vector <TH1F*> hVxaEtaCent_SP, hVyaEtaCent_SP, hVaEtaCent_SP, hVxYaEtaCent_SP, hVyXaEtaCent_SP;
    vector <TH1F*> hVxbEtaCent_SP, hVybEtaCent_SP, hVbEtaCent_SP, hVxYbEtaCent_SP, hVyXbEtaCent_SP;
    vector <TH1F*> hVxcEtaCent_SP, hVycEtaCent_SP, hVcEtaCent_SP, hVxYcEtaCent_SP, hVyXcEtaCent_SP;
    vector <TH1F*> hVxEtaCent_SP, hVyEtaCent_SP, hVEtaCent_SP, hVxYEtaCent_SP, hVyXEtaCent_SP;
    vector <TH1F*> hVxaEtaCent_EP, hVyaEtaCent_EP, hVaEtaCent_EP, hVxYaEtaCent_EP, hVyXaEtaCent_EP;
    vector <TH1F*> hVxbEtaCent_EP, hVybEtaCent_EP, hVbEtaCent_EP, hVxYbEtaCent_EP, hVyXbEtaCent_EP;
    vector <TH1F*> hVxcEtaCent_EP, hVycEtaCent_EP, hVcEtaCent_EP, hVxYcEtaCent_EP, hVyXcEtaCent_EP;
    vector <TH1F*> hVxEtaCent_EP, hVyEtaCent_EP, hVEtaCent_EP, hVxYEtaCent_EP, hVyXEtaCent_EP;

    vector <TH1F*> hVxaEtaMult_SP, hVyaEtaMult_SP, hVaEtaMult_SP, hVxYaEtaMult_SP, hVyXaEtaMult_SP;
    vector <TH1F*> hVxbEtaMult_SP, hVybEtaMult_SP, hVbEtaMult_SP, hVxYbEtaMult_SP, hVyXbEtaMult_SP;
    vector <TH1F*> hVxcEtaMult_SP, hVycEtaMult_SP, hVcEtaMult_SP, hVxYcEtaMult_SP, hVyXcEtaMult_SP;
    vector <TH1F*> hVxEtaMult_SP, hVyEtaMult_SP, hVEtaMult_SP, hVxYEtaMult_SP, hVyXEtaMult_SP;
    vector <TH1F*> hVxaEtaMult_EP, hVyaEtaMult_EP, hVaEtaMult_EP, hVxYaEtaMult_EP, hVyXaEtaMult_EP;
    vector <TH1F*> hVxbEtaMult_EP, hVybEtaMult_EP, hVbEtaMult_EP, hVxYbEtaMult_EP, hVyXbEtaMult_EP;
    vector <TH1F*> hVxcEtaMult_EP, hVycEtaMult_EP, hVcEtaMult_EP, hVxYcEtaMult_EP, hVyXcEtaMult_EP;
    vector <TH1F*> hVxEtaMult_EP, hVyEtaMult_EP, hVEtaMult_EP, hVxYEtaMult_EP, hVyXEtaMult_EP;

    vector <TH1F*> hVxaEtaCentBS_SP, hVyaEtaCentBS_SP, hVaEtaCentBS_SP, hVxYaEtaCentBS_SP, hVyXaEtaCentBS_SP;
    vector <TH1F*> hVxbEtaCentBS_SP, hVybEtaCentBS_SP, hVbEtaCentBS_SP, hVxYbEtaCentBS_SP, hVyXbEtaCentBS_SP;
    vector <TH1F*> hVxcEtaCentBS_SP, hVycEtaCentBS_SP, hVcEtaCentBS_SP, hVxYcEtaCentBS_SP, hVyXcEtaCentBS_SP;
    vector <TH1F*> hVxEtaCentBS_SP, hVyEtaCentBS_SP, hVEtaCentBS_SP, hVxYEtaCentBS_SP, hVyXEtaCentBS_SP;
    vector <TH1F*> hVxaEtaCentBS_EP, hVyaEtaCentBS_EP, hVaEtaCentBS_EP, hVxYaEtaCentBS_EP, hVyXaEtaCentBS_EP;
    vector <TH1F*> hVxbEtaCentBS_EP, hVybEtaCentBS_EP, hVbEtaCentBS_EP, hVxYbEtaCentBS_EP, hVyXbEtaCentBS_EP;
    vector <TH1F*> hVxcEtaCentBS_EP, hVycEtaCentBS_EP, hVcEtaCentBS_EP, hVxYcEtaCentBS_EP, hVyXcEtaCentBS_EP;
    vector <TH1F*> hVxEtaCentBS_EP, hVyEtaCentBS_EP, hVEtaCentBS_EP, hVxYEtaCentBS_EP, hVyXEtaCentBS_EP;

    vector <TH1F*> hVxaEtaMultBS_SP, hVyaEtaMultBS_SP, hVaEtaMultBS_SP, hVxYaEtaMultBS_SP, hVyXaEtaMultBS_SP;
    vector <TH1F*> hVxbEtaMultBS_SP, hVybEtaMultBS_SP, hVbEtaMultBS_SP, hVxYbEtaMultBS_SP, hVyXbEtaMultBS_SP;
    vector <TH1F*> hVxcEtaMultBS_SP, hVycEtaMultBS_SP, hVcEtaMultBS_SP, hVxYcEtaMultBS_SP, hVyXcEtaMultBS_SP;
    vector <TH1F*> hVxEtaMultBS_SP, hVyEtaMultBS_SP, hVEtaMultBS_SP, hVxYEtaMultBS_SP, hVyXEtaMultBS_SP;
    vector <TH1F*> hVxaEtaMultBS_EP, hVyaEtaMultBS_EP, hVaEtaMultBS_EP, hVxYaEtaMultBS_EP, hVyXaEtaMultBS_EP;
    vector <TH1F*> hVxbEtaMultBS_EP, hVybEtaMultBS_EP, hVbEtaMultBS_EP, hVxYbEtaMultBS_EP, hVyXbEtaMultBS_EP;
    vector <TH1F*> hVxcEtaMultBS_EP, hVycEtaMultBS_EP, hVcEtaMultBS_EP, hVxYcEtaMultBS_EP, hVyXcEtaMultBS_EP;
    vector <TH1F*> hVxEtaMultBS_EP, hVyEtaMultBS_EP, hVEtaMultBS_EP, hVxYEtaMultBS_EP, hVyXEtaMultBS_EP;

    vector <TH1F*> hVxaEtaReflCent_SP, hVyaEtaReflCent_SP, hVaEtaReflCent_SP, hVxYaEtaReflCent_SP, hVyXaEtaReflCent_SP;
    vector <TH1F*> hVxbEtaReflCent_SP, hVybEtaReflCent_SP, hVbEtaReflCent_SP, hVxYbEtaReflCent_SP, hVyXbEtaReflCent_SP;
    vector <TH1F*> hVxcEtaReflCent_SP, hVycEtaReflCent_SP, hVcEtaReflCent_SP, hVxYcEtaReflCent_SP, hVyXcEtaReflCent_SP;
    vector <TH1F*> hVxEtaReflCent_SP, hVyEtaReflCent_SP, hVEtaReflCent_SP, hVxYEtaReflCent_SP, hVyXEtaReflCent_SP;
    vector <TH1F*> hVxaEtaReflCent_EP, hVyaEtaReflCent_EP, hVaEtaReflCent_EP, hVxYaEtaReflCent_EP, hVyXaEtaReflCent_EP;
    vector <TH1F*> hVxbEtaReflCent_EP, hVybEtaReflCent_EP, hVbEtaReflCent_EP, hVxYbEtaReflCent_EP, hVyXbEtaReflCent_EP;
    vector <TH1F*> hVxcEtaReflCent_EP, hVycEtaReflCent_EP, hVcEtaReflCent_EP, hVxYcEtaReflCent_EP, hVyXcEtaReflCent_EP;
    vector <TH1F*> hVxEtaReflCent_EP, hVyEtaReflCent_EP, hVEtaReflCent_EP, hVxYEtaReflCent_EP, hVyXEtaReflCent_EP;

    vector <TH1F*> hVxaEtaReflMult_SP, hVyaEtaReflMult_SP, hVaEtaReflMult_SP, hVxYaEtaReflMult_SP, hVyXaEtaReflMult_SP;
    vector <TH1F*> hVxbEtaReflMult_SP, hVybEtaReflMult_SP, hVbEtaReflMult_SP, hVxYbEtaReflMult_SP, hVyXbEtaReflMult_SP;
    vector <TH1F*> hVxcEtaReflMult_SP, hVycEtaReflMult_SP, hVcEtaReflMult_SP, hVxYcEtaReflMult_SP, hVyXcEtaReflMult_SP;
    vector <TH1F*> hVxEtaReflMult_SP, hVyEtaReflMult_SP, hVEtaReflMult_SP, hVxYEtaReflMult_SP, hVyXEtaReflMult_SP;
    vector <TH1F*> hVxaEtaReflMult_EP, hVyaEtaReflMult_EP, hVaEtaReflMult_EP, hVxYaEtaReflMult_EP, hVyXaEtaReflMult_EP;
    vector <TH1F*> hVxbEtaReflMult_EP, hVybEtaReflMult_EP, hVbEtaReflMult_EP, hVxYbEtaReflMult_EP, hVyXbEtaReflMult_EP;
    vector <TH1F*> hVxcEtaReflMult_EP, hVycEtaReflMult_EP, hVcEtaReflMult_EP, hVxYcEtaReflMult_EP, hVyXcEtaReflMult_EP;
    vector <TH1F*> hVxEtaReflMult_EP, hVyEtaReflMult_EP, hVEtaReflMult_EP, hVxYEtaReflMult_EP, hVyXEtaReflMult_EP;

    vector <TH1F*> hVxaEtaReflCentBS_SP, hVyaEtaReflCentBS_SP, hVaEtaReflCentBS_SP, hVxYaEtaReflCentBS_SP, hVyXaEtaReflCentBS_SP;
    vector <TH1F*> hVxbEtaReflCentBS_SP, hVybEtaReflCentBS_SP, hVbEtaReflCentBS_SP, hVxYbEtaReflCentBS_SP, hVyXbEtaReflCentBS_SP;
    vector <TH1F*> hVxcEtaReflCentBS_SP, hVycEtaReflCentBS_SP, hVcEtaReflCentBS_SP, hVxYcEtaReflCentBS_SP, hVyXcEtaReflCentBS_SP;
    vector <TH1F*> hVxEtaReflCentBS_SP, hVyEtaReflCentBS_SP, hVEtaReflCentBS_SP, hVxYEtaReflCentBS_SP, hVyXEtaReflCentBS_SP;
    vector <TH1F*> hVxaEtaReflCentBS_EP, hVyaEtaReflCentBS_EP, hVaEtaReflCentBS_EP, hVxYaEtaReflCentBS_EP, hVyXaEtaReflCentBS_EP;
    vector <TH1F*> hVxbEtaReflCentBS_EP, hVybEtaReflCentBS_EP, hVbEtaReflCentBS_EP, hVxYbEtaReflCentBS_EP, hVyXbEtaReflCentBS_EP;
    vector <TH1F*> hVxcEtaReflCentBS_EP, hVycEtaReflCentBS_EP, hVcEtaReflCentBS_EP, hVxYcEtaReflCentBS_EP, hVyXcEtaReflCentBS_EP;
    vector <TH1F*> hVxEtaReflCentBS_EP, hVyEtaReflCentBS_EP, hVEtaReflCentBS_EP, hVxYEtaReflCentBS_EP, hVyXEtaReflCentBS_EP;

    vector <TH1F*> hVxaEtaReflMultBS_SP, hVyaEtaReflMultBS_SP, hVaEtaReflMultBS_SP, hVxYaEtaReflMultBS_SP, hVyXaEtaReflMultBS_SP;
    vector <TH1F*> hVxbEtaReflMultBS_SP, hVybEtaReflMultBS_SP, hVbEtaReflMultBS_SP, hVxYbEtaReflMultBS_SP, hVyXbEtaReflMultBS_SP;
    vector <TH1F*> hVxcEtaReflMultBS_SP, hVycEtaReflMultBS_SP, hVcEtaReflMultBS_SP, hVxYcEtaReflMultBS_SP, hVyXcEtaReflMultBS_SP;
    vector <TH1F*> hVxEtaReflMultBS_SP, hVyEtaReflMultBS_SP, hVEtaReflMultBS_SP, hVxYEtaReflMultBS_SP, hVyXEtaReflMultBS_SP;
    vector <TH1F*> hVxaEtaReflMultBS_EP, hVyaEtaReflMultBS_EP, hVaEtaReflMultBS_EP, hVxYaEtaReflMultBS_EP, hVyXaEtaReflMultBS_EP;
    vector <TH1F*> hVxbEtaReflMultBS_EP, hVybEtaReflMultBS_EP, hVbEtaReflMultBS_EP, hVxYbEtaReflMultBS_EP, hVyXbEtaReflMultBS_EP;
    vector <TH1F*> hVxcEtaReflMultBS_EP, hVycEtaReflMultBS_EP, hVcEtaReflMultBS_EP, hVxYcEtaReflMultBS_EP, hVyXcEtaReflMultBS_EP;
    vector <TH1F*> hVxEtaReflMultBS_EP, hVyEtaReflMultBS_EP, hVEtaReflMultBS_EP, hVxYEtaReflMultBS_EP, hVyXEtaReflMultBS_EP;

        for (Int_t i = 0; i < nHarmonics; i++) {
            n = harmonicsMap [i];

            pxXaCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%iaCent_SP", n, n)));
            pyYaCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%iaCent_SP", n, n)));
            pyXaCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%iaCent_SP", n, n)));
            pxYaCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%iaCent_SP", n, n)));
            pxXbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%ibCent_SP", n, n)));
            pyYbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%ibCent_SP", n, n)));
            pyXbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%ibCent_SP", n, n)));
            pxYbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%ibCent_SP", n, n)));
            pxXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%icCent_SP", n, n)));
            pyYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%icCent_SP", n, n)));
            pyXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%icCent_SP", n, n)));
            pxYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%icCent_SP", n, n)));
            pxXaCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%iaCent_EP", n, n)));
            pyYaCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%iaCent_EP", n, n)));
            pyXaCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%iaCent_EP", n, n)));
            pxYaCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%iaCent_EP", n, n)));
            pxXbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%ibCent_EP", n, n)));
            pyYbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%ibCent_EP", n, n)));
            pyXbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%ibCent_EP", n, n)));
            pxYbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%ibCent_EP", n, n)));
            pxXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%icCent_EP", n, n)));
            pyYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%icCent_EP", n, n)));
            pyXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%icCent_EP", n, n)));
            pxYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%icCent_EP", n, n)));

            pxXaMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%iaMult_SP", n, n)));
            pyYaMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%iaMult_SP", n, n)));
            pyXaMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%iaMult_SP", n, n)));
            pxYaMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%iaMult_SP", n, n)));
            pxXbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%ibMult_SP", n, n)));
            pyYbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%ibMult_SP", n, n)));
            pyXbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%ibMult_SP", n, n)));
            pxYbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%ibMult_SP", n, n)));
            pxXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%icMult_SP", n, n)));
            pyYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%icMult_SP", n, n)));
            pyXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%icMult_SP", n, n)));
            pxYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%icMult_SP", n, n)));
            pxXaMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%iaMult_EP", n, n)));
            pyYaMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%iaMult_EP", n, n)));
            pyXaMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%iaMult_EP", n, n)));
            pxYaMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%iaMult_EP", n, n)));
            pxXbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%ibMult_EP", n, n)));
            pyYbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%ibMult_EP", n, n)));
            pyXbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%ibMult_EP", n, n)));
            pxYbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%ibMult_EP", n, n)));
            pxXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iX%icMult_EP", n, n)));
            pyYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iY%icMult_EP", n, n)));
            pyXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("py%iX%icMult_EP", n, n)));
            pxYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("px%iY%icMult_EP", n, n)));

            p2xXaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaCent_SP", n, n)));
            p2yYaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaCent_SP", n, n)));
            p2yXaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaCent_SP", n, n)));
            p2xYaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaCent_SP", n, n)));
            p2xXbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibCent_SP", n, n)));
            p2yYbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibCent_SP", n, n)));
            p2yXbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibCent_SP", n, n)));
            p2xYbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibCent_SP", n, n)));
            p2xXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icCent_SP", n, n)));
            p2yYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icCent_SP", n, n)));
            p2yXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icCent_SP", n, n)));
            p2xYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icCent_SP", n, n)));
            p2xXaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaCent_EP", n, n)));
            p2yYaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaCent_EP", n, n)));
            p2yXaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaCent_EP", n, n)));
            p2xYaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaCent_EP", n, n)));
            p2xXbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibCent_EP", n, n)));
            p2yYbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibCent_EP", n, n)));
            p2yXbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibCent_EP", n, n)));
            p2xYbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibCent_EP", n, n)));
            p2xXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icCent_EP", n, n)));
            p2yYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icCent_EP", n, n)));
            p2yXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icCent_EP", n, n)));
            p2xYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icCent_EP", n, n)));

            p2xXaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaMult_SP", n, n)));
            p2yYaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaMult_SP", n, n)));
            p2yXaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaMult_SP", n, n)));
            p2xYaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaMult_SP", n, n)));
            p2xXbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibMult_SP", n, n)));
            p2yYbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibMult_SP", n, n)));
            p2yXbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibMult_SP", n, n)));
            p2xYbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibMult_SP", n, n)));
            p2xXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icMult_SP", n, n)));
            p2yYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icMult_SP", n, n)));
            p2yXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icMult_SP", n, n)));
            p2xYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icMult_SP", n, n)));
            p2xXaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaMult_EP", n, n)));
            p2yYaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaMult_EP", n, n)));
            p2yXaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaMult_EP", n, n)));
            p2xYaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaMult_EP", n, n)));
            p2xXbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibMult_EP", n, n)));
            p2yYbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibMult_EP", n, n)));
            p2yXbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibMult_EP", n, n)));
            p2xYbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibMult_EP", n, n)));
            p2xXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icMult_EP", n, n)));
            p2yYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icMult_EP", n, n)));
            p2yXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icMult_EP", n, n)));
            p2xYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icMult_EP", n, n)));

            p2xXaPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaPtCent_SP", n, n)));
            p2xYaPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaPtCent_SP", n, n)));
            p2yXaPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaPtCent_SP", n, n)));
            p2yYaPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaPtCent_SP", n, n)));
            p2xXbPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibPtCent_SP", n, n)));
            p2xYbPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibPtCent_SP", n, n)));
            p2yXbPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibPtCent_SP", n, n)));
            p2yYbPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibPtCent_SP", n, n)));
            p2xXcPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icPtCent_SP", n, n)));
            p2xYcPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icPtCent_SP", n, n)));
            p2yXcPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icPtCent_SP", n, n)));
            p2yYcPtCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icPtCent_SP", n, n)));
            p2xXaPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaPtCent_EP", n, n)));
            p2xYaPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaPtCent_EP", n, n)));
            p2yXaPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaPtCent_EP", n, n)));
            p2yYaPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaPtCent_EP", n, n)));
            p2xXbPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibPtCent_EP", n, n)));
            p2xYbPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibPtCent_EP", n, n)));
            p2yXbPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibPtCent_EP", n, n)));
            p2yYbPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibPtCent_EP", n, n)));
            p2xXcPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icPtCent_EP", n, n)));
            p2xYcPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icPtCent_EP", n, n)));
            p2yXcPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icPtCent_EP", n, n)));
            p2yYcPtCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icPtCent_EP", n, n)));

            p2xXaEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaEtaCent_SP", n, n)));
            p2xYaEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaEtaCent_SP", n, n)));
            p2yXaEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaEtaCent_SP", n, n)));
            p2yYaEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaEtaCent_SP", n, n)));
            p2xXbEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibEtaCent_SP", n, n)));
            p2xYbEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibEtaCent_SP", n, n)));
            p2yXbEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibEtaCent_SP", n, n)));
            p2yYbEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibEtaCent_SP", n, n)));
            p2xXcEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icEtaCent_SP", n, n)));
            p2xYcEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icEtaCent_SP", n, n)));
            p2yXcEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icEtaCent_SP", n, n)));
            p2yYcEtaCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icEtaCent_SP", n, n)));
            p2xXaEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaEtaCent_EP", n, n)));
            p2xYaEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaEtaCent_EP", n, n)));
            p2yXaEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaEtaCent_EP", n, n)));
            p2yYaEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaEtaCent_EP", n, n)));
            p2xXbEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibEtaCent_EP", n, n)));
            p2xYbEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibEtaCent_EP", n, n)));
            p2yXbEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibEtaCent_EP", n, n)));
            p2yYbEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibEtaCent_EP", n, n)));
            p2xXcEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icEtaCent_EP", n, n)));
            p2xYcEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icEtaCent_EP", n, n)));
            p2yXcEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icEtaCent_EP", n, n)));
            p2yYcEtaCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icEtaCent_EP", n, n)));

            p2xXaPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaPtMult_SP", n, n)));
            p2xYaPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaPtMult_SP", n, n)));
            p2yXaPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaPtMult_SP", n, n)));
            p2yYaPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaPtMult_SP", n, n)));
            p2xXbPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibPtMult_SP", n, n)));
            p2xYbPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibPtMult_SP", n, n)));
            p2yXbPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibPtMult_SP", n, n)));
            p2yYbPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibPtMult_SP", n, n)));
            p2xXcPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icPtMult_SP", n, n)));
            p2xYcPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icPtMult_SP", n, n)));
            p2yXcPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icPtMult_SP", n, n)));
            p2yYcPtMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icPtMult_SP", n, n)));
            p2xXaPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaPtMult_EP", n, n)));
            p2xYaPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaPtMult_EP", n, n)));
            p2yXaPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaPtMult_EP", n, n)));
            p2yYaPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaPtMult_EP", n, n)));
            p2xXbPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibPtMult_EP", n, n)));
            p2xYbPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibPtMult_EP", n, n)));
            p2yXbPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibPtMult_EP", n, n)));
            p2yYbPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibPtMult_EP", n, n)));
            p2xXcPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icPtMult_EP", n, n)));
            p2xYcPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icPtMult_EP", n, n)));
            p2yXcPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icPtMult_EP", n, n)));
            p2yYcPtMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icPtMult_EP", n, n)));

            p2xXaEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaEtaMult_SP", n, n)));
            p2xYaEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaEtaMult_SP", n, n)));
            p2yXaEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaEtaMult_SP", n, n)));
            p2yYaEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaEtaMult_SP", n, n)));
            p2xXbEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibEtaMult_SP", n, n)));
            p2xYbEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibEtaMult_SP", n, n)));
            p2yXbEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibEtaMult_SP", n, n)));
            p2yYbEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibEtaMult_SP", n, n)));
            p2xXcEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icEtaMult_SP", n, n)));
            p2xYcEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icEtaMult_SP", n, n)));
            p2yXcEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icEtaMult_SP", n, n)));
            p2yYcEtaMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icEtaMult_SP", n, n)));
            p2xXaEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%iaEtaMult_EP", n, n)));
            p2xYaEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%iaEtaMult_EP", n, n)));
            p2yXaEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%iaEtaMult_EP", n, n)));
            p2yYaEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%iaEtaMult_EP", n, n)));
            p2xXbEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%ibEtaMult_EP", n, n)));
            p2xYbEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%ibEtaMult_EP", n, n)));
            p2yXbEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%ibEtaMult_EP", n, n)));
            p2yYbEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%ibEtaMult_EP", n, n)));
            p2xXcEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iX%icEtaMult_EP", n, n)));
            p2xYcEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2x%iY%icEtaMult_EP", n, n)));
            p2yXcEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iX%icEtaMult_EP", n, n)));
            p2yYcEtaMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2y%iY%icEtaMult_EP", n, n)));

            p3xXaPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaPtCent_SP", n, n)));
            p3xYaPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaPtCent_SP", n, n)));
            p3yXaPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaPtCent_SP", n, n)));
            p3yYaPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaPtCent_SP", n, n)));
            p3xXbPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibPtCent_SP", n, n)));
            p3xYbPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibPtCent_SP", n, n)));
            p3yXbPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibPtCent_SP", n, n)));
            p3yYbPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibPtCent_SP", n, n)));
            p3xXcPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icPtCent_SP", n, n)));
            p3xYcPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icPtCent_SP", n, n)));
            p3yXcPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icPtCent_SP", n, n)));
            p3yYcPtCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icPtCent_SP", n, n)));
            p3xXaPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaPtCent_EP", n, n)));
            p3xYaPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaPtCent_EP", n, n)));
            p3yXaPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaPtCent_EP", n, n)));
            p3yYaPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaPtCent_EP", n, n)));
            p3xXbPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibPtCent_EP", n, n)));
            p3xYbPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibPtCent_EP", n, n)));
            p3yXbPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibPtCent_EP", n, n)));
            p3yYbPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibPtCent_EP", n, n)));
            p3xXcPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icPtCent_EP", n, n)));
            p3xYcPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icPtCent_EP", n, n)));
            p3yXcPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icPtCent_EP", n, n)));
            p3yYcPtCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icPtCent_EP", n, n)));

            p3xXaEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaEtaCent_SP", n, n)));
            p3xYaEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaEtaCent_SP", n, n)));
            p3yXaEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaEtaCent_SP", n, n)));
            p3yYaEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaEtaCent_SP", n, n)));
            p3xXbEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibEtaCent_SP", n, n)));
            p3xYbEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibEtaCent_SP", n, n)));
            p3yXbEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibEtaCent_SP", n, n)));
            p3yYbEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibEtaCent_SP", n, n)));
            p3xXcEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icEtaCent_SP", n, n)));
            p3xYcEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icEtaCent_SP", n, n)));
            p3yXcEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icEtaCent_SP", n, n)));
            p3yYcEtaCent_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icEtaCent_SP", n, n)));
            p3xXaEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaEtaCent_EP", n, n)));
            p3xYaEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaEtaCent_EP", n, n)));
            p3yXaEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaEtaCent_EP", n, n)));
            p3yYaEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaEtaCent_EP", n, n)));
            p3xXbEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibEtaCent_EP", n, n)));
            p3xYbEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibEtaCent_EP", n, n)));
            p3yXbEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibEtaCent_EP", n, n)));
            p3yYbEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibEtaCent_EP", n, n)));
            p3xXcEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icEtaCent_EP", n, n)));
            p3xYcEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icEtaCent_EP", n, n)));
            p3yXcEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icEtaCent_EP", n, n)));
            p3yYcEtaCent_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icEtaCent_EP", n, n)));

            p3xXaPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaPtMult_SP", n, n)));
            p3xYaPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaPtMult_SP", n, n)));
            p3yXaPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaPtMult_SP", n, n)));
            p3yYaPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaPtMult_SP", n, n)));
            p3xXbPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibPtMult_SP", n, n)));
            p3xYbPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibPtMult_SP", n, n)));
            p3yXbPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibPtMult_SP", n, n)));
            p3yYbPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibPtMult_SP", n, n)));
            p3xXcPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icPtMult_SP", n, n)));
            p3xYcPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icPtMult_SP", n, n)));
            p3yXcPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icPtMult_SP", n, n)));
            p3yYcPtMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icPtMult_SP", n, n)));
            p3xXaPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaPtMult_EP", n, n)));
            p3xYaPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaPtMult_EP", n, n)));
            p3yXaPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaPtMult_EP", n, n)));
            p3yYaPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaPtMult_EP", n, n)));
            p3xXbPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibPtMult_EP", n, n)));
            p3xYbPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibPtMult_EP", n, n)));
            p3yXbPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibPtMult_EP", n, n)));
            p3yYbPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibPtMult_EP", n, n)));
            p3xXcPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icPtMult_EP", n, n)));
            p3xYcPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icPtMult_EP", n, n)));
            p3yXcPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icPtMult_EP", n, n)));
            p3yYcPtMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icPtMult_EP", n, n)));

            p3xXaEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaEtaMult_SP", n, n)));
            p3xYaEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaEtaMult_SP", n, n)));
            p3yXaEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaEtaMult_SP", n, n)));
            p3yYaEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaEtaMult_SP", n, n)));
            p3xXbEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibEtaMult_SP", n, n)));
            p3xYbEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibEtaMult_SP", n, n)));
            p3yXbEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibEtaMult_SP", n, n)));
            p3yYbEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibEtaMult_SP", n, n)));
            p3xXcEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icEtaMult_SP", n, n)));
            p3xYcEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icEtaMult_SP", n, n)));
            p3yXcEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icEtaMult_SP", n, n)));
            p3yYcEtaMult_SP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icEtaMult_SP", n, n)));
            p3xXaEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%iaEtaMult_EP", n, n)));
            p3xYaEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%iaEtaMult_EP", n, n)));
            p3yXaEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%iaEtaMult_EP", n, n)));
            p3yYaEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%iaEtaMult_EP", n, n)));
            p3xXbEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%ibEtaMult_EP", n, n)));
            p3xYbEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%ibEtaMult_EP", n, n)));
            p3yXbEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%ibEtaMult_EP", n, n)));
            p3yYbEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%ibEtaMult_EP", n, n)));
            p3xXcEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iX%icEtaMult_EP", n, n)));
            p3xYcEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3x%iY%icEtaMult_EP", n, n)));
            p3yXcEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iX%icEtaMult_EP", n, n)));
            p3yYcEtaMult_EP.push_back ((TProfile3D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p3y%iY%icEtaMult_EP", n, n)));

            pXaXbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%ibCent_SP", n, n)));
            pXaYbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%ibCent_SP", n, n)));
            pYaXbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%ibCent_SP", n, n)));
            pYaYbCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%ibCent_SP", n, n)));
            pXaXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%icCent_SP", n, n)));
            pXaYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%icCent_SP", n, n)));
            pYaXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%icCent_SP", n, n)));
            pYaYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%icCent_SP", n, n)));
            pXbXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibX%icCent_SP", n, n)));
            pXbYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibY%icCent_SP", n, n)));
            pYbXcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibX%icCent_SP", n, n)));
            pYbYcCent_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibY%icCent_SP", n, n)));
            pXaXbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%ibCent_EP", n, n)));
            pXaYbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%ibCent_EP", n, n)));
            pYaXbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%ibCent_EP", n, n)));
            pYaYbCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%ibCent_EP", n, n)));
            pXaXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%icCent_EP", n, n)));
            pXaYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%icCent_EP", n, n)));
            pYaXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%icCent_EP", n, n)));
            pYaYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%icCent_EP", n, n)));
            pXbXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibX%icCent_EP", n, n)));
            pXbYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibY%icCent_EP", n, n)));
            pYbXcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibX%icCent_EP", n, n)));
            pYbYcCent_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibY%icCent_EP", n, n)));

            pXaXbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%ibMult_SP", n, n)));
            pXaYbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%ibMult_SP", n, n)));
            pYaXbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%ibMult_SP", n, n)));
            pYaYbMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%ibMult_SP", n, n)));
            pXaXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%icMult_SP", n, n)));
            pXaYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%icMult_SP", n, n)));
            pYaXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%icMult_SP", n, n)));
            pYaYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%icMult_SP", n, n)));
            pXbXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibX%icMult_SP", n, n)));
            pXbYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibY%icMult_SP", n, n)));
            pYbXcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibX%icMult_SP", n, n)));
            pYbYcMult_SP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibY%icMult_SP", n, n)));
            pXaXbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%ibMult_EP", n, n)));
            pXaYbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%ibMult_EP", n, n)));
            pYaXbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%ibMult_EP", n, n)));
            pYaYbMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%ibMult_EP", n, n)));
            pXaXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaX%icMult_EP", n, n)));
            pXaYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%iaY%icMult_EP", n, n)));
            pYaXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaX%icMult_EP", n, n)));
            pYaYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%iaY%icMult_EP", n, n)));
            pXbXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibX%icMult_EP", n, n)));
            pXbYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pX%ibY%icMult_EP", n, n)));
            pYbXcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibX%icMult_EP", n, n)));
            pYbYcMult_EP.push_back ((TProfile*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("pY%ibY%icMult_EP", n, n)));

            p2XaXbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%ibCent_SP", n, n)));
            p2XaYbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%ibCent_SP", n, n)));
            p2YaXbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%ibCent_SP", n, n)));
            p2YaYbCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%ibCent_SP", n, n)));
            p2XaXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%icCent_SP", n, n)));
            p2XaYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%icCent_SP", n, n)));
            p2YaXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%icCent_SP", n, n)));
            p2YaYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%icCent_SP", n, n)));
            p2XbXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibX%icCent_SP", n, n)));
            p2XbYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibY%icCent_SP", n, n)));
            p2YbXcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibX%icCent_SP", n, n)));
            p2YbYcCent_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibY%icCent_SP", n, n)));
            p2XaXbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%ibCent_EP", n, n)));
            p2XaYbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%ibCent_EP", n, n)));
            p2YaXbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%ibCent_EP", n, n)));
            p2YaYbCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%ibCent_EP", n, n)));
            p2XaXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%icCent_EP", n, n)));
            p2XaYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%icCent_EP", n, n)));
            p2YaXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%icCent_EP", n, n)));
            p2YaYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%icCent_EP", n, n)));
            p2XbXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibX%icCent_EP", n, n)));
            p2XbYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibY%icCent_EP", n, n)));
            p2YbXcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibX%icCent_EP", n, n)));
            p2YbYcCent_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibY%icCent_EP", n, n)));

            p2XaXbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%ibMult_SP", n, n)));
            p2XaYbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%ibMult_SP", n, n)));
            p2YaXbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%ibMult_SP", n, n)));
            p2YaYbMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%ibMult_SP", n, n)));
            p2XaXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%icMult_SP", n, n)));
            p2XaYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%icMult_SP", n, n)));
            p2YaXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%icMult_SP", n, n)));
            p2YaYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%icMult_SP", n, n)));
            p2XbXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibX%icMult_SP", n, n)));
            p2XbYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibY%icMult_SP", n, n)));
            p2YbXcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibX%icMult_SP", n, n)));
            p2YbYcMult_SP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibY%icMult_SP", n, n)));
            p2XaXbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%ibMult_EP", n, n)));
            p2XaYbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%ibMult_EP", n, n)));
            p2YaXbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%ibMult_EP", n, n)));
            p2YaYbMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%ibMult_EP", n, n)));
            p2XaXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaX%icMult_EP", n, n)));
            p2XaYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%iaY%icMult_EP", n, n)));
            p2YaXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaX%icMult_EP", n, n)));
            p2YaYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%iaY%icMult_EP", n, n)));
            p2XbXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibX%icMult_EP", n, n)));
            p2XbYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2X%ibY%icMult_EP", n, n)));
            p2YbXcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibX%icMult_EP", n, n)));
            p2YbYcMult_EP.push_back ((TProfile2D*) corrFile -> Get (dirName [step] + "/Source Histograms/" + Form ("p2Y%ibY%icMult_EP", n, n)));

            nBinsBS_ = p2XaXbCent_SP [i] -> GetNbinsY ();
			nBinsCent_ = p2XaXbCent_SP [i] -> GetNbinsX ();
			centMin_ = p2XaXbCent_SP [i] -> GetXaxis() -> GetXmin ();
			centMax_ = p2XaXbCent_SP [i] -> GetXaxis() -> GetXmax ();
            nBinsMh_ = p2XaXbMult_SP [i] -> GetNbinsX ();
			mhMin_ = p2XaXbMult_SP [i] -> GetXaxis() -> GetXmin ();
			mhMax_ = p2XaXbMult_SP [i] -> GetXaxis() -> GetXmax ();
			nBinsEta_ = p3xXaEtaCent_SP [i] -> GetNbinsX ();
			etaMin_ = p3xXaEtaCent_SP [i] -> GetXaxis() -> GetXmin ();
			etaMax_ = p3xXaEtaCent_SP [i] -> GetXaxis() -> GetXmax ();
			nBinsPt_ = p3xXaPtCent_SP [i] -> GetNbinsX ();
			ptMin_ = p3xXaPtCent_SP [i] -> GetXaxis() -> GetXmin ();
			ptMax_ = p3xXaPtCent_SP [i] -> GetXaxis() -> GetXmax ();
            mhLowerBin_ = p2XaXbMult_SP [i] -> GetXaxis () -> FindBin (mhLow_ + 0.1);
            mhHigherBin_ = p2XaXbMult_SP [i] -> GetXaxis () -> FindBin (mhHigh_ - 0.1);
            centLowerBin_ = p2XaXbCent_SP [i] -> GetXaxis () -> FindBin (centLow_ + 0.001);
            centHigherBin_ = p2XaXbCent_SP [i] -> GetXaxis () -> FindBin (centHigh_ - 0.001);

            cout << "BS: " << nBinsBS_ << endl;
			cout << "mh: " << nBinsMh_ << " " << mhMin_ << " " << mhMax_ << endl;
			cout << "mhFlowRange" << " " << mhLow_ << " " << mhHigh_ << endl;
			cout << "cent: " << nBinsCent_ << " " << centMin_ << " " << centMax_ << endl;
			cout << "centFlowRange" << " " << centLow_ << " " << centHigh_ << endl;
			cout << "pt: " << nBinsPt_ << " " << ptMin_ << " " << ptMax_ << endl;
			cout << "eta: " << nBinsEta_ << " " << etaMin_ << " " << etaMax_ << endl;
            cout << "mhBins: " << mhLowerBin_ << " to " << mhHigherBin_ << endl;
            cout << "centBins: " << centLowerBin_ << " to " << centHigherBin_ << endl;

			stepDir -> cd ();


			pqQaCent_SP.push_back (new TProfile (Form ("pq%iQ%iaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
			pqQbCent_SP.push_back (new TProfile (Form ("pq%iQ%ibCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
			pqQcCent_SP.push_back (new TProfile (Form ("pq%iQ%icCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
			pqQaCent_EP.push_back (new TProfile (Form ("pq%iQ%iaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
			pqQbCent_EP.push_back (new TProfile (Form ("pq%iQ%ibCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));
			pqQcCent_EP.push_back (new TProfile (Form ("pq%iQ%icCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_));

			pqQaMult_SP.push_back (new TProfile (Form ("pq%iQ%iaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));
			pqQbMult_SP.push_back (new TProfile (Form ("pq%iQ%ibMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));
			pqQcMult_SP.push_back (new TProfile (Form ("pq%iQ%icMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));
			pqQaMult_EP.push_back (new TProfile (Form ("pq%iQ%iaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));
			pqQbMult_EP.push_back (new TProfile (Form ("pq%iQ%ibMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));
			pqQcMult_EP.push_back (new TProfile (Form ("pq%iQ%icMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_));

			p2qQaCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p2qQbCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p2qQcCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%icCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p2qQaCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p2qQbCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p2qQcCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%icCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

			p2qQaMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p2qQbMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p2qQcMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%icMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p2qQaMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p2qQbMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p2qQcMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%icMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

			p2qQaPtCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
			p2qQbPtCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
			p2qQcPtCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%icPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
			p2qQaPtCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
			p2qQbPtCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));
			p2qQcPtCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%icPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_));

			p2qQaEtaCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
			p2qQbEtaCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
			p2qQcEtaCent_SP.push_back (new TProfile2D (Form ("p2q%iQ%icEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
            p2qQaEtaCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
			p2qQbEtaCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));
			p2qQcEtaCent_EP.push_back (new TProfile2D (Form ("p2q%iQ%icEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_));

			p2qQaPtMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQbPtMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQcPtMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%icPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQaPtMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQbPtMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQcPtMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%icPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_));

			p2qQaEtaMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%iaEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQbEtaMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%ibEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQcEtaMult_SP.push_back (new TProfile2D (Form ("p2q%iQ%icEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
            p2qQaEtaMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%iaEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQbEtaMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%ibEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));
			p2qQcEtaMult_EP.push_back (new TProfile2D (Form ("p2q%iQ%icEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_));

			p3qQaPtCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%iaPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbPtCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%ibPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcPtCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%icPtCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQaPtCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%iaPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbPtCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%ibPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcPtCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%icPtCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);P_{T};cent;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

			p3qQaEtaCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%iaEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbEtaCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%ibEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcEtaCent_SP.push_back (new TProfile3D (Form ("p3q%iQ%icEtaCent_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p3qQaEtaCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%iaEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbEtaCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%ibEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcEtaCent_EP.push_back (new TProfile3D (Form ("p3q%iQ%icEtaCent_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);" + varName_ + " ;cent;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

			p3qQaPtMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%iaPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbPtMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%ibPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcPtMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%icPtMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQaPtMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%iaPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbPtMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%ibPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcPtMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%icPtMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);P_{T};mult;sample", n, n), nBinsPt_, ptMin_, ptMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

			p3qQaEtaMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%iaEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbEtaMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%ibEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcEtaMult_SP.push_back (new TProfile3D (Form ("p3q%iQ%icEtaMult_SP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (SP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p3qQaEtaMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%iaEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{a}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQbEtaMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%ibEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{b}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			p3qQcEtaMult_EP.push_back (new TProfile3D (Form ("p3q%iQ%icEtaMult_EP", n, n), Form ("#LTq_{%i}Q_{%i}^{c}#GT (EP);" + varName_ + " ;mult;sample", n, n), nBinsEta_, etaMin_, etaMax_, nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

            pQaQbCent_SP.push_back (new TProfile (Form ("pQ%iaQ%ibCent_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQcCent_SP.push_back (new TProfile (Form ("pQ%iaQ%icCent_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQbQcCent_SP.push_back (new TProfile (Form ("pQ%ibQ%icCent_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQbCent_EP.push_back (new TProfile (Form ("pQ%iaQ%ibCent_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQcCent_EP.push_back (new TProfile (Form ("pQ%iaQ%icCent_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQbQcCent_EP.push_back (new TProfile (Form ("pQ%ibQ%icCent_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));

            pQaQbMult_SP.push_back (new TProfile (Form ("pQ%iaQ%ibMult_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQcMult_SP.push_back (new TProfile (Form ("pQ%iaQ%icMult_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQbQcMult_SP.push_back (new TProfile (Form ("pQ%ibQ%icMult_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQbMult_EP.push_back (new TProfile (Form ("pQ%iaQ%ibMult_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQcMult_EP.push_back (new TProfile (Form ("pQ%iaQ%icMult_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQbQcMult_EP.push_back (new TProfile (Form ("pQ%ibQ%icMult_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));

            p2QaQbCent_SP.push_back (new TProfile2D (Form ("p2Q%iaQ%ibCent_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQcCent_SP.push_back (new TProfile2D (Form ("p2Q%iaQ%icCent_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2QbQcCent_SP.push_back (new TProfile2D (Form ("p2Q%ibQ%icCent_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQbCent_EP.push_back (new TProfile2D (Form ("p2Q%iaQ%ibCent_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQcCent_EP.push_back (new TProfile2D (Form ("p2Q%iaQ%icCent_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2QbQcCent_EP.push_back (new TProfile2D (Form ("p2Q%ibQ%icCent_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);cent;sample", n, n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

            p2QaQbMult_SP.push_back (new TProfile2D (Form ("p2Q%iaQ%ibMult_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQcMult_SP.push_back (new TProfile2D (Form ("p2Q%iaQ%icMult_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2QbQcMult_SP.push_back (new TProfile2D (Form ("p2Q%ibQ%icMult_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQbMult_EP.push_back (new TProfile2D (Form ("p2Q%iaQ%ibMult_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2QaQcMult_EP.push_back (new TProfile2D (Form ("p2Q%iaQ%icMult_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2QbQcMult_EP.push_back (new TProfile2D (Form ("p2Q%ibQ%icMult_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);mult;sample", n, n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

            pQaQbCentBS_SP.push_back (new TProfile (Form ("pQ%iaQ%ibCentBS_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQcCentBS_SP.push_back (new TProfile (Form ("pQ%iaQ%icCentBS_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQbQcCentBS_SP.push_back (new TProfile (Form ("pQ%ibQ%icCentBS_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQbCentBS_EP.push_back (new TProfile (Form ("pQ%iaQ%ibCentBS_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQaQcCentBS_EP.push_back (new TProfile (Form ("pQ%iaQ%icCentBS_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));
            pQbQcCentBS_EP.push_back (new TProfile (Form ("pQ%ibQ%icCentBS_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);cent", n, n), nBinsCent_, centMin_, centMax_));

            pQaQbMultBS_SP.push_back (new TProfile (Form ("pQ%iaQ%ibMultBS_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQcMultBS_SP.push_back (new TProfile (Form ("pQ%iaQ%icMultBS_SP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQbQcMultBS_SP.push_back (new TProfile (Form ("pQ%ibQ%icMultBS_SP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (SP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQbMultBS_EP.push_back (new TProfile (Form ("pQ%iaQ%ibMultBS_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{b}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQaQcMultBS_EP.push_back (new TProfile (Form ("pQ%iaQ%icMultBS_EP", n, n), Form ("#LTQ_{%i}^{a}Q_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));
            pQbQcMultBS_EP.push_back (new TProfile (Form ("pQ%ibQ%icMultBS_EP", n, n), Form ("#LTQ_{%i}^{b}Q_{%i}^{c}#GT (EP);mult", n, n), nBinsMh_, mhMin_, mhMax_));

			pXaXbCentBS_SP.push_back (new TProfile (Form ("pX%iaX%ibCentBS_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaXcCentBS_SP.push_back (new TProfile (Form ("pX%iaX%icCentBS_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXbXcCentBS_SP.push_back (new TProfile (Form ("pX%ibX%icCentBS_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaYbCentBS_SP.push_back (new TProfile (Form ("pY%iaY%ibCentBS_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaYcCentBS_SP.push_back (new TProfile (Form ("pY%iaY%icCentBS_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYbYcCentBS_SP.push_back (new TProfile (Form ("pY%ibY%icCentBS_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaYbCentBS_SP.push_back (new TProfile (Form ("pX%iaY%ibCentBS_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaYcCentBS_SP.push_back (new TProfile (Form ("pX%iaY%icCentBS_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXbYcCentBS_SP.push_back (new TProfile (Form ("pX%ibY%icCentBS_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaXbCentBS_SP.push_back (new TProfile (Form ("pY%iaX%ibCentBS_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaXcCentBS_SP.push_back (new TProfile (Form ("pY%iaX%icCentBS_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYbXcCentBS_SP.push_back (new TProfile (Form ("pY%ibX%icCentBS_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaXbCentBS_EP.push_back (new TProfile (Form ("pX%iaX%ibCentBS_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaXcCentBS_EP.push_back (new TProfile (Form ("pX%iaX%icCentBS_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXbXcCentBS_EP.push_back (new TProfile (Form ("pX%ibX%icCentBS_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaYbCentBS_EP.push_back (new TProfile (Form ("pY%iaY%ibCentBS_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaYcCentBS_EP.push_back (new TProfile (Form ("pY%iaY%icCentBS_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYbYcCentBS_EP.push_back (new TProfile (Form ("pY%ibY%icCentBS_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaYbCentBS_EP.push_back (new TProfile (Form ("pX%iaY%ibCentBS_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXaYcCentBS_EP.push_back (new TProfile (Form ("pX%iaY%icCentBS_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pXbYcCentBS_EP.push_back (new TProfile (Form ("pX%ibY%icCentBS_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaXbCentBS_EP.push_back (new TProfile (Form ("pY%iaX%ibCentBS_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYaXcCentBS_EP.push_back (new TProfile (Form ("pY%iaX%icCentBS_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));
			pYbXcCentBS_EP.push_back (new TProfile (Form ("pY%ibX%icCentBS_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP, sampling);cent", n, n), nBinsCent_, centMin_, centMax_));

			pXaXbMultBS_SP.push_back (new TProfile (Form ("pX%iaX%ibMultBS_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaXcMultBS_SP.push_back (new TProfile (Form ("pX%iaX%icMultBS_SP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXbXcMultBS_SP.push_back (new TProfile (Form ("pX%ibX%icMultBS_SP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaYbMultBS_SP.push_back (new TProfile (Form ("pY%iaY%ibMultBS_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaYcMultBS_SP.push_back (new TProfile (Form ("pY%iaY%icMultBS_SP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYbYcMultBS_SP.push_back (new TProfile (Form ("pY%ibY%icMultBS_SP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaYbMultBS_SP.push_back (new TProfile (Form ("pX%iaY%ibMultBS_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaYcMultBS_SP.push_back (new TProfile (Form ("pX%iaY%icMultBS_SP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXbYcMultBS_SP.push_back (new TProfile (Form ("pX%ibY%icMultBS_SP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaXbMultBS_SP.push_back (new TProfile (Form ("pY%iaX%ibMultBS_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaXcMultBS_SP.push_back (new TProfile (Form ("pY%iaX%icMultBS_SP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYbXcMultBS_SP.push_back (new TProfile (Form ("pY%ibX%icMultBS_SP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (SP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaXbMultBS_EP.push_back (new TProfile (Form ("pX%iaX%ibMultBS_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaXcMultBS_EP.push_back (new TProfile (Form ("pX%iaX%icMultBS_EP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXbXcMultBS_EP.push_back (new TProfile (Form ("pX%ibX%icMultBS_EP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaYbMultBS_EP.push_back (new TProfile (Form ("pY%iaY%ibMultBS_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaYcMultBS_EP.push_back (new TProfile (Form ("pY%iaY%icMultBS_EP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYbYcMultBS_EP.push_back (new TProfile (Form ("pY%ibY%icMultBS_EP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaYbMultBS_EP.push_back (new TProfile (Form ("pX%iaY%ibMultBS_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{b}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXaYcMultBS_EP.push_back (new TProfile (Form ("pX%iaY%icMultBS_EP", n, n), Form ("#LTX_{%i}^{a}Y_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pXbYcMultBS_EP.push_back (new TProfile (Form ("pX%ibY%icMultBS_EP", n, n), Form ("#LTX_{%i}^{b}Y_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaXbMultBS_EP.push_back (new TProfile (Form ("pY%iaX%ibMultBS_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{b}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYaXcMultBS_EP.push_back (new TProfile (Form ("pY%iaX%icMultBS_EP", n, n), Form ("#LTY_{%i}^{a}X_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));
			pYbXcMultBS_EP.push_back (new TProfile (Form ("pY%ibX%icMultBS_EP", n, n), Form ("#LTY_{%i}^{b}X_{%i}^{c}#GT (EP, sampling);mult", n, n), nBinsMh_, mhMin_, mhMax_));

			h2RxaCent_SP.push_back (new TH2F (Form ("h2R%ixaCent_SP", n), Form ("R_{%i, a}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RxbCent_SP.push_back (new TH2F (Form ("h2R%ixbCent_SP", n), Form ("R_{%i, b}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RxcCent_SP.push_back (new TH2F (Form ("h2R%ixcCent_SP", n), Form ("R_{%i, c}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            h2RyaCent_SP.push_back (new TH2F (Form ("h2R%iyaCent_SP", n), Form ("R_{%i, a}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RybCent_SP.push_back (new TH2F (Form ("h2R%iybCent_SP", n), Form ("R_{%i, b}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RycCent_SP.push_back (new TH2F (Form ("h2R%iycCent_SP", n), Form ("R_{%i, c}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            h2RaCent_SP.push_back (new TH2F (Form ("h2R%iaCent_SP", n), Form ("R_{%i, a}^{x+y,SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RbCent_SP.push_back (new TH2F (Form ("h2R%ibCent_SP", n), Form ("R_{%i, b}^{x+y,SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RcCent_SP.push_back (new TH2F (Form ("h2R%icCent_SP", n), Form ("R_{%i, c}^{x+y,SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            h2RxaCent_EP.push_back (new TH2F (Form ("h2R%ixaCent_EP", n), Form ("R_{%i, a}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RxbCent_EP.push_back (new TH2F (Form ("h2R%ixbCent_EP", n), Form ("R_{%i, b}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RxcCent_EP.push_back (new TH2F (Form ("h2R%ixcCent_EP", n), Form ("R_{%i, c}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            h2RyaCent_EP.push_back (new TH2F (Form ("h2R%iyaCent_EP", n), Form ("R_{%i, a}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RybCent_EP.push_back (new TH2F (Form ("h2R%iybCent_EP", n), Form ("R_{%i, b}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RycCent_EP.push_back (new TH2F (Form ("h2R%iycCent_EP", n), Form ("R_{%i, c}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            h2RaCent_EP.push_back (new TH2F (Form ("h2R%iaCent_EP", n), Form ("R_{%i, a}^{x+y,EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RbCent_EP.push_back (new TH2F (Form ("h2R%ibCent_EP", n), Form ("R_{%i, b}^{x+y,EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
			h2RcCent_EP.push_back (new TH2F (Form ("h2R%icCent_EP", n), Form ("R_{%i, c}^{x+y,EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

			h2RxaMult_SP.push_back (new TH2F (Form ("h2R%ixaMult_SP", n), Form ("R_{%i, a}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RxbMult_SP.push_back (new TH2F (Form ("h2R%ixbMult_SP", n), Form ("R_{%i, b}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RxcMult_SP.push_back (new TH2F (Form ("h2R%ixcMult_SP", n), Form ("R_{%i, c}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            h2RyaMult_SP.push_back (new TH2F (Form ("h2R%iyaMult_SP", n), Form ("R_{%i, a}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RybMult_SP.push_back (new TH2F (Form ("h2R%iybMult_SP", n), Form ("R_{%i, b}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RycMult_SP.push_back (new TH2F (Form ("h2R%iycMult_SP", n), Form ("R_{%i, c}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            h2RaMult_SP.push_back (new TH2F (Form ("h2R%iaMult_SP", n), Form ("R_{%i, a}^{x+y,SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RbMult_SP.push_back (new TH2F (Form ("h2R%ibMult_SP", n), Form ("R_{%i, b}^{x+y,SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RcMult_SP.push_back (new TH2F (Form ("h2R%icMult_SP", n), Form ("R_{%i, c}^{x+y,SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            h2RxaMult_EP.push_back (new TH2F (Form ("h2R%ixaMult_EP", n), Form ("R_{%i, a}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RxbMult_EP.push_back (new TH2F (Form ("h2R%ixbMult_EP", n), Form ("R_{%i, b}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RxcMult_EP.push_back (new TH2F (Form ("h2R%ixcMult_EP", n), Form ("R_{%i, c}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            h2RyaMult_EP.push_back (new TH2F (Form ("h2R%iyaMult_EP", n), Form ("R_{%i, a}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RybMult_EP.push_back (new TH2F (Form ("h2R%iybMult_EP", n), Form ("R_{%i, b}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RycMult_EP.push_back (new TH2F (Form ("h2R%iycMult_EP", n), Form ("R_{%i, c}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            h2RaMult_EP.push_back (new TH2F (Form ("h2R%iaMult_EP", n), Form ("R_{%i, a}^{x+y,EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RbMult_EP.push_back (new TH2F (Form ("h2R%ibMult_EP", n), Form ("R_{%i, b}^{x+y,EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
			h2RcMult_EP.push_back (new TH2F (Form ("h2R%icMult_EP", n), Form ("R_{%i, c}^{x+y,EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

            hRxaCent_SP.push_back (new TH1F (Form ("hR%ixaCent_SP", n), Form ("R_{%i, a}^{x, SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRxbCent_SP.push_back (new TH1F (Form ("hR%ixbCent_SP", n), Form ("R_{%i, b}^{x, SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRxcCent_SP.push_back (new TH1F (Form ("hR%ixcCent_SP", n), Form ("R_{%i, c}^{x, SP};cent", n), nBinsCent_, centMin_, centMax_));
            hRyaCent_SP.push_back (new TH1F (Form ("hR%iyaCent_SP", n), Form ("R_{%i, a}^{y, SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRybCent_SP.push_back (new TH1F (Form ("hR%iybCent_SP", n), Form ("R_{%i, b}^{y, SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRycCent_SP.push_back (new TH1F (Form ("hR%iycCent_SP", n), Form ("R_{%i, c}^{y, SP};cent", n), nBinsCent_, centMin_, centMax_));
            hRaCent_SP.push_back (new TH1F (Form ("hR%iaCent_SP", n), Form ("R_{%i, a}^{x+y,SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRbCent_SP.push_back (new TH1F (Form ("hR%ibCent_SP", n), Form ("R_{%i, b}^{x+y,SP};cent", n), nBinsCent_, centMin_, centMax_));
			hRcCent_SP.push_back (new TH1F (Form ("hR%icCent_SP", n), Form ("R_{%i, c}^{x+y,SP};cent", n), nBinsCent_, centMin_, centMax_));
            hRxaCent_EP.push_back (new TH1F (Form ("hR%ixaCent_EP", n), Form ("R_{%i, a}^{x, EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRxbCent_EP.push_back (new TH1F (Form ("hR%ixbCent_EP", n), Form ("R_{%i, b}^{x, EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRxcCent_EP.push_back (new TH1F (Form ("hR%ixcCent_EP", n), Form ("R_{%i, c}^{x, EP};cent", n), nBinsCent_, centMin_, centMax_));
            hRyaCent_EP.push_back (new TH1F (Form ("hR%iyaCent_EP", n), Form ("R_{%i, a}^{y, EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRybCent_EP.push_back (new TH1F (Form ("hR%iybCent_EP", n), Form ("R_{%i, b}^{y, EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRycCent_EP.push_back (new TH1F (Form ("hR%iycCent_EP", n), Form ("R_{%i, c}^{y, EP};cent", n), nBinsCent_, centMin_, centMax_));
            hRaCent_EP.push_back (new TH1F (Form ("hR%iaCent_EP", n), Form ("R_{%i, a}^{x+y,EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRbCent_EP.push_back (new TH1F (Form ("hR%ibCent_EP", n), Form ("R_{%i, b}^{x+y,EP};cent", n), nBinsCent_, centMin_, centMax_));
			hRcCent_EP.push_back (new TH1F (Form ("hR%icCent_EP", n), Form ("R_{%i, c}^{x+y,EP};cent", n), nBinsCent_, centMin_, centMax_));

			hRxaMult_SP.push_back (new TH1F (Form ("hR%ixaMult_SP", n), Form ("R_{%i, a}^{x, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxbMult_SP.push_back (new TH1F (Form ("hR%ixbMult_SP", n), Form ("R_{%i, b}^{x, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxcMult_SP.push_back (new TH1F (Form ("hR%ixcMult_SP", n), Form ("R_{%i, c}^{x, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
            hRyaMult_SP.push_back (new TH1F (Form ("hR%iyaMult_SP", n), Form ("R_{%i, a}^{y, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRybMult_SP.push_back (new TH1F (Form ("hR%iybMult_SP", n), Form ("R_{%i, b}^{y, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRycMult_SP.push_back (new TH1F (Form ("hR%iycMult_SP", n), Form ("R_{%i, c}^{y, SP};mult", n), nBinsMh_, mhMin_, mhMax_));
            hRaMult_SP.push_back (new TH1F (Form ("hR%iaMult_SP", n), Form ("R_{%i, a}^{x+y,SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRbMult_SP.push_back (new TH1F (Form ("hR%ibMult_SP", n), Form ("R_{%i, b}^{x+y,SP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRcMult_SP.push_back (new TH1F (Form ("hR%icMult_SP", n), Form ("R_{%i, c}^{x+y,SP};mult", n), nBinsMh_, mhMin_, mhMax_));
            hRxaMult_EP.push_back (new TH1F (Form ("hR%ixaMult_EP", n), Form ("R_{%i, a}^{x, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxbMult_EP.push_back (new TH1F (Form ("hR%ixbMult_EP", n), Form ("R_{%i, b}^{x, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxcMult_EP.push_back (new TH1F (Form ("hR%ixcMult_EP", n), Form ("R_{%i, c}^{x, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
            hRyaMult_EP.push_back (new TH1F (Form ("hR%iyaMult_EP", n), Form ("R_{%i, a}^{y, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRybMult_EP.push_back (new TH1F (Form ("hR%iybMult_EP", n), Form ("R_{%i, b}^{y, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRycMult_EP.push_back (new TH1F (Form ("hR%iycMult_EP", n), Form ("R_{%i, c}^{y, EP};mult", n), nBinsMh_, mhMin_, mhMax_));
            hRaMult_EP.push_back (new TH1F (Form ("hR%iaMult_EP", n), Form ("R_{%i, a}^{x+y,EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRbMult_EP.push_back (new TH1F (Form ("hR%ibMult_EP", n), Form ("R_{%i, b}^{x+y,EP};mult", n), nBinsMh_, mhMin_, mhMax_));
			hRcMult_EP.push_back (new TH1F (Form ("hR%icMult_EP", n), Form ("R_{%i, c}^{x+y,EP};mult", n), nBinsMh_, mhMin_, mhMax_));

            hRxaCentBS_SP.push_back (new TH1F (Form ("hR%ixaCentBS_SP", n), Form ("R_{%i, a}^{x, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRxbCentBS_SP.push_back (new TH1F (Form ("hR%ixbCentBS_SP", n), Form ("R_{%i, b}^{x, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRxcCentBS_SP.push_back (new TH1F (Form ("hR%ixcCentBS_SP", n), Form ("R_{%i, c}^{x, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
            hRyaCentBS_SP.push_back (new TH1F (Form ("hR%iyaCentBS_SP", n), Form ("R_{%i, a}^{y, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRybCentBS_SP.push_back (new TH1F (Form ("hR%iybCentBS_SP", n), Form ("R_{%i, b}^{y, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRycCentBS_SP.push_back (new TH1F (Form ("hR%iycCentBS_SP", n), Form ("R_{%i, c}^{y, SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
            hRaCentBS_SP.push_back (new TH1F (Form ("hR%iaCentBS_SP", n), Form ("R_{%i, a}^{x+y,SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRbCentBS_SP.push_back (new TH1F (Form ("hR%ibCentBS_SP", n), Form ("R_{%i, b}^{x+y,SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRcCentBS_SP.push_back (new TH1F (Form ("hR%icCentBS_SP", n), Form ("R_{%i, c}^{x+y,SP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
            hRxaCentBS_EP.push_back (new TH1F (Form ("hR%ixaCentBS_EP", n), Form ("R_{%i, a}^{x, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRxbCentBS_EP.push_back (new TH1F (Form ("hR%ixbCentBS_EP", n), Form ("R_{%i, b}^{x, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRxcCentBS_EP.push_back (new TH1F (Form ("hR%ixcCentBS_EP", n), Form ("R_{%i, c}^{x, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
            hRyaCentBS_EP.push_back (new TH1F (Form ("hR%iyaCentBS_EP", n), Form ("R_{%i, a}^{y, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRybCentBS_EP.push_back (new TH1F (Form ("hR%iybCentBS_EP", n), Form ("R_{%i, b}^{y, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRycCentBS_EP.push_back (new TH1F (Form ("hR%iycCentBS_EP", n), Form ("R_{%i, c}^{y, EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
            hRaCentBS_EP.push_back (new TH1F (Form ("hR%iaCentBS_EP", n), Form ("R_{%i, a}^{x+y,EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRbCentBS_EP.push_back (new TH1F (Form ("hR%ibCentBS_EP", n), Form ("R_{%i, b}^{x+y,EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));
			hRcCentBS_EP.push_back (new TH1F (Form ("hR%icCentBS_EP", n), Form ("R_{%i, c}^{x+y,EP}(sampling);cent", n), nBinsCent_, centMin_, centMax_));

			hRxaMultBS_SP.push_back (new TH1F (Form ("hR%ixaMultBS_SP", n), Form ("R_{%i, a}^{x, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxbMultBS_SP.push_back (new TH1F (Form ("hR%ixbMultBS_SP", n), Form ("R_{%i, b}^{x, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxcMultBS_SP.push_back (new TH1F (Form ("hR%ixcMultBS_SP", n), Form ("R_{%i, c}^{x, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
            hRyaMultBS_SP.push_back (new TH1F (Form ("hR%iyaMultBS_SP", n), Form ("R_{%i, a}^{y, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRybMultBS_SP.push_back (new TH1F (Form ("hR%iybMultBS_SP", n), Form ("R_{%i, b}^{y, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRycMultBS_SP.push_back (new TH1F (Form ("hR%iycMultBS_SP", n), Form ("R_{%i, c}^{y, SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
            hRaMultBS_SP.push_back (new TH1F (Form ("hR%iaMultBS_SP", n), Form ("R_{%i, a}^{x+y,SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRbMultBS_SP.push_back (new TH1F (Form ("hR%ibMultBS_SP", n), Form ("R_{%i, b}^{x+y,SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRcMultBS_SP.push_back (new TH1F (Form ("hR%icMultBS_SP", n), Form ("R_{%i, c}^{x+y,SP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
            hRxaMultBS_EP.push_back (new TH1F (Form ("hR%ixaMultBS_EP", n), Form ("R_{%i, a}^{x, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxbMultBS_EP.push_back (new TH1F (Form ("hR%ixbMultBS_EP", n), Form ("R_{%i, b}^{x, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRxcMultBS_EP.push_back (new TH1F (Form ("hR%ixcMultBS_EP", n), Form ("R_{%i, c}^{x, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
            hRyaMultBS_EP.push_back (new TH1F (Form ("hR%iyaMultBS_EP", n), Form ("R_{%i, a}^{y, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRybMultBS_EP.push_back (new TH1F (Form ("hR%iybMultBS_EP", n), Form ("R_{%i, b}^{y, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRycMultBS_EP.push_back (new TH1F (Form ("hR%iycMultBS_EP", n), Form ("R_{%i, c}^{y, EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
            hRaMultBS_EP.push_back (new TH1F (Form ("hR%iaMultBS_EP", n), Form ("R_{%i, a}^{x+y,EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRbMultBS_EP.push_back (new TH1F (Form ("hR%ibMultBS_EP", n), Form ("R_{%i, b}^{x+y,EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));
			hRcMultBS_EP.push_back (new TH1F (Form ("hR%icMultBS_EP", n), Form ("R_{%i, c}^{x+y,EP}(sampling);mult", n), nBinsMh_, mhMin_, mhMax_));

            p2VxaCent_SP.push_back (new TProfile2D (Form ("p2V%ixaCent_SP", n), Form ("V_{%i, a}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbCent_SP.push_back (new TProfile2D (Form ("p2V%ixbCent_SP", n), Form ("V_{%i, b}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcCent_SP.push_back (new TProfile2D (Form ("p2V%ixcCent_SP", n), Form ("V_{%i, c}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxCent_SP.push_back (new TProfile2D (Form ("p2V%ixCent_SP", n), Form ("V_{%i}^{x, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaCent_SP.push_back (new TProfile2D (Form ("p2V%iyaCent_SP", n), Form ("V_{%i, a}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VybCent_SP.push_back (new TProfile2D (Form ("p2V%iybCent_SP", n), Form ("V_{%i, b}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VycCent_SP.push_back (new TProfile2D (Form ("p2V%iycCent_SP", n), Form ("V_{%i, c}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyCent_SP.push_back (new TProfile2D (Form ("p2V%iyCent_SP", n), Form ("V_{%i}^{y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VaCent_SP.push_back (new TProfile2D (Form ("p2V%iaCent_SP", n), Form ("V_{%i, a}^{x+y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VbCent_SP.push_back (new TProfile2D (Form ("p2V%ibCent_SP", n), Form ("V_{%i, b}^{x+y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VcCent_SP.push_back (new TProfile2D (Form ("p2V%icCent_SP", n), Form ("V_{%i, c}^{x+y, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VCent_SP.push_back (new TProfile2D (Form ("p2V%iCent_SP", n), Form ("V_{%i}^{SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaCent_SP.push_back (new TProfile2D (Form ("p2V%ixYaCent_SP", n), Form ("V_{%i, a}^{xY, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbCent_SP.push_back (new TProfile2D (Form ("p2V%ixYbCent_SP", n), Form ("V_{%i, b}^{xY, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcCent_SP.push_back (new TProfile2D (Form ("p2V%ixYcCent_SP", n), Form ("V_{%i, c}^{xY, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaCent_SP.push_back (new TProfile2D (Form ("p2V%iyXaCent_SP", n), Form ("V_{%i, a}^{yX, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbCent_SP.push_back (new TProfile2D (Form ("p2V%iyXbCent_SP", n), Form ("V_{%i, b}^{yX, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcCent_SP.push_back (new TProfile2D (Form ("p2V%iyXcCent_SP", n), Form ("V_{%i, c}^{yX, SP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaCent_EP.push_back (new TProfile2D (Form ("p2V%ixaCent_EP", n), Form ("V_{%i, a}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbCent_EP.push_back (new TProfile2D (Form ("p2V%ixbCent_EP", n), Form ("V_{%i, b}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcCent_EP.push_back (new TProfile2D (Form ("p2V%ixcCent_EP", n), Form ("V_{%i, c}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxCent_EP.push_back (new TProfile2D (Form ("p2V%ixCent_EP", n), Form ("V_{%i}^{x, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaCent_EP.push_back (new TProfile2D (Form ("p2V%iyaCent_EP", n), Form ("V_{%i, a}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VybCent_EP.push_back (new TProfile2D (Form ("p2V%iybCent_EP", n), Form ("V_{%i, b}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VycCent_EP.push_back (new TProfile2D (Form ("p2V%iycCent_EP", n), Form ("V_{%i, c}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyCent_EP.push_back (new TProfile2D (Form ("p2V%iyCent_EP", n), Form ("V_{%i}^{y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VaCent_EP.push_back (new TProfile2D (Form ("p2V%iaCent_EP", n), Form ("V_{%i, a}^{x+y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VbCent_EP.push_back (new TProfile2D (Form ("p2V%ibCent_EP", n), Form ("V_{%i, b}^{x+y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VcCent_EP.push_back (new TProfile2D (Form ("p2V%icCent_EP", n), Form ("V_{%i, c}^{x+y, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VCent_EP.push_back (new TProfile2D (Form ("p2V%iCent_EP", n), Form ("V_{%i}^{EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaCent_EP.push_back (new TProfile2D (Form ("p2V%ixYaCent_EP", n), Form ("V_{%i, a}^{xY, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbCent_EP.push_back (new TProfile2D (Form ("p2V%ixYbCent_EP", n), Form ("V_{%i, b}^{xY, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcCent_EP.push_back (new TProfile2D (Form ("p2V%ixYcCent_EP", n), Form ("V_{%i, c}^{xY, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaCent_EP.push_back (new TProfile2D (Form ("p2V%iyXaCent_EP", n), Form ("V_{%i, a}^{yX, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbCent_EP.push_back (new TProfile2D (Form ("p2V%iyXbCent_EP", n), Form ("V_{%i, b}^{yX, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcCent_EP.push_back (new TProfile2D (Form ("p2V%iyXcCent_EP", n), Form ("V_{%i, c}^{yX, EP};cent;sample", n), nBinsCent_, centMin_, centMax_, nBinsBS_, 0, nBinsBS_));

            p2VxaMult_SP.push_back (new TProfile2D (Form ("p2V%ixaMult_SP", n), Form ("V_{%i, a}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbMult_SP.push_back (new TProfile2D (Form ("p2V%ixbMult_SP", n), Form ("V_{%i, b}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcMult_SP.push_back (new TProfile2D (Form ("p2V%ixcMult_SP", n), Form ("V_{%i, c}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxMult_SP.push_back (new TProfile2D (Form ("p2V%ixMult_SP", n), Form ("V_{%i}^{x, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaMult_SP.push_back (new TProfile2D (Form ("p2V%iyaMult_SP", n), Form ("V_{%i, a}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VybMult_SP.push_back (new TProfile2D (Form ("p2V%iybMult_SP", n), Form ("V_{%i, b}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VycMult_SP.push_back (new TProfile2D (Form ("p2V%iycMult_SP", n), Form ("V_{%i, c}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyMult_SP.push_back (new TProfile2D (Form ("p2V%iyMult_SP", n), Form ("V_{%i}^{y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VaMult_SP.push_back (new TProfile2D (Form ("p2V%iaMult_SP", n), Form ("V_{%i, a}^{x+y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VbMult_SP.push_back (new TProfile2D (Form ("p2V%ibMult_SP", n), Form ("V_{%i, b}^{x+y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VcMult_SP.push_back (new TProfile2D (Form ("p2V%icMult_SP", n), Form ("V_{%i, c}^{x+y, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VMult_SP.push_back (new TProfile2D (Form ("p2V%iMult_SP", n), Form ("V_{%i}^{SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaMult_SP.push_back (new TProfile2D (Form ("p2V%ixYaMult_SP", n), Form ("V_{%i, a}^{xY, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbMult_SP.push_back (new TProfile2D (Form ("p2V%ixYbMult_SP", n), Form ("V_{%i, b}^{xY, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcMult_SP.push_back (new TProfile2D (Form ("p2V%ixYcMult_SP", n), Form ("V_{%i, c}^{xY, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaMult_SP.push_back (new TProfile2D (Form ("p2V%iyXaMult_SP", n), Form ("V_{%i, a}^{yX, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbMult_SP.push_back (new TProfile2D (Form ("p2V%iyXbMult_SP", n), Form ("V_{%i, b}^{yX, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcMult_SP.push_back (new TProfile2D (Form ("p2V%iyXcMult_SP", n), Form ("V_{%i, c}^{yX, SP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaMult_EP.push_back (new TProfile2D (Form ("p2V%ixaMult_EP", n), Form ("V_{%i, a}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbMult_EP.push_back (new TProfile2D (Form ("p2V%ixbMult_EP", n), Form ("V_{%i, b}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcMult_EP.push_back (new TProfile2D (Form ("p2V%ixcMult_EP", n), Form ("V_{%i, c}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxMult_EP.push_back (new TProfile2D (Form ("p2V%ixMult_EP", n), Form ("V_{%i}^{x, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaMult_EP.push_back (new TProfile2D (Form ("p2V%iyaMult_EP", n), Form ("V_{%i, a}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VybMult_EP.push_back (new TProfile2D (Form ("p2V%iybMult_EP", n), Form ("V_{%i, b}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VycMult_EP.push_back (new TProfile2D (Form ("p2V%iycMult_EP", n), Form ("V_{%i, c}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyMult_EP.push_back (new TProfile2D (Form ("p2V%iyMult_EP", n), Form ("V_{%i}^{y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VaMult_EP.push_back (new TProfile2D (Form ("p2V%iaMult_EP", n), Form ("V_{%i, a}^{x+y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VbMult_EP.push_back (new TProfile2D (Form ("p2V%ibMult_EP", n), Form ("V_{%i, b}^{x+y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VcMult_EP.push_back (new TProfile2D (Form ("p2V%icMult_EP", n), Form ("V_{%i, c}^{x+y, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VMult_EP.push_back (new TProfile2D (Form ("p2V%iMult_EP", n), Form ("V_{%i}^{EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaMult_EP.push_back (new TProfile2D (Form ("p2V%ixYaMult_EP", n), Form ("V_{%i, a}^{xY, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbMult_EP.push_back (new TProfile2D (Form ("p2V%ixYbMult_EP", n), Form ("V_{%i, b}^{xY, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcMult_EP.push_back (new TProfile2D (Form ("p2V%ixYcMult_EP", n), Form ("V_{%i, c}^{xY, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaMult_EP.push_back (new TProfile2D (Form ("p2V%iyXaMult_EP", n), Form ("V_{%i, a}^{yX, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbMult_EP.push_back (new TProfile2D (Form ("p2V%iyXbMult_EP", n), Form ("V_{%i, b}^{yX, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcMult_EP.push_back (new TProfile2D (Form ("p2V%iyXcMult_EP", n), Form ("V_{%i, c}^{yX, EP};mult;sample", n), nBinsMh_, mhMin_, mhMax_, nBinsBS_, 0, nBinsBS_));

            hVxaCent_SP.push_back (new TH1F (Form ("hV%ixaCent_SP", n), Form ("V_{%i, a}^{x, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxbCent_SP.push_back (new TH1F (Form ("hV%ixbCent_SP", n), Form ("V_{%i, b}^{x, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxcCent_SP.push_back (new TH1F (Form ("hV%ixcCent_SP", n), Form ("V_{%i, c}^{x, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxCent_SP.push_back (new TH1F (Form ("hV%ixCent_SP", n), Form ("V_{%i}^{x, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyaCent_SP.push_back (new TH1F (Form ("hV%iyaCent_SP", n), Form ("V_{%i, a}^{y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVybCent_SP.push_back (new TH1F (Form ("hV%iybCent_SP", n), Form ("V_{%i, b}^{y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVycCent_SP.push_back (new TH1F (Form ("hV%iycCent_SP", n), Form ("V_{%i, c}^{y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyCent_SP.push_back (new TH1F (Form ("hV%iyCent_SP", n), Form ("V_{%i}^{y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVaCent_SP.push_back (new TH1F (Form ("hV%iaCent_SP", n), Form ("V_{%i, a}^{x+y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVbCent_SP.push_back (new TH1F (Form ("hV%ibCent_SP", n), Form ("V_{%i, b}^{x+y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVcCent_SP.push_back (new TH1F (Form ("hV%icCent_SP", n), Form ("V_{%i, c}^{x+y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVCent_SP.push_back (new TH1F (Form ("hV%iCent_SP", n), Form ("V_{%i}^{x+y, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYaCent_SP.push_back (new TH1F (Form ("hV%ixYaCent_SP", n), Form ("V_{%i, a}^{xY, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYbCent_SP.push_back (new TH1F (Form ("hV%ixYbCent_SP", n), Form ("V_{%i, b}^{xY, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYcCent_SP.push_back (new TH1F (Form ("hV%ixYcCent_SP", n), Form ("V_{%i, c}^{xY, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXaCent_SP.push_back (new TH1F (Form ("hV%iyXaCent_SP", n), Form ("V_{%i, a}^{yX, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXbCent_SP.push_back (new TH1F (Form ("hV%iyXbCent_SP", n), Form ("V_{%i, b}^{yX, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXcCent_SP.push_back (new TH1F (Form ("hV%iyXcCent_SP", n), Form ("V_{%i, c}^{yX, SP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxaCent_EP.push_back (new TH1F (Form ("hV%ixaCent_EP", n), Form ("V_{%i, a}^{x, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxbCent_EP.push_back (new TH1F (Form ("hV%ixbCent_EP", n), Form ("V_{%i, b}^{x, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxcCent_EP.push_back (new TH1F (Form ("hV%ixcCent_EP", n), Form ("V_{%i, c}^{x, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxCent_EP.push_back (new TH1F (Form ("hV%ixCent_EP", n), Form ("V_{%i}^{x, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyaCent_EP.push_back (new TH1F (Form ("hV%iyaCent_EP", n), Form ("V_{%i, a}^{y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVybCent_EP.push_back (new TH1F (Form ("hV%iybCent_EP", n), Form ("V_{%i, b}^{y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVycCent_EP.push_back (new TH1F (Form ("hV%iycCent_EP", n), Form ("V_{%i, c}^{y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyCent_EP.push_back (new TH1F (Form ("hV%iyCent_EP", n), Form ("V_{%i}^{y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVaCent_EP.push_back (new TH1F (Form ("hV%iaCent_EP", n), Form ("V_{%i, a}^{x+y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVbCent_EP.push_back (new TH1F (Form ("hV%ibCent_EP", n), Form ("V_{%i, b}^{x+y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVcCent_EP.push_back (new TH1F (Form ("hV%icCent_EP", n), Form ("V_{%i}^{x+y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVCent_EP.push_back (new TH1F (Form ("hV%iCent_EP", n), Form ("V_{%i, c}^{x+y, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYaCent_EP.push_back (new TH1F (Form ("hV%ixYaCent_EP", n), Form ("V_{%i, a}^{xY, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYbCent_EP.push_back (new TH1F (Form ("hV%ixYbCent_EP", n), Form ("V_{%i, b}^{xY, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYcCent_EP.push_back (new TH1F (Form ("hV%ixYcCent_EP", n), Form ("V_{%i, c}^{xY, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXaCent_EP.push_back (new TH1F (Form ("hV%iyXaCent_EP", n), Form ("V_{%i, a}^{yX, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXbCent_EP.push_back (new TH1F (Form ("hV%iyXbCent_EP", n), Form ("V_{%i, b}^{yX, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXcCent_EP.push_back (new TH1F (Form ("hV%iyXcCent_EP", n), Form ("V_{%i, c}^{yX, EP};cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));

            hVxaMult_SP.push_back (new TH1F (Form ("hV%ixaMult_SP", n), Form ("V_{%i, a}^{x, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxbMult_SP.push_back (new TH1F (Form ("hV%ixbMult_SP", n), Form ("V_{%i, b}^{x, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxcMult_SP.push_back (new TH1F (Form ("hV%ixcMult_SP", n), Form ("V_{%i, c}^{x, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxMult_SP.push_back (new TH1F (Form ("hV%ixMult_SP", n), Form ("V_{%i}^{x, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyaMult_SP.push_back (new TH1F (Form ("hV%iyaMult_SP", n), Form ("V_{%i, a}^{y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVybMult_SP.push_back (new TH1F (Form ("hV%iybMult_SP", n), Form ("V_{%i, b}^{y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVycMult_SP.push_back (new TH1F (Form ("hV%iycMult_SP", n), Form ("V_{%i, c}^{y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyMult_SP.push_back (new TH1F (Form ("hV%iyMult_SP", n), Form ("V_{%i}^{y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVaMult_SP.push_back (new TH1F (Form ("hV%iaMult_SP", n), Form ("V_{%i, a}^{x+y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVbMult_SP.push_back (new TH1F (Form ("hV%ibMult_SP", n), Form ("V_{%i, b}^{x+y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVcMult_SP.push_back (new TH1F (Form ("hV%icMult_SP", n), Form ("V_{%i, c}^{x+y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVMult_SP.push_back (new TH1F (Form ("hV%iMult_SP", n), Form ("V_{%i}^{x+y, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYaMult_SP.push_back (new TH1F (Form ("hV%ixYaMult_SP", n), Form ("V_{%i, a}^{xY, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYbMult_SP.push_back (new TH1F (Form ("hV%ixYbMult_SP", n), Form ("V_{%i, b}^{xY, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYcMult_SP.push_back (new TH1F (Form ("hV%ixYcMult_SP", n), Form ("V_{%i, c}^{xY, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXaMult_SP.push_back (new TH1F (Form ("hV%iyXaMult_SP", n), Form ("V_{%i, a}^{yX, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXbMult_SP.push_back (new TH1F (Form ("hV%iyXbMult_SP", n), Form ("V_{%i, b}^{yX, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXcMult_SP.push_back (new TH1F (Form ("hV%iyXcMult_SP", n), Form ("V_{%i, c}^{yX, SP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxaMult_EP.push_back (new TH1F (Form ("hV%ixaMult_EP", n), Form ("V_{%i, a}^{x, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxbMult_EP.push_back (new TH1F (Form ("hV%ixbMult_EP", n), Form ("V_{%i, b}^{x, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxcMult_EP.push_back (new TH1F (Form ("hV%ixcMult_EP", n), Form ("V_{%i, c}^{x, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxMult_EP.push_back (new TH1F (Form ("hV%ixMult_EP", n), Form ("V_{%i}^{x, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyaMult_EP.push_back (new TH1F (Form ("hV%iyaMult_EP", n), Form ("V_{%i, a}^{y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVybMult_EP.push_back (new TH1F (Form ("hV%iybMult_EP", n), Form ("V_{%i, b}^{y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVycMult_EP.push_back (new TH1F (Form ("hV%iycMult_EP", n), Form ("V_{%i, c}^{y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyMult_EP.push_back (new TH1F (Form ("hV%iyMult_EP", n), Form ("V_{%i}^{y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVaMult_EP.push_back (new TH1F (Form ("hV%iaMult_EP", n), Form ("V_{%i, a}^{x+y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVbMult_EP.push_back (new TH1F (Form ("hV%ibMult_EP", n), Form ("V_{%i, b}^{x+y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVcMult_EP.push_back (new TH1F (Form ("hV%icMult_EP", n), Form ("V_{%i, c}^{x+y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVMult_EP.push_back (new TH1F (Form ("hV%iMult_EP", n), Form ("V_{%i}^{x+y, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYaMult_EP.push_back (new TH1F (Form ("hV%ixYaMult_EP", n), Form ("V_{%i, a}^{xY, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYbMult_EP.push_back (new TH1F (Form ("hV%ixYbMult_EP", n), Form ("V_{%i, b}^{xY, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYcMult_EP.push_back (new TH1F (Form ("hV%ixYcMult_EP", n), Form ("V_{%i, c}^{xY, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXaMult_EP.push_back (new TH1F (Form ("hV%iyXaMult_EP", n), Form ("V_{%i, a}^{yX, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXbMult_EP.push_back (new TH1F (Form ("hV%iyXbMult_EP", n), Form ("V_{%i, b}^{yX, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXcMult_EP.push_back (new TH1F (Form ("hV%iyXcMult_EP", n), Form ("V_{%i, c}^{yX, EP};mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));

            hVxaCentBS_SP.push_back (new TH1F (Form ("hV%ixaCentBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxbCentBS_SP.push_back (new TH1F (Form ("hV%ixbCentBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxcCentBS_SP.push_back (new TH1F (Form ("hV%ixcCentBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxCentBS_SP.push_back (new TH1F (Form ("hV%ixCentBS_SP", n), Form ("V_{%i}^{x, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyaCentBS_SP.push_back (new TH1F (Form ("hV%iyaCentBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVybCentBS_SP.push_back (new TH1F (Form ("hV%iybCentBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVycCentBS_SP.push_back (new TH1F (Form ("hV%iycCentBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyCentBS_SP.push_back (new TH1F (Form ("hV%iyCentBS_SP", n), Form ("V_{%i}^{y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVaCentBS_SP.push_back (new TH1F (Form ("hV%iaCentBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVbCentBS_SP.push_back (new TH1F (Form ("hV%ibCentBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVcCentBS_SP.push_back (new TH1F (Form ("hV%icCentBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVCentBS_SP.push_back (new TH1F (Form ("hV%iCentBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYaCentBS_SP.push_back (new TH1F (Form ("hV%ixYaCentBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYbCentBS_SP.push_back (new TH1F (Form ("hV%ixYbCentBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYcCentBS_SP.push_back (new TH1F (Form ("hV%ixYcCentBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXaCentBS_SP.push_back (new TH1F (Form ("hV%iyXaCentBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXbCentBS_SP.push_back (new TH1F (Form ("hV%iyXbCentBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXcCentBS_SP.push_back (new TH1F (Form ("hV%iyXcCentBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxaCentBS_EP.push_back (new TH1F (Form ("hV%ixaCentBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxbCentBS_EP.push_back (new TH1F (Form ("hV%ixbCentBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxcCentBS_EP.push_back (new TH1F (Form ("hV%ixcCentBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxCentBS_EP.push_back (new TH1F (Form ("hV%ixCentBS_EP", n), Form ("V_{%i}^{x, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyaCentBS_EP.push_back (new TH1F (Form ("hV%iyaCentBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVybCentBS_EP.push_back (new TH1F (Form ("hV%iybCentBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVycCentBS_EP.push_back (new TH1F (Form ("hV%iycCentBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyCentBS_EP.push_back (new TH1F (Form ("hV%iyCentBS_EP", n), Form ("V_{%i}^{y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVaCentBS_EP.push_back (new TH1F (Form ("hV%iaCentBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVbCentBS_EP.push_back (new TH1F (Form ("hV%ibCentBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVcCentBS_EP.push_back (new TH1F (Form ("hV%icCentBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVCentBS_EP.push_back (new TH1F (Form ("hV%iCentBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYaCentBS_EP.push_back (new TH1F (Form ("hV%ixYaCentBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYbCentBS_EP.push_back (new TH1F (Form ("hV%ixYbCentBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVxYcCentBS_EP.push_back (new TH1F (Form ("hV%ixYcCentBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXaCentBS_EP.push_back (new TH1F (Form ("hV%iyXaCentBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXbCentBS_EP.push_back (new TH1F (Form ("hV%iyXbCentBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));
            hVyXcCentBS_EP.push_back (new TH1F (Form ("hV%iyXcCentBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling);cent;V_{%i}", n, n), nBinsCent_, centMin_, centMax_));

            hVxaMultBS_SP.push_back (new TH1F (Form ("hV%ixaMultBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxbMultBS_SP.push_back (new TH1F (Form ("hV%ixbMultBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxcMultBS_SP.push_back (new TH1F (Form ("hV%ixcMultBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxMultBS_SP.push_back (new TH1F (Form ("hV%ixMultBS_SP", n), Form ("V_{%i}^{x, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyaMultBS_SP.push_back (new TH1F (Form ("hV%iyaMultBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVybMultBS_SP.push_back (new TH1F (Form ("hV%iybMultBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVycMultBS_SP.push_back (new TH1F (Form ("hV%iycMultBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyMultBS_SP.push_back (new TH1F (Form ("hV%iyMultBS_SP", n), Form ("V_{%i}^{y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVaMultBS_SP.push_back (new TH1F (Form ("hV%iaMultBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVbMultBS_SP.push_back (new TH1F (Form ("hV%ibMultBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVcMultBS_SP.push_back (new TH1F (Form ("hV%icMultBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVMultBS_SP.push_back (new TH1F (Form ("hV%iMultBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYaMultBS_SP.push_back (new TH1F (Form ("hV%ixYaMultBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYbMultBS_SP.push_back (new TH1F (Form ("hV%ixYbMultBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYcMultBS_SP.push_back (new TH1F (Form ("hV%ixYcMultBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXaMultBS_SP.push_back (new TH1F (Form ("hV%iyXaMultBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXbMultBS_SP.push_back (new TH1F (Form ("hV%iyXbMultBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXcMultBS_SP.push_back (new TH1F (Form ("hV%iyXcMultBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxaMultBS_EP.push_back (new TH1F (Form ("hV%ixaMultBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxbMultBS_EP.push_back (new TH1F (Form ("hV%ixbMultBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxcMultBS_EP.push_back (new TH1F (Form ("hV%ixcMultBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxMultBS_EP.push_back (new TH1F (Form ("hV%ixMultBS_EP", n), Form ("V_{%i}^{x, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyaMultBS_EP.push_back (new TH1F (Form ("hV%iyaMultBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVybMultBS_EP.push_back (new TH1F (Form ("hV%iybMultBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVycMultBS_EP.push_back (new TH1F (Form ("hV%iycMultBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyMultBS_EP.push_back (new TH1F (Form ("hV%iyMultBS_EP", n), Form ("V_{%i}^{y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVaMultBS_EP.push_back (new TH1F (Form ("hV%iaMultBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVbMultBS_EP.push_back (new TH1F (Form ("hV%ibMultBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVcMultBS_EP.push_back (new TH1F (Form ("hV%icMultBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVMultBS_EP.push_back (new TH1F (Form ("hV%iMultBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYaMultBS_EP.push_back (new TH1F (Form ("hV%ixYaMultBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYbMultBS_EP.push_back (new TH1F (Form ("hV%ixYbMultBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVxYcMultBS_EP.push_back (new TH1F (Form ("hV%ixYcMultBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXaMultBS_EP.push_back (new TH1F (Form ("hV%iyXaMultBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXbMultBS_EP.push_back (new TH1F (Form ("hV%iyXbMultBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));
            hVyXcMultBS_EP.push_back (new TH1F (Form ("hV%iyXcMultBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling);mult;V_{%i}", n, n), nBinsMh_, mhMin_, mhMax_));

            p2VxaPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixaPtCent_SP", n), Form ("V_{%i, a}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixbPtCent_SP", n), Form ("V_{%i, b}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixcPtCent_SP", n), Form ("V_{%i, c}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixPtCent_SP", n), Form ("V_{%i}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaPtCent_SP.push_back (new TProfile2D (Form ("p2V%iyaPtCent_SP", n), Form ("V_{%i, a}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VybPtCent_SP.push_back (new TProfile2D (Form ("p2V%iybPtCent_SP", n), Form ("V_{%i, b}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VycPtCent_SP.push_back (new TProfile2D (Form ("p2V%iycPtCent_SP", n), Form ("V_{%i, c}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyPtCent_SP.push_back (new TProfile2D (Form ("p2V%iyPtCent_SP", n), Form ("V_{%i}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VaPtCent_SP.push_back (new TProfile2D (Form ("p2V%iaPtCent_SP", n), Form ("V_{%i, a}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VbPtCent_SP.push_back (new TProfile2D (Form ("p2V%ibPtCent_SP", n), Form ("V_{%i, b}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VcPtCent_SP.push_back (new TProfile2D (Form ("p2V%icPtCent_SP", n), Form ("V_{%i, c}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VPtCent_SP.push_back (new TProfile2D (Form ("p2V%iPtCent_SP", n), Form ("V_{%i}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixYaPtCent_SP", n), Form ("V_{%i, a}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixYbPtCent_SP", n), Form ("V_{%i, b}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcPtCent_SP.push_back (new TProfile2D (Form ("p2V%ixYcPtCent_SP", n), Form ("V_{%i, c}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaPtCent_SP.push_back (new TProfile2D (Form ("p2V%iyXaPtCent_SP", n), Form ("V_{%i, a}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbPtCent_SP.push_back (new TProfile2D (Form ("p2V%iyXbPtCent_SP", n), Form ("V_{%i, b}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcPtCent_SP.push_back (new TProfile2D (Form ("p2V%iyXcPtCent_SP", n), Form ("V_{%i, c}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixaPtCent_EP", n), Form ("V_{%i, a}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixbPtCent_EP", n), Form ("V_{%i, b}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixcPtCent_EP", n), Form ("V_{%i, c}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixPtCent_EP", n), Form ("V_{%i}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaPtCent_EP.push_back (new TProfile2D (Form ("p2V%iyaPtCent_EP", n), Form ("V_{%i, a}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VybPtCent_EP.push_back (new TProfile2D (Form ("p2V%iybPtCent_EP", n), Form ("V_{%i, b}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VycPtCent_EP.push_back (new TProfile2D (Form ("p2V%iycPtCent_EP", n), Form ("V_{%i, c}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyPtCent_EP.push_back (new TProfile2D (Form ("p2V%iyPtCent_EP", n), Form ("V_{%i}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VaPtCent_EP.push_back (new TProfile2D (Form ("p2V%iaPtCent_EP", n), Form ("V_{%i, a}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VbPtCent_EP.push_back (new TProfile2D (Form ("p2V%ibPtCent_EP", n), Form ("V_{%i, b}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VcPtCent_EP.push_back (new TProfile2D (Form ("p2V%icPtCent_EP", n), Form ("V_{%i, c}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VPtCent_EP.push_back (new TProfile2D (Form ("p2V%iPtCent_EP", n), Form ("V_{%i}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixYaPtCent_EP", n), Form ("V_{%i, a}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixYbPtCent_EP", n), Form ("V_{%i, b}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcPtCent_EP.push_back (new TProfile2D (Form ("p2V%ixYcPtCent_EP", n), Form ("V_{%i, c}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaPtCent_EP.push_back (new TProfile2D (Form ("p2V%iyXaPtCent_EP", n), Form ("V_{%i, a}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbPtCent_EP.push_back (new TProfile2D (Form ("p2V%iyXbPtCent_EP", n), Form ("V_{%i, b}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcPtCent_EP.push_back (new TProfile2D (Form ("p2V%iyXcPtCent_EP", n), Form ("V_{%i, c}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));

            p2VxaEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixaEtaCent_SP", n), Form ("V_{%i, a}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixbEtaCent_SP", n), Form ("V_{%i, b}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixcEtaCent_SP", n), Form ("V_{%i, c}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixEtaCent_SP", n), Form ("V_{%i}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iyaEtaCent_SP", n), Form ("V_{%i, a}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VybEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iybEtaCent_SP", n), Form ("V_{%i, b}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VycEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iycEtaCent_SP", n), Form ("V_{%i, c}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iyEtaCent_SP", n), Form ("V_{%i}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VaEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iaEtaCent_SP", n), Form ("V_{%i, a}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VbEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ibEtaCent_SP", n), Form ("V_{%i, b}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VcEtaCent_SP.push_back (new TProfile2D (Form ("p2V%icEtaCent_SP", n), Form ("V_{%i, c}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iEtaCent_SP", n), Form ("V_{%i}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixYaEtaCent_SP", n), Form ("V_{%i, a}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixYbEtaCent_SP", n), Form ("V_{%i, b}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcEtaCent_SP.push_back (new TProfile2D (Form ("p2V%ixYcEtaCent_SP", n), Form ("V_{%i, c}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iyXaEtaCent_SP", n), Form ("V_{%i, a}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iyXbEtaCent_SP", n), Form ("V_{%i, b}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcEtaCent_SP.push_back (new TProfile2D (Form ("p2V%iyXcEtaCent_SP", n), Form ("V_{%i, c}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixaEtaCent_EP", n), Form ("V_{%i, a}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixbEtaCent_EP", n), Form ("V_{%i, b}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixcEtaCent_EP", n), Form ("V_{%i, c}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixEtaCent_EP", n), Form ("V_{%i}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iyaEtaCent_EP", n), Form ("V_{%i, a}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VybEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iybEtaCent_EP", n), Form ("V_{%i, b}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VycEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iycEtaCent_EP", n), Form ("V_{%i, c}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iyEtaCent_EP", n), Form ("V_{%i}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VaEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iaEtaCent_EP", n), Form ("V_{%i, a}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VbEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ibEtaCent_EP", n), Form ("V_{%i, b}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VcEtaCent_EP.push_back (new TProfile2D (Form ("p2V%icEtaCent_EP", n), Form ("V_{%i, c}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iEtaCent_EP", n), Form ("V_{%i}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixYaEtaCent_EP", n), Form ("V_{%i, a}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixYbEtaCent_EP", n), Form ("V_{%i, b}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcEtaCent_EP.push_back (new TProfile2D (Form ("p2V%ixYcEtaCent_EP", n), Form ("V_{%i, c}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iyXaEtaCent_EP", n), Form ("V_{%i, a}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iyXbEtaCent_EP", n), Form ("V_{%i, b}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcEtaCent_EP.push_back (new TProfile2D (Form ("p2V%iyXcEtaCent_EP", n), Form ("V_{%i, c}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));

            p2VxaPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixaPtMult_SP", n), Form ("V_{%i, a}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixbPtMult_SP", n), Form ("V_{%i, b}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixcPtMult_SP", n), Form ("V_{%i, c}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixPtMult_SP", n), Form ("V_{%i}^{x, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaPtMult_SP.push_back (new TProfile2D (Form ("p2V%iyaPtMult_SP", n), Form ("V_{%i, a}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VybPtMult_SP.push_back (new TProfile2D (Form ("p2V%iybPtMult_SP", n), Form ("V_{%i, b}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VycPtMult_SP.push_back (new TProfile2D (Form ("p2V%iycPtMult_SP", n), Form ("V_{%i, c}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyPtMult_SP.push_back (new TProfile2D (Form ("p2V%iyPtMult_SP", n), Form ("V_{%i}^{y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VaPtMult_SP.push_back (new TProfile2D (Form ("p2V%iaPtMult_SP", n), Form ("V_{%i, a}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VbPtMult_SP.push_back (new TProfile2D (Form ("p2V%ibPtMult_SP", n), Form ("V_{%i, b}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VcPtMult_SP.push_back (new TProfile2D (Form ("p2V%icPtMult_SP", n), Form ("V_{%i, c}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VPtMult_SP.push_back (new TProfile2D (Form ("p2V%iPtMult_SP", n), Form ("V_{%i}^{x+y, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixYaPtMult_SP", n), Form ("V_{%i, a}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixYbPtMult_SP", n), Form ("V_{%i, b}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcPtMult_SP.push_back (new TProfile2D (Form ("p2V%ixYcPtMult_SP", n), Form ("V_{%i, c}^{xY, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaPtMult_SP.push_back (new TProfile2D (Form ("p2V%iyXaPtMult_SP", n), Form ("V_{%i, a}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbPtMult_SP.push_back (new TProfile2D (Form ("p2V%iyXbPtMult_SP", n), Form ("V_{%i, b}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcPtMult_SP.push_back (new TProfile2D (Form ("p2V%iyXcPtMult_SP", n), Form ("V_{%i, c}^{yX, SP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixaPtMult_EP", n), Form ("V_{%i, a}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixbPtMult_EP", n), Form ("V_{%i, b}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixcPtMult_EP", n), Form ("V_{%i, c}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixPtMult_EP", n), Form ("V_{%i}^{x, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaPtMult_EP.push_back (new TProfile2D (Form ("p2V%iyaPtMult_EP", n), Form ("V_{%i, a}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VybPtMult_EP.push_back (new TProfile2D (Form ("p2V%iybPtMult_EP", n), Form ("V_{%i, b}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VycPtMult_EP.push_back (new TProfile2D (Form ("p2V%iycPtMult_EP", n), Form ("V_{%i, c}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyPtMult_EP.push_back (new TProfile2D (Form ("p2V%iyPtMult_EP", n), Form ("V_{%i}^{y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VaPtMult_EP.push_back (new TProfile2D (Form ("p2V%iaPtMult_EP", n), Form ("V_{%i, a}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VbPtMult_EP.push_back (new TProfile2D (Form ("p2V%ibPtMult_EP", n), Form ("V_{%i, b}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VcPtMult_EP.push_back (new TProfile2D (Form ("p2V%icPtMult_EP", n), Form ("V_{%i, c}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VPtMult_EP.push_back (new TProfile2D (Form ("p2V%iPtMult_EP", n), Form ("V_{%i}^{x+y, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixYaPtMult_EP", n), Form ("V_{%i, a}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixYbPtMult_EP", n), Form ("V_{%i, b}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcPtMult_EP.push_back (new TProfile2D (Form ("p2V%ixYcPtMult_EP", n), Form ("V_{%i, c}^{xY, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaPtMult_EP.push_back (new TProfile2D (Form ("p2V%iyXaPtMult_EP", n), Form ("V_{%i, a}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbPtMult_EP.push_back (new TProfile2D (Form ("p2V%iyXbPtMult_EP", n), Form ("V_{%i, b}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcPtMult_EP.push_back (new TProfile2D (Form ("p2V%iyXcPtMult_EP", n), Form ("V_{%i, c}^{yX, EP}; P_{T} [GeV/c];sample", n), nBinsPt_, ptMin_, ptMax_, nBinsBS_, 0, nBinsBS_));

            p2VxaEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixaEtaMult_SP", n), Form ("V_{%i, a}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixbEtaMult_SP", n), Form ("V_{%i, b}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixcEtaMult_SP", n), Form ("V_{%i, c}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixEtaMult_SP", n), Form ("V_{%i}^{x, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iyaEtaMult_SP", n), Form ("V_{%i, a}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VybEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iybEtaMult_SP", n), Form ("V_{%i, b}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VycEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iycEtaMult_SP", n), Form ("V_{%i, c}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iyEtaMult_SP", n), Form ("V_{%i}^{y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VaEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iaEtaMult_SP", n), Form ("V_{%i, a}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VbEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ibEtaMult_SP", n), Form ("V_{%i, b}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VcEtaMult_SP.push_back (new TProfile2D (Form ("p2V%icEtaMult_SP", n), Form ("V_{%i, c}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iEtaMult_SP", n), Form ("V_{%i}^{x+y, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixYaEtaMult_SP", n), Form ("V_{%i, a}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixYbEtaMult_SP", n), Form ("V_{%i, b}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcEtaMult_SP.push_back (new TProfile2D (Form ("p2V%ixYcEtaMult_SP", n), Form ("V_{%i, c}^{xY, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iyXaEtaMult_SP", n), Form ("V_{%i, a}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iyXbEtaMult_SP", n), Form ("V_{%i, b}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcEtaMult_SP.push_back (new TProfile2D (Form ("p2V%iyXcEtaMult_SP", n), Form ("V_{%i, c}^{yX, SP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxaEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixaEtaMult_EP", n), Form ("V_{%i, a}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxbEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixbEtaMult_EP", n), Form ("V_{%i, b}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxcEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixcEtaMult_EP", n), Form ("V_{%i, c}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixEtaMult_EP", n), Form ("V_{%i}^{x, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyaEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iyaEtaMult_EP", n), Form ("V_{%i, a}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VybEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iybEtaMult_EP", n), Form ("V_{%i, b}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VycEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iycEtaMult_EP", n), Form ("V_{%i, c}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iyEtaMult_EP", n), Form ("V_{%i}^{y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VaEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iaEtaMult_EP", n), Form ("V_{%i, a}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VbEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ibEtaMult_EP", n), Form ("V_{%i, b}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VcEtaMult_EP.push_back (new TProfile2D (Form ("p2V%icEtaMult_EP", n), Form ("V_{%i, c}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iEtaMult_EP", n), Form ("V_{%i}^{x+y, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYaEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixYaEtaMult_EP", n), Form ("V_{%i, a}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYbEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixYbEtaMult_EP", n), Form ("V_{%i, b}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VxYcEtaMult_EP.push_back (new TProfile2D (Form ("p2V%ixYcEtaMult_EP", n), Form ("V_{%i, c}^{xY, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXaEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iyXaEtaMult_EP", n), Form ("V_{%i, a}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXbEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iyXbEtaMult_EP", n), Form ("V_{%i, b}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));
            p2VyXcEtaMult_EP.push_back (new TProfile2D (Form ("p2V%iyXcEtaMult_EP", n), Form ("V_{%i, c}^{yX, EP};", n) + varName_ + ";sample", nBinsEta_, etaMin_, etaMax_, nBinsBS_, 0, nBinsBS_));

            hVxaPtCent_SP.push_back (new TH1F (Form ("hV%ixaPtCent_SP", n), Form ("V_{%i, a}^{x, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtCent_SP.push_back (new TH1F (Form ("hV%ixbPtCent_SP", n), Form ("V_{%i, b}^{x, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtCent_SP.push_back (new TH1F (Form ("hV%ixcPtCent_SP", n), Form ("V_{%i, c}^{x, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtCent_SP.push_back (new TH1F (Form ("hV%ixPtCent_SP", n), Form ("V_{%i}^{x, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtCent_SP.push_back (new TH1F (Form ("hV%iyaPtCent_SP", n), Form ("V_{%i, a}^{y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtCent_SP.push_back (new TH1F (Form ("hV%iybPtCent_SP", n), Form ("V_{%i, b}^{y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtCent_SP.push_back (new TH1F (Form ("hV%iycPtCent_SP", n), Form ("V_{%i, c}^{y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtCent_SP.push_back (new TH1F (Form ("hV%iyPtCent_SP", n), Form ("V_{%i}^{y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtCent_SP.push_back (new TH1F (Form ("hV%iaPtCent_SP", n), Form ("V_{%i, a}^{x+y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtCent_SP.push_back (new TH1F (Form ("hV%ibPtCent_SP", n), Form ("V_{%i, b}^{x+y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtCent_SP.push_back (new TH1F (Form ("hV%icPtCent_SP", n), Form ("V_{%i, c}^{x+y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtCent_SP.push_back (new TH1F (Form ("hV%iPtCent_SP", n), Form ("V_{%i}^{x+y, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtCent_SP.push_back (new TH1F (Form ("hV%ixYaPtCent_SP", n), Form ("V_{%i, a}^{xY, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtCent_SP.push_back (new TH1F (Form ("hV%ixYbPtCent_SP", n), Form ("V_{%i, b}^{xY, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtCent_SP.push_back (new TH1F (Form ("hV%ixYcPtCent_SP", n), Form ("V_{%i, c}^{xY, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtCent_SP.push_back (new TH1F (Form ("hV%ixYPtCent_SP", n), Form ("V_{%i}^{xY, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtCent_SP.push_back (new TH1F (Form ("hV%iyXaPtCent_SP", n), Form ("V_{%i, a}^{yX, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtCent_SP.push_back (new TH1F (Form ("hV%iyXbPtCent_SP", n), Form ("V_{%i, b}^{yX, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtCent_SP.push_back (new TH1F (Form ("hV%iyXcPtCent_SP", n), Form ("V_{%i, c}^{yX, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtCent_SP.push_back (new TH1F (Form ("hV%iyXPtCent_SP", n), Form ("V_{%i}^{yX, SP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxaPtCent_EP.push_back (new TH1F (Form ("hV%ixaPtCent_EP", n), Form ("V_{%i, a}^{x, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtCent_EP.push_back (new TH1F (Form ("hV%ixbPtCent_EP", n), Form ("V_{%i, b}^{x, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtCent_EP.push_back (new TH1F (Form ("hV%ixcPtCent_EP", n), Form ("V_{%i, c}^{x, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtCent_EP.push_back (new TH1F (Form ("hV%ixPtCent_EP", n), Form ("V_{%i}^{x, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtCent_EP.push_back (new TH1F (Form ("hV%iyaPtCent_EP", n), Form ("V_{%i, a}^{y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtCent_EP.push_back (new TH1F (Form ("hV%iybPtCent_EP", n), Form ("V_{%i, b}^{y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtCent_EP.push_back (new TH1F (Form ("hV%iycPtCent_EP", n), Form ("V_{%i, c}^{y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtCent_EP.push_back (new TH1F (Form ("hV%iyPtCent_EP", n), Form ("V_{%i}^{y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtCent_EP.push_back (new TH1F (Form ("hV%iaPtCent_EP", n), Form ("V_{%i, a}^{x+y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtCent_EP.push_back (new TH1F (Form ("hV%ibPtCent_EP", n), Form ("V_{%i, b}^{x+y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtCent_EP.push_back (new TH1F (Form ("hV%icPtCent_EP", n), Form ("V_{%i, c}^{x+y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtCent_EP.push_back (new TH1F (Form ("hV%iPtCent_EP", n), Form ("V_{%i}^{x+y, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtCent_EP.push_back (new TH1F (Form ("hV%ixYaPtCent_EP", n), Form ("V_{%i, a}^{xY, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtCent_EP.push_back (new TH1F (Form ("hV%ixYbPtCent_EP", n), Form ("V_{%i, b}^{xY, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtCent_EP.push_back (new TH1F (Form ("hV%ixYcPtCent_EP", n), Form ("V_{%i, c}^{xY, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtCent_EP.push_back (new TH1F (Form ("hV%ixYPtCent_EP", n), Form ("V_{%i}^{xY, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtCent_EP.push_back (new TH1F (Form ("hV%iyXaPtCent_EP", n), Form ("V_{%i, a}^{yX, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtCent_EP.push_back (new TH1F (Form ("hV%iyXbPtCent_EP", n), Form ("V_{%i, b}^{yX, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtCent_EP.push_back (new TH1F (Form ("hV%iyXcPtCent_EP", n), Form ("V_{%i, c}^{yX, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtCent_EP.push_back (new TH1F (Form ("hV%iyXPtCent_EP", n), Form ("V_{%i}^{yX, EP} (over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));

            hVxaPtMult_SP.push_back (new TH1F (Form ("hV%ixaPtMult_SP", n), Form ("V_{%i, a}^{x, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtMult_SP.push_back (new TH1F (Form ("hV%ixbPtMult_SP", n), Form ("V_{%i, b}^{x, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtMult_SP.push_back (new TH1F (Form ("hV%ixcPtMult_SP", n), Form ("V_{%i, c}^{x, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtMult_SP.push_back (new TH1F (Form ("hV%ixPtMult_SP", n), Form ("V_{%i}^{x, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtMult_SP.push_back (new TH1F (Form ("hV%iyaPtMult_SP", n), Form ("V_{%i, a}^{y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtMult_SP.push_back (new TH1F (Form ("hV%iybPtMult_SP", n), Form ("V_{%i, b}^{y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtMult_SP.push_back (new TH1F (Form ("hV%iycPtMult_SP", n), Form ("V_{%i, c}^{y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtMult_SP.push_back (new TH1F (Form ("hV%iyPtMult_SP", n), Form ("V_{%i}^{y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtMult_SP.push_back (new TH1F (Form ("hV%iaPtMult_SP", n), Form ("V_{%i, a}^{x+y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtMult_SP.push_back (new TH1F (Form ("hV%ibPtMult_SP", n), Form ("V_{%i, b}^{x+y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtMult_SP.push_back (new TH1F (Form ("hV%icPtMult_SP", n), Form ("V_{%i, c}^{x+y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtMult_SP.push_back (new TH1F (Form ("hV%iPtMult_SP", n), Form ("V_{%i}^{x+y, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtMult_SP.push_back (new TH1F (Form ("hV%ixYaPtMult_SP", n), Form ("V_{%i, a}^{xY, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtMult_SP.push_back (new TH1F (Form ("hV%ixYbPtMult_SP", n), Form ("V_{%i, b}^{xY, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtMult_SP.push_back (new TH1F (Form ("hV%ixYcPtMult_SP", n), Form ("V_{%i, c}^{xY, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtMult_SP.push_back (new TH1F (Form ("hV%ixYPtMult_SP", n), Form ("V_{%i}^{xY, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtMult_SP.push_back (new TH1F (Form ("hV%iyXaPtMult_SP", n), Form ("V_{%i, a}^{yX, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtMult_SP.push_back (new TH1F (Form ("hV%iyXbPtMult_SP", n), Form ("V_{%i, b}^{yX, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtMult_SP.push_back (new TH1F (Form ("hV%iyXcPtMult_SP", n), Form ("V_{%i, c}^{yX, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtMult_SP.push_back (new TH1F (Form ("hV%iyXPtMult_SP", n), Form ("V_{%i}^{yX, SP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxaPtMult_EP.push_back (new TH1F (Form ("hV%ixaPtMult_EP", n), Form ("V_{%i, a}^{x, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtMult_EP.push_back (new TH1F (Form ("hV%ixbPtMult_EP", n), Form ("V_{%i, b}^{x, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtMult_EP.push_back (new TH1F (Form ("hV%ixcPtMult_EP", n), Form ("V_{%i, c}^{x, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtMult_EP.push_back (new TH1F (Form ("hV%ixPtMult_EP", n), Form ("V_{%i}^{x, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtMult_EP.push_back (new TH1F (Form ("hV%iyaPtMult_EP", n), Form ("V_{%i, a}^{y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtMult_EP.push_back (new TH1F (Form ("hV%iybPtMult_EP", n), Form ("V_{%i, b}^{y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtMult_EP.push_back (new TH1F (Form ("hV%iycPtMult_EP", n), Form ("V_{%i, c}^{y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtMult_EP.push_back (new TH1F (Form ("hV%iyPtMult_EP", n), Form ("V_{%i}^{y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtMult_EP.push_back (new TH1F (Form ("hV%iaPtMult_EP", n), Form ("V_{%i, a}^{x+y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtMult_EP.push_back (new TH1F (Form ("hV%ibPtMult_EP", n), Form ("V_{%i, b}^{x+y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtMult_EP.push_back (new TH1F (Form ("hV%icPtMult_EP", n), Form ("V_{%i, c}^{x+y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtMult_EP.push_back (new TH1F (Form ("hV%iPtMult_EP", n), Form ("V_{%i}^{x+y, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtMult_EP.push_back (new TH1F (Form ("hV%ixYaPtMult_EP", n), Form ("V_{%i, a}^{xY, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtMult_EP.push_back (new TH1F (Form ("hV%ixYbPtMult_EP", n), Form ("V_{%i, b}^{xY, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtMult_EP.push_back (new TH1F (Form ("hV%ixYcPtMult_EP", n), Form ("V_{%i, c}^{xY, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtMult_EP.push_back (new TH1F (Form ("hV%ixYPtMult_EP", n), Form ("V_{%i}^{xY, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtMult_EP.push_back (new TH1F (Form ("hV%iyXaPtMult_EP", n), Form ("V_{%i, a}^{yX, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtMult_EP.push_back (new TH1F (Form ("hV%iyXbPtMult_EP", n), Form ("V_{%i, b}^{yX, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtMult_EP.push_back (new TH1F (Form ("hV%iyXcPtMult_EP", n), Form ("V_{%i, c}^{yX, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtMult_EP.push_back (new TH1F (Form ("hV%iyXPtMult_EP", n), Form ("V_{%i}^{yX, EP} (over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));

            hVxaEtaCent_SP.push_back (new TH1F (Form ("hV%ixaEtaCent_SP", n), Form ("V_{%i, a}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaCent_SP.push_back (new TH1F (Form ("hV%ixbEtaCent_SP", n), Form ("V_{%i, b}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaCent_SP.push_back (new TH1F (Form ("hV%ixcEtaCent_SP", n), Form ("V_{%i, c}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaCent_SP.push_back (new TH1F (Form ("hV%ixEtaCent_SP", n), Form ("V_{%i}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaCent_SP.push_back (new TH1F (Form ("hV%iyaEtaCent_SP", n), Form ("V_{%i, a}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaCent_SP.push_back (new TH1F (Form ("hV%iybEtaCent_SP", n), Form ("V_{%i, b}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaCent_SP.push_back (new TH1F (Form ("hV%iycEtaCent_SP", n), Form ("V_{%i, c}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaCent_SP.push_back (new TH1F (Form ("hV%iyEtaCent_SP", n), Form ("V_{%i}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaCent_SP.push_back (new TH1F (Form ("hV%iaEtaCent_SP", n), Form ("V_{%i, a}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaCent_SP.push_back (new TH1F (Form ("hV%ibEtaCent_SP", n), Form ("V_{%i, b}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaCent_SP.push_back (new TH1F (Form ("hV%icEtaCent_SP", n), Form ("V_{%i, c}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaCent_SP.push_back (new TH1F (Form ("hV%iEtaCent_SP", n), Form ("V_{%i}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaCent_SP.push_back (new TH1F (Form ("hV%ixYaEtaCent_SP", n), Form ("V_{%i, a}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaCent_SP.push_back (new TH1F (Form ("hV%ixYbEtaCent_SP", n), Form ("V_{%i, b}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaCent_SP.push_back (new TH1F (Form ("hV%ixYcEtaCent_SP", n), Form ("V_{%i, c}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaCent_SP.push_back (new TH1F (Form ("hV%ixYEtaCent_SP", n), Form ("V_{%i}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaCent_SP.push_back (new TH1F (Form ("hV%iyXaEtaCent_SP", n), Form ("V_{%i, a}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaCent_SP.push_back (new TH1F (Form ("hV%iyXbEtaCent_SP", n), Form ("V_{%i, b}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaCent_SP.push_back (new TH1F (Form ("hV%iyXcEtaCent_SP", n), Form ("V_{%i, }^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaCent_SP.push_back (new TH1F (Form ("hV%iyXEtaCent_SP", n), Form ("V_{%i}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVxaEtaCent_EP.push_back (new TH1F (Form ("hV%ixaEtaCent_EP", n), Form ("V_{%i, a}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaCent_EP.push_back (new TH1F (Form ("hV%ixbEtaCent_EP", n), Form ("V_{%i, b}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaCent_EP.push_back (new TH1F (Form ("hV%ixcEtaCent_EP", n), Form ("V_{%i, c}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaCent_EP.push_back (new TH1F (Form ("hV%ixEtaCent_EP", n), Form ("V_{%i}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaCent_EP.push_back (new TH1F (Form ("hV%iyaEtaCent_EP", n), Form ("V_{%i, a}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaCent_EP.push_back (new TH1F (Form ("hV%iybEtaCent_EP", n), Form ("V_{%i, b}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaCent_EP.push_back (new TH1F (Form ("hV%iycEtaCent_EP", n), Form ("V_{%i, c}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaCent_EP.push_back (new TH1F (Form ("hV%iyEtaCent_EP", n), Form ("V_{%i}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaCent_EP.push_back (new TH1F (Form ("hV%iaEtaCent_EP", n), Form ("V_{%i, a}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaCent_EP.push_back (new TH1F (Form ("hV%ibEtaCent_EP", n), Form ("V_{%i, b}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaCent_EP.push_back (new TH1F (Form ("hV%icEtaCent_EP", n), Form ("V_{%i, c}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaCent_EP.push_back (new TH1F (Form ("hV%iEtaCent_EP", n), Form ("V_{%i}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaCent_EP.push_back (new TH1F (Form ("hV%ixYaEtaCent_EP", n), Form ("V_{%i, a}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaCent_EP.push_back (new TH1F (Form ("hV%ixYbEtaCent_EP", n), Form ("V_{%i, b}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaCent_EP.push_back (new TH1F (Form ("hV%ixYcEtaCent_EP", n), Form ("V_{%i, c}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaCent_EP.push_back (new TH1F (Form ("hV%ixYEtaCent_EP", n), Form ("V_{%i}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaCent_EP.push_back (new TH1F (Form ("hV%iyXaEtaCent_EP", n), Form ("V_{%i, a}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaCent_EP.push_back (new TH1F (Form ("hV%iyXbEtaCent_EP", n), Form ("V_{%i, b}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaCent_EP.push_back (new TH1F (Form ("hV%iyXcEtaCent_EP", n), Form ("V_{%i, }^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaCent_EP.push_back (new TH1F (Form ("hV%iyXEtaCent_EP", n), Form ("V_{%i}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));

            hVxaEtaMult_SP.push_back (new TH1F (Form ("hV%ixaEtaMult_SP", n), Form ("V_{%i, a}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaMult_SP.push_back (new TH1F (Form ("hV%ixbEtaMult_SP", n), Form ("V_{%i, b}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaMult_SP.push_back (new TH1F (Form ("hV%ixcEtaMult_SP", n), Form ("V_{%i, c}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaMult_SP.push_back (new TH1F (Form ("hV%ixEtaMult_SP", n), Form ("V_{%i}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaMult_SP.push_back (new TH1F (Form ("hV%iyaEtaMult_SP", n), Form ("V_{%i, a}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaMult_SP.push_back (new TH1F (Form ("hV%iybEtaMult_SP", n), Form ("V_{%i, b}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaMult_SP.push_back (new TH1F (Form ("hV%iycEtaMult_SP", n), Form ("V_{%i, c}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaMult_SP.push_back (new TH1F (Form ("hV%iyEtaMult_SP", n), Form ("V_{%i}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaMult_SP.push_back (new TH1F (Form ("hV%iaEtaMult_SP", n), Form ("V_{%i, a}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaMult_SP.push_back (new TH1F (Form ("hV%ibEtaMult_SP", n), Form ("V_{%i, b}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaMult_SP.push_back (new TH1F (Form ("hV%icEtaMult_SP", n), Form ("V_{%i, c}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaMult_SP.push_back (new TH1F (Form ("hV%iEtaMult_SP", n), Form ("V_{%i}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaMult_SP.push_back (new TH1F (Form ("hV%ixYaEtaMult_SP", n), Form ("V_{%i, a}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaMult_SP.push_back (new TH1F (Form ("hV%ixYbEtaMult_SP", n), Form ("V_{%i, b}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaMult_SP.push_back (new TH1F (Form ("hV%ixYcEtaMult_SP", n), Form ("V_{%i, c}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaMult_SP.push_back (new TH1F (Form ("hV%ixYEtaMult_SP", n), Form ("V_{%i}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaMult_SP.push_back (new TH1F (Form ("hV%iyXaEtaMult_SP", n), Form ("V_{%i, a}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaMult_SP.push_back (new TH1F (Form ("hV%iyXbEtaMult_SP", n), Form ("V_{%i, b}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaMult_SP.push_back (new TH1F (Form ("hV%iyXcEtaMult_SP", n), Form ("V_{%i, c}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaMult_SP.push_back (new TH1F (Form ("hV%iyXEtaMult_SP", n), Form ("V_{%i}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVxaEtaMult_EP.push_back (new TH1F (Form ("hV%ixaEtaMult_EP", n), Form ("V_{%i, a}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaMult_EP.push_back (new TH1F (Form ("hV%ixbEtaMult_EP", n), Form ("V_{%i, b}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaMult_EP.push_back (new TH1F (Form ("hV%ixcEtaMult_EP", n), Form ("V_{%i, c}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaMult_EP.push_back (new TH1F (Form ("hV%ixEtaMult_EP", n), Form ("V_{%i}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaMult_EP.push_back (new TH1F (Form ("hV%iyaEtaMult_EP", n), Form ("V_{%i, a}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaMult_EP.push_back (new TH1F (Form ("hV%iybEtaMult_EP", n), Form ("V_{%i, b}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaMult_EP.push_back (new TH1F (Form ("hV%iycEtaMult_EP", n), Form ("V_{%i, c}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaMult_EP.push_back (new TH1F (Form ("hV%iyEtaMult_EP", n), Form ("V_{%i}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaMult_EP.push_back (new TH1F (Form ("hV%iaEtaMult_EP", n), Form ("V_{%i, a}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaMult_EP.push_back (new TH1F (Form ("hV%ibEtaMult_EP", n), Form ("V_{%i, b}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaMult_EP.push_back (new TH1F (Form ("hV%icEtaMult_EP", n), Form ("V_{%i, c}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaMult_EP.push_back (new TH1F (Form ("hV%iEtaMult_EP", n), Form ("V_{%i}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaMult_EP.push_back (new TH1F (Form ("hV%ixYaEtaMult_EP", n), Form ("V_{%i, a}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaMult_EP.push_back (new TH1F (Form ("hV%ixYbEtaMult_EP", n), Form ("V_{%i, b}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaMult_EP.push_back (new TH1F (Form ("hV%ixYcEtaMult_EP", n), Form ("V_{%i, c}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaMult_EP.push_back (new TH1F (Form ("hV%ixYEtaMult_EP", n), Form ("V_{%i}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaMult_EP.push_back (new TH1F (Form ("hV%iyXaEtaMult_EP", n), Form ("V_{%i, a}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaMult_EP.push_back (new TH1F (Form ("hV%iyXbEtaMult_EP", n), Form ("V_{%i, b}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaMult_EP.push_back (new TH1F (Form ("hV%iyXcEtaMult_EP", n), Form ("V_{%i, c}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaMult_EP.push_back (new TH1F (Form ("hV%iyXEtaMult_EP", n), Form ("V_{%i}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));

            hVxaPtCentBS_SP.push_back (new TH1F (Form ("hV%ixaPtCentBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtCentBS_SP.push_back (new TH1F (Form ("hV%ixbPtCentBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtCentBS_SP.push_back (new TH1F (Form ("hV%ixcPtCentBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtCentBS_SP.push_back (new TH1F (Form ("hV%ixPtCentBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtCentBS_SP.push_back (new TH1F (Form ("hV%iyaPtCentBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtCentBS_SP.push_back (new TH1F (Form ("hV%iybPtCentBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtCentBS_SP.push_back (new TH1F (Form ("hV%iycPtCentBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtCentBS_SP.push_back (new TH1F (Form ("hV%iyPtCentBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtCentBS_SP.push_back (new TH1F (Form ("hV%iaPtCentBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtCentBS_SP.push_back (new TH1F (Form ("hV%ibPtCentBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtCentBS_SP.push_back (new TH1F (Form ("hV%icPtCentBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtCentBS_SP.push_back (new TH1F (Form ("hV%iPtCentBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtCentBS_SP.push_back (new TH1F (Form ("hV%ixYaPtCentBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtCentBS_SP.push_back (new TH1F (Form ("hV%ixYbPtCentBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtCentBS_SP.push_back (new TH1F (Form ("hV%ixYcPtCentBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtCentBS_SP.push_back (new TH1F (Form ("hV%ixYPtCentBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtCentBS_SP.push_back (new TH1F (Form ("hV%iyXaPtCentBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtCentBS_SP.push_back (new TH1F (Form ("hV%iyXbPtCentBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtCentBS_SP.push_back (new TH1F (Form ("hV%iyXcPtCentBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtCentBS_SP.push_back (new TH1F (Form ("hV%iyXPtCentBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxaPtCentBS_EP.push_back (new TH1F (Form ("hV%ixaPtCentBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtCentBS_EP.push_back (new TH1F (Form ("hV%ixbPtCentBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtCentBS_EP.push_back (new TH1F (Form ("hV%ixcPtCentBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtCentBS_EP.push_back (new TH1F (Form ("hV%ixPtCentBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtCentBS_EP.push_back (new TH1F (Form ("hV%iyaPtCentBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtCentBS_EP.push_back (new TH1F (Form ("hV%iybPtCentBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtCentBS_EP.push_back (new TH1F (Form ("hV%iycPtCentBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtCentBS_EP.push_back (new TH1F (Form ("hV%iyPtCentBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtCentBS_EP.push_back (new TH1F (Form ("hV%iaPtCentBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtCentBS_EP.push_back (new TH1F (Form ("hV%ibPtCentBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtCentBS_EP.push_back (new TH1F (Form ("hV%icPtCentBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtCentBS_EP.push_back (new TH1F (Form ("hV%iPtCentBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtCentBS_EP.push_back (new TH1F (Form ("hV%ixYaPtCentBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtCentBS_EP.push_back (new TH1F (Form ("hV%ixYbPtCentBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtCentBS_EP.push_back (new TH1F (Form ("hV%ixYcPtCentBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtCentBS_EP.push_back (new TH1F (Form ("hV%ixYPtCentBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtCentBS_EP.push_back (new TH1F (Form ("hV%iyXaPtCentBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtCentBS_EP.push_back (new TH1F (Form ("hV%iyXbPtCentBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtCentBS_EP.push_back (new TH1F (Form ("hV%iyXcPtCentBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtCentBS_EP.push_back (new TH1F (Form ("hV%iyXPtCentBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over centrality); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));

            hVxaPtMultBS_SP.push_back (new TH1F (Form ("hV%ixaPtMultBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtMultBS_SP.push_back (new TH1F (Form ("hV%ixbPtMultBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtMultBS_SP.push_back (new TH1F (Form ("hV%ixcPtMultBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtMultBS_SP.push_back (new TH1F (Form ("hV%ixPtMultBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtMultBS_SP.push_back (new TH1F (Form ("hV%iyaPtMultBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtMultBS_SP.push_back (new TH1F (Form ("hV%iybPtMultBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtMultBS_SP.push_back (new TH1F (Form ("hV%iycPtMultBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtMultBS_SP.push_back (new TH1F (Form ("hV%iyPtMultBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtMultBS_SP.push_back (new TH1F (Form ("hV%iaPtMultBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtMultBS_SP.push_back (new TH1F (Form ("hV%ibPtMultBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtMultBS_SP.push_back (new TH1F (Form ("hV%icPtMultBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtMultBS_SP.push_back (new TH1F (Form ("hV%iPtMultBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtMultBS_SP.push_back (new TH1F (Form ("hV%ixYaPtMultBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtMultBS_SP.push_back (new TH1F (Form ("hV%ixYbPtMultBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtMultBS_SP.push_back (new TH1F (Form ("hV%ixYcPtMultBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtMultBS_SP.push_back (new TH1F (Form ("hV%ixYPtMultBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtMultBS_SP.push_back (new TH1F (Form ("hV%iyXaPtMultBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtMultBS_SP.push_back (new TH1F (Form ("hV%iyXbPtMultBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtMultBS_SP.push_back (new TH1F (Form ("hV%iyXcPtMultBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtMultBS_SP.push_back (new TH1F (Form ("hV%iyXPtMultBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxaPtMultBS_EP.push_back (new TH1F (Form ("hV%ixaPtMultBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxbPtMultBS_EP.push_back (new TH1F (Form ("hV%ixbPtMultBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxcPtMultBS_EP.push_back (new TH1F (Form ("hV%ixcPtMultBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxPtMultBS_EP.push_back (new TH1F (Form ("hV%ixPtMultBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyaPtMultBS_EP.push_back (new TH1F (Form ("hV%iyaPtMultBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVybPtMultBS_EP.push_back (new TH1F (Form ("hV%iybPtMultBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVycPtMultBS_EP.push_back (new TH1F (Form ("hV%iycPtMultBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyPtMultBS_EP.push_back (new TH1F (Form ("hV%iyPtMultBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVaPtMultBS_EP.push_back (new TH1F (Form ("hV%iaPtMultBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVbPtMultBS_EP.push_back (new TH1F (Form ("hV%ibPtMultBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVcPtMultBS_EP.push_back (new TH1F (Form ("hV%icPtMultBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVPtMultBS_EP.push_back (new TH1F (Form ("hV%iPtMultBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{x+y}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYaPtMultBS_EP.push_back (new TH1F (Form ("hV%ixYaPtMultBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYbPtMultBS_EP.push_back (new TH1F (Form ("hV%ixYbPtMultBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYcPtMultBS_EP.push_back (new TH1F (Form ("hV%ixYcPtMultBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVxYPtMultBS_EP.push_back (new TH1F (Form ("hV%ixYPtMultBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{xY}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXaPtMultBS_EP.push_back (new TH1F (Form ("hV%iyXaPtMultBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXbPtMultBS_EP.push_back (new TH1F (Form ("hV%iyXbPtMultBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXcPtMultBS_EP.push_back (new TH1F (Form ("hV%iyXcPtMultBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));
            hVyXPtMultBS_EP.push_back (new TH1F (Form ("hV%iyXPtMultBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over multiplicity); P_{T} [GeV/c]; V_{%i}^{yX}", n, n), nBinsPt_, ptMin_, ptMax_));

            hVxaEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixaEtaCentBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixbEtaCentBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixcEtaCentBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixEtaCentBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyaEtaCentBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaCentBS_SP.push_back (new TH1F (Form ("hV%iybEtaCentBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaCentBS_SP.push_back (new TH1F (Form ("hV%iycEtaCentBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyEtaCentBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaCentBS_SP.push_back (new TH1F (Form ("hV%iaEtaCentBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaCentBS_SP.push_back (new TH1F (Form ("hV%ibEtaCentBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaCentBS_SP.push_back (new TH1F (Form ("hV%icEtaCentBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaCentBS_SP.push_back (new TH1F (Form ("hV%iEtaCentBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixYaEtaCentBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixYbEtaCentBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixYcEtaCentBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaCentBS_SP.push_back (new TH1F (Form ("hV%ixYEtaCentBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyXaEtaCentBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyXbEtaCentBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyXcEtaCentBS_SP", n), Form ("V_{%i, }^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaCentBS_SP.push_back (new TH1F (Form ("hV%iyXEtaCentBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVxaEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixaEtaCentBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixbEtaCentBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixcEtaCentBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixEtaCentBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyaEtaCentBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaCentBS_EP.push_back (new TH1F (Form ("hV%iybEtaCentBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaCentBS_EP.push_back (new TH1F (Form ("hV%iycEtaCentBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyEtaCentBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaCentBS_EP.push_back (new TH1F (Form ("hV%iaEtaCentBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaCentBS_EP.push_back (new TH1F (Form ("hV%ibEtaCentBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaCentBS_EP.push_back (new TH1F (Form ("hV%icEtaCentBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaCentBS_EP.push_back (new TH1F (Form ("hV%iEtaCentBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixYaEtaCentBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixYbEtaCentBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixYcEtaCentBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaCentBS_EP.push_back (new TH1F (Form ("hV%ixYEtaCentBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyXaEtaCentBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyXbEtaCentBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyXcEtaCentBS_EP", n), Form ("V_{%i, }^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaCentBS_EP.push_back (new TH1F (Form ("hV%iyXEtaCentBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));

            hVxaEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixaEtaMultBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixbEtaMultBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixcEtaMultBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixEtaMultBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyaEtaMultBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaMultBS_SP.push_back (new TH1F (Form ("hV%iybEtaMultBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaMultBS_SP.push_back (new TH1F (Form ("hV%iycEtaMultBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyEtaMultBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaMultBS_SP.push_back (new TH1F (Form ("hV%iaEtaMultBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaMultBS_SP.push_back (new TH1F (Form ("hV%ibEtaMultBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaMultBS_SP.push_back (new TH1F (Form ("hV%icEtaMultBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaMultBS_SP.push_back (new TH1F (Form ("hV%iEtaMultBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixYaEtaMultBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixYbEtaMultBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixYcEtaMultBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaMultBS_SP.push_back (new TH1F (Form ("hV%ixYEtaMultBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyXaEtaMultBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyXbEtaMultBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyXcEtaMultBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaMultBS_SP.push_back (new TH1F (Form ("hV%iyXEtaMultBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVxaEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixaEtaMultBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxbEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixbEtaMultBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxcEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixcEtaMultBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVxEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixEtaMultBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEta_, etaMin_, etaMax_));
            hVyaEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyaEtaMultBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVybEtaMultBS_EP.push_back (new TH1F (Form ("hV%iybEtaMultBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVycEtaMultBS_EP.push_back (new TH1F (Form ("hV%iycEtaMultBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVyEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyEtaMultBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEta_, etaMin_, etaMax_));
            hVaEtaMultBS_EP.push_back (new TH1F (Form ("hV%iaEtaMultBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVbEtaMultBS_EP.push_back (new TH1F (Form ("hV%ibEtaMultBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVcEtaMultBS_EP.push_back (new TH1F (Form ("hV%icEtaMultBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVEtaMultBS_EP.push_back (new TH1F (Form ("hV%iEtaMultBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYaEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixYaEtaMultBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYbEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixYbEtaMultBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYcEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixYcEtaMultBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVxYEtaMultBS_EP.push_back (new TH1F (Form ("hV%ixYEtaMultBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXaEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyXaEtaMultBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXbEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyXbEtaMultBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXcEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyXcEtaMultBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));
            hVyXEtaMultBS_EP.push_back (new TH1F (Form ("hV%iyXEtaMultBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEta_, etaMin_, etaMax_));

            if (nBinsEtaRefl_ != 0) {
                hVxaEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixaEtaReflCent_SP", n), Form ("V_{%i, a}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixbEtaReflCent_SP", n), Form ("V_{%i, b}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixcEtaReflCent_SP", n), Form ("V_{%i, c}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixEtaReflCent_SP", n), Form ("V_{%i}^{x, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyaEtaReflCent_SP", n), Form ("V_{%i, a}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflCent_SP.push_back (new TH1F (Form ("hV%iybEtaReflCent_SP", n), Form ("V_{%i, b}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflCent_SP.push_back (new TH1F (Form ("hV%iycEtaReflCent_SP", n), Form ("V_{%i, c}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyEtaReflCent_SP", n), Form ("V_{%i}^{y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflCent_SP.push_back (new TH1F (Form ("hV%iaEtaReflCent_SP", n), Form ("V_{%i, a}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflCent_SP.push_back (new TH1F (Form ("hV%ibEtaReflCent_SP", n), Form ("V_{%i, b}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflCent_SP.push_back (new TH1F (Form ("hV%icEtaReflCent_SP", n), Form ("V_{%i, c}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflCent_SP.push_back (new TH1F (Form ("hV%iEtaReflCent_SP", n), Form ("V_{%i}^{x+y, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixYaEtaReflCent_SP", n), Form ("V_{%i, a}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixYbEtaReflCent_SP", n), Form ("V_{%i, b}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixYcEtaReflCent_SP", n), Form ("V_{%i, c}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflCent_SP.push_back (new TH1F (Form ("hV%ixYEtaReflCent_SP", n), Form ("V_{%i}^{xY, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyXaEtaReflCent_SP", n), Form ("V_{%i, a}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyXbEtaReflCent_SP", n), Form ("V_{%i, b}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyXcEtaReflCent_SP", n), Form ("V_{%i, c}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflCent_SP.push_back (new TH1F (Form ("hV%iyXEtaReflCent_SP", n), Form ("V_{%i}^{yX, SP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxaEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixaEtaReflCent_EP", n), Form ("V_{%i, a}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixbEtaReflCent_EP", n), Form ("V_{%i, b}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixcEtaReflCent_EP", n), Form ("V_{%i, c}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixEtaReflCent_EP", n), Form ("V_{%i}^{x, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyaEtaReflCent_EP", n), Form ("V_{%i, a}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflCent_EP.push_back (new TH1F (Form ("hV%iybEtaReflCent_EP", n), Form ("V_{%i, b}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflCent_EP.push_back (new TH1F (Form ("hV%iycEtaReflCent_EP", n), Form ("V_{%i, c}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyEtaReflCent_EP", n), Form ("V_{%i}^{y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflCent_EP.push_back (new TH1F (Form ("hV%iaEtaReflCent_EP", n), Form ("V_{%i, a}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflCent_EP.push_back (new TH1F (Form ("hV%ibEtaReflCent_EP", n), Form ("V_{%i, b}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflCent_EP.push_back (new TH1F (Form ("hV%icEtaReflCent_EP", n), Form ("V_{%i, c}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflCent_EP.push_back (new TH1F (Form ("hV%iEtaReflCent_EP", n), Form ("V_{%i}^{x+y, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixYaEtaReflCent_EP", n), Form ("V_{%i, a}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixYbEtaReflCent_EP", n), Form ("V_{%i, b}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixYcEtaReflCent_EP", n), Form ("V_{%i, c}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflCent_EP.push_back (new TH1F (Form ("hV%ixYEtaReflCent_EP", n), Form ("V_{%i}^{xY, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyXaEtaReflCent_EP", n), Form ("V_{%i, a}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyXbEtaReflCent_EP", n), Form ("V_{%i, b}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyXcEtaReflCent_EP", n), Form ("V_{%i, c}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflCent_EP.push_back (new TH1F (Form ("hV%iyXEtaReflCent_EP", n), Form ("V_{%i}^{yX, EP} (over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));

                hVxaEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixaEtaReflMult_SP", n), Form ("V_{%i, a}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixbEtaReflMult_SP", n), Form ("V_{%i, b}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixcEtaReflMult_SP", n), Form ("V_{%i, c}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixEtaReflMult_SP", n), Form ("V_{%i}^{x, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyaEtaReflMult_SP", n), Form ("V_{%i, a}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflMult_SP.push_back (new TH1F (Form ("hV%iybEtaReflMult_SP", n), Form ("V_{%i, b}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflMult_SP.push_back (new TH1F (Form ("hV%iycEtaReflMult_SP", n), Form ("V_{%i, c}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyEtaReflMult_SP", n), Form ("V_{%i}^{y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflMult_SP.push_back (new TH1F (Form ("hV%iaEtaReflMult_SP", n), Form ("V_{%i, a}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflMult_SP.push_back (new TH1F (Form ("hV%ibEtaReflMult_SP", n), Form ("V_{%i, b}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflMult_SP.push_back (new TH1F (Form ("hV%icEtaReflMult_SP", n), Form ("V_{%i, c}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflMult_SP.push_back (new TH1F (Form ("hV%iEtaReflMult_SP", n), Form ("V_{%i}^{x+y, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixYaEtaReflMult_SP", n), Form ("V_{%i, a}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixYbEtaReflMult_SP", n), Form ("V_{%i, b}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixYcEtaReflMult_SP", n), Form ("V_{%i, c}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflMult_SP.push_back (new TH1F (Form ("hV%ixYEtaReflMult_SP", n), Form ("V_{%i}^{xY, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyXaEtaReflMult_SP", n), Form ("V_{%i, a}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyXbEtaReflMult_SP", n), Form ("V_{%i, b}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyXcEtaReflMult_SP", n), Form ("V_{%i, c}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflMult_SP.push_back (new TH1F (Form ("hV%iyXEtaReflMult_SP", n), Form ("V_{%i}^{yX, SP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxaEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixaEtaReflMult_EP", n), Form ("V_{%i, a}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixbEtaReflMult_EP", n), Form ("V_{%i, b}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixcEtaReflMult_EP", n), Form ("V_{%i, c}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixEtaReflMult_EP", n), Form ("V_{%i}^{x, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyaEtaReflMult_EP", n), Form ("V_{%i, a}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflMult_EP.push_back (new TH1F (Form ("hV%iybEtaReflMult_EP", n), Form ("V_{%i, b}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflMult_EP.push_back (new TH1F (Form ("hV%iycEtaReflMult_EP", n), Form ("V_{%i, c}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyEtaReflMult_EP", n), Form ("V_{%i}^{y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflMult_EP.push_back (new TH1F (Form ("hV%iaEtaReflMult_EP", n), Form ("V_{%i, a}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflMult_EP.push_back (new TH1F (Form ("hV%ibEtaReflMult_EP", n), Form ("V_{%i, b}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflMult_EP.push_back (new TH1F (Form ("hV%icEtaReflMult_EP", n), Form ("V_{%i, c}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflMult_EP.push_back (new TH1F (Form ("hV%iEtaReflMult_EP", n), Form ("V_{%i}^{x+y, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixYaEtaReflMult_EP", n), Form ("V_{%i, a}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixYbEtaReflMult_EP", n), Form ("V_{%i, b}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixYcEtaReflMult_EP", n), Form ("V_{%i, c}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflMult_EP.push_back (new TH1F (Form ("hV%ixYEtaReflMult_EP", n), Form ("V_{%i}^{xY, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyXaEtaReflMult_EP", n), Form ("V_{%i, a}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyXbEtaReflMult_EP", n), Form ("V_{%i, b}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyXcEtaReflMult_EP", n), Form ("V_{%i, c}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflMult_EP.push_back (new TH1F (Form ("hV%iyXEtaReflMult_EP", n), Form ("V_{%i}^{yX, EP} (over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));

                hVxaEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixaEtaReflCentBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixbEtaReflCentBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixcEtaReflCentBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixEtaReflCentBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyaEtaReflCentBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iybEtaReflCentBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iycEtaReflCentBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyEtaReflCentBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iaEtaReflCentBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ibEtaReflCentBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%icEtaReflCentBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iEtaReflCentBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixYaEtaReflCentBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixYbEtaReflCentBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixYcEtaReflCentBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%ixYEtaReflCentBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyXaEtaReflCentBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyXbEtaReflCentBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyXcEtaReflCentBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflCentBS_SP.push_back (new TH1F (Form ("hV%iyXEtaReflCentBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxaEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixaEtaReflCentBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixbEtaReflCentBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixcEtaReflCentBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixEtaReflCentBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyaEtaReflCentBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iybEtaReflCentBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iycEtaReflCentBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyEtaReflCentBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iaEtaReflCentBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ibEtaReflCentBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%icEtaReflCentBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iEtaReflCentBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixYaEtaReflCentBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixYbEtaReflCentBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixYcEtaReflCentBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%ixYEtaReflCentBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyXaEtaReflCentBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyXbEtaReflCentBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyXcEtaReflCentBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflCentBS_EP.push_back (new TH1F (Form ("hV%iyXEtaReflCentBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over centrality);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));

                hVxaEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixaEtaReflMultBS_SP", n), Form ("V_{%i, a}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixbEtaReflMultBS_SP", n), Form ("V_{%i, b}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixcEtaReflMultBS_SP", n), Form ("V_{%i, c}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixEtaReflMultBS_SP", n), Form ("V_{%i}^{x, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyaEtaReflMultBS_SP", n), Form ("V_{%i, a}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iybEtaReflMultBS_SP", n), Form ("V_{%i, b}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iycEtaReflMultBS_SP", n), Form ("V_{%i, c}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyEtaReflMultBS_SP", n), Form ("V_{%i}^{y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iaEtaReflMultBS_SP", n), Form ("V_{%i, a}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ibEtaReflMultBS_SP", n), Form ("V_{%i, b}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%icEtaReflMultBS_SP", n), Form ("V_{%i, c}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iEtaReflMultBS_SP", n), Form ("V_{%i}^{x+y, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixYaEtaReflMultBS_SP", n), Form ("V_{%i, a}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixYbEtaReflMultBS_SP", n), Form ("V_{%i, b}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixYcEtaReflMultBS_SP", n), Form ("V_{%i, c}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%ixYEtaReflMultBS_SP", n), Form ("V_{%i}^{xY, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyXaEtaReflMultBS_SP", n), Form ("V_{%i, a}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyXbEtaReflMultBS_SP", n), Form ("V_{%i, b}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyXcEtaReflMultBS_SP", n), Form ("V_{%i, c}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflMultBS_SP.push_back (new TH1F (Form ("hV%iyXEtaReflMultBS_SP", n), Form ("V_{%i}^{yX, SP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxaEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixaEtaReflMultBS_EP", n), Form ("V_{%i, a}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxbEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixbEtaReflMultBS_EP", n), Form ("V_{%i, b}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxcEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixcEtaReflMultBS_EP", n), Form ("V_{%i, c}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixEtaReflMultBS_EP", n), Form ("V_{%i}^{x, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyaEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyaEtaReflMultBS_EP", n), Form ("V_{%i, a}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVybEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iybEtaReflMultBS_EP", n), Form ("V_{%i, b}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVycEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iycEtaReflMultBS_EP", n), Form ("V_{%i, c}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyEtaReflMultBS_EP", n), Form ("V_{%i}^{y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVaEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iaEtaReflMultBS_EP", n), Form ("V_{%i, a}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVbEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ibEtaReflMultBS_EP", n), Form ("V_{%i, b}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVcEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%icEtaReflMultBS_EP", n), Form ("V_{%i, c}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iEtaReflMultBS_EP", n), Form ("V_{%i}^{x+y, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{x+y}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYaEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixYaEtaReflMultBS_EP", n), Form ("V_{%i, a}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYbEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixYbEtaReflMultBS_EP", n), Form ("V_{%i, b}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYcEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixYcEtaReflMultBS_EP", n), Form ("V_{%i, c}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVxYEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%ixYEtaReflMultBS_EP", n), Form ("V_{%i}^{xY, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{xY}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXaEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyXaEtaReflMultBS_EP", n), Form ("V_{%i, a}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXbEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyXbEtaReflMultBS_EP", n), Form ("V_{%i, b}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXcEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyXcEtaReflMultBS_EP", n), Form ("V_{%i, c}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
                hVyXEtaReflMultBS_EP.push_back (new TH1F (Form ("hV%iyXEtaReflMultBS_EP", n), Form ("V_{%i}^{yX, EP} (sampling, over multiplicity);" ,n) + varName_ + Form (";V_{%i}^{yX}", n), nBinsEtaRefl_, (-1) * etaMax_, (-1) * etaMax_ + (etaMax_ - etaMin_) / nBinsEta_ * nBinsEtaRefl_));
            }

// test
//            if (step == 1) {
//                for (Int_t j = 1; j <= nBinsCent_; j++) {
//                    for (Int_t k = 1; k <= nBinsBS_; k++) {
//                        cent = h2RxaCent_SP [i] -> GetXaxis () -> GetBinCenter (j);
//                        //QQ [AB, AC, BC][x, y][x, y]
//                        QQ [0][0][0] = p2XaXbCent_SP [i] -> GetBinContent (j, k);
//                        QQ [1][0][0] = p2XaXcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [2][0][0] = p2XbXcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [0][1][1] = p2YaYbCent_SP [i] -> GetBinContent (j, k);
//                        QQ [1][1][1] = p2YaYcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [2][1][1] = p2YbYcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [0][0][1] = p2XaYbCent_SP [i] -> GetBinContent (j, k);
//                        QQ [1][0][1] = p2XaYcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [2][0][1] = p2XbYcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [0][1][0] = p2YaXbCent_SP [i] -> GetBinContent (j, k);
//                        QQ [1][1][0] = p2YaXcCent_SP [i] -> GetBinContent (j, k);
//                        QQ [2][1][0] = p2YbXcCent_SP [i] -> GetBinContent (j, k);
//                        testTree -> Fill();
//                    }
//                }
//            }
// test

            pqQaCent_SP [i] -> Add (pxXaCent_SP [i], pyYaCent_SP [i], 0.5, 0.5);
            pqQbCent_SP [i] -> Add (pxXbCent_SP [i], pyYbCent_SP [i], 0.5, 0.5);
            pqQcCent_SP [i] -> Add (pxXcCent_SP [i], pyYcCent_SP [i], 0.5, 0.5);
            pqQaCent_EP [i] -> Add (pxXaCent_EP [i], pyYaCent_EP [i], 0.5, 0.5);
            pqQbCent_EP [i] -> Add (pxXbCent_EP [i], pyYbCent_EP [i], 0.5, 0.5);
            pqQcCent_EP [i] -> Add (pxXcCent_EP [i], pyYcCent_EP [i], 0.5, 0.5);

            pqQaMult_SP [i] -> Add (pxXaMult_SP [i], pyYaMult_SP [i], 0.5, 0.5);
            pqQbMult_SP [i] -> Add (pxXbMult_SP [i], pyYbMult_SP [i], 0.5, 0.5);
            pqQcMult_SP [i] -> Add (pxXcMult_SP [i], pyYcMult_SP [i], 0.5, 0.5);
            pqQaMult_EP [i] -> Add (pxXaMult_EP [i], pyYaMult_EP [i], 0.5, 0.5);
            pqQbMult_EP [i] -> Add (pxXbMult_EP [i], pyYbMult_EP [i], 0.5, 0.5);
            pqQcMult_EP [i] -> Add (pxXcMult_EP [i], pyYcMult_EP [i], 0.5, 0.5);

            p2qQaCent_SP [i] -> Add (p2xXaCent_SP [i], p2yYaCent_SP [i], 0.5, 0.5);
            p2qQbCent_SP [i] -> Add (p2xXbCent_SP [i], p2yYbCent_SP [i], 0.5, 0.5);
            p2qQcCent_SP [i] -> Add (p2xXcCent_SP [i], p2yYcCent_SP [i], 0.5, 0.5);
            p2qQaCent_EP [i] -> Add (p2xXaCent_EP [i], p2yYaCent_EP [i], 0.5, 0.5);
            p2qQbCent_EP [i] -> Add (p2xXbCent_EP [i], p2yYbCent_EP [i], 0.5, 0.5);
            p2qQcCent_EP [i] -> Add (p2xXcCent_EP [i], p2yYcCent_EP [i], 0.5, 0.5);

            p2qQaMult_SP [i] -> Add (p2xXaMult_SP [i], p2yYaMult_SP [i], 0.5, 0.5);
            p2qQbMult_SP [i] -> Add (p2xXbMult_SP [i], p2yYbMult_SP [i], 0.5, 0.5);
            p2qQcMult_SP [i] -> Add (p2xXcMult_SP [i], p2yYcMult_SP [i], 0.5, 0.5);
            p2qQaMult_EP [i] -> Add (p2xXaMult_EP [i], p2yYaMult_EP [i], 0.5, 0.5);
            p2qQbMult_EP [i] -> Add (p2xXbMult_EP [i], p2yYbMult_EP [i], 0.5, 0.5);
            p2qQcMult_EP [i] -> Add (p2xXcMult_EP [i], p2yYcMult_EP [i], 0.5, 0.5);

            p2qQaPtCent_SP [i] -> Add (p2xXaPtCent_SP [i], p2yYaPtCent_SP [i], 0.5, 0.5);
            p2qQbPtCent_SP [i] -> Add (p2xXbPtCent_SP [i], p2yYbPtCent_SP [i], 0.5, 0.5);
            p2qQcPtCent_SP [i] -> Add (p2xXcPtCent_SP [i], p2yYcPtCent_SP [i], 0.5, 0.5);
            p2qQaPtCent_EP [i] -> Add (p2xXaPtCent_EP [i], p2yYaPtCent_EP [i], 0.5, 0.5);
            p2qQbPtCent_EP [i] -> Add (p2xXbPtCent_EP [i], p2yYbPtCent_EP [i], 0.5, 0.5);
            p2qQcPtCent_EP [i] -> Add (p2xXcPtCent_EP [i], p2yYcPtCent_EP [i], 0.5, 0.5);

            p2qQaEtaCent_SP [i] -> Add (p2xXaEtaCent_SP [i], p2yYaEtaCent_SP [i], 0.5, 0.5);
            p2qQbEtaCent_SP [i] -> Add (p2xXbEtaCent_SP [i], p2yYbEtaCent_SP [i], 0.5, 0.5);
            p2qQcEtaCent_SP [i] -> Add (p2xXcEtaCent_SP [i], p2yYcEtaCent_SP [i], 0.5, 0.5);
            p2qQaEtaCent_EP [i] -> Add (p2xXaEtaCent_EP [i], p2yYaEtaCent_EP [i], 0.5, 0.5);
            p2qQbEtaCent_EP [i] -> Add (p2xXbEtaCent_EP [i], p2yYbEtaCent_EP [i], 0.5, 0.5);
            p2qQcEtaCent_EP [i] -> Add (p2xXcEtaCent_EP [i], p2yYcEtaCent_EP [i], 0.5, 0.5);

            p2qQaPtMult_SP [i] -> Add (p2xXaPtMult_SP [i], p2yYaPtMult_SP [i], 0.5, 0.5);
            p2qQbPtMult_SP [i] -> Add (p2xXbPtMult_SP [i], p2yYbPtMult_SP [i], 0.5, 0.5);
            p2qQcPtMult_SP [i] -> Add (p2xXcPtMult_SP [i], p2yYcPtMult_SP [i], 0.5, 0.5);
            p2qQaPtMult_EP [i] -> Add (p2xXaPtMult_EP [i], p2yYaPtMult_EP [i], 0.5, 0.5);
            p2qQbPtMult_EP [i] -> Add (p2xXbPtMult_EP [i], p2yYbPtMult_EP [i], 0.5, 0.5);
            p2qQcPtMult_EP [i] -> Add (p2xXcPtMult_EP [i], p2yYcPtMult_EP [i], 0.5, 0.5);

            p2qQaEtaMult_SP [i] -> Add (p2xXaEtaMult_SP [i], p2yYaEtaMult_SP [i], 0.5, 0.5);
            p2qQbEtaMult_SP [i] -> Add (p2xXbEtaMult_SP [i], p2yYbEtaMult_SP [i], 0.5, 0.5);
            p2qQcEtaMult_SP [i] -> Add (p2xXcEtaMult_SP [i], p2yYcEtaMult_SP [i], 0.5, 0.5);
            p2qQaEtaMult_EP [i] -> Add (p2xXaEtaMult_EP [i], p2yYaEtaMult_EP [i], 0.5, 0.5);
            p2qQbEtaMult_EP [i] -> Add (p2xXbEtaMult_EP [i], p2yYbEtaMult_EP [i], 0.5, 0.5);
            p2qQcEtaMult_EP [i] -> Add (p2xXcEtaMult_EP [i], p2yYcEtaMult_EP [i], 0.5, 0.5);

            p3qQaPtCent_SP [i] -> Add (p3xXaPtCent_SP [i], p3yYaPtCent_SP [i], 0.5, 0.5);
            p3qQbPtCent_SP [i] -> Add (p3xXbPtCent_SP [i], p3yYbPtCent_SP [i], 0.5, 0.5);
            p3qQcPtCent_SP [i] -> Add (p3xXcPtCent_SP [i], p3yYcPtCent_SP [i], 0.5, 0.5);
            p3qQaPtCent_EP [i] -> Add (p3xXaPtCent_EP [i], p3yYaPtCent_EP [i], 0.5, 0.5);
            p3qQbPtCent_EP [i] -> Add (p3xXbPtCent_EP [i], p3yYbPtCent_EP [i], 0.5, 0.5);
            p3qQcPtCent_EP [i] -> Add (p3xXcPtCent_EP [i], p3yYcPtCent_EP [i], 0.5, 0.5);

            p3qQaEtaCent_SP [i] -> Add (p3xXaEtaCent_SP [i], p3yYaEtaCent_SP [i], 0.5, 0.5);
            p3qQbEtaCent_SP [i] -> Add (p3xXbEtaCent_SP [i], p3yYbEtaCent_SP [i], 0.5, 0.5);
            p3qQcEtaCent_SP [i] -> Add (p3xXcEtaCent_SP [i], p3yYcEtaCent_SP [i], 0.5, 0.5);
            p3qQaEtaCent_EP [i] -> Add (p3xXaEtaCent_EP [i], p3yYaEtaCent_EP [i], 0.5, 0.5);
            p3qQbEtaCent_EP [i] -> Add (p3xXbEtaCent_EP [i], p3yYbEtaCent_EP [i], 0.5, 0.5);
            p3qQcEtaCent_EP [i] -> Add (p3xXcEtaCent_EP [i], p3yYcEtaCent_EP [i], 0.5, 0.5);

            p3qQaPtMult_SP [i] -> Add (p3xXaPtMult_SP [i], p3yYaPtMult_SP [i], 0.5, 0.5);
            p3qQbPtMult_SP [i] -> Add (p3xXbPtMult_SP [i], p3yYbPtMult_SP [i], 0.5, 0.5);
            p3qQcPtMult_SP [i] -> Add (p3xXcPtMult_SP [i], p3yYcPtMult_SP [i], 0.5, 0.5);
            p3qQaPtMult_EP [i] -> Add (p3xXaPtMult_EP [i], p3yYaPtMult_EP [i], 0.5, 0.5);
            p3qQbPtMult_EP [i] -> Add (p3xXbPtMult_EP [i], p3yYbPtMult_EP [i], 0.5, 0.5);
            p3qQcPtMult_EP [i] -> Add (p3xXcPtMult_EP [i], p3yYcPtMult_EP [i], 0.5, 0.5);

            p3qQaEtaMult_SP [i] -> Add (p3xXaEtaMult_SP [i], p3yYaEtaMult_SP [i], 0.5, 0.5);
            p3qQbEtaMult_SP [i] -> Add (p3xXbEtaMult_SP [i], p3yYbEtaMult_SP [i], 0.5, 0.5);
            p3qQcEtaMult_SP [i] -> Add (p3xXcEtaMult_SP [i], p3yYcEtaMult_SP [i], 0.5, 0.5);
            p3qQaEtaMult_EP [i] -> Add (p3xXaEtaMult_EP [i], p3yYaEtaMult_EP [i], 0.5, 0.5);
            p3qQbEtaMult_EP [i] -> Add (p3xXbEtaMult_EP [i], p3yYbEtaMult_EP [i], 0.5, 0.5);
            p3qQcEtaMult_EP [i] -> Add (p3xXcEtaMult_EP [i], p3yYcEtaMult_EP [i], 0.5, 0.5);

            pQaQbCent_SP [i] -> Add (pXaXbCent_SP [i], pYaYbCent_SP [i], 0.5, 0.5);
            pQaQcCent_SP [i] -> Add (pXaXcCent_SP [i], pYaYcCent_SP [i], 0.5, 0.5);
            pQbQcCent_SP [i] -> Add (pXbXcCent_SP [i], pYbYcCent_SP [i], 0.5, 0.5);
            pQaQbCent_EP [i] -> Add (pXaXbCent_EP [i], pYaYbCent_EP [i], 0.5, 0.5);
            pQaQcCent_EP [i] -> Add (pXaXcCent_EP [i], pYaYcCent_EP [i], 0.5, 0.5);
            pQbQcCent_EP [i] -> Add (pXbXcCent_EP [i], pYbYcCent_EP [i], 0.5, 0.5);

            pQaQbMult_SP [i] -> Add (pXaXbMult_SP [i], pYaYbMult_SP [i], 0.5, 0.5);
            pQaQcMult_SP [i] -> Add (pXaXcMult_SP [i], pYaYcMult_SP [i], 0.5, 0.5);
            pQbQcMult_SP [i] -> Add (pXbXcMult_SP [i], pYbYcMult_SP [i], 0.5, 0.5);
            pQaQbMult_EP [i] -> Add (pXaXbMult_EP [i], pYaYbMult_EP [i], 0.5, 0.5);
            pQaQcMult_EP [i] -> Add (pXaXcMult_EP [i], pYaYcMult_EP [i], 0.5, 0.5);
            pQbQcMult_EP [i] -> Add (pXbXcMult_EP [i], pYbYcMult_EP [i], 0.5, 0.5);

            p2QaQbCent_SP [i] -> Add (p2XaXbCent_SP [i], p2YaYbCent_SP [i], 0.5, 0.5);
            p2QaQcCent_SP [i] -> Add (p2XaXcCent_SP [i], p2YaYcCent_SP [i], 0.5, 0.5);
            p2QbQcCent_SP [i] -> Add (p2XbXcCent_SP [i], p2YbYcCent_SP [i], 0.5, 0.5);
            p2QaQbCent_EP [i] -> Add (p2XaXbCent_EP [i], p2YaYbCent_EP [i], 0.5, 0.5);
            p2QaQcCent_EP [i] -> Add (p2XaXcCent_EP [i], p2YaYcCent_EP [i], 0.5, 0.5);
            p2QbQcCent_EP [i] -> Add (p2XbXcCent_EP [i], p2YbYcCent_EP [i], 0.5, 0.5);

            p2QaQbMult_SP [i] -> Add (p2XaXbMult_SP [i], p2YaYbMult_SP [i], 0.5, 0.5);
            p2QaQcMult_SP [i] -> Add (p2XaXcMult_SP [i], p2YaYcMult_SP [i], 0.5, 0.5);
            p2QbQcMult_SP [i] -> Add (p2XbXcMult_SP [i], p2YbYcMult_SP [i], 0.5, 0.5);
            p2QaQbMult_EP [i] -> Add (p2XaXbMult_EP [i], p2YaYbMult_EP [i], 0.5, 0.5);
            p2QaQcMult_EP [i] -> Add (p2XaXcMult_EP [i], p2YaYcMult_EP [i], 0.5, 0.5);
            p2QbQcMult_EP [i] -> Add (p2XbXcMult_EP [i], p2YbYcMult_EP [i], 0.5, 0.5);

            CalculateCorrelationsWithSampling (p2XaXbCent_SP [i], pXaXbCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaXcCent_SP [i], pXaXcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XbXcCent_SP [i], pXbXcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaYbCent_SP [i], pYaYbCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaYcCent_SP [i], pYaYcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YbYcCent_SP [i], pYbYcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaYbCent_SP [i], pXaYbCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaYcCent_SP [i], pXaYcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XbYcCent_SP [i], pXbYcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaXbCent_SP [i], pYaXbCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaXcCent_SP [i], pYaXcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YbXcCent_SP [i], pYbXcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QaQbCent_SP [i], pQaQbCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QaQcCent_SP [i], pQaQcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QbQcCent_SP [i], pQbQcCentBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaXbCent_EP [i], pXaXbCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaXcCent_EP [i], pXaXcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XbXcCent_EP [i], pXbXcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaYbCent_EP [i], pYaYbCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaYcCent_EP [i], pYaYcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YbYcCent_EP [i], pYbYcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaYbCent_EP [i], pXaYbCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaYcCent_EP [i], pXaYcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XbYcCent_EP [i], pXbYcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaXbCent_EP [i], pYaXbCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaXcCent_EP [i], pYaXcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YbXcCent_EP [i], pYbXcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QaQbCent_EP [i], pQaQbCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QaQcCent_EP [i], pQaQcCentBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QbQcCent_EP [i], pQbQcCentBS_EP [i]);

            CalculateCorrelationsWithSampling (p2XaXbMult_SP [i], pXaXbMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaXcMult_SP [i], pXaXcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XbXcMult_SP [i], pXbXcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaYbMult_SP [i], pYaYbMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaYcMult_SP [i], pYaYcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YbYcMult_SP [i], pYbYcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaYbMult_SP [i], pXaYbMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaYcMult_SP [i], pXaYcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XbYcMult_SP [i], pXbYcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaXbMult_SP [i], pYaXbMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YaXcMult_SP [i], pYaXcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2YbXcMult_SP [i], pYbXcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QaQbMult_SP [i], pQaQbMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QaQcMult_SP [i], pQaQcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2QbQcMult_SP [i], pQbQcMultBS_SP [i]);
            CalculateCorrelationsWithSampling (p2XaXbMult_EP [i], pXaXbMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaXcMult_EP [i], pXaXcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XbXcMult_EP [i], pXbXcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaYbMult_EP [i], pYaYbMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaYcMult_EP [i], pYaYcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YbYcMult_EP [i], pYbYcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaYbMult_EP [i], pXaYbMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XaYcMult_EP [i], pXaYcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2XbYcMult_EP [i], pXbYcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaXbMult_EP [i], pYaXbMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YaXcMult_EP [i], pYaXcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2YbXcMult_EP [i], pYbXcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QaQbMult_EP [i], pQaQbMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QaQcMult_EP [i], pQaQcMultBS_EP [i]);
            CalculateCorrelationsWithSampling (p2QbQcMult_EP [i], pQbQcMultBS_EP [i]);

            CalculateResolutionNoSampling (pXaXbCent_SP [i], pXaXcCent_SP [i], pXbXcCent_SP [i], hRxaCent_SP [i], hRxbCent_SP [i], hRxcCent_SP [i]);
            CalculateResolutionNoSampling (pYaYbCent_SP [i], pYaYcCent_SP [i], pYbYcCent_SP [i], hRyaCent_SP [i], hRybCent_SP [i], hRycCent_SP [i]);
            CalculateResolutionNoSampling (pQaQbCent_SP [i], pQaQcCent_SP [i], pQbQcCent_SP [i], hRaCent_SP [i], hRbCent_SP [i], hRcCent_SP [i]);
            CalculateResolutionNoSampling (pXaXbCent_EP [i], pXaXcCent_EP [i], pXbXcCent_EP [i], hRxaCent_EP [i], hRxbCent_EP [i], hRxcCent_EP [i]);
            CalculateResolutionNoSampling (pYaYbCent_EP [i], pYaYcCent_EP [i], pYbYcCent_EP [i], hRyaCent_EP [i], hRybCent_EP [i], hRycCent_EP [i]);
            CalculateResolutionNoSampling (pQaQbCent_EP [i], pQaQcCent_EP [i], pQbQcCent_EP [i], hRaCent_EP [i], hRbCent_EP [i], hRcCent_EP [i]);

            CalculateResolutionNoSampling (pXaXbMult_SP [i], pXaXcMult_SP [i], pXbXcMult_SP [i], hRxaMult_SP [i], hRxbMult_SP [i], hRxcMult_SP [i]);
            CalculateResolutionNoSampling (pYaYbMult_SP [i], pYaYcMult_SP [i], pYbYcMult_SP [i], hRyaMult_SP [i], hRybMult_SP [i], hRycMult_SP [i]);
            CalculateResolutionNoSampling (pQaQbMult_SP [i], pQaQcMult_SP [i], pQbQcMult_SP [i], hRaMult_SP [i], hRbMult_SP [i], hRcMult_SP [i]);
            CalculateResolutionNoSampling (pXaXbMult_EP [i], pXaXcMult_EP [i], pXbXcMult_EP [i], hRxaMult_EP [i], hRxbMult_EP [i], hRxcMult_EP [i]);
            CalculateResolutionNoSampling (pYaYbMult_EP [i], pYaYcMult_EP [i], pYbYcMult_EP [i], hRyaMult_EP [i], hRybMult_EP [i], hRycMult_EP [i]);
            CalculateResolutionNoSampling (pQaQbMult_EP [i], pQaQcMult_EP [i], pQbQcMult_EP [i], hRaMult_EP [i], hRbMult_EP [i], hRcMult_EP [i]);

            CalculateResolutionWithSampling (p2XaXbCent_SP [i], p2XaXcCent_SP [i], p2XbXcCent_SP [i], h2RxaCent_SP [i], h2RxbCent_SP [i], h2RxcCent_SP [i]);
            CalculateResolutionWithSampling (p2YaYbCent_SP [i], p2YaYcCent_SP [i], p2YbYcCent_SP [i], h2RyaCent_SP [i], h2RybCent_SP [i], h2RycCent_SP [i]);
            CalculateResolutionWithSampling (p2QaQbCent_SP [i], p2QaQcCent_SP [i], p2QbQcCent_SP [i], h2RaCent_SP [i], h2RbCent_SP [i], h2RcCent_SP [i]);
            CalculateResolutionWithSampling (p2XaXbCent_EP [i], p2XaXcCent_EP [i], p2XbXcCent_EP [i], h2RxaCent_EP [i], h2RxbCent_EP [i], h2RxcCent_EP [i]);
            CalculateResolutionWithSampling (p2YaYbCent_EP [i], p2YaYcCent_EP [i], p2YbYcCent_EP [i], h2RyaCent_EP [i], h2RybCent_EP [i], h2RycCent_EP [i]);
            CalculateResolutionWithSampling (p2QaQbCent_EP [i], p2QaQcCent_EP [i], p2QbQcCent_EP [i], h2RaCent_EP [i], h2RbCent_EP [i], h2RcCent_EP [i]);

            TH2toTH1withSampling (h2RxaCent_SP [i], hRxaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxbCent_SP [i], hRxbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxcCent_SP [i], hRxcCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RyaCent_SP [i], hRyaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RybCent_SP [i], hRybCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RycCent_SP [i], hRycCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RaCent_SP [i], hRaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RbCent_SP [i], hRbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RcCent_SP [i], hRcCentBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxaCent_EP [i], hRxaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RxbCent_EP [i], hRxbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RxcCent_EP [i], hRxcCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RyaCent_EP [i], hRyaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RybCent_EP [i], hRybCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RycCent_EP [i], hRycCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RaCent_EP [i], hRaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RbCent_EP [i], hRbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RcCent_EP [i], hRcCentBS_EP [i], distrDir);

            TH2toTH1withSampling (h2RxaMult_SP [i], hRxaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxbMult_SP [i], hRxbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxcMult_SP [i], hRxcMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RyaMult_SP [i], hRyaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RybMult_SP [i], hRybMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RycMult_SP [i], hRycMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RaMult_SP [i], hRaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RbMult_SP [i], hRbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RcMult_SP [i], hRcMultBS_SP [i], distrDir);
            TH2toTH1withSampling (h2RxaMult_EP [i], hRxaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RxbMult_EP [i], hRxbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RxcMult_EP [i], hRxcMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RyaMult_EP [i], hRyaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RybMult_EP [i], hRybMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RycMult_EP [i], hRycMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RaMult_EP [i], hRaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RbMult_EP [i], hRbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (h2RcMult_EP [i], hRcMultBS_EP [i], distrDir);

            Sqrt (hRxaCent_SP [i]);
            Sqrt (hRxbCent_SP [i]);
            Sqrt (hRxcCent_SP [i]);
            Sqrt (hRyaCent_SP [i]);
            Sqrt (hRybCent_SP [i]);
            Sqrt (hRycCent_SP [i]);
            Sqrt (hRaCent_SP [i]);
            Sqrt (hRbCent_SP [i]);
            Sqrt (hRcCent_SP [i]);
            Sqrt (hRxaCent_EP [i]);
            Sqrt (hRxbCent_EP [i]);
            Sqrt (hRxcCent_EP [i]);
            Sqrt (hRyaCent_EP [i]);
            Sqrt (hRybCent_EP [i]);
            Sqrt (hRycCent_EP [i]);
            Sqrt (hRaCent_EP [i]);
            Sqrt (hRbCent_EP [i]);
            Sqrt (hRcCent_EP [i]);
            Sqrt (hRxaCentBS_SP [i]);
            Sqrt (hRxbCentBS_SP [i]);
            Sqrt (hRxcCentBS_SP [i]);
            Sqrt (hRyaCentBS_SP [i]);
            Sqrt (hRybCentBS_SP [i]);
            Sqrt (hRycCentBS_SP [i]);
            Sqrt (hRaCentBS_SP [i]);
            Sqrt (hRbCentBS_SP [i]);
            Sqrt (hRcCentBS_SP [i]);
            Sqrt (hRxaCentBS_EP [i]);
            Sqrt (hRxbCentBS_EP [i]);
            Sqrt (hRxcCentBS_EP [i]);
            Sqrt (hRyaCentBS_EP [i]);
            Sqrt (hRybCentBS_EP [i]);
            Sqrt (hRycCentBS_EP [i]);
            Sqrt (hRaCentBS_EP [i]);
            Sqrt (hRbCentBS_EP [i]);
            Sqrt (hRcCentBS_EP [i]);

            Sqrt (hRxaMult_SP [i]);
            Sqrt (hRxbMult_SP [i]);
            Sqrt (hRxcMult_SP [i]);
            Sqrt (hRyaMult_SP [i]);
            Sqrt (hRybMult_SP [i]);
            Sqrt (hRycMult_SP [i]);
            Sqrt (hRaMult_SP [i]);
            Sqrt (hRbMult_SP [i]);
            Sqrt (hRcMult_SP [i]);
            Sqrt (hRxaMult_EP [i]);
            Sqrt (hRxbMult_EP [i]);
            Sqrt (hRxcMult_EP [i]);
            Sqrt (hRyaMult_EP [i]);
            Sqrt (hRybMult_EP [i]);
            Sqrt (hRycMult_EP [i]);
            Sqrt (hRaMult_EP [i]);
            Sqrt (hRbMult_EP [i]);
            Sqrt (hRcMult_EP [i]);
            Sqrt (hRxaMultBS_SP [i]);
            Sqrt (hRxbMultBS_SP [i]);
            Sqrt (hRxcMultBS_SP [i]);
            Sqrt (hRyaMultBS_SP [i]);
            Sqrt (hRybMultBS_SP [i]);
            Sqrt (hRycMultBS_SP [i]);
            Sqrt (hRaMultBS_SP [i]);
            Sqrt (hRbMultBS_SP [i]);
            Sqrt (hRcMultBS_SP [i]);
            Sqrt (hRxaMultBS_EP [i]);
            Sqrt (hRxbMultBS_EP [i]);
            Sqrt (hRxcMultBS_EP [i]);
            Sqrt (hRyaMultBS_EP [i]);
            Sqrt (hRybMultBS_EP [i]);
            Sqrt (hRycMultBS_EP [i]);
            Sqrt (hRaMultBS_EP [i]);
            Sqrt (hRbMultBS_EP [i]);
            Sqrt (hRcMultBS_EP [i]);

            Sqrt (h2RxaCent_SP [i]);
            Sqrt (h2RxbCent_SP [i]);
            Sqrt (h2RxcCent_SP [i]);
            Sqrt (h2RyaCent_SP [i]);
            Sqrt (h2RybCent_SP [i]);
            Sqrt (h2RycCent_SP [i]);
            Sqrt (h2RaCent_SP [i]);
            Sqrt (h2RbCent_SP [i]);
            Sqrt (h2RcCent_SP [i]);
            Sqrt (h2RxaCent_EP [i]);
            Sqrt (h2RxbCent_EP [i]);
            Sqrt (h2RxcCent_EP [i]);
            Sqrt (h2RyaCent_EP [i]);
            Sqrt (h2RybCent_EP [i]);
            Sqrt (h2RycCent_EP [i]);
            Sqrt (h2RaCent_EP [i]);
            Sqrt (h2RbCent_EP [i]);
            Sqrt (h2RcCent_EP [i]);

            Sqrt (h2RxaMult_SP [i]);
            Sqrt (h2RxbMult_SP [i]);
            Sqrt (h2RxcMult_SP [i]);
            Sqrt (h2RyaMult_SP [i]);
            Sqrt (h2RybMult_SP [i]);
            Sqrt (h2RycMult_SP [i]);
            Sqrt (h2RaMult_SP [i]);
            Sqrt (h2RbMult_SP [i]);
            Sqrt (h2RcMult_SP [i]);
            Sqrt (h2RxaMult_EP [i]);
            Sqrt (h2RxbMult_EP [i]);
            Sqrt (h2RxcMult_EP [i]);
            Sqrt (h2RyaMult_EP [i]);
            Sqrt (h2RybMult_EP [i]);
            Sqrt (h2RycMult_EP [i]);
            Sqrt (h2RaMult_EP [i]);
            Sqrt (h2RbMult_EP [i]);
            Sqrt (h2RcMult_EP [i]);

            CalculateFlowNoSampling (pxXaCent_SP [i], hRxaCent_SP [i], hVxaCent_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxXbCent_SP [i], hRxbCent_SP [i], hVxbCent_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxXcCent_SP [i], hRxcCent_SP [i], hVxcCent_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyYaCent_SP [i], hRyaCent_SP [i], hVyaCent_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyYbCent_SP [i], hRybCent_SP [i], hVybCent_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyYcCent_SP [i], hRycCent_SP [i], hVycCent_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pqQaCent_SP [i], hRaCent_SP [i], hVaCent_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pqQbCent_SP [i], hRbCent_SP [i], hVbCent_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pqQcCent_SP [i], hRcCent_SP [i], hVcCent_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyXaCent_SP [i], hRyaCent_SP [i], hVxYaCent_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyXbCent_SP [i], hRybCent_SP [i], hVxYbCent_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyXcCent_SP [i], hRycCent_SP [i], hVxYcCent_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxYaCent_SP [i], hRxaCent_SP [i], hVyXaCent_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxYbCent_SP [i], hRxbCent_SP [i], hVyXbCent_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxYcCent_SP [i], hRxcCent_SP [i], hVyXcCent_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxXaCent_EP [i], hRxaCent_EP [i], hVxaCent_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxXbCent_EP [i], hRxbCent_EP [i], hVxbCent_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxXcCent_EP [i], hRxcCent_EP [i], hVxcCent_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyYaCent_EP [i], hRyaCent_EP [i], hVyaCent_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyYbCent_EP [i], hRybCent_EP [i], hVybCent_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyYcCent_EP [i], hRycCent_EP [i], hVycCent_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pqQaCent_EP [i], hRaCent_EP [i], hVaCent_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pqQbCent_EP [i], hRbCent_EP [i], hVbCent_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pqQcCent_EP [i], hRcCent_EP [i], hVcCent_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyXaCent_EP [i], hRyaCent_EP [i], hVxYaCent_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyXbCent_EP [i], hRybCent_EP [i], hVxYbCent_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyXcCent_EP [i], hRycCent_EP [i], hVxYcCent_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxYaCent_EP [i], hRxaCent_EP [i], hVyXaCent_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxYbCent_EP [i], hRxbCent_EP [i], hVyXbCent_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxYcCent_EP [i], hRxcCent_EP [i], hVyXcCent_EP [i], resSign_[i][2]);

            CalculateFlowNoSampling (pxXaMult_SP [i], hRxaMult_SP [i], hVxaMult_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxXbMult_SP [i], hRxbMult_SP [i], hVxbMult_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxXcMult_SP [i], hRxcMult_SP [i], hVxcMult_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyYaMult_SP [i], hRyaMult_SP [i], hVyaMult_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyYbMult_SP [i], hRybMult_SP [i], hVybMult_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyYcMult_SP [i], hRycMult_SP [i], hVycMult_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pqQaMult_SP [i], hRaMult_SP [i], hVaMult_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pqQbMult_SP [i], hRbMult_SP [i], hVbMult_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pqQcMult_SP [i], hRcMult_SP [i], hVcMult_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyXaMult_SP [i], hRyaMult_SP [i], hVxYaMult_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyXbMult_SP [i], hRybMult_SP [i], hVxYbMult_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyXcMult_SP [i], hRycMult_SP [i], hVxYcMult_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxYaMult_SP [i], hRxaMult_SP [i], hVyXaMult_SP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxYbMult_SP [i], hRxbMult_SP [i], hVyXbMult_SP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxYcMult_SP [i], hRxcMult_SP [i], hVyXcMult_SP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxXaMult_EP [i], hRxaMult_EP [i], hVxaMult_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxXbMult_EP [i], hRxbMult_EP [i], hVxbMult_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxXcMult_EP [i], hRxcMult_EP [i], hVxcMult_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyYaMult_EP [i], hRyaMult_EP [i], hVyaMult_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyYbMult_EP [i], hRybMult_EP [i], hVybMult_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyYcMult_EP [i], hRycMult_EP [i], hVycMult_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pqQaMult_EP [i], hRaMult_EP [i], hVaMult_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pqQbMult_EP [i], hRbMult_EP [i], hVbMult_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pqQcMult_EP [i], hRcMult_EP [i], hVcMult_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pyXaMult_EP [i], hRyaMult_EP [i], hVxYaMult_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pyXbMult_EP [i], hRybMult_EP [i], hVxYbMult_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pyXcMult_EP [i], hRycMult_EP [i], hVxYcMult_EP [i], resSign_[i][2]);
            CalculateFlowNoSampling (pxYaMult_EP [i], hRxaMult_EP [i], hVyXaMult_EP [i], resSign_[i][0]);
            CalculateFlowNoSampling (pxYbMult_EP [i], hRxbMult_EP [i], hVyXbMult_EP [i], resSign_[i][1]);
            CalculateFlowNoSampling (pxYcMult_EP [i], hRxcMult_EP [i], hVyXcMult_EP [i], resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaPtCent_SP [i], hRxaCent_SP [i], hVxaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbPtCent_SP [i], hRxbCent_SP [i], hVxbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcPtCent_SP [i], hRxcCent_SP [i], hVxcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaPtCent_SP [i], hRyaCent_SP [i], hVyaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbPtCent_SP [i], hRybCent_SP [i], hVybPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcPtCent_SP [i], hRycCent_SP [i], hVycPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaPtCent_SP [i], hRaCent_SP [i], hVaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbPtCent_SP [i], hRbCent_SP [i], hVbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcPtCent_SP [i], hRcCent_SP [i], hVcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaPtCent_SP [i], hRyaCent_SP [i], hVxYaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbPtCent_SP [i], hRybCent_SP [i], hVxYbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcPtCent_SP [i], hRycCent_SP [i], hVxYcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaPtCent_SP [i], hRxaCent_SP [i], hVyXaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbPtCent_SP [i], hRxbCent_SP [i], hVyXbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcPtCent_SP [i], hRxcCent_SP [i], hVyXcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaPtCent_EP [i], hRxaCent_EP [i], hVxaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbPtCent_EP [i], hRxbCent_EP [i], hVxbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcPtCent_EP [i], hRxcCent_EP [i], hVxcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaPtCent_EP [i], hRyaCent_EP [i], hVyaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbPtCent_EP [i], hRybCent_EP [i], hVybPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcPtCent_EP [i], hRycCent_EP [i], hVycPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaPtCent_EP [i], hRaCent_EP [i], hVaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbPtCent_EP [i], hRbCent_EP [i], hVbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcPtCent_EP [i], hRcCent_EP [i], hVcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaPtCent_EP [i], hRyaCent_EP [i], hVxYaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbPtCent_EP [i], hRybCent_EP [i], hVxYbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcPtCent_EP [i], hRycCent_EP [i], hVxYcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaPtCent_EP [i], hRxaCent_EP [i], hVyXaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbPtCent_EP [i], hRxbCent_EP [i], hVyXbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcPtCent_EP [i], hRxcCent_EP [i], hVyXcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaEtaCent_SP [i], hRxaCent_SP [i], hVxaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbEtaCent_SP [i], hRxbCent_SP [i], hVxbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcEtaCent_SP [i], hRxcCent_SP [i], hVxcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaEtaCent_SP [i], hRyaCent_SP [i], hVyaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbEtaCent_SP [i], hRybCent_SP [i], hVybEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcEtaCent_SP [i], hRycCent_SP [i], hVycEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaEtaCent_SP [i], hRaCent_SP [i], hVaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbEtaCent_SP [i], hRbCent_SP [i], hVbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcEtaCent_SP [i], hRcCent_SP [i], hVcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaEtaCent_SP [i], hRyaCent_SP [i], hVxYaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbEtaCent_SP [i], hRybCent_SP [i], hVxYbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcEtaCent_SP [i], hRycCent_SP [i], hVxYcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaEtaCent_SP [i], hRxaCent_SP [i], hVyXaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbEtaCent_SP [i], hRxbCent_SP [i], hVyXbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcEtaCent_SP [i], hRxcCent_SP [i], hVyXcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaEtaCent_EP [i], hRxaCent_EP [i], hVxaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbEtaCent_EP [i], hRxbCent_EP [i], hVxbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcEtaCent_EP [i], hRxcCent_EP [i], hVxcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaEtaCent_EP [i], hRyaCent_EP [i], hVyaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbEtaCent_EP [i], hRybCent_EP [i], hVybEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcEtaCent_EP [i], hRycCent_EP [i], hVycEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaEtaCent_EP [i], hRaCent_EP [i], hVaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbEtaCent_EP [i], hRbCent_EP [i], hVbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcEtaCent_EP [i], hRcCent_EP [i], hVcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaEtaCent_EP [i], hRyaCent_EP [i], hVxYaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbEtaCent_EP [i], hRybCent_EP [i], hVxYbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcEtaCent_EP [i], hRycCent_EP [i], hVxYcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaEtaCent_EP [i], hRxaCent_EP [i], hVyXaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbEtaCent_EP [i], hRxbCent_EP [i], hVyXbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcEtaCent_EP [i], hRxcCent_EP [i], hVyXcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaPtMult_SP [i], hRxaMult_SP [i], hVxaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbPtMult_SP [i], hRxbMult_SP [i], hVxbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcPtMult_SP [i], hRxcMult_SP [i], hVxcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaPtMult_SP [i], hRyaMult_SP [i], hVyaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbPtMult_SP [i], hRybMult_SP [i], hVybPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcPtMult_SP [i], hRycMult_SP [i], hVycPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaPtMult_SP [i], hRaMult_SP [i], hVaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbPtMult_SP [i], hRbMult_SP [i], hVbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcPtMult_SP [i], hRcMult_SP [i], hVcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaPtMult_SP [i], hRyaMult_SP [i], hVxYaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbPtMult_SP [i], hRybMult_SP [i], hVxYbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcPtMult_SP [i], hRycMult_SP [i], hVxYcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaPtMult_SP [i], hRxaMult_SP [i], hVyXaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbPtMult_SP [i], hRxbMult_SP [i], hVyXbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcPtMult_SP [i], hRxcMult_SP [i], hVyXcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaPtMult_EP [i], hRxaMult_EP [i], hVxaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbPtMult_EP [i], hRxbMult_EP [i], hVxbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcPtMult_EP [i], hRxcMult_EP [i], hVxcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaPtMult_EP [i], hRyaMult_EP [i], hVyaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbPtMult_EP [i], hRybMult_EP [i], hVybPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcPtMult_EP [i], hRycMult_EP [i], hVycPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaPtMult_EP [i], hRaMult_EP [i], hVaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbPtMult_EP [i], hRbMult_EP [i], hVbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcPtMult_EP [i], hRcMult_EP [i], hVcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaPtMult_EP [i], hRyaMult_EP [i], hVxYaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbPtMult_EP [i], hRybMult_EP [i], hVxYbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcPtMult_EP [i], hRycMult_EP [i], hVxYcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaPtMult_EP [i], hRxaMult_EP [i], hVyXaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbPtMult_EP [i], hRxbMult_EP [i], hVyXbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcPtMult_EP [i], hRxcMult_EP [i], hVyXcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaEtaMult_SP [i], hRxaMult_SP [i], hVxaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbEtaMult_SP [i], hRxbMult_SP [i], hVxbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcEtaMult_SP [i], hRxcMult_SP [i], hVxcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaEtaMult_SP [i], hRyaMult_SP [i], hVyaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbEtaMult_SP [i], hRybMult_SP [i], hVybEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcEtaMult_SP [i], hRycMult_SP [i], hVycEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaEtaMult_SP [i], hRaMult_SP [i], hVaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbEtaMult_SP [i], hRbMult_SP [i], hVbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcEtaMult_SP [i], hRcMult_SP [i], hVcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaEtaMult_SP [i], hRyaMult_SP [i], hVxYaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbEtaMult_SP [i], hRybMult_SP [i], hVxYbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcEtaMult_SP [i], hRycMult_SP [i], hVxYcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaEtaMult_SP [i], hRxaMult_SP [i], hVyXaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbEtaMult_SP [i], hRxbMult_SP [i], hVyXbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcEtaMult_SP [i], hRxcMult_SP [i], hVyXcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowNoSampling (p2xXaEtaMult_EP [i], hRxaMult_EP [i], hVxaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xXbEtaMult_EP [i], hRxbMult_EP [i], hVxbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xXcEtaMult_EP [i], hRxcMult_EP [i], hVxcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yYaEtaMult_EP [i], hRyaMult_EP [i], hVyaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yYbEtaMult_EP [i], hRybMult_EP [i], hVybEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yYcEtaMult_EP [i], hRycMult_EP [i], hVycEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2qQaEtaMult_EP [i], hRaMult_EP [i], hVaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2qQbEtaMult_EP [i], hRbMult_EP [i], hVbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2qQcEtaMult_EP [i], hRcMult_EP [i], hVcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2yXaEtaMult_EP [i], hRyaMult_EP [i], hVxYaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2yXbEtaMult_EP [i], hRybMult_EP [i], hVxYbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2yXcEtaMult_EP [i], hRycMult_EP [i], hVxYcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowNoSampling (p2xYaEtaMult_EP [i], hRxaMult_EP [i], hVyXaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowNoSampling (p2xYbEtaMult_EP [i], hRxbMult_EP [i], hVyXbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowNoSampling (p2xYcEtaMult_EP [i], hRxcMult_EP [i], hVyXcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p2xXaCent_SP [i], h2RxaCent_SP [i], p2VxaCent_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xXbCent_SP [i], h2RxbCent_SP [i], p2VxbCent_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xXcCent_SP [i], h2RxcCent_SP [i], p2VxcCent_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yYaCent_SP [i], h2RyaCent_SP [i], p2VyaCent_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yYbCent_SP [i], h2RybCent_SP [i], p2VybCent_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yYcCent_SP [i], h2RycCent_SP [i], p2VycCent_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2qQaCent_SP [i], h2RaCent_SP [i], p2VaCent_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2qQbCent_SP [i], h2RbCent_SP [i], p2VbCent_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2qQcCent_SP [i], h2RcCent_SP [i], p2VcCent_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yXaCent_SP [i], h2RyaCent_SP [i], p2VxYaCent_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yXbCent_SP [i], h2RybCent_SP [i], p2VxYbCent_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yXcCent_SP [i], h2RycCent_SP [i], p2VxYcCent_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xYaCent_SP [i], h2RxaCent_SP [i], p2VyXaCent_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xYbCent_SP [i], h2RxbCent_SP [i], p2VyXbCent_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xYcCent_SP [i], h2RxcCent_SP [i], p2VyXcCent_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xXaCent_EP [i], h2RxaCent_EP [i], p2VxaCent_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xXbCent_EP [i], h2RxbCent_EP [i], p2VxbCent_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xXcCent_EP [i], h2RxcCent_EP [i], p2VxcCent_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yYaCent_EP [i], h2RyaCent_EP [i], p2VyaCent_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yYbCent_EP [i], h2RybCent_EP [i], p2VybCent_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yYcCent_EP [i], h2RycCent_EP [i], p2VycCent_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2qQaCent_EP [i], h2RaCent_EP [i], p2VaCent_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2qQbCent_EP [i], h2RbCent_EP [i], p2VbCent_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2qQcCent_EP [i], h2RcCent_EP [i], p2VcCent_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yXaCent_EP [i], h2RyaCent_EP [i], p2VxYaCent_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yXbCent_EP [i], h2RybCent_EP [i], p2VxYbCent_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yXcCent_EP [i], h2RycCent_EP [i], p2VxYcCent_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xYaCent_EP [i], h2RxaCent_EP [i], p2VyXaCent_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xYbCent_EP [i], h2RxbCent_EP [i], p2VyXbCent_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xYcCent_EP [i], h2RxcCent_EP [i], p2VyXcCent_EP [i], resSign_[i][2]);

            CalculateFlowWithSampling (p2xXaMult_SP [i], h2RxaMult_SP [i], p2VxaMult_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xXbMult_SP [i], h2RxbMult_SP [i], p2VxbMult_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xXcMult_SP [i], h2RxcMult_SP [i], p2VxcMult_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yYaMult_SP [i], h2RyaMult_SP [i], p2VyaMult_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yYbMult_SP [i], h2RybMult_SP [i], p2VybMult_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yYcMult_SP [i], h2RycMult_SP [i], p2VycMult_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2qQaMult_SP [i], h2RaMult_SP [i], p2VaMult_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2qQbMult_SP [i], h2RbMult_SP [i], p2VbMult_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2qQcMult_SP [i], h2RcMult_SP [i], p2VcMult_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yXaMult_SP [i], h2RyaMult_SP [i], p2VxYaMult_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yXbMult_SP [i], h2RybMult_SP [i], p2VxYbMult_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yXcMult_SP [i], h2RycMult_SP [i], p2VxYcMult_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xYaMult_SP [i], h2RxaMult_SP [i], p2VyXaMult_SP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xYbMult_SP [i], h2RxbMult_SP [i], p2VyXbMult_SP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xYcMult_SP [i], h2RxcMult_SP [i], p2VyXcMult_SP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xXaMult_EP [i], h2RxaMult_EP [i], p2VxaMult_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xXbMult_EP [i], h2RxbMult_EP [i], p2VxbMult_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xXcMult_EP [i], h2RxcMult_EP [i], p2VxcMult_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yYaMult_EP [i], h2RyaMult_EP [i], p2VyaMult_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yYbMult_EP [i], h2RybMult_EP [i], p2VybMult_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yYcMult_EP [i], h2RycMult_EP [i], p2VycMult_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2qQaMult_EP [i], h2RaMult_EP [i], p2VaMult_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2qQbMult_EP [i], h2RbMult_EP [i], p2VbMult_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2qQcMult_EP [i], h2RcMult_EP [i], p2VcMult_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2yXaMult_EP [i], h2RyaMult_EP [i], p2VxYaMult_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2yXbMult_EP [i], h2RybMult_EP [i], p2VxYbMult_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2yXcMult_EP [i], h2RycMult_EP [i], p2VxYcMult_EP [i], resSign_[i][2]);
            CalculateFlowWithSampling (p2xYaMult_EP [i], h2RxaMult_EP [i], p2VyXaMult_EP [i], resSign_[i][0]);
            CalculateFlowWithSampling (p2xYbMult_EP [i], h2RxbMult_EP [i], p2VyXbMult_EP [i], resSign_[i][1]);
            CalculateFlowWithSampling (p2xYcMult_EP [i], h2RxcMult_EP [i], p2VyXcMult_EP [i], resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaPtCent_SP [i], h2RxaCent_SP [i], p2VxaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbPtCent_SP [i], h2RxbCent_SP [i], p2VxbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcPtCent_SP [i], h2RxcCent_SP [i], p2VxcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaPtCent_SP [i], h2RyaCent_SP [i], p2VyaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbPtCent_SP [i], h2RybCent_SP [i], p2VybPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcPtCent_SP [i], h2RycCent_SP [i], p2VycPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaPtCent_SP [i], h2RaCent_SP [i], p2VaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbPtCent_SP [i], h2RbCent_SP [i], p2VbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcPtCent_SP [i], h2RcCent_SP [i], p2VcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaPtCent_SP [i], h2RyaCent_SP [i], p2VxYaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbPtCent_SP [i], h2RybCent_SP [i], p2VxYbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcPtCent_SP [i], h2RycCent_SP [i], p2VxYcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaPtCent_SP [i], h2RxaCent_SP [i], p2VyXaPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbPtCent_SP [i], h2RxbCent_SP [i], p2VyXbPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcPtCent_SP [i], h2RxcCent_SP [i], p2VyXcPtCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaPtCent_EP [i], h2RxaCent_EP [i], p2VxaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbPtCent_EP [i], h2RxbCent_EP [i], p2VxbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcPtCent_EP [i], h2RxcCent_EP [i], p2VxcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaPtCent_EP [i], h2RyaCent_EP [i], p2VyaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbPtCent_EP [i], h2RybCent_EP [i], p2VybPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcPtCent_EP [i], h2RycCent_EP [i], p2VycPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaPtCent_EP [i], h2RaCent_EP [i], p2VaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbPtCent_EP [i], h2RbCent_EP [i], p2VbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcPtCent_EP [i], h2RcCent_EP [i], p2VcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaPtCent_EP [i], h2RyaCent_EP [i], p2VxYaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbPtCent_EP [i], h2RybCent_EP [i], p2VxYbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcPtCent_EP [i], h2RycCent_EP [i], p2VxYcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaPtCent_EP [i], h2RxaCent_EP [i], p2VyXaPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbPtCent_EP [i], h2RxbCent_EP [i], p2VyXbPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcPtCent_EP [i], h2RxcCent_EP [i], p2VyXcPtCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaEtaCent_SP [i], h2RxaCent_SP [i], p2VxaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbEtaCent_SP [i], h2RxbCent_SP [i], p2VxbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcEtaCent_SP [i], h2RxcCent_SP [i], p2VxcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaEtaCent_SP [i], h2RyaCent_SP [i], p2VyaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbEtaCent_SP [i], h2RybCent_SP [i], p2VybEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcEtaCent_SP [i], h2RycCent_SP [i], p2VycEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaEtaCent_SP [i], h2RaCent_SP [i], p2VaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbEtaCent_SP [i], h2RbCent_SP [i], p2VbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcEtaCent_SP [i], h2RcCent_SP [i], p2VcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaEtaCent_SP [i], h2RyaCent_SP [i], p2VxYaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbEtaCent_SP [i], h2RybCent_SP [i], p2VxYbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcEtaCent_SP [i], h2RycCent_SP [i], p2VxYcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaEtaCent_SP [i], h2RxaCent_SP [i], p2VyXaEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbEtaCent_SP [i], h2RxbCent_SP [i], p2VyXbEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcEtaCent_SP [i], h2RxcCent_SP [i], p2VyXcEtaCent_SP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaEtaCent_EP [i], h2RxaCent_EP [i], p2VxaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbEtaCent_EP [i], h2RxbCent_EP [i], p2VxbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcEtaCent_EP [i], h2RxcCent_EP [i], p2VxcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaEtaCent_EP [i], h2RyaCent_EP [i], p2VyaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbEtaCent_EP [i], h2RybCent_EP [i], p2VybEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcEtaCent_EP [i], h2RycCent_EP [i], p2VycEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaEtaCent_EP [i], h2RaCent_EP [i], p2VaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbEtaCent_EP [i], h2RbCent_EP [i], p2VbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcEtaCent_EP [i], h2RcCent_EP [i], p2VcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaEtaCent_EP [i], h2RyaCent_EP [i], p2VxYaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbEtaCent_EP [i], h2RybCent_EP [i], p2VxYbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcEtaCent_EP [i], h2RycCent_EP [i], p2VxYcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaEtaCent_EP [i], h2RxaCent_EP [i], p2VyXaEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbEtaCent_EP [i], h2RxbCent_EP [i], p2VyXbEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcEtaCent_EP [i], h2RxcCent_EP [i], p2VyXcEtaCent_EP [i], centLowerBin_, centHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaPtMult_SP [i], h2RxaMult_SP [i], p2VxaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbPtMult_SP [i], h2RxbMult_SP [i], p2VxbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcPtMult_SP [i], h2RxcMult_SP [i], p2VxcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaPtMult_SP [i], h2RyaMult_SP [i], p2VyaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbPtMult_SP [i], h2RybMult_SP [i], p2VybPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcPtMult_SP [i], h2RycMult_SP [i], p2VycPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaPtMult_SP [i], h2RaMult_SP [i], p2VaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbPtMult_SP [i], h2RbMult_SP [i], p2VbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcPtMult_SP [i], h2RcMult_SP [i], p2VcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaPtMult_SP [i], h2RyaMult_SP [i], p2VxYaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbPtMult_SP [i], h2RybMult_SP [i], p2VxYbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcPtMult_SP [i], h2RycMult_SP [i], p2VxYcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaPtMult_SP [i], h2RxaMult_SP [i], p2VyXaPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbPtMult_SP [i], h2RxbMult_SP [i], p2VyXbPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcPtMult_SP [i], h2RxcMult_SP [i], p2VyXcPtMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaPtMult_EP [i], h2RxaMult_EP [i], p2VxaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbPtMult_EP [i], h2RxbMult_EP [i], p2VxbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcPtMult_EP [i], h2RxcMult_EP [i], p2VxcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaPtMult_EP [i], h2RyaMult_EP [i], p2VyaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbPtMult_EP [i], h2RybMult_EP [i], p2VybPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcPtMult_EP [i], h2RycMult_EP [i], p2VycPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaPtMult_EP [i], h2RaMult_EP [i], p2VaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbPtMult_EP [i], h2RbMult_EP [i], p2VbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcPtMult_EP [i], h2RcMult_EP [i], p2VcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaPtMult_EP [i], h2RyaMult_EP [i], p2VxYaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbPtMult_EP [i], h2RybMult_EP [i], p2VxYbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcPtMult_EP [i], h2RycMult_EP [i], p2VxYcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaPtMult_EP [i], h2RxaMult_EP [i], p2VyXaPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbPtMult_EP [i], h2RxbMult_EP [i], p2VyXbPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcPtMult_EP [i], h2RxcMult_EP [i], p2VyXcPtMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaEtaMult_SP [i], h2RxaMult_SP [i], p2VxaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbEtaMult_SP [i], h2RxbMult_SP [i], p2VxbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcEtaMult_SP [i], h2RxcMult_SP [i], p2VxcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaEtaMult_SP [i], h2RyaMult_SP [i], p2VyaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbEtaMult_SP [i], h2RybMult_SP [i], p2VybEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcEtaMult_SP [i], h2RycMult_SP [i], p2VycEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaEtaMult_SP [i], h2RaMult_SP [i], p2VaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbEtaMult_SP [i], h2RbMult_SP [i], p2VbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcEtaMult_SP [i], h2RcMult_SP [i], p2VcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaEtaMult_SP [i], h2RyaMult_SP [i], p2VxYaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbEtaMult_SP [i], h2RybMult_SP [i], p2VxYbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcEtaMult_SP [i], h2RycMult_SP [i], p2VxYcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaEtaMult_SP [i], h2RxaMult_SP [i], p2VyXaEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbEtaMult_SP [i], h2RxbMult_SP [i], p2VyXbEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcEtaMult_SP [i], h2RxcMult_SP [i], p2VyXcEtaMult_SP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CalculateFlowWithSampling (p3xXaEtaMult_EP [i], h2RxaMult_EP [i], p2VxaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xXbEtaMult_EP [i], h2RxbMult_EP [i], p2VxbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xXcEtaMult_EP [i], h2RxcMult_EP [i], p2VxcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yYaEtaMult_EP [i], h2RyaMult_EP [i], p2VyaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yYbEtaMult_EP [i], h2RybMult_EP [i], p2VybEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yYcEtaMult_EP [i], h2RycMult_EP [i], p2VycEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3qQaEtaMult_EP [i], h2RaMult_EP [i], p2VaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3qQbEtaMult_EP [i], h2RbMult_EP [i], p2VbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3qQcEtaMult_EP [i], h2RcMult_EP [i], p2VcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3yXaEtaMult_EP [i], h2RyaMult_EP [i], p2VxYaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3yXbEtaMult_EP [i], h2RybMult_EP [i], p2VxYbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3yXcEtaMult_EP [i], h2RycMult_EP [i], p2VxYcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);
            CalculateFlowWithSampling (p3xYaEtaMult_EP [i], h2RxaMult_EP [i], p2VyXaEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][0]);
            CalculateFlowWithSampling (p3xYbEtaMult_EP [i], h2RxbMult_EP [i], p2VyXbEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][1]);
            CalculateFlowWithSampling (p3xYcEtaMult_EP [i], h2RxcMult_EP [i], p2VyXcEtaMult_EP [i], mhLowerBin_, mhHigherBin_, resSign_[i][2]);

            CombineSubevents (hVaCent_SP [i], hVbCent_SP [i], hVcCent_SP [i], hVCent_SP [i]);
            CombineSubevents (hVxaCent_SP [i], hVxbCent_SP [i], hVxcCent_SP [i], hVxCent_SP [i]);
            CombineSubevents (hVyaCent_SP [i], hVybCent_SP [i], hVycCent_SP [i], hVyCent_SP [i]);
            CombineSubevents (hVaCent_EP [i], hVbCent_EP [i], hVcCent_EP [i], hVCent_EP [i]);
            CombineSubevents (hVxaCent_EP [i], hVxbCent_EP [i], hVxcCent_EP [i], hVxCent_EP [i]);
            CombineSubevents (hVyaCent_EP [i], hVybCent_EP [i], hVycCent_EP [i], hVyCent_EP [i]);

            CombineSubevents (hVaMult_SP [i], hVbMult_SP [i], hVcMult_SP [i], hVMult_SP [i]);
            CombineSubevents (hVxaMult_SP [i], hVxbMult_SP [i], hVxcMult_SP [i], hVxMult_SP [i]);
            CombineSubevents (hVyaMult_SP [i], hVybMult_SP [i], hVycMult_SP [i], hVyMult_SP [i]);
            CombineSubevents (hVaMult_EP [i], hVbMult_EP [i], hVcMult_EP [i], hVMult_EP [i]);
            CombineSubevents (hVxaMult_EP [i], hVxbMult_EP [i], hVxcMult_EP [i], hVxMult_EP [i]);
            CombineSubevents (hVyaMult_EP [i], hVybMult_EP [i], hVycMult_EP [i], hVyMult_EP [i]);

            CombineSubevents (hVaPtCent_SP [i], hVbPtCent_SP [i], hVcPtCent_SP [i], hVPtCent_SP [i]);
            CombineSubevents (hVxaPtCent_SP [i], hVxbPtCent_SP [i], hVxcPtCent_SP [i], hVxPtCent_SP [i]);
            CombineSubevents (hVyaPtCent_SP [i], hVybPtCent_SP [i], hVycPtCent_SP [i], hVyPtCent_SP [i]);
            CombineSubevents (hVaPtCent_EP [i], hVbPtCent_EP [i], hVcPtCent_EP [i], hVPtCent_EP [i]);
            CombineSubevents (hVxaPtCent_EP [i], hVxbPtCent_EP [i], hVxcPtCent_EP [i], hVxPtCent_EP [i]);
            CombineSubevents (hVyaPtCent_EP [i], hVybPtCent_EP [i], hVycPtCent_EP [i], hVyPtCent_EP [i]);

            CombineSubevents (hVaEtaCent_SP [i], hVbEtaCent_SP [i], hVcEtaCent_SP [i], hVEtaCent_SP [i]);
            CombineSubevents (hVxaEtaCent_SP [i], hVxbEtaCent_SP [i], hVxcEtaCent_SP [i], hVxEtaCent_SP [i]);
            CombineSubevents (hVyaEtaCent_SP [i], hVybEtaCent_SP [i], hVycEtaCent_SP [i], hVyEtaCent_SP [i]);
            CombineSubevents (hVaEtaCent_EP [i], hVbEtaCent_EP [i], hVcEtaCent_EP [i], hVEtaCent_EP [i]);
            CombineSubevents (hVxaEtaCent_EP [i], hVxbEtaCent_EP [i], hVxcEtaCent_EP [i], hVxEtaCent_EP [i]);
            CombineSubevents (hVyaEtaCent_EP [i], hVybEtaCent_EP [i], hVycEtaCent_EP [i], hVyEtaCent_EP [i]);

            CombineSubevents (hVaPtMult_SP [i], hVbPtMult_SP [i], hVcPtMult_SP [i], hVPtMult_SP [i]);
            CombineSubevents (hVxaPtMult_SP [i], hVxbPtMult_SP [i], hVxcPtMult_SP [i], hVxPtMult_SP [i]);
            CombineSubevents (hVyaPtMult_SP [i], hVybPtMult_SP [i], hVycPtMult_SP [i], hVyPtMult_SP [i]);
            CombineSubevents (hVaPtMult_EP [i], hVbPtMult_EP [i], hVcPtMult_EP [i], hVPtMult_EP [i]);
            CombineSubevents (hVxaPtMult_EP [i], hVxbPtMult_EP [i], hVxcPtMult_EP [i], hVxPtMult_EP [i]);
            CombineSubevents (hVyaPtMult_EP [i], hVybPtMult_EP [i], hVycPtMult_EP [i], hVyPtMult_EP [i]);

            CombineSubevents (hVaEtaMult_SP [i], hVbEtaMult_SP [i], hVcEtaMult_SP [i], hVEtaMult_SP [i]);
            CombineSubevents (hVxaEtaMult_SP [i], hVxbEtaMult_SP [i], hVxcEtaMult_SP [i], hVxEtaMult_SP [i]);
            CombineSubevents (hVyaEtaMult_SP [i], hVybEtaMult_SP [i], hVycEtaMult_SP [i], hVyEtaMult_SP [i]);
            CombineSubevents (hVaEtaMult_EP [i], hVbEtaMult_EP [i], hVcEtaMult_EP [i], hVEtaMult_EP [i]);
            CombineSubevents (hVxaEtaMult_EP [i], hVxbEtaMult_EP [i], hVxcEtaMult_EP [i], hVxEtaMult_EP [i]);
            CombineSubevents (hVyaEtaMult_EP [i], hVybEtaMult_EP [i], hVycEtaMult_EP [i], hVyEtaMult_EP [i]);

            CombineSubevents (p2VaCent_SP [i], p2VbCent_SP [i], p2VcCent_SP [i], p2VCent_SP [i]);
            CombineSubevents (p2VxaCent_SP [i], p2VxbCent_SP [i], p2VxcCent_SP [i], p2VxCent_SP [i]);
            CombineSubevents (p2VyaCent_SP [i], p2VybCent_SP [i], p2VycCent_SP [i], p2VyCent_SP [i]);
            CombineSubevents (p2VaCent_EP [i], p2VbCent_EP [i], p2VcCent_EP [i], p2VCent_EP [i]);
            CombineSubevents (p2VxaCent_EP [i], p2VxbCent_EP [i], p2VxcCent_EP [i], p2VxCent_EP [i]);
            CombineSubevents (p2VyaCent_EP [i], p2VybCent_EP [i], p2VycCent_EP [i], p2VyCent_EP [i]);

            CombineSubevents (p2VaMult_SP [i], p2VbMult_SP [i], p2VcMult_SP [i], p2VMult_SP [i]);
            CombineSubevents (p2VxaMult_SP [i], p2VxbMult_SP [i], p2VxcMult_SP [i], p2VxMult_SP [i]);
            CombineSubevents (p2VyaMult_SP [i], p2VybMult_SP [i], p2VycMult_SP [i], p2VyMult_SP [i]);
            CombineSubevents (p2VaMult_EP [i], p2VbMult_EP [i], p2VcMult_EP [i], p2VMult_EP [i]);
            CombineSubevents (p2VxaMult_EP [i], p2VxbMult_EP [i], p2VxcMult_EP [i], p2VxMult_EP [i]);
            CombineSubevents (p2VyaMult_EP [i], p2VybMult_EP [i], p2VycMult_EP [i], p2VyMult_EP [i]);

            CombineSubevents (p2VaPtCent_SP [i], p2VbPtCent_SP [i], p2VcPtCent_SP [i], p2VPtCent_SP [i]);
            CombineSubevents (p2VxaPtCent_SP [i], p2VxbPtCent_SP [i], p2VxcPtCent_SP [i], p2VxPtCent_SP [i]);
            CombineSubevents (p2VyaPtCent_SP [i], p2VybPtCent_SP [i], p2VycPtCent_SP [i], p2VyPtCent_SP [i]);
            CombineSubevents (p2VaPtCent_EP [i], p2VbPtCent_EP [i], p2VcPtCent_EP [i], p2VPtCent_EP [i]);
            CombineSubevents (p2VxaPtCent_EP [i], p2VxbPtCent_EP [i], p2VxcPtCent_EP [i], p2VxPtCent_EP [i]);
            CombineSubevents (p2VyaPtCent_EP [i], p2VybPtCent_EP [i], p2VycPtCent_EP [i], p2VyPtCent_EP [i]);

            CombineSubevents (p2VaEtaCent_SP [i], p2VbEtaCent_SP [i], p2VcEtaCent_SP [i], p2VEtaCent_SP [i]);
            CombineSubevents (p2VxaEtaCent_SP [i], p2VxbEtaCent_SP [i], p2VxcEtaCent_SP [i], p2VxEtaCent_SP [i]);
            CombineSubevents (p2VyaEtaCent_SP [i], p2VybEtaCent_SP [i], p2VycEtaCent_SP [i], p2VyEtaCent_SP [i]);
            CombineSubevents (p2VaEtaCent_EP [i], p2VbEtaCent_EP [i], p2VcEtaCent_EP [i], p2VEtaCent_EP [i]);
            CombineSubevents (p2VxaEtaCent_EP [i], p2VxbEtaCent_EP [i], p2VxcEtaCent_EP [i], p2VxEtaCent_EP [i]);
            CombineSubevents (p2VyaEtaCent_EP [i], p2VybEtaCent_EP [i], p2VycEtaCent_EP [i], p2VyEtaCent_EP [i]);

            CombineSubevents (p2VaPtMult_SP [i], p2VbPtMult_SP [i], p2VcPtMult_SP [i], p2VPtMult_SP [i]);
            CombineSubevents (p2VxaPtMult_SP [i], p2VxbPtMult_SP [i], p2VxcPtMult_SP [i], p2VxPtMult_SP [i]);
            CombineSubevents (p2VyaPtMult_SP [i], p2VybPtMult_SP [i], p2VycPtMult_SP [i], p2VyPtMult_SP [i]);
            CombineSubevents (p2VaPtMult_EP [i], p2VbPtMult_EP [i], p2VcPtMult_EP [i], p2VPtMult_EP [i]);
            CombineSubevents (p2VxaPtMult_EP [i], p2VxbPtMult_EP [i], p2VxcPtMult_EP [i], p2VxPtMult_EP [i]);
            CombineSubevents (p2VyaPtMult_EP [i], p2VybPtMult_EP [i], p2VycPtMult_EP [i], p2VyPtMult_EP [i]);

            CombineSubevents (p2VaEtaMult_SP [i], p2VbEtaMult_SP [i], p2VcEtaMult_SP [i], p2VEtaMult_SP [i]);
            CombineSubevents (p2VxaEtaMult_SP [i], p2VxbEtaMult_SP [i], p2VxcEtaMult_SP [i], p2VxEtaMult_SP [i]);
            CombineSubevents (p2VyaEtaMult_SP [i], p2VybEtaMult_SP [i], p2VycEtaMult_SP [i], p2VyEtaMult_SP [i]);
            CombineSubevents (p2VaEtaMult_EP [i], p2VbEtaMult_EP [i], p2VcEtaMult_EP [i], p2VEtaMult_EP [i]);
            CombineSubevents (p2VxaEtaMult_EP [i], p2VxbEtaMult_EP [i], p2VxcEtaMult_EP [i], p2VxEtaMult_EP [i]);
            CombineSubevents (p2VyaEtaMult_EP [i], p2VybEtaMult_EP [i], p2VycEtaMult_EP [i], p2VyEtaMult_EP [i]);

            TH2toTH1withSampling (p2VxaCent_SP [i], hVxaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbCent_SP [i], hVxbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcCent_SP [i], hVxcCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaCent_SP [i], hVyaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybCent_SP [i], hVybCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycCent_SP [i], hVycCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaCent_SP [i], hVaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbCent_SP [i], hVbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcCent_SP [i], hVcCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VCent_SP [i], hVCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxCent_SP [i], hVxCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyCent_SP [i], hVyCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaCent_SP [i], hVxYaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbCent_SP [i], hVxYbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcCent_SP [i], hVxYcCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaCent_SP [i], hVyXaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbCent_SP [i], hVyXbCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcCent_SP [i], hVyXcCentBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaCent_EP [i], hVxaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbCent_EP [i], hVxbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcCent_EP [i], hVxcCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaCent_EP [i], hVyaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybCent_EP [i], hVybCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycCent_EP [i], hVycCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaCent_EP [i], hVaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbCent_EP [i], hVbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcCent_EP [i], hVcCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VCent_EP [i], hVCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxCent_EP [i], hVxCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyCent_EP [i], hVyCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaCent_EP [i], hVxYaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbCent_EP [i], hVxYbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcCent_EP [i], hVxYcCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaCent_EP [i], hVyXaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbCent_EP [i], hVyXbCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcCent_EP [i], hVyXcCentBS_EP [i], distrDir);

            TH2toTH1withSampling (p2VxaMult_SP [i], hVxaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbMult_SP [i], hVxbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcMult_SP [i], hVxcMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaMult_SP [i], hVyaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybMult_SP [i], hVybMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycMult_SP [i], hVycMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaMult_SP [i], hVaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbMult_SP [i], hVbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcMult_SP [i], hVcMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VMult_SP [i], hVMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxMult_SP [i], hVxMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyMult_SP [i], hVyMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaMult_SP [i], hVxYaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbMult_SP [i], hVxYbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcMult_SP [i], hVxYcMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaMult_SP [i], hVyXaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbMult_SP [i], hVyXbMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcMult_SP [i], hVyXcMultBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaMult_EP [i], hVxaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbMult_EP [i], hVxbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcMult_EP [i], hVxcMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaMult_EP [i], hVyaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybMult_EP [i], hVybMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycMult_EP [i], hVycMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaMult_EP [i], hVaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbMult_EP [i], hVbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcMult_EP [i], hVcMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VMult_EP [i], hVMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxMult_EP [i], hVxMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyMult_EP [i], hVyMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaMult_EP [i], hVxYaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbMult_EP [i], hVxYbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcMult_EP [i], hVxYcMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaMult_EP [i], hVyXaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbMult_EP [i], hVyXbMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcMult_EP [i], hVyXcMultBS_EP [i], distrDir);

            TH2toTH1withSampling (p2VxaPtCent_SP [i], hVxaPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbPtCent_SP [i], hVxbPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcPtCent_SP [i], hVxcPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaPtCent_SP [i], hVyaPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybPtCent_SP [i], hVybPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycPtCent_SP [i], hVycPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaPtCent_SP [i], hVaPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbPtCent_SP [i], hVbPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcPtCent_SP [i], hVcPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VPtCent_SP [i], hVPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxPtCent_SP [i], hVxPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyPtCent_SP [i], hVyPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaPtCent_SP [i], hVxYaPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbPtCent_SP [i], hVxYbPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcPtCent_SP [i], hVxYcPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaPtCent_SP [i], hVyXaPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbPtCent_SP [i], hVyXbPtCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcPtCent_SP [i], hVyXcPtCentBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaPtCent_EP [i], hVxaPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbPtCent_EP [i], hVxbPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcPtCent_EP [i], hVxcPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaPtCent_EP [i], hVyaPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybPtCent_EP [i], hVybPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycPtCent_EP [i], hVycPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaPtCent_EP [i], hVaPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbPtCent_EP [i], hVbPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcPtCent_EP [i], hVcPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VPtCent_EP [i], hVPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxPtCent_EP [i], hVxPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyPtCent_EP [i], hVyPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaPtCent_EP [i], hVxYaPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbPtCent_EP [i], hVxYbPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcPtCent_EP [i], hVxYcPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaPtCent_EP [i], hVyXaPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbPtCent_EP [i], hVyXbPtCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcPtCent_EP [i], hVyXcPtCentBS_EP [i], distrDir);

            TH2toTH1withSampling (p2VxaEtaCent_SP [i], hVxaEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbEtaCent_SP [i], hVxbEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcEtaCent_SP [i], hVxcEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaEtaCent_SP [i], hVyaEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybEtaCent_SP [i], hVybEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycEtaCent_SP [i], hVycEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaEtaCent_SP [i], hVaEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbEtaCent_SP [i], hVbEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcEtaCent_SP [i], hVcEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VEtaCent_SP [i], hVEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxEtaCent_SP [i], hVxEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyEtaCent_SP [i], hVyEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaEtaCent_SP [i], hVxYaEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbEtaCent_SP [i], hVxYbEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcEtaCent_SP [i], hVxYcEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaEtaCent_SP [i], hVyXaEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbEtaCent_SP [i], hVyXbEtaCentBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcEtaCent_SP [i], hVyXcEtaCentBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaEtaCent_EP [i], hVxaEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbEtaCent_EP [i], hVxbEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcEtaCent_EP [i], hVxcEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaEtaCent_EP [i], hVyaEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybEtaCent_EP [i], hVybEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycEtaCent_EP [i], hVycEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaEtaCent_EP [i], hVaEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbEtaCent_EP [i], hVbEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcEtaCent_EP [i], hVcEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VEtaCent_EP [i], hVEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxEtaCent_EP [i], hVxEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyEtaCent_EP [i], hVyEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaEtaCent_EP [i], hVxYaEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbEtaCent_EP [i], hVxYbEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcEtaCent_EP [i], hVxYcEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaEtaCent_EP [i], hVyXaEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbEtaCent_EP [i], hVyXbEtaCentBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcEtaCent_EP [i], hVyXcEtaCentBS_EP [i], distrDir);

            TH2toTH1withSampling (p2VxaPtMult_SP [i], hVxaPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbPtMult_SP [i], hVxbPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcPtMult_SP [i], hVxcPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaPtMult_SP [i], hVyaPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybPtMult_SP [i], hVybPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycPtMult_SP [i], hVycPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaPtMult_SP [i], hVaPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbPtMult_SP [i], hVbPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcPtMult_SP [i], hVcPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VPtMult_SP [i], hVPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxPtMult_SP [i], hVxPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyPtMult_SP [i], hVyPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaPtMult_SP [i], hVxYaPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbPtMult_SP [i], hVxYbPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcPtMult_SP [i], hVxYcPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaPtMult_SP [i], hVyXaPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbPtMult_SP [i], hVyXbPtMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcPtMult_SP [i], hVyXcPtMultBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaPtMult_EP [i], hVxaPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbPtMult_EP [i], hVxbPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcPtMult_EP [i], hVxcPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaPtMult_EP [i], hVyaPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybPtMult_EP [i], hVybPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycPtMult_EP [i], hVycPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaPtMult_EP [i], hVaPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbPtMult_EP [i], hVbPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcPtMult_EP [i], hVcPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VPtMult_EP [i], hVPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxPtMult_EP [i], hVxPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyPtMult_EP [i], hVyPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaPtMult_EP [i], hVxYaPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbPtMult_EP [i], hVxYbPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcPtMult_EP [i], hVxYcPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaPtMult_EP [i], hVyXaPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbPtMult_EP [i], hVyXbPtMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcPtMult_EP [i], hVyXcPtMultBS_EP [i], distrDir);

            TH2toTH1withSampling (p2VxaEtaMult_SP [i], hVxaEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxbEtaMult_SP [i], hVxbEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxcEtaMult_SP [i], hVxcEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyaEtaMult_SP [i], hVyaEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VybEtaMult_SP [i], hVybEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VycEtaMult_SP [i], hVycEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VaEtaMult_SP [i], hVaEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VbEtaMult_SP [i], hVbEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VcEtaMult_SP [i], hVcEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VEtaMult_SP [i], hVEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxEtaMult_SP [i], hVxEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyEtaMult_SP [i], hVyEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYaEtaMult_SP [i], hVxYaEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYbEtaMult_SP [i], hVxYbEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VxYcEtaMult_SP [i], hVxYcEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXaEtaMult_SP [i], hVyXaEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXbEtaMult_SP [i], hVyXbEtaMultBS_SP [i], distrDir);
            TH2toTH1withSampling (p2VyXcEtaMult_SP [i], hVyXcEtaMultBS_SP [i], distrDir);

            TH2toTH1withSampling (p2VxaEtaMult_EP [i], hVxaEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxbEtaMult_EP [i], hVxbEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxcEtaMult_EP [i], hVxcEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyaEtaMult_EP [i], hVyaEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VybEtaMult_EP [i], hVybEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VycEtaMult_EP [i], hVycEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VaEtaMult_EP [i], hVaEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VbEtaMult_EP [i], hVbEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VcEtaMult_EP [i], hVcEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VEtaMult_EP [i], hVEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxEtaMult_EP [i], hVxEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyEtaMult_EP [i], hVyEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYaEtaMult_EP [i], hVxYaEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYbEtaMult_EP [i], hVxYbEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VxYcEtaMult_EP [i], hVxYcEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXaEtaMult_EP [i], hVyXaEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXbEtaMult_EP [i], hVyXbEtaMultBS_EP [i], distrDir);
            TH2toTH1withSampling (p2VyXcEtaMult_EP [i], hVyXcEtaMultBS_EP [i], distrDir);

            // reflected rapidity bins
            ReflectRapidity (hVxaEtaCent_SP [i], hVxaEtaReflCent_SP [i], n);
            ReflectRapidity (hVxbEtaCent_SP [i], hVxbEtaReflCent_SP [i], n);
            ReflectRapidity (hVxcEtaCent_SP [i], hVxcEtaReflCent_SP [i], n);
            ReflectRapidity (hVxEtaCent_SP [i], hVxEtaReflCent_SP [i], n);
            ReflectRapidity (hVyaEtaCent_SP [i], hVyaEtaReflCent_SP [i], n);
            ReflectRapidity (hVybEtaCent_SP [i], hVybEtaReflCent_SP [i], n);
            ReflectRapidity (hVycEtaCent_SP [i], hVycEtaReflCent_SP [i], n);
            ReflectRapidity (hVyEtaCent_SP [i], hVyEtaReflCent_SP [i], n);
            ReflectRapidity (hVaEtaCent_SP [i], hVaEtaReflCent_SP [i], n);
            ReflectRapidity (hVbEtaCent_SP [i], hVbEtaReflCent_SP [i], n);
            ReflectRapidity (hVcEtaCent_SP [i], hVcEtaReflCent_SP [i], n);
            ReflectRapidity (hVEtaCent_SP [i], hVEtaReflCent_SP [i], n);
            ReflectRapidity (hVxYaEtaCent_SP [i], hVxYaEtaReflCent_SP [i], n);
            ReflectRapidity (hVxYbEtaCent_SP [i], hVxYbEtaReflCent_SP [i], n);
            ReflectRapidity (hVxYcEtaCent_SP [i], hVxYcEtaReflCent_SP [i], n);
            ReflectRapidity (hVyXaEtaCent_SP [i], hVyXaEtaReflCent_SP [i], n);
            ReflectRapidity (hVyXbEtaCent_SP [i], hVyXbEtaReflCent_SP [i], n);
            ReflectRapidity (hVyXcEtaCent_SP [i], hVyXcEtaReflCent_SP [i], n);

            ReflectRapidity (hVxaEtaCent_EP [i], hVxaEtaReflCent_EP [i], n);
            ReflectRapidity (hVxbEtaCent_EP [i], hVxbEtaReflCent_EP [i], n);
            ReflectRapidity (hVxcEtaCent_EP [i], hVxcEtaReflCent_EP [i], n);
            ReflectRapidity (hVxEtaCent_EP [i], hVxEtaReflCent_EP [i], n);
            ReflectRapidity (hVyaEtaCent_EP [i], hVyaEtaReflCent_EP [i], n);
            ReflectRapidity (hVybEtaCent_EP [i], hVybEtaReflCent_EP [i], n);
            ReflectRapidity (hVycEtaCent_EP [i], hVycEtaReflCent_EP [i], n);
            ReflectRapidity (hVyEtaCent_EP [i], hVyEtaReflCent_EP [i], n);
            ReflectRapidity (hVaEtaCent_EP [i], hVaEtaReflCent_EP [i], n);
            ReflectRapidity (hVbEtaCent_EP [i], hVbEtaReflCent_EP [i], n);
            ReflectRapidity (hVcEtaCent_EP [i], hVcEtaReflCent_EP [i], n);
            ReflectRapidity (hVEtaCent_EP [i], hVEtaReflCent_EP [i], n);
            ReflectRapidity (hVxYaEtaCent_EP [i], hVxYaEtaReflCent_EP [i], n);
            ReflectRapidity (hVxYbEtaCent_EP [i], hVxYbEtaReflCent_EP [i], n);
            ReflectRapidity (hVxYcEtaCent_EP [i], hVxYcEtaReflCent_EP [i], n);
            ReflectRapidity (hVyXaEtaCent_EP [i], hVyXaEtaReflCent_EP [i], n);
            ReflectRapidity (hVyXbEtaCent_EP [i], hVyXbEtaReflCent_EP [i], n);
            ReflectRapidity (hVyXcEtaCent_EP [i], hVyXcEtaReflCent_EP [i], n);

            ReflectRapidity (hVxaEtaMult_SP [i], hVxaEtaReflMult_SP [i], n);
            ReflectRapidity (hVxbEtaMult_SP [i], hVxbEtaReflMult_SP [i], n);
            ReflectRapidity (hVxcEtaMult_SP [i], hVxcEtaReflMult_SP [i], n);
            ReflectRapidity (hVxEtaMult_SP [i], hVxEtaReflMult_SP [i], n);
            ReflectRapidity (hVyaEtaMult_SP [i], hVyaEtaReflMult_SP [i], n);
            ReflectRapidity (hVybEtaMult_SP [i], hVybEtaReflMult_SP [i], n);
            ReflectRapidity (hVycEtaMult_SP [i], hVycEtaReflMult_SP [i], n);
            ReflectRapidity (hVyEtaMult_SP [i], hVyEtaReflMult_SP [i], n);
            ReflectRapidity (hVaEtaMult_SP [i], hVaEtaReflMult_SP [i], n);
            ReflectRapidity (hVbEtaMult_SP [i], hVbEtaReflMult_SP [i], n);
            ReflectRapidity (hVcEtaMult_SP [i], hVcEtaReflMult_SP [i], n);
            ReflectRapidity (hVEtaMult_SP [i], hVEtaReflMult_SP [i], n);
            ReflectRapidity (hVxYaEtaMult_SP [i], hVxYaEtaReflMult_SP [i], n);
            ReflectRapidity (hVxYbEtaMult_SP [i], hVxYbEtaReflMult_SP [i], n);
            ReflectRapidity (hVxYcEtaMult_SP [i], hVxYcEtaReflMult_SP [i], n);
            ReflectRapidity (hVyXaEtaMult_SP [i], hVyXaEtaReflMult_SP [i], n);
            ReflectRapidity (hVyXbEtaMult_SP [i], hVyXbEtaReflMult_SP [i], n);
            ReflectRapidity (hVyXcEtaMult_SP [i], hVyXcEtaReflMult_SP [i], n);

            ReflectRapidity (hVxaEtaMult_EP [i], hVxaEtaReflMult_EP [i], n);
            ReflectRapidity (hVxbEtaMult_EP [i], hVxbEtaReflMult_EP [i], n);
            ReflectRapidity (hVxcEtaMult_EP [i], hVxcEtaReflMult_EP [i], n);
            ReflectRapidity (hVxEtaMult_EP [i], hVxEtaReflMult_EP [i], n);
            ReflectRapidity (hVyaEtaMult_EP [i], hVyaEtaReflMult_EP [i], n);
            ReflectRapidity (hVybEtaMult_EP [i], hVybEtaReflMult_EP [i], n);
            ReflectRapidity (hVycEtaMult_EP [i], hVycEtaReflMult_EP [i], n);
            ReflectRapidity (hVyEtaMult_EP [i], hVyEtaReflMult_EP [i], n);
            ReflectRapidity (hVaEtaMult_EP [i], hVaEtaReflMult_EP [i], n);
            ReflectRapidity (hVbEtaMult_EP [i], hVbEtaReflMult_EP [i], n);
            ReflectRapidity (hVcEtaMult_EP [i], hVcEtaReflMult_EP [i], n);
            ReflectRapidity (hVEtaMult_EP [i], hVEtaReflMult_EP [i], n);
            ReflectRapidity (hVxYaEtaMult_EP [i], hVxYaEtaReflMult_EP [i], n);
            ReflectRapidity (hVxYbEtaMult_EP [i], hVxYbEtaReflMult_EP [i], n);
            ReflectRapidity (hVxYcEtaMult_EP [i], hVxYcEtaReflMult_EP [i], n);
            ReflectRapidity (hVyXaEtaMult_EP [i], hVyXaEtaReflMult_EP [i], n);
            ReflectRapidity (hVyXbEtaMult_EP [i], hVyXbEtaReflMult_EP [i], n);
            ReflectRapidity (hVyXcEtaMult_EP [i], hVyXcEtaReflMult_EP [i], n);

            ReflectRapidity (hVxaEtaCentBS_SP [i], hVxaEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxbEtaCentBS_SP [i], hVxbEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxcEtaCentBS_SP [i], hVxcEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxEtaCentBS_SP [i], hVxEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVyaEtaCentBS_SP [i], hVyaEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVybEtaCentBS_SP [i], hVybEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVycEtaCentBS_SP [i], hVycEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVyEtaCentBS_SP [i], hVyEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVaEtaCentBS_SP [i], hVaEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVbEtaCentBS_SP [i], hVbEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVcEtaCentBS_SP [i], hVcEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVEtaCentBS_SP [i], hVEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxYaEtaCentBS_SP [i], hVxYaEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxYbEtaCentBS_SP [i], hVxYbEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVxYcEtaCentBS_SP [i], hVxYcEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVyXaEtaCentBS_SP [i], hVyXaEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVyXbEtaCentBS_SP [i], hVyXbEtaReflCentBS_SP [i], n);
            ReflectRapidity (hVyXcEtaCentBS_SP [i], hVyXcEtaReflCentBS_SP [i], n);

            ReflectRapidity (hVxaEtaCentBS_EP [i], hVxaEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxbEtaCentBS_EP [i], hVxbEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxcEtaCentBS_EP [i], hVxcEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxEtaCentBS_EP [i], hVxEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVyaEtaCentBS_EP [i], hVyaEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVybEtaCentBS_EP [i], hVybEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVycEtaCentBS_EP [i], hVycEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVyEtaCentBS_EP [i], hVyEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVaEtaCentBS_EP [i], hVaEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVbEtaCentBS_EP [i], hVbEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVcEtaCentBS_EP [i], hVcEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVEtaCentBS_EP [i], hVEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxYaEtaCentBS_EP [i], hVxYaEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxYbEtaCentBS_EP [i], hVxYbEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVxYcEtaCentBS_EP [i], hVxYcEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVyXaEtaCentBS_EP [i], hVyXaEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVyXbEtaCentBS_EP [i], hVyXbEtaReflCentBS_EP [i], n);
            ReflectRapidity (hVyXcEtaCentBS_EP [i], hVyXcEtaReflCentBS_EP [i], n);

            ReflectRapidity (hVxaEtaMultBS_SP [i], hVxaEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxbEtaMultBS_SP [i], hVxbEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxcEtaMultBS_SP [i], hVxcEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxEtaMultBS_SP [i], hVxEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVyaEtaMultBS_SP [i], hVyaEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVybEtaMultBS_SP [i], hVybEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVycEtaMultBS_SP [i], hVycEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVyEtaMultBS_SP [i], hVyEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVaEtaMultBS_SP [i], hVaEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVbEtaMultBS_SP [i], hVbEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVcEtaMultBS_SP [i], hVcEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVEtaMultBS_SP [i], hVEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxYaEtaMultBS_SP [i], hVxYaEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxYbEtaMultBS_SP [i], hVxYbEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVxYcEtaMultBS_SP [i], hVxYcEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVyXaEtaMultBS_SP [i], hVyXaEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVyXbEtaMultBS_SP [i], hVyXbEtaReflMultBS_SP [i], n);
            ReflectRapidity (hVyXcEtaMultBS_SP [i], hVyXcEtaReflMultBS_SP [i], n);

            ReflectRapidity (hVxaEtaMultBS_EP [i], hVxaEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxbEtaMultBS_EP [i], hVxbEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxcEtaMultBS_EP [i], hVxcEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxEtaMultBS_EP [i], hVxEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVyaEtaMultBS_EP [i], hVyaEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVybEtaMultBS_EP [i], hVybEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVycEtaMultBS_EP [i], hVycEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVyEtaMultBS_EP [i], hVyEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVaEtaMultBS_EP [i], hVaEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVbEtaMultBS_EP [i], hVbEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVcEtaMultBS_EP [i], hVcEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVEtaMultBS_EP [i], hVEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxYaEtaMultBS_EP [i], hVxYaEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxYbEtaMultBS_EP [i], hVxYbEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVxYcEtaMultBS_EP [i], hVxYcEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVyXaEtaMultBS_EP [i], hVyXaEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVyXbEtaMultBS_EP [i], hVyXbEtaReflMultBS_EP [i], n);
            ReflectRapidity (hVyXcEtaMultBS_EP [i], hVyXcEtaReflMultBS_EP [i], n);

            if (uniformSet) { // make something with Monte-Carlo
//                pCosnPhi_PsiRPPt [i] -> SetTitle (Form ("V_{%i} versus P_{T} (reaction plane method)", n));
//                pCosnPhi_PsiRPPt [i] -> ProjectionX () -> Write (Form ("hV%iPt_RP", n));
//                pCosnPhi_PsiRPEta [i] -> SetTitle (Form ("V_{%i} versus " + varName_ + " (reaction plane method)", n));
//                pCosnPhi_PsiRPEta [i] -> ProjectionX () -> Write (Form ("hV%iEta_RP", n));
            }

            PlotKinematics (corrFile, stepDir, n, step); // hard coded

            // PLOT CORRELATIONS AND RESOLUTION//

            hList1 [0] = pXaXbCent_SP [i];
            hList1 [1] = pYaYbCent_SP [i];
            hList1 [2] = pXaYbCent_SP [i];
            hList1 [3] = pYaXbCent_SP [i];
            hList1 [4] = pXaXcCent_SP [i];
            hList1 [5] = pYaYcCent_SP [i];
            hList1 [6] = pXaYcCent_SP [i];
            hList1 [7] = pYaXcCent_SP [i];
            hList1 [8] = pXbXcCent_SP [i];
            hList1 [9] = pYbYcCent_SP [i];
            hList1 [10] = pXbYcCent_SP [i];
            hList1 [11] = pYbXcCent_SP [i];

            hList2 [0] = pXaXbCentBS_SP [i];
            hList2 [1] = pYaYbCentBS_SP [i];
            hList2 [2] = pXaYbCentBS_SP [i];
            hList2 [3] = pYaXbCentBS_SP [i];
            hList2 [4] = pXaXcCentBS_SP [i];
            hList2 [5] = pYaYcCentBS_SP [i];
            hList2 [6] = pXaYcCentBS_SP [i];
            hList2 [7] = pYaXcCentBS_SP [i];
            hList2 [8] = pXbXcCentBS_SP [i];
            hList2 [9] = pYbYcCentBS_SP [i];
            hList2 [10] = pXbYcCentBS_SP [i];
            hList2 [11] = pYbXcCentBS_SP [i];

            hList3 [0] = hRxaCent_SP [i];
            hList3 [1] = hRyaCent_SP [i];
            hList3 [2] = hRaCent_SP [i];
            hList3 [3] = hRxbCent_SP [i];
            hList3 [4] = hRybCent_SP [i];
            hList3 [5] = hRbCent_SP [i];
            hList3 [6] = hRxcCent_SP [i];
            hList3 [7] = hRycCent_SP [i];
            hList3 [8] = hRcCent_SP [i];

            hList4 [0] = hRxaCentBS_SP [i];
            hList4 [1] = hRyaCentBS_SP [i];
            hList4 [2] = hRaCentBS_SP [i];
            hList4 [3] = hRxbCentBS_SP [i];
            hList4 [4] = hRybCentBS_SP [i];
            hList4 [5] = hRbCentBS_SP [i];
            hList4 [6] = hRxcCentBS_SP [i];
            hList4 [7] = hRycCentBS_SP [i];
            hList4 [8] = hRcCentBS_SP [i];

            PlotResolution (hList1, hList2, hList3, hList4, nProfs, stepDir);

            hList1 [0] = pXaXbCent_EP [i];
            hList1 [1] = pYaYbCent_EP [i];
            hList1 [2] = pXaYbCent_EP [i];
            hList1 [3] = pYaXbCent_EP [i];
            hList1 [4] = pXaXcCent_EP [i];
            hList1 [5] = pYaYcCent_EP [i];
            hList1 [6] = pXaYcCent_EP [i];
            hList1 [7] = pYaXcCent_EP [i];
            hList1 [8] = pXbXcCent_EP [i];
            hList1 [9] = pYbYcCent_EP [i];
            hList1 [10] = pXbYcCent_EP [i];
            hList1 [11] = pYbXcCent_EP [i];

            hList2 [0] = pXaXbCentBS_EP [i];
            hList2 [1] = pYaYbCentBS_EP [i];
            hList2 [2] = pXaYbCentBS_EP [i];
            hList2 [3] = pYaXbCentBS_EP [i];
            hList2 [4] = pXaXcCentBS_EP [i];
            hList2 [5] = pYaYcCentBS_EP [i];
            hList2 [6] = pXaYcCentBS_EP [i];
            hList2 [7] = pYaXcCentBS_EP [i];
            hList2 [8] = pXbXcCentBS_EP [i];
            hList2 [9] = pYbYcCentBS_EP [i];
            hList2 [10] = pXbYcCentBS_EP [i];
            hList2 [11] = pYbXcCentBS_EP [i];

            hList3 [0] = hRxaCent_EP [i];
            hList3 [1] = hRyaCent_EP [i];
            hList3 [2] = hRaCent_EP [i];
            hList3 [3] = hRxbCent_EP [i];
            hList3 [4] = hRybCent_EP [i];
            hList3 [5] = hRbCent_EP [i];
            hList3 [6] = hRxcCent_EP [i];
            hList3 [7] = hRycCent_EP [i];
            hList3 [8] = hRcCent_EP [i];

            hList4 [0] = hRxaCentBS_EP [i];
            hList4 [1] = hRyaCentBS_EP [i];
            hList4 [2] = hRaCentBS_EP [i];
            hList4 [3] = hRxbCentBS_EP [i];
            hList4 [4] = hRybCentBS_EP [i];
            hList4 [5] = hRbCentBS_EP [i];
            hList4 [6] = hRxcCentBS_EP [i];
            hList4 [7] = hRycCentBS_EP [i];
            hList4 [8] = hRcCentBS_EP [i];

            PlotResolution (hList1, hList2, hList3, hList4, nProfs, stepDir);

            hList1 [0] = hVxCent_SP [i];
            hList1 [1] = hVxaCent_SP [i];
            hList1 [2] = hVxbCent_SP [i];
            hList1 [3] = hVxcCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxCentBS_SP [i];
            hList2 [1] = hVxaCentBS_SP [i];
            hList2 [2] = hVxbCentBS_SP [i];
            hList2 [3] = hVxcCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyCent_SP [i];
            hList1 [1] = hVyaCent_SP [i];
            hList1 [2] = hVybCent_SP [i];
            hList1 [3] = hVycCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyCentBS_SP [i];
            hList2 [1] = hVyaCentBS_SP [i];
            hList2 [2] = hVybCentBS_SP [i];
            hList2 [3] = hVycCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVCent_SP [i];
            hList1 [1] = hVaCent_SP [i];
            hList1 [2] = hVbCent_SP [i];
            hList1 [3] = hVcCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVCentBS_SP [i];
            hList2 [1] = hVaCentBS_SP [i];
            hList2 [2] = hVbCentBS_SP [i];
            hList2 [3] = hVcCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaCent_SP [i];
            hList1 [1] = hVxaCent_SP [i];
            hList1 [2] = hVyaCent_SP [i];
            hList1 [3] = hVxYaCent_SP [i];
            hList1 [4] = hVyXaCent_SP [i];

            hList2 [0] = hVaCentBS_SP [i];
            hList2 [1] = hVxaCentBS_SP [i];
            hList2 [2] = hVyaCentBS_SP [i];
            hList2 [3] = hVxYaCentBS_SP [i];
            hList2 [4] = hVyXaCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbCent_SP [i];
            hList1 [1] = hVxbCent_SP [i];
            hList1 [2] = hVybCent_SP [i];
            hList1 [3] = hVxYbCent_SP [i];
            hList1 [4] = hVyXbCent_SP [i];

            hList2 [0] = hVbCentBS_SP [i];
            hList2 [1] = hVxbCentBS_SP [i];
            hList2 [2] = hVybCentBS_SP [i];
            hList2 [3] = hVxYbCentBS_SP [i];
            hList2 [4] = hVyXbCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcCent_SP [i];
            hList1 [1] = hVxcCent_SP [i];
            hList1 [2] = hVycCent_SP [i];
            hList1 [3] = hVxYcCent_SP [i];
            hList1 [4] = hVyXcCent_SP [i];

            hList2 [0] = hVcCentBS_SP [i];
            hList2 [1] = hVxcCentBS_SP [i];
            hList2 [2] = hVycCentBS_SP [i];
            hList2 [3] = hVxYcCentBS_SP [i];
            hList2 [4] = hVyXcCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxCent_EP [i];
            hList1 [1] = hVxaCent_EP [i];
            hList1 [2] = hVxbCent_EP [i];
            hList1 [3] = hVxcCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxCentBS_EP [i];
            hList2 [1] = hVxaCentBS_EP [i];
            hList2 [2] = hVxbCentBS_EP [i];
            hList2 [3] = hVxcCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyCent_EP [i];
            hList1 [1] = hVyaCent_EP [i];
            hList1 [2] = hVybCent_EP [i];
            hList1 [3] = hVycCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyCentBS_EP [i];
            hList2 [1] = hVyaCentBS_EP [i];
            hList2 [2] = hVybCentBS_EP [i];
            hList2 [3] = hVycCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVCent_EP [i];
            hList1 [1] = hVaCent_EP [i];
            hList1 [2] = hVbCent_EP [i];
            hList1 [3] = hVcCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVCentBS_EP [i];
            hList2 [1] = hVaCentBS_EP [i];
            hList2 [2] = hVbCentBS_EP [i];
            hList2 [3] = hVcCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaCent_EP [i];
            hList1 [1] = hVxaCent_EP [i];
            hList1 [2] = hVyaCent_EP [i];
            hList1 [3] = hVxYaCent_EP [i];
            hList1 [4] = hVyXaCent_EP [i];

            hList2 [0] = hVaCentBS_EP [i];
            hList2 [1] = hVxaCentBS_EP [i];
            hList2 [2] = hVyaCentBS_EP [i];
            hList2 [3] = hVxYaCentBS_EP [i];
            hList2 [4] = hVyXaCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbCent_EP [i];
            hList1 [1] = hVxbCent_EP [i];
            hList1 [2] = hVybCent_EP [i];
            hList1 [3] = hVxYbCent_EP [i];
            hList1 [4] = hVyXbCent_EP [i];

            hList2 [0] = hVbCentBS_EP [i];
            hList2 [1] = hVxbCentBS_EP [i];
            hList2 [2] = hVybCentBS_EP [i];
            hList2 [3] = hVxYbCentBS_EP [i];
            hList2 [4] = hVyXbCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcCent_EP [i];
            hList1 [1] = hVxcCent_EP [i];
            hList1 [2] = hVycCent_EP [i];
            hList1 [3] = hVxYcCent_EP [i];
            hList1 [4] = hVyXcCent_EP [i];

            hList2 [0] = hVcCentBS_EP [i];
            hList2 [1] = hVxcCentBS_EP [i];
            hList2 [2] = hVycCentBS_EP [i];
            hList2 [3] = hVxYcCentBS_EP [i];
            hList2 [4] = hVyXcCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxMult_SP [i];
            hList1 [1] = hVxaMult_SP [i];
            hList1 [2] = hVxbMult_SP [i];
            hList1 [3] = hVxcMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxMultBS_SP [i];
            hList2 [1] = hVxaMultBS_SP [i];
            hList2 [2] = hVxbMultBS_SP [i];
            hList2 [3] = hVxcMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyMult_SP [i];
            hList1 [1] = hVyaMult_SP [i];
            hList1 [2] = hVybMult_SP [i];
            hList1 [3] = hVycMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyMultBS_SP [i];
            hList2 [1] = hVyaMultBS_SP [i];
            hList2 [2] = hVybMultBS_SP [i];
            hList2 [3] = hVycMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVMult_SP [i];
            hList1 [1] = hVaMult_SP [i];
            hList1 [2] = hVbMult_SP [i];
            hList1 [3] = hVcMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVMultBS_SP [i];
            hList2 [1] = hVaMultBS_SP [i];
            hList2 [2] = hVbMultBS_SP [i];
            hList2 [3] = hVcMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaMult_SP [i];
            hList1 [1] = hVxaMult_SP [i];
            hList1 [2] = hVyaMult_SP [i];
            hList1 [3] = hVxYaMult_SP [i];
            hList1 [4] = hVyXaMult_SP [i];

            hList2 [0] = hVaMultBS_SP [i];
            hList2 [1] = hVxaMultBS_SP [i];
            hList2 [2] = hVyaMultBS_SP [i];
            hList2 [3] = hVxYaMultBS_SP [i];
            hList2 [4] = hVyXaMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbMult_SP [i];
            hList1 [1] = hVxbMult_SP [i];
            hList1 [2] = hVybMult_SP [i];
            hList1 [3] = hVxYbMult_SP [i];
            hList1 [4] = hVyXbMult_SP [i];

            hList2 [0] = hVbMultBS_SP [i];
            hList2 [1] = hVxbMultBS_SP [i];
            hList2 [2] = hVybMultBS_SP [i];
            hList2 [3] = hVxYbMultBS_SP [i];
            hList2 [4] = hVyXbMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcMult_SP [i];
            hList1 [1] = hVxcMult_SP [i];
            hList1 [2] = hVycMult_SP [i];
            hList1 [3] = hVxYcMult_SP [i];
            hList1 [4] = hVyXcMult_SP [i];

            hList2 [0] = hVcMultBS_SP [i];
            hList2 [1] = hVxcMultBS_SP [i];
            hList2 [2] = hVycMultBS_SP [i];
            hList2 [3] = hVxYcMultBS_SP [i];
            hList2 [4] = hVyXcMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxMult_EP [i];
            hList1 [1] = hVxaMult_EP [i];
            hList1 [2] = hVxbMult_EP [i];
            hList1 [3] = hVxcMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxMultBS_EP [i];
            hList2 [1] = hVxaMultBS_EP [i];
            hList2 [2] = hVxbMultBS_EP [i];
            hList2 [3] = hVxcMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyMult_EP [i];
            hList1 [1] = hVyaMult_EP [i];
            hList1 [2] = hVybMult_EP [i];
            hList1 [3] = hVycMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyMultBS_EP [i];
            hList2 [1] = hVyaMultBS_EP [i];
            hList2 [2] = hVybMultBS_EP [i];
            hList2 [3] = hVycMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVMult_EP [i];
            hList1 [1] = hVaMult_EP [i];
            hList1 [2] = hVbMult_EP [i];
            hList1 [3] = hVcMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVMultBS_EP [i];
            hList2 [1] = hVaMultBS_EP [i];
            hList2 [2] = hVbMultBS_EP [i];
            hList2 [3] = hVcMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaMult_EP [i];
            hList1 [1] = hVxaMult_EP [i];
            hList1 [2] = hVyaMult_EP [i];
            hList1 [3] = hVxYaMult_EP [i];
            hList1 [4] = hVyXaMult_EP [i];

            hList2 [0] = hVaMultBS_EP [i];
            hList2 [1] = hVxaMultBS_EP [i];
            hList2 [2] = hVyaMultBS_EP [i];
            hList2 [3] = hVxYaMultBS_EP [i];
            hList2 [4] = hVyXaMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbMult_EP [i];
            hList1 [1] = hVxbMult_EP [i];
            hList1 [2] = hVybMult_EP [i];
            hList1 [3] = hVxYbMult_EP [i];
            hList1 [4] = hVyXbMult_EP [i];

            hList2 [0] = hVbMultBS_EP [i];
            hList2 [1] = hVxbMultBS_EP [i];
            hList2 [2] = hVybMultBS_EP [i];
            hList2 [3] = hVxYbMultBS_EP [i];
            hList2 [4] = hVyXbMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcMult_EP [i];
            hList1 [1] = hVxcMult_EP [i];
            hList1 [2] = hVycMult_EP [i];
            hList1 [3] = hVxYcMult_EP [i];
            hList1 [4] = hVyXcMult_EP [i];

            hList2 [0] = hVcMultBS_EP [i];
            hList2 [1] = hVxcMultBS_EP [i];
            hList2 [2] = hVycMultBS_EP [i];
            hList2 [3] = hVxYcMultBS_EP [i];
            hList2 [4] = hVyXcMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxPtCent_SP [i];
            hList1 [1] = hVxaPtCent_SP [i];
            hList1 [2] = hVxbPtCent_SP [i];
            hList1 [3] = hVxcPtCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxPtCentBS_SP [i];
            hList2 [1] = hVxaPtCentBS_SP [i];
            hList2 [2] = hVxbPtCentBS_SP [i];
            hList2 [3] = hVxcPtCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyPtCent_SP [i];
            hList1 [1] = hVyaPtCent_SP [i];
            hList1 [2] = hVybPtCent_SP [i];
            hList1 [3] = hVycPtCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyPtCentBS_SP [i];
            hList2 [1] = hVyaPtCentBS_SP [i];
            hList2 [2] = hVybPtCentBS_SP [i];
            hList2 [3] = hVycPtCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVPtCent_SP [i];
            hList1 [1] = hVaPtCent_SP [i];
            hList1 [2] = hVbPtCent_SP [i];
            hList1 [3] = hVcPtCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVPtCentBS_SP [i];
            hList2 [1] = hVaPtCentBS_SP [i];
            hList2 [2] = hVbPtCentBS_SP [i];
            hList2 [3] = hVcPtCentBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaPtCent_SP [i];
            hList1 [1] = hVxaPtCent_SP [i];
            hList1 [2] = hVyaPtCent_SP [i];
            hList1 [3] = hVxYaPtCent_SP [i];
            hList1 [4] = hVyXaPtCent_SP [i];

            hList2 [0] = hVaPtCentBS_SP [i];
            hList2 [1] = hVxaPtCentBS_SP [i];
            hList2 [2] = hVyaPtCentBS_SP [i];
            hList2 [3] = hVxYaPtCentBS_SP [i];
            hList2 [4] = hVyXaPtCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbPtCent_SP [i];
            hList1 [1] = hVxbPtCent_SP [i];
            hList1 [2] = hVybPtCent_SP [i];
            hList1 [3] = hVxYbPtCent_SP [i];
            hList1 [4] = hVyXbPtCent_SP [i];

            hList2 [0] = hVbPtCentBS_SP [i];
            hList2 [1] = hVxbPtCentBS_SP [i];
            hList2 [2] = hVybPtCentBS_SP [i];
            hList2 [3] = hVxYbPtCentBS_SP [i];
            hList2 [4] = hVyXbPtCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcPtCent_SP [i];
            hList1 [1] = hVxcPtCent_SP [i];
            hList1 [2] = hVycPtCent_SP [i];
            hList1 [3] = hVxYcPtCent_SP [i];
            hList1 [4] = hVyXcPtCent_SP [i];

            hList2 [0] = hVcPtCentBS_SP [i];
            hList2 [1] = hVxcPtCentBS_SP [i];
            hList2 [2] = hVycPtCentBS_SP [i];
            hList2 [3] = hVxYcPtCentBS_SP [i];
            hList2 [4] = hVyXcPtCentBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxPtCent_EP [i];
            hList1 [1] = hVxaPtCent_EP [i];
            hList1 [2] = hVxbPtCent_EP [i];
            hList1 [3] = hVxcPtCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxPtCentBS_EP [i];
            hList2 [1] = hVxaPtCentBS_EP [i];
            hList2 [2] = hVxbPtCentBS_EP [i];
            hList2 [3] = hVxcPtCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyPtCent_EP [i];
            hList1 [1] = hVyaPtCent_EP [i];
            hList1 [2] = hVybPtCent_EP [i];
            hList1 [3] = hVycPtCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyPtCentBS_EP [i];
            hList2 [1] = hVyaPtCentBS_EP [i];
            hList2 [2] = hVybPtCentBS_EP [i];
            hList2 [3] = hVycPtCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVPtCent_EP [i];
            hList1 [1] = hVaPtCent_EP [i];
            hList1 [2] = hVbPtCent_EP [i];
            hList1 [3] = hVcPtCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVPtCentBS_EP [i];
            hList2 [1] = hVaPtCentBS_EP [i];
            hList2 [2] = hVbPtCentBS_EP [i];
            hList2 [3] = hVcPtCentBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaPtCent_EP [i];
            hList1 [1] = hVxaPtCent_EP [i];
            hList1 [2] = hVyaPtCent_EP [i];
            hList1 [3] = hVxYaPtCent_EP [i];
            hList1 [4] = hVyXaPtCent_EP [i];

            hList2 [0] = hVaPtCentBS_EP [i];
            hList2 [1] = hVxaPtCentBS_EP [i];
            hList2 [2] = hVyaPtCentBS_EP [i];
            hList2 [3] = hVxYaPtCentBS_EP [i];
            hList2 [4] = hVyXaPtCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbPtCent_EP [i];
            hList1 [1] = hVxbPtCent_EP [i];
            hList1 [2] = hVybPtCent_EP [i];
            hList1 [3] = hVxYbPtCent_EP [i];
            hList1 [4] = hVyXbPtCent_EP [i];

            hList2 [0] = hVbPtCentBS_EP [i];
            hList2 [1] = hVxbPtCentBS_EP [i];
            hList2 [2] = hVybPtCentBS_EP [i];
            hList2 [3] = hVxYbPtCentBS_EP [i];
            hList2 [4] = hVyXbPtCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcPtCent_EP [i];
            hList1 [1] = hVxcPtCent_EP [i];
            hList1 [2] = hVycPtCent_EP [i];
            hList1 [3] = hVxYcPtCent_EP [i];
            hList1 [4] = hVyXcPtCent_EP [i];

            hList2 [0] = hVcPtCentBS_EP [i];
            hList2 [1] = hVxcPtCentBS_EP [i];
            hList2 [2] = hVycPtCentBS_EP [i];
            hList2 [3] = hVxYcPtCentBS_EP [i];
            hList2 [4] = hVyXcPtCentBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxPtMult_SP [i];
            hList1 [1] = hVxaPtMult_SP [i];
            hList1 [2] = hVxbPtMult_SP [i];
            hList1 [3] = hVxcPtMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxPtMultBS_SP [i];
            hList2 [1] = hVxaPtMultBS_SP [i];
            hList2 [2] = hVxbPtMultBS_SP [i];
            hList2 [3] = hVxcPtMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyPtMult_SP [i];
            hList1 [1] = hVyaPtMult_SP [i];
            hList1 [2] = hVybPtMult_SP [i];
            hList1 [3] = hVycPtMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyPtMultBS_SP [i];
            hList2 [1] = hVyaPtMultBS_SP [i];
            hList2 [2] = hVybPtMultBS_SP [i];
            hList2 [3] = hVycPtMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVPtMult_SP [i];
            hList1 [1] = hVaPtMult_SP [i];
            hList1 [2] = hVbPtMult_SP [i];
            hList1 [3] = hVcPtMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVPtMultBS_SP [i];
            hList2 [1] = hVaPtMultBS_SP [i];
            hList2 [2] = hVbPtMultBS_SP [i];
            hList2 [3] = hVcPtMultBS_SP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaPtMult_SP [i];
            hList1 [1] = hVxaPtMult_SP [i];
            hList1 [2] = hVyaPtMult_SP [i];
            hList1 [3] = hVxYaPtMult_SP [i];
            hList1 [4] = hVyXaPtMult_SP [i];

            hList2 [0] = hVaPtMultBS_SP [i];
            hList2 [1] = hVxaPtMultBS_SP [i];
            hList2 [2] = hVyaPtMultBS_SP [i];
            hList2 [3] = hVxYaPtMultBS_SP [i];
            hList2 [4] = hVyXaPtMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbPtMult_SP [i];
            hList1 [1] = hVxbPtMult_SP [i];
            hList1 [2] = hVybPtMult_SP [i];
            hList1 [3] = hVxYbPtMult_SP [i];
            hList1 [4] = hVyXbPtMult_SP [i];

            hList2 [0] = hVbPtMultBS_SP [i];
            hList2 [1] = hVxbPtMultBS_SP [i];
            hList2 [2] = hVybPtMultBS_SP [i];
            hList2 [3] = hVxYbPtMultBS_SP [i];
            hList2 [4] = hVyXbPtMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcPtMult_SP [i];
            hList1 [1] = hVxcPtMult_SP [i];
            hList1 [2] = hVycPtMult_SP [i];
            hList1 [3] = hVxYcPtMult_SP [i];
            hList1 [4] = hVyXcPtMult_SP [i];

            hList2 [0] = hVcPtMultBS_SP [i];
            hList2 [1] = hVxcPtMultBS_SP [i];
            hList2 [2] = hVycPtMultBS_SP [i];
            hList2 [3] = hVxYcPtMultBS_SP [i];
            hList2 [4] = hVyXcPtMultBS_SP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxPtMult_EP [i];
            hList1 [1] = hVxaPtMult_EP [i];
            hList1 [2] = hVxbPtMult_EP [i];
            hList1 [3] = hVxcPtMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxPtMultBS_EP [i];
            hList2 [1] = hVxaPtMultBS_EP [i];
            hList2 [2] = hVxbPtMultBS_EP [i];
            hList2 [3] = hVxcPtMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVyPtMult_EP [i];
            hList1 [1] = hVyaPtMult_EP [i];
            hList1 [2] = hVybPtMult_EP [i];
            hList1 [3] = hVycPtMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyPtMultBS_EP [i];
            hList2 [1] = hVyaPtMultBS_EP [i];
            hList2 [2] = hVybPtMultBS_EP [i];
            hList2 [3] = hVycPtMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVPtMult_EP [i];
            hList1 [1] = hVaPtMult_EP [i];
            hList1 [2] = hVbPtMult_EP [i];
            hList1 [3] = hVcPtMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVPtMultBS_EP [i];
            hList2 [1] = hVaPtMultBS_EP [i];
            hList2 [2] = hVbPtMultBS_EP [i];
            hList2 [3] = hVcPtMultBS_EP [i];
            hList2 [4] = 0;

            PlotFlow (hList1, hList2);

            hList1 [0] = hVaPtMult_EP [i];
            hList1 [1] = hVxaPtMult_EP [i];
            hList1 [2] = hVyaPtMult_EP [i];
            hList1 [3] = hVxYaPtMult_EP [i];
            hList1 [4] = hVyXaPtMult_EP [i];

            hList2 [0] = hVaPtMultBS_EP [i];
            hList2 [1] = hVxaPtMultBS_EP [i];
            hList2 [2] = hVyaPtMultBS_EP [i];
            hList2 [3] = hVxYaPtMultBS_EP [i];
            hList2 [4] = hVyXaPtMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVbPtMult_EP [i];
            hList1 [1] = hVxbPtMult_EP [i];
            hList1 [2] = hVybPtMult_EP [i];
            hList1 [3] = hVxYbPtMult_EP [i];
            hList1 [4] = hVyXbPtMult_EP [i];

            hList2 [0] = hVbPtMultBS_EP [i];
            hList2 [1] = hVxbPtMultBS_EP [i];
            hList2 [2] = hVybPtMultBS_EP [i];
            hList2 [3] = hVxYbPtMultBS_EP [i];
            hList2 [4] = hVyXbPtMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVcPtMult_EP [i];
            hList1 [1] = hVxcPtMult_EP [i];
            hList1 [2] = hVycPtMult_EP [i];
            hList1 [3] = hVxYcPtMult_EP [i];
            hList1 [4] = hVyXcPtMult_EP [i];

            hList2 [0] = hVcPtMultBS_EP [i];
            hList2 [1] = hVxcPtMultBS_EP [i];
            hList2 [2] = hVycPtMultBS_EP [i];
            hList2 [3] = hVxYcPtMultBS_EP [i];
            hList2 [4] = hVyXcPtMultBS_EP [i];

            PlotFlow (hList1, hList2);

            hList1 [0] = hVxEtaCent_SP [i];
            hList1 [1] = hVxaEtaCent_SP [i];
            hList1 [2] = hVxbEtaCent_SP [i];
            hList1 [3] = hVxcEtaCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxEtaCentBS_SP [i];
            hList2 [1] = hVxaEtaCentBS_SP [i];
            hList2 [2] = hVxbEtaCentBS_SP [i];
            hList2 [3] = hVxcEtaCentBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVxEtaReflCent_SP [i];
            hList3 [1] = hVxaEtaReflCent_SP [i];
            hList3 [2] = hVxbEtaReflCent_SP [i];
            hList3 [3] = hVxcEtaReflCent_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVxEtaReflCentBS_SP [i];
            hList4 [1] = hVxaEtaReflCentBS_SP [i];
            hList4 [2] = hVxbEtaReflCentBS_SP [i];
            hList4 [3] = hVxcEtaReflCentBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VxEtaCent_SP

            hList1 [0] = hVyEtaCent_SP [i];
            hList1 [1] = hVyaEtaCent_SP [i];
            hList1 [2] = hVybEtaCent_SP [i];
            hList1 [3] = hVycEtaCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyEtaCentBS_SP [i];
            hList2 [1] = hVyaEtaCentBS_SP [i];
            hList2 [2] = hVybEtaCentBS_SP [i];
            hList2 [3] = hVycEtaCentBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVyEtaReflCent_SP [i];
            hList3 [1] = hVyaEtaReflCent_SP [i];
            hList3 [2] = hVybEtaReflCent_SP [i];
            hList3 [3] = hVycEtaReflCent_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVyEtaReflCentBS_SP [i];
            hList4 [1] = hVyaEtaReflCentBS_SP [i];
            hList4 [2] = hVybEtaReflCentBS_SP [i];
            hList4 [3] = hVycEtaReflCentBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VyEtaCent_SP

            hList1 [0] = hVEtaCent_SP [i];
            hList1 [1] = hVaEtaCent_SP [i];
            hList1 [2] = hVbEtaCent_SP [i];
            hList1 [3] = hVcEtaCent_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVEtaCentBS_SP [i];
            hList2 [1] = hVaEtaCentBS_SP [i];
            hList2 [2] = hVbEtaCentBS_SP [i];
            hList2 [3] = hVcEtaCentBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVEtaReflCent_SP [i];
            hList3 [1] = hVaEtaReflCent_SP [i];
            hList3 [2] = hVbEtaReflCent_SP [i];
            hList3 [3] = hVcEtaReflCent_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVEtaReflCentBS_SP [i];
            hList4 [1] = hVaEtaReflCentBS_SP [i];
            hList4 [2] = hVbEtaReflCentBS_SP [i];
            hList4 [3] = hVcEtaReflCentBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VEtaCent_SP

            hList1 [0] = hVaEtaCent_SP [i];
            hList1 [1] = hVxaEtaCent_SP [i];
            hList1 [2] = hVyaEtaCent_SP [i];
            hList1 [3] = hVxYaEtaCent_SP [i];
            hList1 [4] = hVyXaEtaCent_SP [i];

            hList2 [0] = hVaEtaCentBS_SP [i];
            hList2 [1] = hVxaEtaCentBS_SP [i];
            hList2 [2] = hVyaEtaCentBS_SP [i];
            hList2 [3] = hVxYaEtaCentBS_SP [i];
            hList2 [4] = hVyXaEtaCentBS_SP [i];

            hList3 [0] = hVaEtaReflCent_SP [i];
            hList3 [1] = hVxaEtaReflCent_SP [i];
            hList3 [2] = hVyaEtaReflCent_SP [i];
            hList3 [3] = hVxYaEtaReflCent_SP [i];
            hList3 [4] = hVyXaEtaReflCent_SP [i];

            hList4 [0] = hVaEtaReflCentBS_SP [i];
            hList4 [1] = hVxaEtaReflCentBS_SP [i];
            hList4 [2] = hVyaEtaReflCentBS_SP [i];
            hList4 [3] = hVxYaEtaReflCentBS_SP [i];
            hList4 [4] = hVyXaEtaReflCentBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VaEtaCent_SP

            hList1 [0] = hVbEtaCent_SP [i];
            hList1 [1] = hVxbEtaCent_SP [i];
            hList1 [2] = hVybEtaCent_SP [i];
            hList1 [3] = hVxYbEtaCent_SP [i];
            hList1 [4] = hVyXbEtaCent_SP [i];

            hList2 [0] = hVbEtaCentBS_SP [i];
            hList2 [1] = hVxbEtaCentBS_SP [i];
            hList2 [2] = hVybEtaCentBS_SP [i];
            hList2 [3] = hVxYbEtaCentBS_SP [i];
            hList2 [4] = hVyXbEtaCentBS_SP [i];

            hList3 [0] = hVbEtaReflCent_SP [i];
            hList3 [1] = hVxbEtaReflCent_SP [i];
            hList3 [2] = hVybEtaReflCent_SP [i];
            hList3 [3] = hVxYbEtaReflCent_SP [i];
            hList3 [4] = hVyXbEtaReflCent_SP [i];

            hList4 [0] = hVbEtaReflCentBS_SP [i];
            hList4 [1] = hVxbEtaReflCentBS_SP [i];
            hList4 [2] = hVybEtaReflCentBS_SP [i];
            hList4 [3] = hVxYbEtaReflCentBS_SP [i];
            hList4 [4] = hVyXbEtaReflCentBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VbEtaCent_SP

            hList1 [0] = hVcEtaCent_SP [i];
            hList1 [1] = hVxcEtaCent_SP [i];
            hList1 [2] = hVycEtaCent_SP [i];
            hList1 [3] = hVxYcEtaCent_SP [i];
            hList1 [4] = hVyXcEtaCent_SP [i];

            hList2 [0] = hVcEtaCentBS_SP [i];
            hList2 [1] = hVxcEtaCentBS_SP [i];
            hList2 [2] = hVycEtaCentBS_SP [i];
            hList2 [3] = hVxYcEtaCentBS_SP [i];
            hList2 [4] = hVyXcEtaCentBS_SP [i];

            hList3 [0] = hVcEtaReflCent_SP [i];
            hList3 [1] = hVxcEtaReflCent_SP [i];
            hList3 [2] = hVycEtaReflCent_SP [i];
            hList3 [3] = hVxYcEtaReflCent_SP [i];
            hList3 [4] = hVyXcEtaReflCent_SP [i];

            hList4 [0] = hVcEtaReflCentBS_SP [i];
            hList4 [1] = hVxcEtaReflCentBS_SP [i];
            hList4 [2] = hVycEtaReflCentBS_SP [i];
            hList4 [3] = hVxYcEtaReflCentBS_SP [i];
            hList4 [4] = hVyXcEtaReflCentBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VcEtaCent_SP

            hList1 [0] = hVxEtaCent_EP [i];
            hList1 [1] = hVxaEtaCent_EP [i];
            hList1 [2] = hVxbEtaCent_EP [i];
            hList1 [3] = hVxcEtaCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxEtaCentBS_EP [i];
            hList2 [1] = hVxaEtaCentBS_EP [i];
            hList2 [2] = hVxbEtaCentBS_EP [i];
            hList2 [3] = hVxcEtaCentBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVxEtaReflCent_EP [i];
            hList3 [1] = hVxaEtaReflCent_EP [i];
            hList3 [2] = hVxbEtaReflCent_EP [i];
            hList3 [3] = hVxcEtaReflCent_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVxEtaReflCentBS_EP [i];
            hList4 [1] = hVxaEtaReflCentBS_EP [i];
            hList4 [2] = hVxbEtaReflCentBS_EP [i];
            hList4 [3] = hVxcEtaReflCentBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VxEtaCent_EP

            hList1 [0] = hVyEtaCent_EP [i];
            hList1 [1] = hVyaEtaCent_EP [i];
            hList1 [2] = hVybEtaCent_EP [i];
            hList1 [3] = hVycEtaCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyEtaCentBS_EP [i];
            hList2 [1] = hVyaEtaCentBS_EP [i];
            hList2 [2] = hVybEtaCentBS_EP [i];
            hList2 [3] = hVycEtaCentBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVyEtaReflCent_EP [i];
            hList3 [1] = hVyaEtaReflCent_EP [i];
            hList3 [2] = hVybEtaReflCent_EP [i];
            hList3 [3] = hVycEtaReflCent_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVyEtaReflCentBS_EP [i];
            hList4 [1] = hVyaEtaReflCentBS_EP [i];
            hList4 [2] = hVybEtaReflCentBS_EP [i];
            hList4 [3] = hVycEtaReflCentBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VyEtaCent_EP

            hList1 [0] = hVEtaCent_EP [i];
            hList1 [1] = hVaEtaCent_EP [i];
            hList1 [2] = hVbEtaCent_EP [i];
            hList1 [3] = hVcEtaCent_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVEtaCentBS_EP [i];
            hList2 [1] = hVaEtaCentBS_EP [i];
            hList2 [2] = hVbEtaCentBS_EP [i];
            hList2 [3] = hVcEtaCentBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVEtaReflCent_EP [i];
            hList3 [1] = hVaEtaReflCent_EP [i];
            hList3 [2] = hVbEtaReflCent_EP [i];
            hList3 [3] = hVcEtaReflCent_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVEtaReflCentBS_EP [i];
            hList4 [1] = hVaEtaReflCentBS_EP [i];
            hList4 [2] = hVbEtaReflCentBS_EP [i];
            hList4 [3] = hVcEtaReflCentBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VEtaCent_EP

            hList1 [0] = hVaEtaCent_EP [i];
            hList1 [1] = hVxaEtaCent_EP [i];
            hList1 [2] = hVyaEtaCent_EP [i];
            hList1 [3] = hVxYaEtaCent_EP [i];
            hList1 [4] = hVyXaEtaCent_EP [i];

            hList2 [0] = hVaEtaCentBS_EP [i];
            hList2 [1] = hVxaEtaCentBS_EP [i];
            hList2 [2] = hVyaEtaCentBS_EP [i];
            hList2 [3] = hVxYaEtaCentBS_EP [i];
            hList2 [4] = hVyXaEtaCentBS_EP [i];

            hList3 [0] = hVaEtaReflCent_EP [i];
            hList3 [1] = hVxaEtaReflCent_EP [i];
            hList3 [2] = hVyaEtaReflCent_EP [i];
            hList3 [3] = hVxYaEtaReflCent_EP [i];
            hList3 [4] = hVyXaEtaReflCent_EP [i];

            hList4 [0] = hVaEtaReflCentBS_EP [i];
            hList4 [1] = hVxaEtaReflCentBS_EP [i];
            hList4 [2] = hVyaEtaReflCentBS_EP [i];
            hList4 [3] = hVxYaEtaReflCentBS_EP [i];
            hList4 [4] = hVyXaEtaReflCentBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VaEtaCent_EP

            hList1 [0] = hVbEtaCent_EP [i];
            hList1 [1] = hVxbEtaCent_EP [i];
            hList1 [2] = hVybEtaCent_EP [i];
            hList1 [3] = hVxYbEtaCent_EP [i];
            hList1 [4] = hVyXbEtaCent_EP [i];

            hList2 [0] = hVbEtaCentBS_EP [i];
            hList2 [1] = hVxbEtaCentBS_EP [i];
            hList2 [2] = hVybEtaCentBS_EP [i];
            hList2 [3] = hVxYbEtaCentBS_EP [i];
            hList2 [4] = hVyXbEtaCentBS_EP [i];

            hList3 [0] = hVbEtaReflCent_EP [i];
            hList3 [1] = hVxbEtaReflCent_EP [i];
            hList3 [2] = hVybEtaReflCent_EP [i];
            hList3 [3] = hVxYbEtaReflCent_EP [i];
            hList3 [4] = hVyXbEtaReflCent_EP [i];

            hList4 [0] = hVbEtaReflCentBS_EP [i];
            hList4 [1] = hVxbEtaReflCentBS_EP [i];
            hList4 [2] = hVybEtaReflCentBS_EP [i];
            hList4 [3] = hVxYbEtaReflCentBS_EP [i];
            hList4 [4] = hVyXbEtaReflCentBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VbEtaCent_EP

            hList1 [0] = hVcEtaCent_EP [i];
            hList1 [1] = hVxcEtaCent_EP [i];
            hList1 [2] = hVycEtaCent_EP [i];
            hList1 [3] = hVxYcEtaCent_EP [i];
            hList1 [4] = hVyXcEtaCent_EP [i];

            hList2 [0] = hVcEtaCentBS_EP [i];
            hList2 [1] = hVxcEtaCentBS_EP [i];
            hList2 [2] = hVycEtaCentBS_EP [i];
            hList2 [3] = hVxYcEtaCentBS_EP [i];
            hList2 [4] = hVyXcEtaCentBS_EP [i];

            hList3 [0] = hVcEtaReflCent_EP [i];
            hList3 [1] = hVxcEtaReflCent_EP [i];
            hList3 [2] = hVycEtaReflCent_EP [i];
            hList3 [3] = hVxYcEtaReflCent_EP [i];
            hList3 [4] = hVyXcEtaReflCent_EP [i];

            hList4 [0] = hVcEtaReflCentBS_EP [i];
            hList4 [1] = hVxcEtaReflCentBS_EP [i];
            hList4 [2] = hVycEtaReflCentBS_EP [i];
            hList4 [3] = hVxYcEtaReflCentBS_EP [i];
            hList4 [4] = hVyXcEtaReflCentBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VcEtaCent_EP

            hList1 [0] = hVxEtaMult_SP [i];
            hList1 [1] = hVxaEtaMult_SP [i];
            hList1 [2] = hVxbEtaMult_SP [i];
            hList1 [3] = hVxcEtaMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxEtaMultBS_SP [i];
            hList2 [1] = hVxaEtaMultBS_SP [i];
            hList2 [2] = hVxbEtaMultBS_SP [i];
            hList2 [3] = hVxcEtaMultBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVxEtaReflMult_SP [i];
            hList3 [1] = hVxaEtaReflMult_SP [i];
            hList3 [2] = hVxbEtaReflMult_SP [i];
            hList3 [3] = hVxcEtaReflMult_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVxEtaReflMultBS_SP [i];
            hList4 [1] = hVxaEtaReflMultBS_SP [i];
            hList4 [2] = hVxbEtaReflMultBS_SP [i];
            hList4 [3] = hVxcEtaReflMultBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VxEtaMult_SP

            hList1 [0] = hVyEtaMult_SP [i];
            hList1 [1] = hVyaEtaMult_SP [i];
            hList1 [2] = hVybEtaMult_SP [i];
            hList1 [3] = hVycEtaMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyEtaMultBS_SP [i];
            hList2 [1] = hVyaEtaMultBS_SP [i];
            hList2 [2] = hVybEtaMultBS_SP [i];
            hList2 [3] = hVycEtaMultBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVyEtaReflMult_SP [i];
            hList3 [1] = hVyaEtaReflMult_SP [i];
            hList3 [2] = hVybEtaReflMult_SP [i];
            hList3 [3] = hVycEtaReflMult_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVyEtaReflMultBS_SP [i];
            hList4 [1] = hVyaEtaReflMultBS_SP [i];
            hList4 [2] = hVybEtaReflMultBS_SP [i];
            hList4 [3] = hVycEtaReflMultBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VyEtaMult_SP

            hList1 [0] = hVEtaMult_SP [i];
            hList1 [1] = hVaEtaMult_SP [i];
            hList1 [2] = hVbEtaMult_SP [i];
            hList1 [3] = hVcEtaMult_SP [i];
            hList1 [4] = 0;

            hList2 [0] = hVEtaMultBS_SP [i];
            hList2 [1] = hVaEtaMultBS_SP [i];
            hList2 [2] = hVbEtaMultBS_SP [i];
            hList2 [3] = hVcEtaMultBS_SP [i];
            hList2 [4] = 0;

            hList3 [0] = hVEtaReflMult_SP [i];
            hList3 [1] = hVaEtaReflMult_SP [i];
            hList3 [2] = hVbEtaReflMult_SP [i];
            hList3 [3] = hVcEtaReflMult_SP [i];
            hList3 [4] = 0;

            hList4 [0] = hVEtaReflMultBS_SP [i];
            hList4 [1] = hVaEtaReflMultBS_SP [i];
            hList4 [2] = hVbEtaReflMultBS_SP [i];
            hList4 [3] = hVcEtaReflMultBS_SP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VEtaMult_SP

            hList1 [0] = hVaEtaMult_SP [i];
            hList1 [1] = hVxaEtaMult_SP [i];
            hList1 [2] = hVyaEtaMult_SP [i];
            hList1 [3] = hVxYaEtaMult_SP [i];
            hList1 [4] = hVyXaEtaMult_SP [i];

            hList2 [0] = hVaEtaMultBS_SP [i];
            hList2 [1] = hVxaEtaMultBS_SP [i];
            hList2 [2] = hVyaEtaMultBS_SP [i];
            hList2 [3] = hVxYaEtaMultBS_SP [i];
            hList2 [4] = hVyXaEtaMultBS_SP [i];

            hList3 [0] = hVaEtaReflMult_SP [i];
            hList3 [1] = hVxaEtaReflMult_SP [i];
            hList3 [2] = hVyaEtaReflMult_SP [i];
            hList3 [3] = hVxYaEtaReflMult_SP [i];
            hList3 [4] = hVyXaEtaReflMult_SP [i];

            hList4 [0] = hVaEtaReflMultBS_SP [i];
            hList4 [1] = hVxaEtaReflMultBS_SP [i];
            hList4 [2] = hVyaEtaReflMultBS_SP [i];
            hList4 [3] = hVxYaEtaReflMultBS_SP [i];
            hList4 [4] = hVyXaEtaReflMultBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VaEtaMult_SP

            hList1 [0] = hVbEtaMult_SP [i];
            hList1 [1] = hVxbEtaMult_SP [i];
            hList1 [2] = hVybEtaMult_SP [i];
            hList1 [3] = hVxYbEtaMult_SP [i];
            hList1 [4] = hVyXbEtaMult_SP [i];

            hList2 [0] = hVbEtaMultBS_SP [i];
            hList2 [1] = hVxbEtaMultBS_SP [i];
            hList2 [2] = hVybEtaMultBS_SP [i];
            hList2 [3] = hVxYbEtaMultBS_SP [i];
            hList2 [4] = hVyXbEtaMultBS_SP [i];

            hList3 [0] = hVbEtaReflMult_SP [i];
            hList3 [1] = hVxbEtaReflMult_SP [i];
            hList3 [2] = hVybEtaReflMult_SP [i];
            hList3 [3] = hVxYbEtaReflMult_SP [i];
            hList3 [4] = hVyXbEtaReflMult_SP [i];

            hList4 [0] = hVbEtaReflMultBS_SP [i];
            hList4 [1] = hVxbEtaReflMultBS_SP [i];
            hList4 [2] = hVybEtaReflMultBS_SP [i];
            hList4 [3] = hVxYbEtaReflMultBS_SP [i];
            hList4 [4] = hVyXbEtaReflMultBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VbEtaMult_SP

            hList1 [0] = hVcEtaMult_SP [i];
            hList1 [1] = hVxcEtaMult_SP [i];
            hList1 [2] = hVycEtaMult_SP [i];
            hList1 [3] = hVxYcEtaMult_SP [i];
            hList1 [4] = hVyXcEtaMult_SP [i];

            hList2 [0] = hVcEtaMultBS_SP [i];
            hList2 [1] = hVxcEtaMultBS_SP [i];
            hList2 [2] = hVycEtaMultBS_SP [i];
            hList2 [3] = hVxYcEtaMultBS_SP [i];
            hList2 [4] = hVyXcEtaMultBS_SP [i];

            hList3 [0] = hVcEtaReflMult_SP [i];
            hList3 [1] = hVxcEtaReflMult_SP [i];
            hList3 [2] = hVycEtaReflMult_SP [i];
            hList3 [3] = hVxYcEtaReflMult_SP [i];
            hList3 [4] = hVyXcEtaReflMult_SP [i];

            hList4 [0] = hVcEtaReflMultBS_SP [i];
            hList4 [1] = hVxcEtaReflMultBS_SP [i];
            hList4 [2] = hVycEtaReflMultBS_SP [i];
            hList4 [3] = hVxYcEtaReflMultBS_SP [i];
            hList4 [4] = hVyXcEtaReflMultBS_SP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VcEtaMult_SP

            hList1 [0] = hVxEtaMult_EP [i];
            hList1 [1] = hVxaEtaMult_EP [i];
            hList1 [2] = hVxbEtaMult_EP [i];
            hList1 [3] = hVxcEtaMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVxEtaMultBS_EP [i];
            hList2 [1] = hVxaEtaMultBS_EP [i];
            hList2 [2] = hVxbEtaMultBS_EP [i];
            hList2 [3] = hVxcEtaMultBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVxEtaReflMult_EP [i];
            hList3 [1] = hVxaEtaReflMult_EP [i];
            hList3 [2] = hVxbEtaReflMult_EP [i];
            hList3 [3] = hVxcEtaReflMult_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVxEtaReflMultBS_EP [i];
            hList4 [1] = hVxaEtaReflMultBS_EP [i];
            hList4 [2] = hVxbEtaReflMultBS_EP [i];
            hList4 [3] = hVxcEtaReflMultBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VxEtaMult_EP

            hList1 [0] = hVyEtaMult_EP [i];
            hList1 [1] = hVyaEtaMult_EP [i];
            hList1 [2] = hVybEtaMult_EP [i];
            hList1 [3] = hVycEtaMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVyEtaMultBS_EP [i];
            hList2 [1] = hVyaEtaMultBS_EP [i];
            hList2 [2] = hVybEtaMultBS_EP [i];
            hList2 [3] = hVycEtaMultBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVyEtaReflMult_EP [i];
            hList3 [1] = hVyaEtaReflMult_EP [i];
            hList3 [2] = hVybEtaReflMult_EP [i];
            hList3 [3] = hVycEtaReflMult_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVyEtaReflMultBS_EP [i];
            hList4 [1] = hVyaEtaReflMultBS_EP [i];
            hList4 [2] = hVybEtaReflMultBS_EP [i];
            hList4 [3] = hVycEtaReflMultBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VyEtaMult_EP

            hList1 [0] = hVEtaMult_EP [i];
            hList1 [1] = hVaEtaMult_EP [i];
            hList1 [2] = hVbEtaMult_EP [i];
            hList1 [3] = hVcEtaMult_EP [i];
            hList1 [4] = 0;

            hList2 [0] = hVEtaMultBS_EP [i];
            hList2 [1] = hVaEtaMultBS_EP [i];
            hList2 [2] = hVbEtaMultBS_EP [i];
            hList2 [3] = hVcEtaMultBS_EP [i];
            hList2 [4] = 0;

            hList3 [0] = hVEtaReflMult_EP [i];
            hList3 [1] = hVaEtaReflMult_EP [i];
            hList3 [2] = hVbEtaReflMult_EP [i];
            hList3 [3] = hVcEtaReflMult_EP [i];
            hList3 [4] = 0;

            hList4 [0] = hVEtaReflMultBS_EP [i];
            hList4 [1] = hVaEtaReflMultBS_EP [i];
            hList4 [2] = hVbEtaReflMultBS_EP [i];
            hList4 [3] = hVcEtaReflMultBS_EP [i];
            hList4 [4] = 0;

            PlotFlow (hList1, hList2, hList3, hList4); // VEtaMult_EP

            hList1 [0] = hVaEtaMult_EP [i];
            hList1 [1] = hVxaEtaMult_EP [i];
            hList1 [2] = hVyaEtaMult_EP [i];
            hList1 [3] = hVxYaEtaMult_EP [i];
            hList1 [4] = hVyXaEtaMult_EP [i];

            hList2 [0] = hVaEtaMultBS_EP [i];
            hList2 [1] = hVxaEtaMultBS_EP [i];
            hList2 [2] = hVyaEtaMultBS_EP [i];
            hList2 [3] = hVxYaEtaMultBS_EP [i];
            hList2 [4] = hVyXaEtaMultBS_EP [i];

            hList3 [0] = hVaEtaReflMult_EP [i];
            hList3 [1] = hVxaEtaReflMult_EP [i];
            hList3 [2] = hVyaEtaReflMult_EP [i];
            hList3 [3] = hVxYaEtaReflMult_EP [i];
            hList3 [4] = hVyXaEtaReflMult_EP [i];

            hList4 [0] = hVaEtaReflMultBS_EP [i];
            hList4 [1] = hVxaEtaReflMultBS_EP [i];
            hList4 [2] = hVyaEtaReflMultBS_EP [i];
            hList4 [3] = hVxYaEtaReflMultBS_EP [i];
            hList4 [4] = hVyXaEtaReflMultBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VaEtaMult_EP

            hList1 [0] = hVbEtaMult_EP [i];
            hList1 [1] = hVxbEtaMult_EP [i];
            hList1 [2] = hVybEtaMult_EP [i];
            hList1 [3] = hVxYbEtaMult_EP [i];
            hList1 [4] = hVyXbEtaMult_EP [i];

            hList2 [0] = hVbEtaMultBS_EP [i];
            hList2 [1] = hVxbEtaMultBS_EP [i];
            hList2 [2] = hVybEtaMultBS_EP [i];
            hList2 [3] = hVxYbEtaMultBS_EP [i];
            hList2 [4] = hVyXbEtaMultBS_EP [i];

            hList3 [0] = hVbEtaReflMult_EP [i];
            hList3 [1] = hVxbEtaReflMult_EP [i];
            hList3 [2] = hVybEtaReflMult_EP [i];
            hList3 [3] = hVxYbEtaReflMult_EP [i];
            hList3 [4] = hVyXbEtaReflMult_EP [i];

            hList4 [0] = hVbEtaReflMultBS_EP [i];
            hList4 [1] = hVxbEtaReflMultBS_EP [i];
            hList4 [2] = hVybEtaReflMultBS_EP [i];
            hList4 [3] = hVxYbEtaReflMultBS_EP [i];
            hList4 [4] = hVyXbEtaReflMultBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VbEtaMult_EP

            hList1 [0] = hVcEtaMult_EP [i];
            hList1 [1] = hVxcEtaMult_EP [i];
            hList1 [2] = hVycEtaMult_EP [i];
            hList1 [3] = hVxYcEtaMult_EP [i];
            hList1 [4] = hVyXcEtaMult_EP [i];

            hList2 [0] = hVcEtaMultBS_EP [i];
            hList2 [1] = hVxcEtaMultBS_EP [i];
            hList2 [2] = hVycEtaMultBS_EP [i];
            hList2 [3] = hVxYcEtaMultBS_EP [i];
            hList2 [4] = hVyXcEtaMultBS_EP [i];

            hList3 [0] = hVcEtaReflMult_EP [i];
            hList3 [1] = hVxcEtaReflMult_EP [i];
            hList3 [2] = hVycEtaReflMult_EP [i];
            hList3 [3] = hVxYcEtaReflMult_EP [i];
            hList3 [4] = hVyXcEtaReflMult_EP [i];

            hList4 [0] = hVcEtaReflMultBS_EP [i];
            hList4 [1] = hVxcEtaReflMultBS_EP [i];
            hList4 [2] = hVycEtaReflMultBS_EP [i];
            hList4 [3] = hVxYcEtaReflMultBS_EP [i];
            hList4 [4] = hVyXcEtaReflMultBS_EP [i];

            PlotFlow (hList1, hList2, hList3, hList4); // VcEtaMult_EP

        }
        stepDir -> Write ();

    testFile -> cd (); // test
    testTree -> Write (); // test
    testFile -> Close (); // test
	corrFile -> Close ();
	flowFile -> Close ();
}

void CFlowReconstructor::FillReferenceHist () {
    //histFile = new TFile (histFileName + ".root", "update");
	Int_t n, nDiv = 100;
	map <TString, Float_t> variables;
	Float_t h1, h2;
	Float_t pt, eta, v;
	TH2D* h2VnPtEta;
	TProfile2D* p2VnPtEta;
	TProfile *pVnPtRef, *pVnEtaRef;
	TH1D *hVnPtRef, *hVnEtaRef;
	variables ["ptMin_"] = ptMin_;
	variables ["ptMax_"] = ptMax_;
	variables ["etaMin_"] = etaMin_;
	variables ["etaMax_"] = etaMax_;
	TFile *histFile = new TFile (histFileName_ + ".root", "update");
	TDirectory *histDir = histFile -> mkdir (dirName [4]);
	histDir -> cd ();
	h1 = (ptMax_ - ptMin_) / nBinsPt_ / nDiv;
	h2 = (etaMax_ - etaMin_) / nBinsEta_ / nDiv;

	for (Int_t i = 0; i < nHarmonics; i++) {
		n = harmonicsMap [i];
		p2VnPtEta = new TProfile2D (Form ("p2V%iPtEtaRef", n), Form ("Reference V_{%i} vs P_{T} and ", n) + varName_, nBinsPt_, ptMin_, ptMax_, nBinsEta_, etaMin_, etaMax_);
		pVnPtRef = new TProfile (Form ("pV%iPtRef", n), Form ("Reference V_{%i} vs P_{T}", n), nBinsPt_, ptMin_, ptMax_);
		pVnEtaRef = new TProfile (Form ("pV%iEtaRef", n), Form ("Reference V_{%i} vs ", n) + varName_, nBinsEta_, etaMin_, etaMax_);
		pt = ptMin_ + 0.5 * h1;
		for (int j = 1; j <= nBinsPt_ * nDiv; j++, pt += h1) {
			eta = etaMin_ + 0.5 * h2;
			for (int k = 1; k <= nBinsEta_ * nDiv; k++, eta += h2) {
				variables ["pt"] = pt;
				variables ["eta"] = eta;
				v = harmonicFunctions [n - 1] (variables);
				p2VnPtEta -> Fill (pt, eta, v);
				pVnPtRef -> Fill (pt, v);
				pVnEtaRef -> Fill (eta, v);
				//cout << "Pt = " << pt << "\tEta = " << eta << "\tV" << n << " = " << v << endl;
			}
		}

        h2VnPtEta = p2VnPtEta -> ProjectionXY (Form ("h2V%iPtEtaRef", n));
		hVnPtRef = pVnPtRef -> ProjectionX (Form ("hV%iPtRef", n));
		hVnEtaRef = pVnEtaRef -> ProjectionX (Form ("hV%iEtaRef", n));
		for (int j = 1; j <= nBinsPt_; j++) {
			hVnPtRef -> SetBinError (j, 0.0);
            for (int k = 1; k <= nBinsEta_; k++) {
                hVnEtaRef -> SetBinError (j, 0.0);
                h2VnPtEta -> SetBinError (j, k, 0.0);
            }
		}

		hVnPtRef -> SetTitle (Form ("Reference V_{%i} averaged over ", n) + varName_ + " vs P_{T}");
		hVnEtaRef -> SetTitle (Form ("Reference V_{%i} averaged over P_{T} vs ", n) + varName_);
		hVnPtRef -> Write ();
		hVnEtaRef -> Write ();
		h2VnPtEta -> Write ();
	}
	//histFile -> Close ();
	cout << "Reference histograms built!" << endl;
}


void CFlowReconstructor::SetReferenceOption (Int_t harmonic, TString option) {
    refOptions [harmonic] = option;
}


void CFlowReconstructor::Reference (Float_t ptLow, Float_t ptHigh, Float_t etaLow, Float_t etaHigh) {
	Float_t shift1 [5] = {0.4, 0.3, 0.2, 0.1, 0.15}; // SP, x & SP, y
	Float_t shift2 [5] = {0.35, 0.25, 0.15, 0.05, 0.15}; // EP & RP
	Int_t markerStyle [5] = {24, 25, 26, 27, 20};
	Float_t markerSize [5] = {1.0, 0.5, 1.0, 1.0, 0.75};
	Color_t yColor = kGreen;
	Color_t xColor = kBlue;
	Color_t refColor = kRed;
	Color_t RPColor = kRed;
	Color_t EPColor = kViolet;

	TCanvas* c1;
	TDirectory* dir [4];

	Int_t n, size, minPtBin, minEtaBin, maxPtBin, maxEtaBin, nSteps = 3;
	TH2D *h2x, *h2y, *h2RP, *h2EP, *h2Ref;
	TH1D *hRefPt, *hRefEta, *hx, *hy, *hEP, *hRP, *newHist;
	vector <TH1D*> *allStepsPt, *allStepsEta, *allStepsPtDiv, *allStepsEtaDiv;
	TLegend *leg1, *leg2;
	THStack *histStack;
	TH1F* servHist [4];
	TF1* unityPt = new TF1 ("unityPt", "1", ptMin_, ptMax_);
	TF1* unityEta = new TF1 ("unityEta", "1", etaMin_, etaMax_);
	if (uniformSet) nSteps = 4;

	allStepsPt = new vector <TH1D*> [nHarmonics];
	allStepsEta = new vector <TH1D*> [nHarmonics];
	allStepsPtDiv = new vector <TH1D*> [nHarmonics];
	allStepsEtaDiv = new vector <TH1D*> [nHarmonics];
//	c1 = new TCanvas ("c1","c1", 800, 600);
	c1 = new TCanvas ("c1","c1", 640, 480);
	gStyle -> SetLegendBorderSize (0);


	TFile *histFile = new TFile (histFileName_ + "_flow.root", "READ");
	TFile *outputFile = new TFile (histFileName_ + "_ref.root", "RECREATE");

	for (Int_t i = 0; i < nSteps; i++) {
		dir [i] = outputFile -> mkdir (dirName [i]);
		servHist [i] = new TH1F ();
		servHist [i] -> SetMarkerStyle (markerStyle [i]);
		servHist [i] -> SetMarkerSize (markerSize [i]);
		servHist [i] -> SetMarkerColor (kBlack);
	}

	for (Int_t i = 0; i < nHarmonics; i++) {
		n = harmonicsMap [i];
		TString option = refOptions [n];
		cout << "Reference option " << n << ": " << option << endl;
		if (harmonicFunctionSet) {
//			h2Ref = (TH2D*) histFile -> Get (dirName [4] + Form ("/h2V%iPtEtaRef", n));
//			minPtBin = h2Ref -> GetXaxis () -> FindBin (ptLow);
//			maxPtBin = h2Ref -> GetXaxis () -> FindBin (ptHigh);
//			minEtaBin = h2Ref -> GetYaxis () -> FindBin (etaLow);
//			maxEtaBin = h2Ref -> GetYaxis () -> FindBin (etaHigh);
//			h2Ref -> SetFillColor (refColor);
//			h2Ref -> SetMarkerStyle (markerStyle [4]);
//			h2Ref -> SetMarkerSize (markerSize [4]);
//			h2Ref -> SetMarkerColor (refColor);
			hRefPt = (TH1D*) histFile -> Get (dirName [4] + Form ("/hV%iPtRef", n));
			//hRefPt = h2Ref -> ProjectionX ("hRefPt", minEtaBin, maxEtaBin);
			//hRefPt -> Scale (1.0 / (maxEtaBin - minEtaBin));
			hRefPt -> Sumw2 ();
			hRefPt -> SetMarkerStyle (markerStyle [4]);
			hRefPt -> SetMarkerSize (markerSize [4]);
			hRefPt -> SetMarkerColor (refColor);
			hRefPt -> SetLineColor (refColor);
			allStepsPt [i].push_back (hRefPt);
			hRefEta = (TH1D*) histFile -> Get (dirName [4] + Form ("/hV%iEtaRef", n));
			//hRefEta = h2Ref -> ProjectionY ("hRefEta", minPtBin, maxPtBin);
			//hRefEta -> Scale (1.0 / (maxPtBin - minPtBin));
			hRefEta -> Sumw2 ();
			hRefEta -> SetMarkerStyle (markerStyle [4]);
			hRefEta -> SetMarkerSize (markerSize [4]);
			hRefEta -> SetMarkerColor (refColor);
			hRefEta -> SetLineColor (refColor);
			allStepsEta [i].push_back (hRefEta);
		}

		else {
//			h2x = (TH2D*) histFile -> Get (dirName [0] + Form ("/h2V%ixPtEta_SP", n));
//			minPtBin = h2x -> GetXaxis () -> FindBin (ptLow);
//			maxPtBin = h2x -> GetXaxis () -> FindBin (ptHigh);
//			minEtaBin = h2x -> GetYaxis () -> FindBin (etaLow);
//			maxEtaBin = h2x -> GetYaxis () -> FindBin (etaHigh);
		}
		//printf ("%i\t%i\t%i\t%i\n", minPtBin, maxPtBin, minEtaBin, maxEtaBin); // test
		//Loop over correction steps
		for (int j = 0; j < nSteps; j++) {
			dir [j] -> cd ();
//			h2x = (TH2D*) histFile -> Get (dirName [j] + Form ("/h2V%ixPtEta_SP", n));
//			h2y = (TH2D*) histFile -> Get (dirName [j] + Form ("/h2V%iyPtEta_SP", n));
//			h2EP = (TH2D*) histFile -> Get (dirName [j] + Form ("/h2V%iPtEta_EP", n));

			hx = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%ixPt_SP", n));
//			hx = h2x -> ProjectionX ("hx", minEtaBin, maxEtaBin, "e");
//			hx -> Scale (1.0 / (maxEtaBin - minEtaBin));
			hy = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iyPt_SP", n));
//			hy = h2y -> ProjectionX ("hy", minEtaBin, maxEtaBin, "e");
//			hy -> Scale (1.0 / (maxEtaBin - minEtaBin));
			hEP = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iPt_EP", n));
//			hEP = h2EP -> ProjectionX ("hEP", minEtaBin, maxEtaBin, "e");
//			hEP -> Scale (1.0 / (maxEtaBin - minEtaBin));

            if (option == "average") hx -> Add (hx, hy, 0.5, 0.5);

			if (uniformSet) {
				hRP = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iPt_RP", n));
//				h2RP = (TH2D*) histFile -> Get (dirName [j] + Form ("/h2V%iPtEta_RP", n));
//				hRP = h2RP -> ProjectionX ("hRP", minEtaBin, maxEtaBin, "e");
//				hRP -> Scale (1.0 / (maxEtaBin - minEtaBin));
				hRP -> SetMarkerColor (RPColor);
				hRP -> SetLineColor (RPColor);
				hRP -> SetMarkerStyle (markerStyle [j]);
				hRP -> SetMarkerSize (markerSize [j]);
			}

			gStyle -> SetOptStat (0);
			hx -> SetLineColor (xColor);
			hx -> SetMarkerColor (xColor);
			hx -> SetMarkerStyle (markerStyle [j]);
			hx -> SetMarkerSize (markerSize [j]);
			hy -> SetLineColor (yColor);
			hy -> SetMarkerColor (yColor);
			hy -> SetMarkerStyle (markerStyle [j]);
			hy -> SetMarkerSize (markerSize [j]);
			hEP -> SetMarkerColor (EPColor);
			hEP -> SetLineColor (EPColor);
			hEP -> SetMarkerStyle (markerStyle [j]);
			hEP -> SetMarkerSize (markerSize [j]);

			leg1 = new TLegend (0.7, 0.2, 0.85, 0.5);
			leg1 -> SetFillColor (0);
			leg1 -> SetTextSize(0.04);
            if (option == "" || option == "x") leg1 -> AddEntry (hx, methodName [1], "p");
            if (option == "" || option == "y") leg1 -> AddEntry (hy, methodName [2], "p");
			if (option == "average") leg1 -> AddEntry (hx, methodName [5], "p");
			leg1 -> AddEntry (hEP, methodName [3], "p");
			if (uniformSet) leg1 -> AddEntry (hRP, methodName [4], "p");
			if (harmonicFunctionSet) leg1 -> AddEntry (hRefPt, methodName [0], "p");

			newHist = (TH1D*) hx -> Clone("hV%ixPt_");
			HistShift (newHist, shift1 [j]);
			if (option == "" || option == "average" || option == "x") allStepsPt [i].push_back (newHist);
			newHist = (TH1D*) hy -> Clone("hV%iyPt_");
			HistShift (newHist, -shift1 [j]);
			if (option == "" || option == "y") allStepsPt [i].push_back (newHist);
			newHist = (TH1D*) hEP -> Clone("hV%iPt_EP_");
			HistShift (newHist, -shift2 [j]);
			allStepsPt [i].push_back (newHist);
			if (uniformSet) {
				if (j == 0 || j == 3) {
					newHist = (TH1D*) hRP -> Clone("hV%iPt_RP_");
					HistShift (newHist, shift2 [j]);
					allStepsPt [i].push_back (newHist);
				}
				HistShift (hRP, 0.5 * shift2 [4]);
			}
			HistShift (hx, shift1 [4]);
			HistShift (hy, -shift1 [4]);
			HistShift (hEP, -0.5 * shift2 [4]);

			histStack = new THStack("histStack", Form ("V_{%i} averaged over ", n) + varName_ + Form (" #in [%.1f,%.1f] versus P_{T} (", etaLow, etaHigh) + stepName [j] + Form (");P_{T} [GeV/c];V_{%i}", n));
			if (option == "" || option == "average" || option == "x") histStack -> Add (hx);
			if (option == "" || option == "y") histStack -> Add (hy);
			histStack -> Add (hEP);
			if (uniformSet) histStack -> Add (hRP);
			if (harmonicFunctionSet) histStack -> Add (hRefPt);
			histStack -> Draw ("nostack p e1X0");
			leg1 -> Draw ();
			c1 -> Write (Form ("hV%iPt", n));
			delete histStack;

			if (harmonicFunctionSet) {
				HistShift (hx, -shift1 [4]);
				HistShift (hy, shift1 [4]);
				HistShift (hRP, -0.5 * shift2 [4]);
				HistShift (hEP, 0.5 * shift2 [4]);

				hx -> Divide (hRefPt);
				hy -> Divide (hRefPt);
				hEP -> Divide (hRefPt);
				hRP -> Divide (hRefPt);

				leg2 = new TLegend (0.7, 0.2, 0.85, 0.5);
				leg2 -> SetTextSize(0.04);
				if (option == "" || option == "x") leg2 -> AddEntry (hx, methodName [1], "p");
				if (option == "" || option == "y") leg2 -> AddEntry (hy, methodName [2], "p");
                if (option == "average") leg2 -> AddEntry (hx, methodName [5], "p");
				leg2 -> AddEntry (hEP, methodName [3], "p");
				leg2 -> AddEntry (hRP, methodName [4], "p");
				leg2 -> SetFillColor (0);

				newHist = (TH1D*) hx -> Clone("hV%ixPt/Ref_");
				HistShift (newHist, shift1 [j]);
				if (option == "" || option == "average" || option == "x") allStepsPtDiv [i].push_back (newHist);
				newHist = (TH1D*) hy -> Clone("hV%iyPt/Ref_");
				HistShift (newHist, -shift1 [j]);
				if (option == "" || option == "y") allStepsPtDiv [i].push_back (newHist);
				newHist = (TH1D*) hEP -> Clone("hV%iPt_EP/Ref_");
				HistShift (newHist, -shift2 [j]);
				allStepsPtDiv [i].push_back (newHist);
				if (j == 0 || j == 3) {
					newHist = (TH1D*) hRP -> Clone("hV%iPt_RP/Ref_");
					HistShift (newHist, shift2 [j]);
					allStepsPtDiv [i].push_back (newHist);
				}
				HistShift (hx, shift1 [4]);
				HistShift (hy, -shift1 [4]);
				HistShift (hRP, 0.5 * shift2 [4]);
				HistShift (hEP, -0.5 * shift2 [4]);

				histStack = new THStack ("histStack", Form ("V_{%i} to V_{%i}^{ref} ratio versus P_{T} (", n, n) + stepName [j] + ")");
				if (option == "" || option == "average" || option == "x") histStack -> Add (hx);
				if (option == "" || option == "y") histStack -> Add (hy);
				histStack -> Add (hRP);
				histStack -> Add (hEP);
				histStack -> Draw ("nostack p e1X0");
				histStack -> GetXaxis () -> SetTitle ("P_{T} [GeV/c]");
				histStack -> GetYaxis () -> SetTitle (Form ("V_{%i} / V_{%i}^{ref}", n, n));
				leg2 -> Draw ();
				unityPt -> Draw ("same");
				c1 -> Write (Form ("hV%iPt/Ref", n));
				delete histStack;
			}

			hx = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%ixEta_SP", n));
			hy = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iyEta_SP", n));
			hEP = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iEta_EP", n));

//			hx = h2x -> ProjectionY ("hx", minPtBin, maxPtBin, "e");
//			hx -> Scale (1.0 / (maxPtBin - minPtBin));
//			hy = h2y -> ProjectionY ("hy", minPtBin, maxPtBin, "e");
//			hy -> Scale (1.0 / (maxPtBin - minPtBin));
//			hEP = h2EP -> ProjectionY ("hEP", minPtBin, maxPtBin, "e");
//			hEP -> Scale (1.0 / (maxPtBin - minPtBin));

			if (option == "average") hx -> Add (hx, hy, 0.5, 0.5);

			if (uniformSet) {
//				hRP = h2RP -> ProjectionY ("hRP", minPtBin, maxPtBin, "e");
//				hRP -> Scale (1.0 / (maxPtBin - minPtBin));
//				hRP -> SetMarkerColor (RPColor);
//				hRP -> SetLineColor (RPColor);
//				hRP -> SetMarkerStyle (markerStyle [j]);
//				hRP -> SetMarkerSize (markerSize [j]);

				hRP = (TH1D*) histFile -> Get (dirName [j] + Form ("/hV%iEta_RP", n));
				hRP -> SetMarkerColor (RPColor);
				hRP -> SetLineColor (RPColor);
				hRP -> SetMarkerStyle (markerStyle [j]);
				hRP -> SetMarkerSize (markerSize [j]);
			}

			gStyle -> SetOptStat (0);
			hx -> SetLineColor (xColor);
			hx -> SetMarkerColor (xColor);
			hx -> SetMarkerStyle (markerStyle [j]);
			hx -> SetMarkerSize (markerSize [j]);
			hy -> SetLineColor (yColor);
			hy -> SetMarkerColor (yColor);
			hy -> SetMarkerStyle (markerStyle [j]);
			hy -> SetMarkerSize (markerSize [j]);
			hEP -> SetMarkerColor (EPColor);
			hEP -> SetLineColor (EPColor);
			hEP -> SetMarkerStyle (markerStyle [j]);
			hEP -> SetMarkerSize (markerSize [j]);

			newHist = (TH1D*) hx -> Clone("hV%ixEta_");
			HistShift (newHist, shift1 [j]);
			if (option == "" || option == "average" || option == "x") allStepsEta [i].push_back (newHist);
			newHist = (TH1D*) hy -> Clone("hV%iyEta_");
			HistShift (newHist, -shift1 [j]);
			if (option == "" || option == "y") allStepsEta [i].push_back (newHist);
			newHist = (TH1D*) hEP -> Clone("hV%iEta_EP_");
			HistShift (newHist, -shift2 [j]);
			allStepsEta [i].push_back (newHist);
			if (uniformSet) {
				if (j == 0 || j == 3) {
					newHist = (TH1D*) hRP -> Clone("hV%iEta_RP_");
					HistShift (newHist, shift2 [j]);
					allStepsEta [i].push_back (newHist);
				}
				HistShift (hRP, 0.5 * shift2 [4]);
			}
			HistShift (hx, shift1 [4]);
			HistShift (hy, -shift1 [4]);
			HistShift (hEP, -0.5 * shift2 [4]);

			histStack = new THStack("histStack", Form ("V_{%i} averaged over P_{T} #in [%.1f,%.1f] versus ", n, ptLow, ptHigh) + varName_ + " (" + stepName [j] + ");" + varName_ + Form (";V_{%i}", n));
			if (option == "" || option == "average" || option == "x") histStack -> Add (hx);
			if (option == "" || option == "y") histStack -> Add (hy);
			if (uniformSet) histStack -> Add (hRP);
			histStack -> Add (hEP);
			if (harmonicFunctionSet) histStack -> Add (hRefEta);
			histStack -> Draw ("nostack p e1X0");
			leg1 -> Draw ();
			c1 -> Write (Form ("hV%iEta", n));
			delete histStack;
			delete leg1;

			if (harmonicFunctionSet) {
				HistShift (hx, -shift1 [4]);
				HistShift (hy, shift1 [4]);
				HistShift (hRP, -0.5 * shift2 [4]);
				HistShift (hEP, 0.5 * shift2 [4]);

				hx -> Divide (hRefEta);
				hy -> Divide (hRefEta);
				hEP -> Divide (hRefEta);
				hRP -> Divide (hRefEta);

				newHist = (TH1D*) hx -> Clone("hV%ixEta_");
				HistShift (newHist, shift1 [j]);
				if (option == "" || option == "average" || option == "x") allStepsEtaDiv [i].push_back (newHist);
				newHist = (TH1D*) hy -> Clone("hV%iyEta_");
				HistShift (newHist, -shift1 [j]);
				if (option == "" || option == "y") allStepsEtaDiv [i].push_back (newHist);
				newHist = (TH1D*) hEP -> Clone("hV%iEta_EP_");
				HistShift (newHist, -shift2 [j]);
				allStepsEtaDiv [i].push_back (newHist);
				if (j == 0 || j == 3) {
					newHist = (TH1D*) hRP -> Clone("hV%iEta_RP_");
					HistShift (newHist, shift2 [j]);
					allStepsEtaDiv [i].push_back (newHist);
				}
				HistShift (hx, shift1 [4]);
				HistShift (hy, -shift1 [4]);
				HistShift (hRP, 0.5 * shift2 [4]);
				HistShift (hEP, -0.5 * shift2 [4]);

				histStack = new THStack("histStack", Form ("V_{%i} to V_{%i}^{ref} ratio versus " + varName_ + " (", n, n) + stepName [j] + Form (");P_{T} [GeV/c];V_{%i} / V_{%i}^{ref}", n, n));
				if (option == "" || option == "average" || option == "x") histStack -> Add (hx);
				if (option == "" || option == "y") histStack -> Add (hy);
				histStack -> Add (hRP);
				histStack -> Add (hEP);
				histStack -> Draw ("nostack p e1X0");
				leg2 -> Draw ();
				unityEta -> Draw ("same");
				c1 -> Write (Form ("hV%iEta/Ref", n));
				delete histStack;
				delete leg2;
			}

			/*
			histStack = new THStack ("histStack", Form ("V_{%i} versus P_{T} and " + varName_ + " (", n) + stepName [j] + ")");
			h2x -> SetFillColor (xColor);
			h2x -> SetLineColor (xColor);
			h2x -> SetMarkerColor (xColor);
			h2x -> SetMarkerStyle (markerStyle [j]);
			h2y -> SetFillColor (yColor);
			h2y -> SetLineColor (yColor);
			h2y -> SetMarkerColor (yColor);
			h2y -> SetMarkerStyle (markerStyle [j]);

			histStack -> Add (h2x);
			histStack -> Add (h2y);
			histStack -> Add (h2Ref);
			histStack -> Draw ("p e1X0");
			histStack -> GetXaxis () -> SetTitle ("P_{T}");
			histStack -> GetYaxis () -> SetTitle (varName_);

			leg1 = new TLegend (0.7, 0.2, 0.85, 0.4);
			leg1 -> SetTextSize (0.035);
			leg1 -> SetFillColor (0);
			leg1 -> AddEntry (h2x, Form ("V_{%i}^{x}", n), "pe");
			leg1 -> AddEntry (h2y, Form ("V_{%i}^{y}", n), "pe");
			leg1 -> AddEntry (h2Ref, Form ("V_{%i}^{ref}", n), "pe");
			leg1 -> Draw ();
			c1 -> Write (Form ("h2V%ixyPtEta_p", n));
			delete histStack;
			delete leg1;

			histStack = new THStack ("histStack", Form ("V_{%i} versus P_{T} and " + varName_ + " (", n) + stepName [j] + ")");
			histStack -> Add (h2x);
			histStack -> Add (h2y);
			histStack -> Add (h2Ref);
			histStack -> Draw ();
			histStack -> GetXaxis () -> SetTitle ("P_{T}");
			histStack -> GetYaxis () -> SetTitle (varName_);
			histStack -> Draw ();

			leg1 = new TLegend (0.7, 0.2, 0.85, 0.4);
			leg1 -> SetTextSize (0.035);
			leg1 -> SetFillColor (0);
			leg1 -> AddEntry (h2x, Form ("V_{%i}^{x}", n), "f");
			leg1 -> AddEntry (h2y, Form ("V_{%i}^{y}", n), "f");
			leg1 -> AddEntry (h2Ref, Form ("V_{%i}^{ref}", n), "f");
			leg1 -> Draw ();

			c1 -> Write (Form ("h2V%ixyPtEta_lego", n));
			delete histStack;
			delete leg1;

			h2x -> Divide (h2Ref);
			h2x -> SetTitle (Form ("V_{%i}^{x} to V_{%i}^{ref} ratio versus P_{T} and " + varName_ + " (", n, n) + stepName [j] + ")");
			h2x -> GetYaxis () -> SetTitle (varName_);
			h2x -> GetXaxis () -> SetTitle ("P_{T}");
			//h2x -> Fit (func);
			h2x -> Draw ("COLZ");
			//c1 -> Update();
			//stats_ = (TPaveStats*) h2x -> GetListOfFunctions () -> FindObject ("stats");
			//stats_ -> SetOptFit (1);
			c1 -> Write (Form ("h2V%ixPtEta/Ref", n));

			h2y -> Divide (h2Ref);
			h2y -> Draw ("COLZ");
			c1 -> Write (Form ("h2V%iyPtEta/Ref", n));
			*/
		}

		outputFile -> cd ();

		leg1 = new TLegend (0.6, 0.2, 0.8, 0.6);
		leg1 -> SetFillColor (0);
		leg1 -> SetTextSize (0.04);
        if (option == "" || option == "x") leg1 -> AddEntry (hx, methodName [1], "p");
        if (option == "" || option == "y") leg1 -> AddEntry (hy, methodName [2], "p");
        if (option == "average") leg1 -> AddEntry (hx, methodName [5], "p");
		leg1 -> AddEntry (hEP, methodName [3], "l");
		if (uniformSet) leg1 -> AddEntry (hRP, methodName [4], "l");
		if (harmonicFunctionSet) leg1 -> AddEntry (hRefPt, methodName [0], "p");
		for (int j = 0; j < nSteps; j++) {
			leg1 -> AddEntry (servHist [j], stepName [j], "p");
		}

		histStack = new THStack ("histStack", Form ("V_{%i} averaged over ", n) + varName_ + Form (" #in [%.1f,%.1f] versus P_{T} for all correction steps;P_{T} [GeV/c];V_{%i}", etaLow, etaHigh, n));
		size = allStepsPt [i].size ();
		for (int j = 0; j < size; j++) {
			histStack -> Add (allStepsPt [i][j]);
		}
		histStack -> Draw ("nostack p e1X0");
		leg1 -> Draw ();
		c1 -> Write (Form ("allStepsV%iPt", n));
		delete histStack;

		histStack = new THStack ("histStack", Form ("V_{%i} averaged over P_{T} #in [%.1f,%.1f] versus ", n, ptLow, ptHigh) + varName_ + " for all correction steps;" + varName_ + Form (";V_{%i}", n));
		size = allStepsEta [i].size ();
		for (int j = 0; j < size; j++) {
			histStack -> Add (allStepsEta [i][j]);
		}
		histStack -> Draw ("nostack p e1X0");
		leg1 -> Draw ();
		c1 -> Write (Form ("allStepsV%iEta", n));
		delete histStack;
		delete leg1;

		if (harmonicFunctionSet) {
			leg2 = new TLegend (0.6, 0.2, 0.8, 0.6);
			leg2 -> SetFillColor (0);
			leg2 -> SetTextSize (0.04);
			if (option == "" || option == "x") leg2 -> AddEntry (hx, methodName [1], "l");
			if (option == "" || option == "y") leg2 -> AddEntry (hy, methodName [2], "l");
            if (option == "average") leg2 -> AddEntry (hx, methodName [5], "p");
			leg2 -> AddEntry (hEP, methodName [3], "l");
			leg2 -> AddEntry (hRP, methodName [4], "l");
			for (int j = 0; j < nSteps; j++) {
				leg2 -> AddEntry (servHist [j], stepName [j], "p");
			}

			histStack = new THStack ("histStack", Form ("V_{%i} to V_{%i}^{ref} ratio versus P_{T} for all correction steps;P_{T} [GeV/c];V_{%i} / V_{%i}^{ref}", n, n, n, n));
			size = allStepsPtDiv [i].size ();
			for (int j = 0; j < size; j++) {
				histStack -> Add (allStepsPtDiv [i][j]);
			}
			histStack -> Draw ("nostack p e1X0");
			leg2 -> Draw ();
			unityPt -> Draw ("same");
			c1 -> Write (Form ("allStepsV%iPt/Ref", n));
			delete histStack;

			histStack = new THStack ("histStack", Form ("V_{%i} to V_{%i}^{ref} ratio versus ", n, n) + varName_ + " for all correction steps;" + varName_ + Form (";V_{%i} / V_{%i}^{ref}", n, n));
			size = allStepsEtaDiv [i].size ();
			for (int j = 0; j < size; j++) {
				histStack -> Add (allStepsEtaDiv [i][j]);
			}
			histStack -> Draw ("nostack p e1X0");
			leg2 -> Draw ();
			unityEta -> Draw ("same");
			c1 -> Write (Form ("allStepsV%iEta/Ref", n));
			delete histStack;
			delete leg2;
		}
	}
	gStyle -> SetOptStat (1);
	histFile -> Close ();
	outputFile -> Close ();
	cout << "Reconstructed and reference values compared!" << endl;
}

#endif // CFLOWRECONSTRUCTOR_CXX
