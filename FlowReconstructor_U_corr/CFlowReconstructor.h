#ifndef CFLOWRECONSTRUCTOR_H
#define CFLOWRECONSTRUCTOR_H

#include <map>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "../QnCorrections/QnCorrectionsManager.h"
#include "defines.h"

using namespace std;

static const Int_t MAXNSAMPLES = 1000;

enum Detectors {
	kDetector1 = 0,
	kDetector1A,
	kDetector1B,
	kDetector1C,
	kNDetectors
};

enum Variables {
    kNrun = 0,
	kMult,
	kCent,
	kPt,
	kEta,
	kNVars
};

enum ResolutionMethods {
    kThreeSubevents = 0,
    kRandomSubevent,
    kNResMethods
};

enum SamplingMethods {
    kNoSampling = 0,
    kSubsampling,
    kBootStrapping,
    kNSamplingMethods
};

const TString VarNames [kNVars] = {"Run number", "Multiplicity", "Centrality", "Pt", "Eta"};
const TString DetectorNames [kNDetectors] = {"Detector one", "Detector oneA", "Detector oneB", "Detector oneC"};
const TString dirName [5] = {"Not_Corrected", "Recentered", "Diagonalized", "Uniform_Acceptance", "Analytic"};
const TString stepName [4] = {"not corrected", "recentered", "diagonalized", "uniform"};
const TString methodName [6] = {"analytic", "scalar product, x", "scalar product, y", "event plane", "reaction plane", "scalar product"};

class CFlowReconstructor {
public:
	CFlowReconstructor ();
	//virtual ~CFlowReconstructor();
	void SetUniformInputFileName (TString name);
	void SetNonUniformInputFileName (TString name);
	void SetComment (TString comment);
	void SetHistFileName (TString histFileName);

	bool AddHarmonic (Int_t n);
	void UseAutoHistRanges (Bool_t useAutoHistRanges = 1);
	void UseZeroSubevents (Bool_t useZeroSubevents = 1);
	void SetSamplingMethod (Int_t samplingMethod);
    void PropagateResolutionSign (Bool_t propagateResolutionSign = 1);
	void SetNrunRange (Int_t nRunMin, Int_t nRunMax);
	void SetMhRange (Int_t mhMin, Int_t mhMax);
	void SetNbinsMh (Int_t nBinsMh);
	void SetCentRange (Float_t centMin, Float_t centMax);
	void SetNbinsCent (Int_t nBinsCent);
	void SetPtRange (Float_t ptMin, Float_t ptMax);
	bool SetPtAveragingRange (Int_t harmonic, Float_t ptMin, Float_t ptMax);
	void SetNbinsPt (Int_t nBinsPt);
	void SetEtaRange (Float_t etaMin, Float_t etaMax);
	bool SetEtaAveragingRange (Int_t harmonic, Float_t etaMin, Float_t etaMax);
	void SetNbinsEta (Int_t nBinsEta);
	void SetNbinsEtaRefl (Int_t nBinsEtaRefl);
	void AnalyzeTree ();
	void GetCorrelations ();
	void GetFlow ();
	void FillReferenceHist ();
	void Reference (Float_t ptLow, Float_t ptHigh, Float_t etaLow, Float_t etaHigh);
	void SetHarmonicFunction (Int_t n, floatFunction func);
	bool SetEtaSubeventsLimits (Int_t harmonic, Float_t lim1, Float_t lim2, Float_t lim3, Float_t lim4, Float_t lim5 = 0.0, Float_t lim6 = 0.0);
	bool SetPtSubeventsLimits (Int_t harmonic, Float_t lim1, Float_t lim2, Float_t lim3, Float_t lim4, Float_t lim5 = 0.0, Float_t lim6 = 0.0);
	bool AddResolutionParticle (Int_t harmonic, Int_t subevent, Int_t pid);
	void SetResolutionCharge (Int_t charge);
    void AddFlowParticle (Int_t pid);
    bool SetResolutionSigns (Int_t harmonic, Int_t signA, Int_t signB, Int_t signC);
    void SetReferenceOption (Int_t harmonic, TString option = "");
    void SetVariable (TString var = "eta");
    void SetNbinsBS (Int_t nBinsBS);
    void SetMhRangeForFlow (Int_t mhLow, Int_t mhHigh);
    void SetCentRangeForFlow (Float_t centLow, Float_t centHigh);
    void ExcludeRun (Int_t nRun);
    void SetResolutionMethod (Int_t resMethod = 0);
    void SetNsteps (Int_t nSteps = 3);

private:
	Int_t nHarmonics;
	Int_t *harmonicsMap;
	Int_t nSteps_;
	Int_t nRuns_, nRunMin_, nRunMax_;
	Int_t nBinsMh_, nBinsCent_, nBinsPt_, nBinsEta_, nBinsEtaRefl_;
	Float_t ptMin_, ptMax_, etaMin_, etaMax_, centMin_, centMax_, centLow_, centHigh_;
	Int_t mhMin_, mhMax_, mhLow_, mhHigh_;
    Int_t mhLowerBin_, mhHigherBin_, centLowerBin_, centHigherBin_; // flow is analyzed in these bin ranges
	Int_t resMethod_;
	Int_t nBinsBS_;
	vector <Float_t*> etaLim_, ptLim_;
	vector <Float_t*> etaAveragingRange_, ptAveragingRange_;
	Int_t resCharge;
	vector <Int_t> flowParticles;
	vector <vector <Int_t>*> resParticles;
	vector <Int_t*> resSign_;
	vector <Int_t> excludedRuns_;
	map <Int_t, TString> refOptions;

	TString uniformInputFileName;
	TString nonUniformInputFileName;
	TString qnInputFileName;
	TString qnOutputFileName;
	TString histFileName_;
	TString comment_;
	TString varName_;
	TFile* outputFile;
	floatFunction harmonicFunctions [5];

	Bool_t useAutoHistRanges_;
	Bool_t useZeroSubevents_;
	Int_t samplingMethod_;
	Bool_t propagateResolutionSign_;
	Bool_t uniformSet;
	Bool_t harmonicFunctionSet;
	Bool_t resChargeSet;
	Bool_t mhRangeForFlowSet_;
	Bool_t centRangeForFlowSet_;
	void GetVariableRanges (TTree *inputTree);
	void SetupQnCorrectionsManager (TFile *qnInputFile, TFile *qnPtInputFile, TFile *qnEtaInputFile, QnCorrectionsManager *QnMan, QnCorrectionsManager *QnManPt, QnCorrectionsManager *QnManEta);
	void GetCorrelationsLoop (Int_t step);
	void GetFlowLoop (Int_t step);
	void FinalizeQnCorrectionsManager (TFile *qnInputFile, TFile *qnOutputFile, QnCorrectionsManager *QnMan);
	void HistShift (TH1* h, Float_t shiftX = 0.05, Float_t shiftY = 0.0);
	void ShiftArray (Float_t *arr, Int_t n, Float_t shift);
	void Sqrt (TH1 *hR);
	void Sqrt (TH2 *h2R);
	void CalculateCorrelationsWithSampling (TProfile2D *p2XaXb, TProfile *pXaXb);
	void CalculateResolutionNoSampling (TProfile *pXaXb, TProfile *pXaXc, TProfile *pXbXc, TH1F *hRa, TH1F *hRb, TH1F *hRc);
	void CalculateResolutionWithSampling (TProfile2D *p2XaXb, TProfile2D *p2XaXc, TProfile2D *p2XbXc, TH2F *h2Ra, TH2F *h2Rb, TH2F *h2Rc);
	void CalculateFlow (TProfile3D *p3xX, TH2F *h2R, TProfile2D *h2V, Int_t lowerBin, Int_t higherBin, Int_t sign);
	void CalculateFlow (TProfile2D *p2xX, TH2 *h2R, TProfile2D *h2V, Int_t sign);
	void TH2toTH1withSampling (TH2 *h2In, TH1 *hOut, TDirectory *dir = 0);
	void ReflectRapidity (TH1F *hVEta, TH1F *hVEtaRefl, Int_t nHarmonic);
	void PlotKinematics (TFile *corrFile, TDirectory *outputDir, Int_t nHarmonic, Int_t step);
	void PlotResolution (TH1 *hList1 [12], TH1 *hList2 [12], TH1 *hList3 [9], TH1 *hList4 [9], Int_t nHist, TDirectory *dir);
	void CombineSubevents (TH2 *pa, TH2 *pb, TH2 *pc, TH2 *p);
	void PlotFlow (TH1 *h1, TH1 *h2, TH1 *h3, TH1 *h4, TH1 *h5 = 0, TH1 *h6 = 0, TH1 *h7 = 0, TH1 *h8 = 0, TH1 *h9 = 0, TH1 *h10 = 0);
	void BuildSampleTree (TTree *inputTree);
	void WritePreviousResults (TDirectory *dir);
};

#endif // CFLOWRECONSTRUCTOR_H
