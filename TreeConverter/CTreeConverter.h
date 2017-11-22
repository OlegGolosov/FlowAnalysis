#ifndef CTREECONVERTER_H
#define CTREECONVERTER_H

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "../MyDataTree/CEvent.h"
#include "../defines.h"
#include "../ManualFunctions.h"

class CTreeConverter {

public:
	CTreeConverter ();
	CTreeConverter (TString inputFileName, TString outputFileName, TString inputTreeName);
	void SetInputFileName (TString inputFileName = "../Source.root");
	void SetOutputFileName (TString outputFileName = "../Converted.root");
	void SetInputTreeName (TString inputTreeName);
	void SetMhRange (Int_t mhMin, Int_t mhMax);
	void SetNbinsMh (Int_t nBinsMh);
	void SetCentRange (Float_t centMin, Float_t centMax);
	void SetNbinsCent (Int_t nBinsCent);
	void SetdEdxSource (Int_t TPCid = 3);
	void SetCentralityMethod (Int_t centMethod = 1); // 1 - multiplicity, 2 - Eveto
    void SetSNN (Float_t sNN);
	Bool_t ConvertTree ();

private:
	TFile* inputFile_;
	TFile* outputFile_;
	TTree* inputTree_;
	TTree* outputTree_;
	TString inputFileName_;
	TString outputFileName_;
	TString inputTreeName_;
	Bool_t initFlag_;
	Int_t nBinsMh_, nBinsCent_;
	Int_t mhMin_, mhMax_;
	Float_t centMin_, centMax_;
	Int_t dEdxSource_;
    Int_t centMethod_;
    Float_t sNN_;


	static const Int_t mhMax = 5000;

	// Declaration of leaf types
   Float_t         fVertex_fPchi2;
   Float_t         fVeto_fAdcHadron[4];
   Int_t           fNRun;
   Int_t           fNEvent;
   Int_t           fTriggerMask;
   Int_t           fDate;
   Int_t           fTime;
   Float_t         fEveto;
   Float_t         fVertexX;
   Float_t         fVertexY;
   Float_t         fVertexZ;
   Int_t           fWfaNbeam;
   Int_t           fWfaNinter;
   Int_t           fWfaBeamTime;
   Int_t           fWfaInterTime;
   Int_t           fNPrimaryParticles;
   Char_t          fPrimaryParticles_fIdDet[mhMax];   //[fNPrimaryParticles]
   Char_t          fPrimaryParticles_fCharge[mhMax];   //[fNPrimaryParticles]
   UChar_t         fPrimaryParticles_fNPoint[mhMax][4];   //[fNPrimaryParticles]
   UChar_t         fPrimaryParticles_fNFitPoint[mhMax][4];   //[fNPrimaryParticles]
   UChar_t         fPrimaryParticles_fNDedxPoint[mhMax][4];   //[fNPrimaryParticles]
   UChar_t         fPrimaryParticles_fNMaxPoint[mhMax][4];   //[fNPrimaryParticles]
   UShort_t        fPrimaryParticles_fTmeanCharge[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fPz[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fEta[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fPhi[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fPt[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fSigPx[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fSigPy[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fBx[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fBy[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fPchi2[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fXFirst[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fYFirst[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fZFirst[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fXLast[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fYLast[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fZLast[mhMax][4];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fLabel [mhMax];   //[fNPrimaryParticles]
   Int_t           fPrimaryParticles_fTofIflag[mhMax];   //[fNPrimaryParticles]
   Int_t           fPrimaryParticles_fTofId[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofX[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofY[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofPathl[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofCharge[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofMass2[mhMax];   //[fNPrimaryParticles]
   Float_t         fPrimaryParticles_fTofSigMass2[mhMax];   //[fNPrimaryParticles]
   Float_t         fRing_fADChadron[240];

   // List of branches
   TBranch        *b_fVertex_fPchi2;   //!
   TBranch        *b_Veto_fAdcHadron;   //!
   TBranch        *b_fNRun;   //!
   TBranch        *b_fNEvent;   //!
   TBranch        *b_fTriggerMask;   //!
   TBranch        *b_fDate;   //!
   TBranch        *b_fTime;   //!
   TBranch        *b_fEveto;   //!
   TBranch        *b_fVertexX;   //!
   TBranch        *b_fVertexY;   //!
   TBranch        *b_fVertexZ;   //!
   TBranch        *b_fWfaNbeam;   //!
   TBranch        *b_fWfaNinter;   //!
   TBranch        *b_fWfaBeamTime;   //!
   TBranch        *b_fWfaInterTime;   //!
   TBranch        *b_fNPrimaryParticles;   //!
   TBranch        *b_fPrimaryParticles_fIdDet;   //!
   TBranch        *b_fPrimaryParticles_fCharge;   //!
   TBranch        *b_fPrimaryParticles_fNPoint;   //!
   TBranch        *b_fPrimaryParticles_fNFitPoint;   //!
   TBranch        *b_fPrimaryParticles_fNDedxPoint;   //!
   TBranch        *b_fPrimaryParticles_fNMaxPoint;   //!
   TBranch        *b_fPrimaryParticles_fTmeanCharge;   //!
   TBranch        *b_fPrimaryParticles_fPz;   //!
   TBranch        *b_fPrimaryParticles_fEta;   //!
   TBranch        *b_fPrimaryParticles_fPhi;   //!
   TBranch        *b_fPrimaryParticles_fPt;   //!
   TBranch        *b_fPrimaryParticles_fSigPx;   //!
   TBranch        *b_fPrimaryParticles_fSigPy;   //!
   TBranch        *b_fPrimaryParticles_fBx;   //!
   TBranch        *b_fPrimaryParticles_fBy;   //!
   TBranch        *b_fPrimaryParticles_fPchi2;   //!
   TBranch        *b_fPrimaryParticles_fXFirst;   //!
   TBranch        *b_fPrimaryParticles_fYFirst;   //!
   TBranch        *b_fPrimaryParticles_fZFirst;   //!
   TBranch        *b_fPrimaryParticles_fXLast;   //!
   TBranch        *b_fPrimaryParticles_fYLast;   //!
   TBranch        *b_fPrimaryParticles_fZLast;   //!
   TBranch        *b_fPrimaryParticles_fLabel ;   //!
   TBranch        *b_fPrimaryParticles_fTofIflag;   //!
   TBranch        *b_fPrimaryParticles_fTofId;   //!
   TBranch        *b_fPrimaryParticles_fTofX;   //!
   TBranch        *b_fPrimaryParticles_fTofY;   //!
   TBranch        *b_fPrimaryParticles_fTofPathl;   //!
   TBranch        *b_fPrimaryParticles_fTofCharge;   //!
   TBranch        *b_fPrimaryParticles_fTofMass2;   //!
   TBranch        *b_fPrimaryParticles_fTofSigMass2;   //!
   TBranch        *b_fRing_fADChadron;   //!

	Bool_t SetInputFile ();
	Bool_t SetOutputFile ();
	Bool_t SetInputTree ();
	Bool_t Init ();
	Bool_t CheckEventCuts ();
    Float_t GetCentralityClass (Int_t mh);
	Float_t GetCentralityClass (Float_t Eveto);
	Bool_t CheckTrackCuts (Int_t itrack);
	Int_t GetTrackPid (Int_t itrack);
	Float_t GetRapidity (Float_t pt, Float_t eta, Int_t pid);
	void Finish ();

};

#endif // CTREECONVERTER_H
