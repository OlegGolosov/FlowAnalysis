#ifndef CTREECONVERTERHADES_H
#define CTREECONVERTERHADES_H

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "../MyDataTree/CEvent.h"
#include "../config.h"
#include "../ManualFunctions.h"

class CTreeConverterHADES {

public:
	CTreeConverterHADES ();
	CTreeConverterHADES (TString inputFileName, TString outputFileName, TString inputTreeName);
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
	void  GetAverages ();

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
   vector <TH1*> averages_;


	static const Int_t mhMax = 5000;
  static const Short_t nCuts = 8;
  static const Short_t nTrigger = 4;
  static  const Short_t maxNWallHits = 200;
  static const Short_t maxNTracks = 200;

	// Declaration of leaf types
   Int_t           nRpcClust;
   Int_t           nRpcClustCut;
   Int_t           nRpcHits;
   Int_t           nRpcHitsCut;
   Int_t           nTofHits;
   Int_t           nTofHitsCut;
   Int_t           primaryTracks;
   Int_t           selectedTracks;
   Float_t         psiEPa;
   Float_t         psiEPb;
   Int_t           nA;
   Int_t           nB;
   Bool_t          trigInd[4];
   Short_t         runId;
   Short_t         nWallHitsTot;
   Float_t         runTime;
   Bool_t          cuts[8];
   Short_t         wallModuleIndex[200];   //[nWallHitsTot]
   Float_t         wallHitTime[200];   //[nWallHitsTot]
   Float_t         wallHitCharge[200];   //[nWallHitsTot]
   Float_t         wallHitDistance[200];   //[nWallHitsTot]
   Short_t         wallHitRing[200];   //[nWallHitsTot]
   Float_t         wallHitPhi[200];   //[nWallHitsTot]
   Float_t         wallHitX[200];   //[nWallHitsTot]
   Float_t         wallHitY[200];   //[nWallHitsTot]
   Float_t         wallHitZ[200];   //[nWallHitsTot]
   Float_t         wallHitEta[200];   //[nWallHitsTot]
   Bool_t          isWallHitOk[200];   //[nWallHitsTot]
   Float_t         wallChargeTot_mod[400];
   Float_t         wallChargeTot;
   Float_t         wallChargeTot_ring[9];
   Short_t         nTracks;
   Short_t         nProtons;
   Short_t         nPiP;
   Short_t         nPiM;
   Float_t         vX;
   Float_t         vY;
   Float_t         vZ;
   Float_t         vChi2;
   Int_t           time;
   Short_t         pid_[200];   //[nTracks]
   Float_t         p_[200];   //[nTracks]
   Float_t         phi_[200];   //[nTracks]
   Float_t         theta_[200];   //[nTracks]
   Float_t         pt_[200];   //[nTracks]
   Float_t         rapidity[200];   //[nTracks]
   Float_t         eta_[200];   //[nTracks]
   Float_t         metaBeta[200];   //[nTracks]
   Float_t         metaMass[200];   //[nTracks]
   Short_t         charge_[200];   //[nTracks]
   Float_t         mdcdEdx[200];   //[nTracks]
   Float_t         tofdEdx[200];   //[nTracks]
   Float_t         DCAxy[200];   //[nTracks]
   Float_t         DCAz[200];   //[nTracks]
   Short_t         mdcNhits[200];   //[nTracks]
   Short_t         mdcNhitsInner[200];   //[nTracks]
   Short_t         mdcNhitsOuter[200];   //[nTracks]
   Float_t         chi2all[200];   //[nTracks]
   Float_t         chi2inner[200];   //[nTracks]
   Float_t         chi2outer[200];   //[nTracks]
   Float_t         metaQ[200];   //[nTracks]
   Float_t         metaMatchRadius[200];   //[nTracks]
   Float_t         pCorr[200];   //[nTracks]
   Float_t         pt_corr[200];   //[nTracks]
   Float_t         rapidity_corr[200];   //[nTracks]
   Float_t         metaDx[200];   //[nTracks]
   Float_t         metaDy[200];   //[nTracks]
   Float_t         mdcSecId[200];   //[nTracks]

   // List of branches
   TBranch        *b_nRpcClust;   //!
   TBranch        *b_nRpcClustCut;   //!
   TBranch        *b_nRpcHits;   //!
   TBranch        *b_nRpcHitsCut;   //!
   TBranch        *b_nTofHits;   //!
   TBranch        *b_nTofHitsCut;   //!
   TBranch        *b_primaryTracks;   //!
   TBranch        *b_selectedTracks;   //!
   TBranch        *b_psiEPa;   //!
   TBranch        *b_psiEPb;   //!
   TBranch        *b_nA;   //!
   TBranch        *b_nB;   //!
   TBranch        *b_trigInd;   //!
   TBranch        *b_runId;   //!
   TBranch        *b_nWallHitsTot;   //!
   TBranch        *b_runTime;   //!
   TBranch        *b_cuts;   //!
   TBranch        *b_wallModuleIndex;   //!
   TBranch        *b_wallHitTime;   //!
   TBranch        *b_wallHitCharge;   //!
   TBranch        *b_wallHitDistance;   //!
   TBranch        *b_wallHitRing;   //!
   TBranch        *b_wallHitPhi;   //!
   TBranch        *b_wallHitX;   //!
   TBranch        *b_wallHitY;   //!
   TBranch        *b_wallHitZ;   //!
   TBranch        *b_wallHitEta;   //!
   TBranch        *b_isWallHitOk;   //!
   TBranch        *b_wallChargeTot_mod;   //!
   TBranch        *b_wallChargeTot;   //!
   TBranch        *b_wallChargeTot_ring;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nProtons;   //!
   TBranch        *b_nPiP;   //!
   TBranch        *b_nPiM;   //!
   TBranch        *b_vX;   //!
   TBranch        *b_vY;   //!
   TBranch        *b_vZ;   //!
   TBranch        *b_vChi2;   //!
   TBranch        *b_time;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_p;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_rapidity;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_metaBeta;   //!
   TBranch        *b_metaMass;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_mdcdEdx;   //!
   TBranch        *b_tofdEdx;   //!
   TBranch        *b_DCAxy;   //!
   TBranch        *b_DCAz;   //!
   TBranch        *b_mdcNhits;   //!
   TBranch        *b_mdcNhitsInner;   //!
   TBranch        *b_mdcNhitsOuter;   //!
   TBranch        *b_chi2all;   //!
   TBranch        *b_chi2inner;   //!
   TBranch        *b_chi2outer;   //!
   TBranch        *b_metaQ;   //!
   TBranch        *b_metaMatchRadius;   //!
   TBranch        *b_pCorr;   //!
   TBranch        *b_pt_corr;   //!
   TBranch        *b_rapidity_corr;   //!
   TBranch        *b_metaDx;   //!
   TBranch        *b_metaDy;   //!
   TBranch        *b_mdcSecId;   //!

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

#endif // CTREECONVERTERHADES_H
