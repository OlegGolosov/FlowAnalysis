#include <iostream>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <THStack.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>

#include "../FlowReconstructor/CFlowReconstructor.h"

static const int nFiles = 6, nStacks = 24;
//Int_t nBinsBS [8] = {10, 20, 30, 40, 50, 60, 70, 80};
//Int_t nBinsBS [8] = {100, 125, 150, 175, 200};
Int_t nBinsBS [8] = {10, 30, 50, 100, 150, 200};

void Scan (TString location);
void HistShift (TH1F* h, Float_t shiftX, Float_t shiftY = 0.0);
void HistShift (TProfile* h, Float_t shiftX, Float_t shiftY = 0.0);

void HistShift (TH1F* h, Float_t shiftX, Float_t shiftY) {
	Float_t xMin = h -> GetXaxis() -> GetXmin (), xMax = h -> GetXaxis() -> GetXmax ();
	Float_t deltaX = shiftX * h -> GetBinWidth (1);
	TF1 *f = new TF1 ("f", Form ("1 + %f * x", shiftY), xMin, xMax);
	h -> Multiply (f, 1);
	h -> GetXaxis () -> SetLimits (xMin + deltaX, xMax + deltaX);
}

void HistShift (TProfile* h, Float_t shiftX, Float_t shiftY) {
	Float_t xMin = h -> GetXaxis() -> GetXmin (), xMax = h -> GetXaxis() -> GetXmax ();
	Float_t deltaX = shiftX * h -> GetBinWidth (1);
	TF1 *f = new TF1 ("f", Form ("1 + %f * x", shiftY), xMin, xMax);
	h -> Multiply (f, 1);
	h -> GetXaxis () -> SetLimits (xMin + deltaX, xMax + deltaX);
}


void BSscan (Int_t nFile = 0) {
//    TString location = "/mnt/pool/rhic/2/ovgol/NA49_flow/01d/BSscan/";
//    TString nonUniformInputFileName = "/mnt/pool/rhic/2/ovgol/NA49_conv/01d/01d_eta";
    TString location = "/mnt/pool/rhic/2/ovgol/NA49_flow/BSscan_10/";
    TString nonUniformInputFileName = "/mnt/pool/rhic/2/ovgol/NA49_conv/NA49_MTPC_mh";
//    TString location = "/mnt/pool/rhic/2/ovgol/NA49_flow/02c/BSscan/";
//    TString nonUniformInputFileName = "/mnt/pool/rhic/2/ovgol/NA49_conv/02c/02c_eta";
    TString histFileName;

	CFlowReconstructor flowReconstructor;

    if (nFile == 0) {
        Scan (location);
    }

    else {
        histFileName = location + Form ("NA49_y_pion_%iSS", nBinsBS [nFile - 1]);
        flowReconstructor.SetNonUniformInputFileName (nonUniformInputFileName);
        flowReconstructor.SetHistFileName (histFileName);
        flowReconstructor.AddHarmonic (1);
        flowReconstructor.AddHarmonic (2);
        flowReconstructor.UseZeroSubevents ();
        flowReconstructor.UseSubsampling ();
        flowReconstructor.PropagateResolutionSign ();
        flowReconstructor.SetVariable("y");
        flowReconstructor.SetNrunRange (3003, 3166);
        flowReconstructor.ExcludeRun (3141);
        flowReconstructor.ExcludeRun (3010);
        flowReconstructor.ExcludeRun (3030);
        flowReconstructor.ExcludeRun (3134);
        flowReconstructor.ExcludeRun (3159);
        flowReconstructor.ExcludeRun (3014);
        flowReconstructor.ExcludeRun (3161);
        flowReconstructor.SetMhRange (10, 530); // all
        flowReconstructor.SetNbinsMh (10); // all
        flowReconstructor.SetCentRange (0.0, 0.8); // all
        flowReconstructor.SetNbinsCent (8); // all
        flowReconstructor.SetCentRangeForFlow (0.2, 0.4);
        flowReconstructor.SetMhRangeForFlow (10, 530);
        flowReconstructor.SetPtRange (0.0, 2.5);
        flowReconstructor.SetNbinsPt (10);
        flowReconstructor.SetEtaRange (-1.0, 3.0); //symm
        flowReconstructor.SetNbinsEta (8); // symm
        flowReconstructor.SetNbinsEtaRefl (4); // symm
        flowReconstructor.SetPtAveragingRange (1, 0.1, 3.0); // mix pion
        flowReconstructor.SetPtAveragingRange (2, 0.1, 2.0); // mix pion
        flowReconstructor.SetEtaAveragingRange (1, 0.0, 1.0); // mix pion
        flowReconstructor.SetEtaAveragingRange (2, 0.0, 0.5); // mix pion
        flowReconstructor.SetPtSubeventsLimits (1, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0); // mix pion
        flowReconstructor.SetPtSubeventsLimits (2, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0); // mix pion
        flowReconstructor.SetEtaSubeventsLimits (1, 0.5, 3.0, 1.0, 3.0, 0.5, 3.0); // test
        flowReconstructor.SetEtaSubeventsLimits (2, -3.0, 0.0, 0.5, 1.2, 1.2, 3.0); // test
        flowReconstructor.AddFlowParticle (kPionMinus); // mix pion
        flowReconstructor.AddResolutionParticle (1, 1, kProton); // mix pion
        flowReconstructor.AddResolutionParticle (1, 2, kPionMinus); // mix pion
        flowReconstructor.AddResolutionParticle (1, 3, kPionPlus); // mix pion
        flowReconstructor.AddResolutionParticle (2, 1, kPionMinus); // test
        flowReconstructor.AddResolutionParticle (2, 1, kPionPlus); // mix pion
        flowReconstructor.AddResolutionParticle (2, 1, kProton); // test
        flowReconstructor.AddResolutionParticle (2, 2, kPionMinus); // test
        flowReconstructor.AddResolutionParticle (2, 2, kPionPlus); // test
        flowReconstructor.AddResolutionParticle (2, 2, kProton); // mix pion
        flowReconstructor.AddResolutionParticle (2, 3, kPionMinus); // mix pion
        flowReconstructor.AddResolutionParticle (2, 3, kPionPlus); // test
        flowReconstructor.AddResolutionParticle (2, 3, kProton); // test
        flowReconstructor.SetResolutionSigns(1, -1, -1, -1); // mix pion
        flowReconstructor.SetResolutionSigns(2, 1, 1, 1); // mix pion

        flowReconstructor.SetNbinsBS (nBinsBS [nFile - 1]);
        flowReconstructor.GetCorrelations ();
//        flowReconstructor.CalculateFlow (histFileName);
    }
}

void Scan (TString location) {
    static const Int_t nBinsCent = 12, nBinsMult = 10;
    Int_t nBinsR = 100, nBinsCorr = 50;
    Float_t Rmin = -0.25, Rmax = 0.25;
//    Float_t Rmin = -0.008, Rmax = 0.008;
    Float_t corrMin = -0.01, corrMax = 0.01;
    Int_t n, nSteps = 3;
    static const Int_t nHarmonics = 2;
    Int_t harmonicsMap [2] = {1, 2};
    TString	dirName [5] = {"Not_Corrected/", "Recentered/", "Diagonalized/", "Uniform_Acceptance/", "Analytic/"};
    TString	stepName [4] = {"not corrected", "recentered", "diagonalized", "uniform"};
    TString canvasName [4] = {"X_SP_", "Y_SP_", "X_EP_", "Y_EP_"};
    Int_t markerColors [8] = {1, 2, 3, 4, 6, 49, 41, 9};
    Int_t markerStyles [8] = {24, 25, 26, 27, 28, 30, 32, 5};
    TDirectory *outputDir;

    Float_t shifts [nFiles];
   for (Int_t i = 0; i < nFiles; i++) {
        shifts [i] = 0.1 * (i - nFiles / 2);
    }

    static const int nStacks = 24;

    TString outputFileName = "NA49_pion_BSscan.root";
    TFile *inputFile, *inputFile2;
    TFile *outputFile = new TFile (location + outputFileName, "recreate");

    TCanvas *c1 = new TCanvas ("c1", "c1", 800, 600);
    gStyle -> SetLegendBorderSize (0);
    c1 -> Divide (3, 2);

    TLegend *leg = new TLegend (0.7, 0.15, 0.85, 0.5);
    leg -> SetFillColor (0);
    leg -> SetTextSize(0.04);

    TLegend *leg2 = new TLegend (0.7, 0.15, 0.85, 0.5);
    leg2 -> SetFillColor (0);
    leg2 -> SetTextSize(0.04);

    TH1F *fakeHist [nFiles];
   for (Int_t f = 0; f < nFiles; f++) {
       fakeHist [f] = new TH1F ();
       fakeHist [f] -> SetMarkerStyle (markerStyles [f]);
       fakeHist [f] -> SetMarkerColor (markerColors [f]);
       fakeHist [f] -> SetLineColor (markerColors [f]);
        leg -> AddEntry (fakeHist [f], Form ("%i", nBinsBS [f]), "p");
        leg2 -> AddEntry (fakeHist [f], Form ("%i", nBinsBS [f]), "l");
    }

    TProfile *pListCent [nStacks], *pListMult [nStacks];
    TH1F *hListCent [nBinsCent * nStacks], *hListMult [nBinsMult * nStacks];
    TProfile2D *p2ListCent [nStacks], *p2ListMult [nStacks];
    THStack *hsListCent [nStacks], *hsListCent2 [nBinsCent * nStacks], *hsListMult [nStacks], *hsListMult2 [nBinsMult * nStacks];

   for (Int_t step = 0; step < nSteps; step++) {
        outputDir = outputFile -> mkdir (dirName [step]);

       for (Int_t i = 0; i < nHarmonics; i++) {
            n = harmonicsMap [i];

            hsListCent [0] = new THStack (Form ("X%iaX%ibCentSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT ;centrality", n, n));
            hsListCent [1] = new THStack (Form ("X%iaX%icCentSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT ;centrality", n, n));
            hsListCent [2] = new THStack (Form ("X%ibX%icCentSP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT ;centrality", n, n));
            hsListCent [3] = new THStack (Form ("R%ixCentSPa", n, n), Form ("R_{%i,a}^{x,SP} ;centrality", n));
            hsListCent [4] = new THStack (Form ("R%ixCentSPb", n, n), Form ("R_{%i,b}^{x,SP} ;centrality", n));
            hsListCent [5] = new THStack (Form ("R%ixCentSPc", n, n), Form ("R_{%i,c}^{x,SP} ;centrality", n));
            hsListCent [6] = new THStack (Form ("Y%iaY%ibCentSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT ;centrality", n, n));
            hsListCent [7] = new THStack (Form ("Y%iaY%icCentSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT ;centrality", n, n));
            hsListCent [8] = new THStack (Form ("Y%ibY%icCentSP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT ;centrality", n, n));
            hsListCent [9] = new THStack (Form ("R%iyCentSPa", n), Form ("R_{%i,a}^{y,SP} ;centrality", n));
            hsListCent [10] = new THStack (Form ("R%iyCentSPb", n), Form ("R_{%i,b}^{y,SP} ;centrality", n));
            hsListCent [11] = new THStack (Form ("R%iyCentSPc", n), Form ("R_{%i,c}^{y,SP} ;centrality", n));

            hsListCent [12] = new THStack (Form ("X%iaX%ibCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT ;centrality", n, n, n));
            hsListCent [13] = new THStack (Form ("X%iaX%icCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT ;centrality", n, n, n));
            hsListCent [14] = new THStack (Form ("X%ibX%icCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT ;centrality", n, n, n));
            hsListCent [15] = new THStack (Form ("R%ixCentEPa", n, n), Form ("R_{%i,a}^{x,EP} ;centrality", n));
            hsListCent [16] = new THStack (Form ("R%ixCentEPb", n, n), Form ("R_{%i,b}^{x,EP} ;centrality", n));
            hsListCent [17] = new THStack (Form ("R%ixCentEPc", n, n), Form ("R_{%i,c}^{x,EP} ;centrality", n));
            hsListCent [18] = new THStack (Form ("Y%iaY%ibCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT ;centrality", n, n, n, n));
            hsListCent [19] = new THStack (Form ("Y%iaY%icCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT ;centrality", n, n, n, n));
            hsListCent [20] = new THStack (Form ("Y%ibY%icCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT ;centrality", n, n, n, n));
            hsListCent [21] = new THStack (Form ("R%iyCentEPa", n, n), Form ("R_{%i,a}^{y,EP} ;centrality", n));
            hsListCent [22] = new THStack (Form ("R%iyCentEPb", n, n), Form ("R_{%i,b}^{y,EP} ;centrality", n));
            hsListCent [23] = new THStack (Form ("R%iyCentEPc", n, n), Form ("R_{%i,c}^{y,EP} ;centrality", n));

            for (Int_t j = 0; j < nBinsCent; j++) {
                hsListCent2 [j * nStacks + 0] = new THStack (Form ("X%iaX%ibCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 1] = new THStack (Form ("X%iaX%icCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 2] = new THStack (Form ("X%ibX%icCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 3] = new THStack (Form ("R%ixCentSPa_%i", n, j + 1), Form ("R_{%i,a}^{x,SP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 4] = new THStack (Form ("R%ixCentSPb_%i", n, j + 1), Form ("R_{%i,b}^{x,SP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 5] = new THStack (Form ("R%ixCentSPc_%i", n, j + 1), Form ("R_{%i,c}^{x,SP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 6] = new THStack (Form ("Y%iaY%ibCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 7] = new THStack (Form ("Y%iaY%icCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 8] = new THStack (Form ("Y%ibY%icCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT, centrality class %i", n, n, j + 1));
                hsListCent2 [j * nStacks + 9] = new THStack (Form ("R%iyCentSPa_%i", n, j + 1), Form ("R_{%i,a}^{y,SP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 10] = new THStack (Form ("R%iyCentSPb_%i", n, j + 1), Form ("R_{%i,b}^{y,SP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 11] = new THStack (Form ("R%iyCentSPc_%i", n, j + 1), Form ("R_{%i,c}^{y,SP}, centrality class %i", n, j + 1));

                hsListCent2 [j * nStacks + 12] = new THStack (Form ("X%iaX%ibCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT, centrality class %i", n, n, n, j + 1));
                hsListCent2 [j * nStacks + 13] = new THStack (Form ("X%iaX%icCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT, centrality class %i", n, n, n, j + 1));
                hsListCent2 [j * nStacks + 14] = new THStack (Form ("X%ibX%icCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT, centrality class %i", n, n, n, j + 1));
                hsListCent2 [j * nStacks + 15] = new THStack (Form ("R%ixCentEPa_%i", n, j + 1), Form ("R_{%i,a}^{x,EP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 16] = new THStack (Form ("R%ixCentEPb_%i", n, j + 1), Form ("R_{%i,b}^{x,EP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 17] = new THStack (Form ("R%ixCentEPc_%i", n, j + 1), Form ("R_{%i,c}^{x,EP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 18] = new THStack (Form ("Y%iaY%ibCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT, centrality class %i", n, n, n, n, j + 1));
                hsListCent2 [j * nStacks + 19] = new THStack (Form ("Y%iaY%icCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT, centrality class %i", n, n, n, n, j + 1));
                hsListCent2 [j * nStacks + 20] = new THStack (Form ("Y%ibY%icCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT, centrality class %i", n, n, n, n, j + 1));
                hsListCent2 [j * nStacks + 21] = new THStack (Form ("R%iyCentEPa_%i", n, j + 1), Form ("R_{%i,a}^{y,EP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 22] = new THStack (Form ("R%iyCentEPb_%i", n, j + 1), Form ("R_{%i,b}^{y,EP}, centrality class %i", n, j + 1));
                hsListCent2 [j * nStacks + 23] = new THStack (Form ("R%iyCentEPc_%i", n, j + 1), Form ("R_{%i,c}^{y,EP}, centrality class %i", n, j + 1));
            }

            hsListMult [0] = new THStack (Form ("X%iaX%ibMultSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT ;multiplicity", n, n));
            hsListMult [1] = new THStack (Form ("X%iaX%icMultSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT ;multiplicity", n, n));
            hsListMult [2] = new THStack (Form ("X%ibX%icMultSP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT ;multiplicity", n, n));
            hsListMult [3] = new THStack (Form ("R%ixMultSPa", n, n), Form ("R_{%i,a}^{x,SP} ;multiplicity", n));
            hsListMult [4] = new THStack (Form ("R%ixMultSPb", n, n), Form ("R_{%i,b}^{x,SP} ;multiplicity", n));
            hsListMult [5] = new THStack (Form ("R%ixMultSPc", n, n), Form ("R_{%i,c}^{x,SP} ;multiplicity", n));
            hsListMult [6] = new THStack (Form ("Y%iaY%ibMultSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT ;multiplicity", n, n));
            hsListMult [7] = new THStack (Form ("Y%iaY%icMultSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT ;multiplicity", n, n));
            hsListMult [8] = new THStack (Form ("Y%ibY%icMultSP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT ;multiplicity", n, n));
            hsListMult [9] = new THStack (Form ("R%iyMultSPa", n, n), Form ("R_{%i,a}^{y,SP} ;multiplicity", n));
            hsListMult [10] = new THStack (Form ("R%iyMultSPb", n, n), Form ("R_{%i,b}^{y,SP} ;multiplicity", n));
            hsListMult [11] = new THStack (Form ("R%iyMultSPc", n, n), Form ("R_{%i,c}^{y,SP} ;multiplicity", n));

            hsListMult [12] = new THStack (Form ("X%iaX%ibMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT ;multiplicity", n, n, n));
            hsListMult [13] = new THStack (Form ("X%iaX%icMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT ;multiplicity", n, n, n));
            hsListMult [14] = new THStack (Form ("X%ibX%icMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT ;multiplicity", n, n, n));
            hsListMult [15] = new THStack (Form ("R%ixMultEPa", n, n), Form ("R_{%i,a}^{x,EP} ;multiplicity", n));
            hsListMult [16] = new THStack (Form ("R%ixMultEPb", n, n), Form ("R_{%i,b}^{x,EP} ;multiplicity", n));
            hsListMult [17] = new THStack (Form ("R%ixMultEPc", n, n), Form ("R_{%i,c}^{x,EP} ;multiplicity", n));
            hsListMult [18] = new THStack (Form ("Y%iaY%ibMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT ;multiplicity", n, n, n, n));
            hsListMult [19] = new THStack (Form ("Y%iaY%icMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT ;multiplicity", n, n, n, n));
            hsListMult [20] = new THStack (Form ("Y%ibY%icMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT ;multiplicity", n, n, n, n));
            hsListMult [21] = new THStack (Form ("R%iyMultEPa", n, n), Form ("R_{%i,a}^{y,EP} ;multiplicity", n));
            hsListMult [22] = new THStack (Form ("R%iyMultEPb", n, n), Form ("R_{%i,b}^{y,EP} ;multiplicity", n));
            hsListMult [23] = new THStack (Form ("R%iyMultEPc", n, n), Form ("R_{%i,c}^{y,EP} ;multiplicity", n));

            for (Int_t j = 0; j < nBinsMult; j++) {
                hsListMult2 [j * nStacks + 0] = new THStack (Form ("X%iaX%ibMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 1] = new THStack (Form ("X%iaX%icMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 2] = new THStack (Form ("X%ibX%icMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 3] = new THStack (Form ("R%ixMultSPa_%i", n, j + 1), Form ("R_{%i,a}^{x,SP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 4] = new THStack (Form ("R%ixMultSPb_%i", n, j + 1), Form ("R_{%i,b}^{x,SP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 5] = new THStack (Form ("R%ixMultSPc_%i", n, j + 1), Form ("R_{%i,c}^{x,SP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 6] = new THStack (Form ("Y%iaY%ibMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 7] = new THStack (Form ("Y%iaY%icMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 8] = new THStack (Form ("Y%ibY%icMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT, multiplicity bin %i", n, n, j + 1));
                hsListMult2 [j * nStacks + 9] = new THStack (Form ("R%iyMultSPa_%i", n, j + 1), Form ("R_{%i,a}^{y,SP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 10] = new THStack (Form ("R%iyMultSPb_%i", n, j + 1), Form ("R_{%i,b}^{y,SP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 11] = new THStack (Form ("R%iyMultSPc_%i", n, j + 1), Form ("R_{%i,c}^{y,SP}, multiplicity bin %i", n, j + 1));

                hsListMult2 [j * nStacks + 12] = new THStack (Form ("X%iaX%ibMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT, multiplicity bin %i", n, n, n, j + 1));
                hsListMult2 [j * nStacks + 13] = new THStack (Form ("X%iaX%icMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT, multiplicity bin %i", n, n, n, j + 1));
                hsListMult2 [j * nStacks + 14] = new THStack (Form ("X%ibX%icMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT, multiplicity bin %i", n, n, n, j + 1));
                hsListMult2 [j * nStacks + 15] = new THStack (Form ("R%ixMultEPa_%i", n, j + 1), Form ("R_{%i,a}^{x,EP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 16] = new THStack (Form ("R%ixMultEPb_%i", n, j + 1), Form ("R_{%i,b}^{x,EP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 17] = new THStack (Form ("R%ixMultEPc_%i", n, j + 1), Form ("R_{%i,c}^{x,EP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 18] = new THStack (Form ("Y%iaY%ibMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT, multiplicity bin %i", n, n, n, n, j + 1));
                hsListMult2 [j * nStacks + 19] = new THStack (Form ("Y%iaY%icMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT, multiplicity bin %i", n, n, n, n, j + 1));
                hsListMult2 [j * nStacks + 20] = new THStack (Form ("Y%ibY%icMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT, multiplicity bin %i", n, n, n, n, j + 1));
                hsListMult2 [j * nStacks + 21] = new THStack (Form ("R%iyMultEPa_%i", n, j + 1), Form ("R_{%i,a}^{y,EP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 22] = new THStack (Form ("R%iyMultEPb_%i", n, j + 1), Form ("R_{%i,b}^{y,EP}, multiplicity bin %i", n, j + 1));
                hsListMult2 [j * nStacks + 23] = new THStack (Form ("R%iyMultEPc_%i", n, j + 1), Form ("R_{%i,c}^{y,EP}, multiplicity bin %i", n, j + 1));
           }

           for (Int_t f = 0; f < nFiles; f++) {
                if (f == 0) cout << stepName [step] << endl;
                cout << location + Form ("NA49_y_pion_%iSS_flow.root", nBinsBS [f]) << endl;
                inputFile = new TFile (location + Form ("NA49_y_pion_%iSS_flow.root", nBinsBS [f]), "read");
                inputFile2 = new TFile (location + Form ("NA49_y_pion_%iSS_corr.root", nBinsBS [f]), "read");

                pListCent [0] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%ibCentBS", n, n));
                pListCent [1] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%icCentBS", n, n));
                pListCent [2] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%ibX%icCentBS", n, n));
                pListCent [3] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPaCentBS", n));
                pListCent [4] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPbCentBS", n));
                pListCent [5] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPcCentBS", n));
                pListCent [6] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%ibCentBS", n, n));
                pListCent [7] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%icCentBS", n, n));
                pListCent [8] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%ibY%icCentBS", n, n));
                pListCent [9] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPaCentBS", n));
                pListCent [10] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPbCentBS", n));
                pListCent [11] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPcCentBS", n));

                pListCent [12] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsibCentBS", n));
                pListCent [13] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsicCentBS", n));
                pListCent [14] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsib_PsicCentBS", n));
                pListCent [15] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPaCentBS", n));
                pListCent [16] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPbCentBS", n));
                pListCent [17] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPcCentBS", n));
                pListCent [18] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsibCentBS", n, n));
                pListCent [19] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsicCentBS", n, n));
                pListCent [20] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsibSin%iPsicCentBS", n, n));
                pListCent [21] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPaCentBS", n));
                pListCent [22] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPbCentBS", n));
                pListCent [23] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPcCentBS", n));

                p2ListCent [0] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%iaX%ibCent", n, n));
                p2ListCent [1] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%iaX%icCent", n, n));
                p2ListCent [2] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%ibX%icCent", n, n));
                p2ListCent [3] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPaCent", n));
                p2ListCent [4] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPbCent", n));
                p2ListCent [5] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPcCent", n));
                p2ListCent [6] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%iaY%ibCent", n, n));
                p2ListCent [7] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%iaY%icCent", n, n));
                p2ListCent [8] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%ibY%icCent", n, n));
                p2ListCent [9] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPaCent", n));
                p2ListCent [10] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPbCent", n));
                p2ListCent [11] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPcCent", n));

                p2ListCent [12] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsia_PsibCent", n));
                p2ListCent [13] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsia_PsicCent", n));
                p2ListCent [14] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsib_PsicCent", n));
                p2ListCent [15] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPaCent", n));
                p2ListCent [16] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPbCent", n));
                p2ListCent [17] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPcCent", n));
                p2ListCent [18] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsiaSin%iPsibCent", n, n));
                p2ListCent [19] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsiaSin%iPsicCent", n, n));
                p2ListCent [20] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsibSin%iPsicCent", n, n));
                p2ListCent [21] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPaCent", n));
                p2ListCent [22] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPbCent", n));
                p2ListCent [23] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPcCent", n));

                pListMult [0] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%ibMultBS", n, n));
                pListMult [1] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%icMultBS", n, n));
                pListMult [2] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%ibX%icMultBS", n, n));
                pListMult [3] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPaMultBS", n));
                pListMult [4] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPbMultBS", n));
                pListMult [5] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPcMultBS", n));
                pListMult [6] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%ibMultBS", n, n));
                pListMult [7] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%icMultBS", n, n));
                pListMult [8] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%ibY%icMultBS", n, n));
                pListMult [9] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPaMultBS", n));
                pListMult [10] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPbMultBS", n));
                pListMult [11] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPcMultBS", n));

                pListMult [12] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsibMultBS", n));
                pListMult [13] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsicMultBS", n));
                pListMult [14] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsib_PsicMultBS", n));
                pListMult [15] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPaMultBS", n));
                pListMult [16] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPbMultBS", n));
                pListMult [17] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPcMultBS", n));
                pListMult [18] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsibMultBS", n, n));
                pListMult [19] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsicMultBS", n, n));
                pListMult [20] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsibSin%iPsicMultBS", n, n));
                pListMult [21] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPaMultBS", n));
                pListMult [22] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPbMultBS", n));
                pListMult [23] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPcMultBS", n));

                p2ListMult [0] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%iaX%ibMult", n, n));
                p2ListMult [1] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%iaX%icMult", n, n));
                p2ListMult [2] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2X%ibX%icMult", n, n));
                p2ListMult [3] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPaMult", n));
                p2ListMult [4] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPbMult", n));
                p2ListMult [5] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixSPcMult", n));
                p2ListMult [6] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%iaY%ibMult", n, n));
                p2ListMult [7] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%iaY%icMult", n, n));
                p2ListMult [8] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Y%ibY%icMult", n, n));
                p2ListMult [9] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPaMult", n));
                p2ListMult [10] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPbMult", n));
                p2ListMult [11] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iySPcMult", n));

                p2ListMult [12] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsia_PsibMult", n));
                p2ListMult [13] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsia_PsicMult", n));
                p2ListMult [14] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Cos%iPsib_PsicMult", n));
                p2ListMult [15] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPaMult", n));
                p2ListMult [16] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPbMult", n));
                p2ListMult [17] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%ixEPcMult", n));
                p2ListMult [18] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsiaSin%iPsibMult", n, n));
                p2ListMult [19] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsiaSin%iPsicMult", n, n));
                p2ListMult [20] = (TProfile2D*) inputFile2 -> Get (dirName [step] + Form ("Source Histograms/p2Sin%iPsibSin%iPsicMult", n, n));
                p2ListMult [21] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPaMult", n));
                p2ListMult [22] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPbMult", n));
                p2ListMult [23] = (TProfile2D*) inputFile -> Get (dirName [step] + Form ("p2R%iyEPcMult", n));

                for (Int_t j = 0; j < nBinsCent; j++) {
                    hListCent [j * nStacks + 0] = new TH1F (Form ("X%iaX%ibCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 1] = new TH1F (Form ("X%iaX%icCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 2] = new TH1F (Form ("X%ibX%icCentSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 3] = new TH1F (Form ("R%ixCentSPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{x,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 4] = new TH1F (Form ("R%ixCentSPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{x,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 5] = new TH1F (Form ("R%ixCentSPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{x,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 6] = new TH1F (Form ("Y%iaY%ibCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 7] = new TH1F (Form ("Y%iaY%icCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 8] = new TH1F (Form ("Y%ibY%icCentSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT, centrality class %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 9] = new TH1F (Form ("R%iyCentSPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{y,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 10] = new TH1F (Form ("R%iyCentSPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{y,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 11] = new TH1F (Form ("R%iyCentSPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{y,SP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);

                    hListCent [j * nStacks + 12] = new TH1F (Form ("X%iaX%ibCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT, centrality class %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 13] = new TH1F (Form ("X%iaX%icCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT, centrality class %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 14] = new TH1F (Form ("X%ibX%icCentEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT, centrality class %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 15] = new TH1F (Form ("R%ixCentEPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{x,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 16] = new TH1F (Form ("R%ixCentEPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{x,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 17] = new TH1F (Form ("R%ixCentEPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{x,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 18] = new TH1F (Form ("Y%iaY%ibCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT, centrality class %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 19] = new TH1F (Form ("Y%iaY%icCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT, centrality class %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 20] = new TH1F (Form ("Y%ibY%icCentEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT, centrality class %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListCent [j * nStacks + 21] = new TH1F (Form ("R%iyCentEPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{y,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 22] = new TH1F (Form ("R%iyCentEPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{y,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListCent [j * nStacks + 23] = new TH1F (Form ("R%iyCentEPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{y,EP}, centrality class %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                }

                for (Int_t j = 0; j < nBinsMult; j++) {
                    hListMult [j * nStacks + 0] = new TH1F (Form ("X%iaX%ibMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 1] = new TH1F (Form ("X%iaX%icMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 2] = new TH1F (Form ("X%ibX%icMultSP_%i", n, n, j + 1), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 3] = new TH1F (Form ("R%ixMultSPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{x,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 4] = new TH1F (Form ("R%ixMultSPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{x,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 5] = new TH1F (Form ("R%ixMultSPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{x,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 6] = new TH1F (Form ("Y%iaY%ibMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 7] = new TH1F (Form ("Y%iaY%icMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 8] = new TH1F (Form ("Y%ibY%icMultSP_%i", n, n, j + 1), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT, multiplicity bin %i, file %i", n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 9] = new TH1F (Form ("R%iyMultSPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{y,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 10] = new TH1F (Form ("R%iyMultSPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{y,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 11] = new TH1F (Form ("R%iyMultSPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{y,SP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);

                    hListMult [j * nStacks + 12] = new TH1F (Form ("X%iaX%ibMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT, multiplicity bin %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 13] = new TH1F (Form ("X%iaX%icMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT, multiplicity bin %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 14] = new TH1F (Form ("X%ibX%icMultEP_%i", n, n, j + 1), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT, multiplicity bin %i, file %i", n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 15] = new TH1F (Form ("R%ixMultEPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{x,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 16] = new TH1F (Form ("R%ixMultEPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{x,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 17] = new TH1F (Form ("R%ixMultEPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{x,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 18] = new TH1F (Form ("Y%iaY%ibMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT, multiplicity bin %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 19] = new TH1F (Form ("Y%iaY%icMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT, multiplicity bin %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 20] = new TH1F (Form ("Y%ibY%icMultEP_%i", n, n, j + 1), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT, multiplicity bin %i, file %i", n, n, n, n, j + 1, f), nBinsCorr, corrMin, corrMax);
                    hListMult [j * nStacks + 21] = new TH1F (Form ("R%iyMultEPa_%i_%i", n, j + 1, f), Form ("R_{%i,a}^{y,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 22] = new TH1F (Form ("R%iyMultEPb_%i_%i", n, j + 1, f), Form ("R_{%i,b}^{y,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                    hListMult [j * nStacks + 23] = new TH1F (Form ("R%iyMultEPc_%i_%i", n, j + 1, f), Form ("R_{%i,c}^{y,EP}, multiplicity bin %i, file %i", n, j + 1, f), nBinsR, Rmin, Rmax);
                }

                for (Int_t j = 0; j < nStacks; j++) {
                    pListCent [j] -> SetMarkerStyle (markerStyles [f]);
                    pListCent [j] -> SetMarkerColor (markerColors [f]);
                    pListCent [j] -> SetLineColor (markerColors [f]);
                    HistShift (pListCent [j], shifts [f]);
                    hsListCent [j] -> Add (pListCent [j]);

                    pListMult [j] -> SetMarkerStyle (markerStyles [f]);
                    pListMult [j] -> SetMarkerColor (markerColors [f]);
                    pListMult [j] -> SetLineColor (markerColors [f]);
                    HistShift (pListMult [j], shifts [f]);
                    hsListMult [j] -> Add (pListMult [j]);

                    for (Int_t k = 0; k < nBinsCent; k++) {
                       for (Int_t l = 0; l < nBinsBS [f]; l++) {
                            hListCent [k * nStacks + j] -> Fill (p2ListCent [j] -> GetBinContent (k + 1, l + 1));
                        }
                        hListCent [k * nStacks + j] -> SetLineColor (markerColors [f]);
                        hsListCent2 [k * nStacks + j] -> Add (hListCent [k * nStacks + j]);
                    }

                    for (Int_t k = 0; k < nBinsMult; k++) {
                       for (Int_t l = 0; l < nBinsBS [f]; l++) {
                            hListMult [k * nStacks + j] -> Fill (p2ListMult [j] -> GetBinContent (k + 1, l + 1));
                        }
                        hListMult [k * nStacks + j] -> SetLineColor (markerColors [f]);
                        hsListMult2 [k * nStacks + j] -> Add (hListMult [k * nStacks + j]);
                    }
                }
            }

            outputDir -> cd ();

            for (Int_t j = 0; j < 4; j++) {
                for (Int_t k = 0; k < 6; k++) {
                    c1 -> cd (k + 1);
                    if (k == 3 || k == 4 || k == 5) {
                        hsListCent [j * 6 + k] -> SetMaximum (0.2);
                        hsListMult [j * 6 + k] -> SetMinimum (-0.05);
                    }
                    hsListCent [j * 6 + k] -> Draw ("nostack p e1X0");
                    if (k == 0) leg -> Draw ();
                }
                c1 -> Write (canvasName [j] + Form ("Cent_%i", n));

                for (Int_t k = 0; k < 6; k++) {
                    c1 -> cd (k + 1);
                    if (k == 3 || k == 4 || k == 5) {
                        hsListMult [j * 6 + k] -> SetMaximum (0.2);
                        hsListMult [j * 6 + k] -> SetMinimum (-0.05);
                    }
                    hsListCent [j * 6 + k] -> Draw ("nostack p e1X0");
                    if (k == 0) leg -> Draw ();
                }
                c1 -> Write (canvasName [j] + Form ("Mult_%i", n));
            }

            for (Int_t l = 0; l < nBinsCent; l++) {
                for (Int_t j = 0; j < 4; j++) {
                    for (Int_t k = 0; k < 6; k++) {
                        c1 -> cd (k + 1);
                        hsListCent2 [l * nStacks + j * 6 + k] -> Draw ("nostack");
                        if (k == 0) leg2 -> Draw ();
                    }
                    c1 -> Write (canvasName [j] + Form ("distr_Cent_%i_%i", n, l + 1));
                }
            }

            for (Int_t l = 0; l < nBinsMult; l++) {
                for (Int_t j = 0; j < 4; j++) {
                    for (Int_t k = 0; k < 6; k++) {
                        c1 -> cd (k + 1);
                        hsListMult2 [l * nStacks + j * 6 + k] -> Draw ("nostack");
                        if (k == 0) leg2 -> Draw ();
                    }
                    c1 -> Write (canvasName [j] + Form ("distr_Mult_%i_%i", n, l + 1));
                }
            }
        }
    }
    outputFile -> Close ();
}
