#include <iostream>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TProfile.h>
#include <THStack.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>

#include "../FlowReconstructor/CFlowReconstructor.h"

static const int nFiles = 8, nStacks = 24;

void RapidityScan (TString location, TString subevent);
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


void RapScan (Int_t nFile = 0, Int_t se = 1) {
    TString location = "/mnt/pool/rhic/2/ovgol/NA49_flow/01d/RapidityScan/";
    TString nonUniformInputFileName = "/mnt/pool/rhic/2/ovgol/NA49_conv/01d/01d_eta"; // NA49
    TString histFileName;
    TString subevent [3] = {"a", "b", "c"};

	CFlowReconstructor flowReconstructor;
	flowReconstructor.SetNonUniformInputFileName (nonUniformInputFileName);
    histFileName = location + "NA49_y_" + subevent [se - 1] + Form ("%i_pion", nFile);

	flowReconstructor.AddHarmonic (1);
	flowReconstructor.AddHarmonic (2);
	flowReconstructor.UseZeroSubevents ();
	flowReconstructor.SetVariable("y");
	flowReconstructor.SetNbinsBS (10);
	flowReconstructor.SetNrunRange (3133, 3166);
	flowReconstructor.SetMhRange (10, 510); // all
	flowReconstructor.SetNbinsMh (10); // all
	flowReconstructor.SetCentRange (0.0, 0.6); // all
	flowReconstructor.SetNbinsCent (6); // all

	flowReconstructor.SetPtRange (0.0, 2.5);
	flowReconstructor.SetNbinsPt (10);
	flowReconstructor.SetEtaRange (-1.0, 3.0); //symm
	flowReconstructor.SetNbinsEta (8); // symm


	flowReconstructor.SetPtAveragingRange (1, 0.1, 2.0); // mix pion
	flowReconstructor.SetPtAveragingRange (2, 0.1, 2.0); // mix pion
	flowReconstructor.SetEtaAveragingRange (1, 0.0, 1.0); // mix pion
	flowReconstructor.SetEtaAveragingRange (2, -1.0, 1.0); // mix pion
	flowReconstructor.SetPtSubeventsLimits (1, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0); // mix pion
	flowReconstructor.SetPtSubeventsLimits (2, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0); // mix pion
	flowReconstructor.SetEtaSubeventsLimits (1, 1.0, 3.0, 1.0, 3.0, 1.0, 3.0); // mix pion
    flowReconstructor.SetEtaSubeventsLimits (2, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0); // mix pion
	flowReconstructor.AddFlowParticle (kPionMinus); // mix pion

    flowReconstructor.AddResolutionParticle (1, 1, kPionPlus); // mix pion
    flowReconstructor.AddResolutionParticle (1, 2, kProton); // mix pion
    flowReconstructor.AddResolutionParticle (1, 3, kPionMinus); // mix pion
    flowReconstructor.AddResolutionParticle (2, 1, kPionPlus); // mix pion
    flowReconstructor.AddResolutionParticle (2, 2, kProton); // mix pion
    flowReconstructor.AddResolutionParticle (2, 3, kPionMinus); // mix pion

    flowReconstructor.SetResolutionSigns(1, -1, 1, -1); // mix pion
    flowReconstructor.SetResolutionSigns(2, 1, 1, 1); // mix pion

    if (nFile == 0) {
        RapidityScan (location, subevent [se - 1]);
    }

    else {
        flowReconstructor.SetPtAveragingRange (1, 0.1, 2.0);
        flowReconstructor.SetPtAveragingRange (2, 0.1, 2.0);
        flowReconstructor.SetEtaAveragingRange (1, 0.0, 1.0);
        flowReconstructor.SetEtaAveragingRange (2, -1.0, 1.0);
        flowReconstructor.SetPtSubeventsLimits (1, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0);
        flowReconstructor.SetPtSubeventsLimits (2, 0.1, 2.0, 0.1, 2.0, 0.1, 2.0);
        flowReconstructor.AddFlowParticle (kPionMinus);
        flowReconstructor.AddResolutionParticle (1, 1, kPionPlus);
        flowReconstructor.AddResolutionParticle (1, 2, kProton);
        flowReconstructor.AddResolutionParticle (1, 3, kPionMinus);
        flowReconstructor.AddResolutionParticle (2, 1, kPionPlus);
        flowReconstructor.AddResolutionParticle (2, 2, kProton);
        flowReconstructor.AddResolutionParticle (2, 3, kPionMinus);
        flowReconstructor.SetResolutionSigns (1, -1, 1, -1);
        flowReconstructor.SetResolutionSigns (2, 1, 1, 1);


        Float_t etaLimits [nFiles + 1] = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

        if (se == 1) {
            flowReconstructor.SetEtaSubeventsLimits (1, etaLimits [nFile - 1], etaLimits [nFile], 1.0, 3.0, 1.0, 3.0);
            flowReconstructor.SetEtaSubeventsLimits (2, etaLimits [nFile - 1], etaLimits [nFile], -1.0, 1.0, -1.0, 0.0);
        }

        if (se == 2) {
            flowReconstructor.SetEtaSubeventsLimits (1, 1.0, 3.0, etaLimits [nFile - 1], etaLimits [nFile], 1.0, 3.0);
            flowReconstructor.SetEtaSubeventsLimits (2, -1.0, 1.0, etaLimits [nFile - 1], etaLimits [nFile], -1.0, 0.0);
        }

        if (se == 3) {
            flowReconstructor.SetEtaSubeventsLimits (1, 1.0, 3.0, 2.0, 2.5, etaLimits [nFile - 1], etaLimits [nFile]);
            flowReconstructor.SetEtaSubeventsLimits (2, -1.0, 1.0, 2.0, 2.5, etaLimits [nFile - 1], etaLimits [nFile]);
        }


        flowReconstructor.SetHistFileName (histFileName);
        flowReconstructor.GetCorrelations ();
        flowReconstructor.CalculateFlow (histFileName);
    }
}

void RapidityScan (TString location, TString subevent) {
    Int_t n, nSteps = 3;
    static const Int_t nHarmonics = 2;
    Int_t harmonicsMap [2] = {1, 2};
    TString	dirName [5] = {"Not_Corrected/", "Recentered/", "Diagonalized/", "Uniform_Acceptance/", "Analytic/"};
    TString	stepName [4] = {"not corrected", "recentered", "diagonalized", "uniform"};
    TString canvasName [4] = {"X_SP_", "Y_SP_", "X_EP_", "Y_EP_"};
    Int_t markerColors [8] = {1, 2, 3, 4, 6, 8, 9, 46};
    Int_t markerStyles [8] = {24, 25, 26, 27, 28, 30, 32, 5};
    Float_t shifts [8] = {0.0, 0.1, -0.1, 0.2, -0.2, 0.3, -0.3, 0.4};
    TString yIntervals [8] = {"#it{y}#in[-1.0, -0.5]", "#it{y}#in[-0.5, 0.0]", "#it{y}#in[0.0, 0.5]", "#it{y}#in[0.5, 1.0]", "#it{y}#in[1.0, 1.5]", "#it{y}#in[1.5, 2.0]", "#it{y}#in[2.0, 2.5]", "#it{y}#in[2.5, 3.0]"};
    TDirectory *outputDir;


    static const int nStacks = 24;
    TString outputFileName = "NA49_pion_RS_" + subevent + ".root";
    TFile *inputFile;
    TFile *outputFile = new TFile (location + outputFileName, "recreate");

    TCanvas *c1 = new TCanvas ("c1", "c1", 800, 600);
    gStyle -> SetLegendBorderSize (0);
    c1 -> Divide (3, 2);

    TLegend *leg  = new TLegend (0.7, 0.15, 0.85, 0.5);
    leg -> SetFillColor (0);
    leg -> SetTextSize(0.04);

    TH1F *fakeHist [nFiles];
    for (Int_t f = 0; f < nFiles; f++) {
        fakeHist [f] = new TH1F ();
        fakeHist [f] -> SetMarkerStyle (markerStyles [f]);
        fakeHist [f] -> SetMarkerColor (markerColors [f]);
        leg -> AddEntry (fakeHist [f], yIntervals [f], "p");
    }

    TProfile *histListCent [nStacks], *histListMult [nStacks];
    THStack *hstackListCent [nStacks], *hstackListMult [nStacks];

    for (Int_t step = 0; step < nSteps; step++) {
        outputDir = outputFile -> mkdir (dirName [step]);

        for (Int_t i = 0; i < nHarmonics; i++) {
            n = harmonicsMap [i];

            hstackListCent [0] = new THStack (Form ("X%iaX%ibCentSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT ;centrality", n, n));
            hstackListCent [1] = new THStack (Form ("X%iaX%icCentSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT ;centrality", n, n));
            hstackListCent [2] = new THStack (Form ("X%ibX%icCentSP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT ;centrality", n, n));
            hstackListCent [3]  = new THStack (Form ("R%ixCentSPa", n, n), Form ("R_{%i,a}^{x,SP} ;centrality", n));
            hstackListCent [4]  = new THStack (Form ("R%ixCentSPb", n, n), Form ("R_{%i,b}^{x,SP} ;centrality", n));
            hstackListCent [5]  = new THStack (Form ("R%ixCentSPc", n, n), Form ("R_{%i,c}^{x,SP} ;centrality", n));
            hstackListCent [6]  = new THStack (Form ("Y%iaY%ibCentSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT ;centrality", n, n));
            hstackListCent [7]  = new THStack (Form ("Y%iaY%icCentSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT ;centrality", n, n));
            hstackListCent [8]  = new THStack (Form ("Y%ibY%icCentSP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT ;centrality", n, n));
            hstackListCent [9]  = new THStack (Form ("R%iyCentSPa", n), Form ("R_{%i,a}^{y,SP} ;centrality", n));
            hstackListCent [10]  = new THStack (Form ("R%iyCentSPb", n), Form ("R_{%i,b}^{y,SP} ;centrality", n));
            hstackListCent [11]  = new THStack (Form ("R%iyCentSPc", n), Form ("R_{%i,c}^{y,SP} ;centrality", n));

            hstackListCent [12]  = new THStack (Form ("X%iaX%ibCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT ;centrality", n, n, n));
            hstackListCent [13]  = new THStack (Form ("X%iaX%icCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT ;centrality", n, n, n));
            hstackListCent [14]  = new THStack (Form ("X%ibX%icCentEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT ;centrality", n, n, n));
            hstackListCent [15]  = new THStack (Form ("R%ixCentEPa", n, n), Form ("R_{%i,a}^{x,EP} ;centrality", n));
            hstackListCent [16]  = new THStack (Form ("R%ixCentEPb", n, n), Form ("R_{%i,b}^{x,EP} ;centrality", n));
            hstackListCent [17]  = new THStack (Form ("R%ixCentEPc", n, n), Form ("R_{%i,c}^{x,EP} ;centrality", n));
            hstackListCent [18]  = new THStack (Form ("Y%iaY%ibCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT ;centrality", n, n, n, n));
            hstackListCent [19]  = new THStack (Form ("Y%iaY%icCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT ;centrality", n, n, n, n));
            hstackListCent [20]  = new THStack (Form ("Y%ibY%icCentEP", n, n), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT ;centrality", n, n, n, n));
            hstackListCent [21]  = new THStack (Form ("R%iyCentEPa", n, n), Form ("R_{%i,a}^{y,EP} ;centrality", n));
            hstackListCent [22]  = new THStack (Form ("R%iyCentEPb", n, n), Form ("R_{%i,b}^{y,EP} ;centrality", n));
            hstackListCent [23]  = new THStack (Form ("R%iyCentEPc", n, n), Form ("R_{%i,c}^{y,EP} ;centrality", n));

            hstackListMult [0] = new THStack (Form ("X%iaX%ibMultSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{b}#GT ;multiplicity", n, n));
            hstackListMult [1] = new THStack (Form ("X%iaX%icMultSP", n, n), Form ("#LTX_{%i}^{a}X_{%i}^{c}#GT ;multiplicity", n, n));
            hstackListMult [2] = new THStack (Form ("X%ibX%icMultSP", n, n), Form ("#LTX_{%i}^{b}X_{%i}^{c}#GT ;multiplicity", n, n));
            hstackListMult [3] = new THStack (Form ("R%ixMultSPa", n, n), Form ("R_{%i,a}^{x,SP} ;multiplicity", n));
            hstackListMult [4] = new THStack (Form ("R%ixMultSPb", n, n), Form ("R_{%i,b}^{x,SP} ;multiplicity", n));
            hstackListMult [5] = new THStack (Form ("R%ixMultSPc", n, n), Form ("R_{%i,c}^{x,SP} ;multiplicity", n));
            hstackListMult [6] = new THStack (Form ("Y%iaY%ibMultSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{b}#GT ;multiplicity", n, n));
            hstackListMult [7] = new THStack (Form ("Y%iaY%icMultSP", n, n), Form ("#LTY_{%i}^{a}Y_{%i}^{c}#GT ;multiplicity", n, n));
            hstackListMult [8] = new THStack (Form ("Y%ibY%icMultSP", n, n), Form ("#LTY_{%i}^{b}Y_{%i}^{c}#GT ;multiplicity", n, n));
            hstackListMult [9] = new THStack (Form ("R%iyMultSPa", n, n), Form ("R_{%i,a}^{y,SP} ;multiplicity", n));
            hstackListMult [10] = new THStack (Form ("R%iyMultSPb", n, n), Form ("R_{%i,b}^{y,SP} ;multiplicity", n));
            hstackListMult [11] = new THStack (Form ("R%iyMultSPc", n, n), Form ("R_{%i,c}^{y,SP} ;multiplicity", n));

            hstackListMult [12] = new THStack (Form ("X%iaX%ibMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{b}))#GT ;multiplicity", n, n, n));
            hstackListMult [13] = new THStack (Form ("X%iaX%icMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{a}-#psi_{%i}^{c}))#GT ;multiplicity", n, n, n));
            hstackListMult [14] = new THStack (Form ("X%ibX%icMultEP", n, n), Form ("#LTcos(%i(#psi_{%i}^{b}-#psi_{%i}^{c}))#GT ;multiplicity", n, n, n));
            hstackListMult [15] = new THStack (Form ("R%ixMultEPa", n, n), Form ("R_{%i,a}^{x,EP} ;multiplicity", n));
            hstackListMult [16] = new THStack (Form ("R%ixMultEPb", n, n), Form ("R_{%i,b}^{x,EP} ;multiplicity", n));
            hstackListMult [17] = new THStack (Form ("R%ixMultEPc", n, n), Form ("R_{%i,c}^{x,EP} ;multiplicity", n));
            hstackListMult [18] = new THStack (Form ("Y%iaY%ibMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{b})#GT ;multiplicity", n, n, n, n));
            hstackListMult [19] = new THStack (Form ("Y%iaY%icMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{a}sin%i#psi_{%i}^{c})#GT ;multiplicity", n, n, n, n));
            hstackListMult [20] = new THStack (Form ("Y%ibY%icMultEP", n, n), Form ("#LTsin(%i#psi_{%i}^{b}sin%i#psi_{%i}^{c})#GT ;multiplicity", n, n, n, n));
            hstackListMult [21] = new THStack (Form ("R%iyMultEPa", n, n), Form ("R_{%i,a}^{y,EP} ;multiplicity", n));
            hstackListMult [22] = new THStack (Form ("R%iyMultEPb", n, n), Form ("R_{%i,b}^{y,EP} ;multiplicity", n));
            hstackListMult [23] = new THStack (Form ("R%iyMultEPc", n, n), Form ("R_{%i,c}^{y,EP} ;multiplicity", n));

            for (Int_t f = 0, shift = 0.0; f < nFiles; f++, shift += 0.05) {
                shift = f * 0.05;
                shift *= -1.0;
                if (f == 0) cout << stepName [step] << endl;
                cout << location + "NA49_y_" + subevent + Form ("%i_pion_flow.root", f + 1) << endl;
                inputFile = new TFile (location + "NA49_y_" + subevent + Form ("%i_pion_flow.root", f + 1), "read");

                histListCent [0] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%ibCentBS", n, n));
                histListCent [1] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%icCentBS", n, n));
                histListCent [2] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%ibX%icCentBS", n, n));
                histListCent [3] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPaCentBS", n));
                histListCent [4] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPbCentBS", n));
                histListCent [5] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPcCentBS", n));
                histListCent [6] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%ibCentBS", n, n));
                histListCent [7] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%icCentBS", n, n));
                histListCent [8] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%ibY%icCentBS", n, n));
                histListCent [9] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPaCentBS", n));
                histListCent [10] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPbCentBS", n));
                histListCent [11] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPcCentBS", n));
                histListCent [12] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsibCentBS", n));
                histListCent [13] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsicCentBS", n));
                histListCent [14] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsib_PsicCentBS", n));
                histListCent [15] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPaCentBS", n));
                histListCent [16] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPbCentBS", n));
                histListCent [17] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPcCentBS", n));
                histListCent [18] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsibCentBS", n, n));
                histListCent [19] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsicCentBS", n, n));
                histListCent [20] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsibSin%iPsicCentBS", n, n));
                histListCent [21] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPaCentBS", n));
                histListCent [22] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPbCentBS", n));
                histListCent [23] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPcCentBS", n));

                histListMult [0] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%ibMultBS", n, n));
                histListMult [1] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%iaX%icMultBS", n, n));
                histListMult [2] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pX%ibX%icMultBS", n, n));
                histListMult [3] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPaMultBS", n));
                histListMult [4] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPbMultBS", n));
                histListMult [5] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixSPcMultBS", n));
                histListMult [6] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%ibMultBS", n, n));
                histListMult [7] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%iaY%icMultBS", n, n));
                histListMult [8] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pY%ibY%icMultBS", n, n));
                histListMult [9] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPaMultBS", n));
                histListMult [10] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPbMultBS", n));
                histListMult [11] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iySPcMultBS", n));
                histListMult [12] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsibMultBS", n));
                histListMult [13] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsia_PsicMultBS", n));
                histListMult [14] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pCos%iPsib_PsicMultBS", n));
                histListMult [15] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPaMultBS", n));
                histListMult [16] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPbMultBS", n));
                histListMult [17] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%ixEPcMultBS", n));
                histListMult [18] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsibMultBS", n, n));
                histListMult [19] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsiaSin%iPsicMultBS", n, n));
                histListMult [20] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pSin%iPsibSin%iPsicMultBS", n, n));
                histListMult [21] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPaMultBS", n));
                histListMult [22] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPbMultBS", n));
                histListMult [23] = (TProfile*) inputFile -> Get (dirName [step] + Form ("pR%iyEPcMultBS", n));

                for (Int_t j = 0; j < nStacks; j++) {
                    histListCent [j] -> SetMarkerStyle (markerStyles [f]);
                    histListCent [j] -> SetMarkerColor (markerColors [f]);
                    histListCent [j] -> SetLineColor (markerColors [f]);
                    HistShift (histListCent [j], shifts [f]);
                    hstackListCent [j] -> Add (histListCent [j]);
                    histListMult [j] -> SetMarkerStyle (markerStyles [f]);
                    histListMult [j] -> SetMarkerColor (markerColors [f]);
                    histListMult [j] -> SetLineColor (markerColors [f]);
                    HistShift (histListMult [j], shifts [f]);
                    hstackListMult [j] -> Add (histListMult [j]);
                }
                //inputFile -> Close ();
            }

            outputDir -> cd ();

            for (Int_t j = 0; j < 4; j++) {
                for (Int_t k = 0; k < 6; k++) {
                    c1 -> cd (k + 1);
                    hstackListCent [j * 6 + k] -> Draw ("nostack p e1X0");
                    if (k == 0) leg -> Draw ();
                }
                c1 -> Write (canvasName [j] + Form ("Cent_%i", n));

                for (Int_t k = 0; k < 6; k++) {
                    c1 -> cd (k + 1);
                    hstackListMult [j * 6 + k] -> Draw ("nostack p e1X0");
                    if (k == 0) leg -> Draw ();
                }
                c1 -> Write (canvasName [j] + Form ("Mult_%i", n));
            }
        }
    }
    outputFile -> Close ();
}
