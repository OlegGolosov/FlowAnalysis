#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>

using namespace std;

void results ()
{
    static const int nSets = 3;
    static const int nMethods = 2;
    const int step = 1;
    int colors [6] = {kMagenta + 2, kBlue + 2, kGreen +2, kRed + 2, kViolet, kCyan + 2};
    int styles [6] = {20, 21 ,22 ,23, 33, 34};



    TH1 *R_EP [nSets] [nMethods], *Vpt [nSets] [nMethods], *Vy [nSets] [nMethods], *R_OG [nSets], *R_Oli [nSets], *R_AS;
    TGraphErrors *R_BK, *Vpt_BK, *Vy_BK;
    TString stepNames [] = {"Not corrected", "Recentered", "Diagonalized"};
    TString setNames [] = {"protons", "weight", "no_weight"};
    TString methods [] = {"FWRS", "MDC3S", "FW3S"};
    TString methodNames [] = {"RS", "3Sub [protons]", "3Sub [FW]"};
    TString path;
    TFile *flowFile, *resFile;

//    GET HISTOGRAMS FROM FILES
    for (int set = 0; set < nSets; set++)
    {
        path = "../../HADES_flow_" + setNames [set] + "/";
        resFile = new TFile (path + "resolution.root", "read");
        if (!resFile) continue;
        R_Oli [set] = (TH1F*) resFile -> Get ("hR_Oli_" + setNames [set]);
        R_Oli [set] -> SetTitle ("RS extrapolation (OG)");
        R_OG [set] = (TH1D*) resFile -> Get ("hR_OG_" + setNames [set]);
        R_OG [set] -> SetTitle ("RS Tree (OG)");
        if (set == 2)
        {
            R_AS = (TH1F*) resFile -> Get ("hR_AS_" + setNames [set]);
            R_AS -> SetTitle ("RS HYDRA (AS)");
        }

        for (int meth = 0; meth < nMethods; meth++)
        {
            flowFile = new TFile (path + "HADES_" + methods [meth] + "_flow.root", "read");
            if (!flowFile)
            {
                cout << "No flow file! \n";
                continue;
            }
            if (set == 0 && meth == 0)
            {
                R_BK = (TGraphErrors*) flowFile -> Get ("gR1cent");
                R_BK -> SetName ("R_BK");
                R_BK -> SetTitle ("RS extrapolation (BK)");
                Vpt_BK = (TGraphErrors*) flowFile -> Get ("gV1ptPos");
                Vpt_BK -> SetTitle ("V1_pt_BK");
                Vy_BK = (TGraphErrors*) flowFile -> Get ("gV1y");
                Vy_BK -> SetTitle ("V1_y_BK");
            }
            if (meth == 0) R_EP [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Resolution/hR1aCent_EP");
            if (meth == 1) R_EP [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Resolution/hR1cCent_EP");
            R_EP [set] [meth] -> SetTitle (methodNames [meth] + " Code (OG)");
            if (1)
            {
                Vpt [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Flow/hV1PtCent_EP");
                Vy [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Flow/hV1EtaCent_EP");
//                Vpt [set] [meth] = (TH1F*) flowFile -> Get ("Recentered/Flow/hV1PtCent_EP");
//                Vy [set] [meth] = (TH1F*) flowFile -> Get ("Recentered/Flow/hV1EtaCent_EP");
                cout << stepNames [step] + "/Flow/hV1PtCent_EP\n";
                cout << stepNames [step] + "/Flow/hV1EtaCent_EP\n";
            }
            if (meth == 1)
            {
//                Vpt [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Flow/hV1cPtCent_EP");
//                Vy [set] [meth] = (TH1F*) flowFile -> Get (stepNames [step] + "/Flow/hV1cEtaCent_EP");
                Vpt [set] [meth] = (TH1F*) flowFile -> Get ("Recentered/Flow/hV1cPtCent_EP");
                Vy [set] [meth] = (TH1F*) flowFile -> Get ("Recentered/Flow/hV1cEtaCent_EP");
                cout << stepNames [step] + "/Flow/hV1cPtCent_EP\n";
                cout << stepNames [step] + "/Flow/hV1cEtaCent_EP\n";
            }
                if (Vy [set] [meth] == 0) cout << "y NULL!!!\n";
                if (Vpt [set] [meth] == 0) cout << "pt NULL!!!\n";
            Vpt [set] [meth] -> SetTitle ("V_{1} (p_{T}) (" + setNames [set] + ", " + methods [meth] + ")");
            Vy [set] [meth] -> SetTitle ("V_{1} (#it{y}) (" + setNames [set] + ", " + methods [meth] + ")");
        }
    }
    TFile *resultsFile = new TFile ("results.root", "recreate");
    TDirectory *setDir [nSets], *methDir [nMethods] [nSets];

//    DRAW RESOLUTION
    gStyle -> SetOptStat (0);
    gROOT -> ForceStyle ();
    TCanvas *c1 = new TCanvas ();

cout << "HERE0!\n";
    for (int set = 0; set < nSets; set++)
    {
        setDir [set] = resultsFile -> mkdir (setNames [set]);
        R_Oli [set] -> Draw ("p");
        R_Oli [set] -> GetYaxis () -> SetRangeUser (0., 1.);
        R_Oli [set] -> SetLineColor (colors [0]);
        R_Oli [set] -> SetMarkerColor (colors [0]);
        R_Oli [set] -> SetMarkerStyle (styles [0]);
        R_Oli [set] -> SetMarkerSize (1.5);
//        R_OG [set] -> Draw ("same e1x0");
        R_OG [set] -> SetLineColor (colors [1]);
        R_OG [set] -> SetMarkerColor (colors [1]);
        R_OG [set] -> SetMarkerStyle (styles [1]);
        R_OG [set] -> SetMarkerSize (1.5);
        R_BK -> Draw ("same p");
        R_BK -> SetLineColor (colors [2]);
        R_BK -> SetMarkerColor (colors [2]);
        R_BK -> SetMarkerStyle (styles [2]);
        R_BK -> SetMarkerSize (1.5);

cout << "HERE1.5!\n";
        if (set == 2)
        {
            R_AS -> Draw ("same e1x0");
            R_AS -> SetLineColor (colors [3]);
            R_AS -> SetMarkerColor (colors [3]);
            R_AS -> SetMarkerStyle (styles [3]);
            R_AS -> SetMarkerSize (1.5);
        }
        cout << "HERE1!\n";
        for (int meth = 0; meth < nMethods; meth++)
        {
            methDir [set] [meth] = setDir [set] -> mkdir (methods [meth]);
            R_EP [set] [meth] -> Draw ("same e1x0");
            R_EP [set] [meth] -> SetLineColor (colors [meth + 4]);
            R_EP [set] [meth] -> SetMarkerColor (colors [meth + 4]);
            R_EP [set] [meth] -> SetMarkerStyle (styles [meth + 4]);
            R_EP [set] [meth] -> SetMarkerSize (1.5);
        }

cout << "HERE2!\n";
        c1 -> BuildLegend ();
        setDir [set] -> cd ();
        c1 -> Write ("cR_" + setNames [set]);
    }

cout << "HERE3!\n";
//    DRAW FLOW PT
    for (int set = 0; set < nSets; set++)
    {
        Vpt_BK -> Draw ("AP");
        Vpt_BK -> SetLineColor (colors [0]);
        Vpt_BK -> SetMarkerColor (colors [0]);
        Vpt_BK -> SetMarkerStyle (styles [0]);
        Vpt_BK -> SetMarkerSize (1.5);

        for (int meth = 0; meth < nMethods; meth++)
        {
            Vpt [set] [meth] -> Draw ("same e1x0");
            Vpt [set] [meth] -> SetLineColor (colors [meth + 1]);
            Vpt [set] [meth] -> SetMarkerColor (colors [meth + 1]);
            Vpt [set] [meth] -> SetMarkerStyle (styles [meth + 1]);
            Vpt [set] [meth] -> SetMarkerSize (1.5);
        }

        c1 -> BuildLegend ();
        setDir [set] -> cd ();
        c1 -> Write ("cVpt_" + setNames [set]);
    }

//    DRAW FLOW Y
    for (int set = 0; set < nSets; set++)
    {
        Vy_BK -> Draw ("AP");
        Vy_BK -> SetLineColor (colors [0]);
        Vy_BK -> SetMarkerColor (colors [0]);
        Vy_BK -> SetMarkerStyle (styles [0]);
        Vy_BK -> SetMarkerSize (1.5);

        for (int meth = 0; meth < nMethods; meth++)
        {
            Vy [set] [meth] -> Draw ("same e1x0");
            Vy [set] [meth] -> SetLineColor (colors [meth + 1]);
            Vy [set] [meth] -> SetMarkerColor (colors [meth + 1]);
            Vy [set] [meth] -> SetMarkerStyle (styles [meth + 1]);
            Vy [set] [meth] -> SetMarkerSize (1.5);
        }

        c1 -> BuildLegend ();
        setDir [set] -> cd ();
        c1 -> Write ("cVy_" + setNames [set]);
    }

    cout << "Writing results to file...\n";
    resultsFile -> cd ();
    R_BK -> Write ("gR1");
    Vpt_BK -> Write ("gV1ptPos");
    Vy_BK -> Write ("gVy");


    for (int set = 0; set < nSets; set++)
    {
        setDir [set] -> cd ();
        R_Oli [set] -> Write ();
        R_OG [set] -> Write ();
        if (set == 3) R_AS -> Write ();
        for (int meth = 0; meth < nMethods; meth++)
        {
            methDir [set] [meth] -> cd ();
            R_EP [set] [meth] -> Write ();
            Vpt [set] [meth] -> Write ();
            Vy [set] [meth] -> Write ();
        }
    }

    resultsFile -> Close ();
    cout << "ALIVE!!!\n";
}
