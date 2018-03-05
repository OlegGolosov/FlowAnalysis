#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>

static const int nSets = 3;
const int step = 1;
const int set = 1;

void getHist () {
    TString setNames [nSets] = {"protons", "weight", "no_weight"};
    TChain *ch = new TChain ("treeQ");
    TString path = "../../HADES_corr_no_weight/";
    ch -> Draw (Form ("acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) : cent >> deltaPsi2bins (10, 0, 50, 2, 0, %f)", pi), "", "colz");
    TString filename = "*FWRS" + "_Q_1.root";
    ch -> Add (path + filename);
    Float_t pi = TMath::Pi();
//    TCanvas *c1 = new TCanvas ();
//    ch -> Draw ("TMath::Abs ((atan2 (Ya [0], Xa[0])) - (atan2 (Yb[0], Xb[0]))) >> habs", "", "", 100000);
//    ch -> Draw ("acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) >> hcos", "", "");
//    ch -> Draw ("TMath::Abs (TMath::Abs (atan2 (Ya [0], Xa[0])) - TMath::Abs (atan2 (Yb[0], Xb[0])))", "", "", 100000);
//    c1 -> Print ("subeventMult.C");
//    TCanvas *c2 = new TCanvas ();
//    ch -> Draw (Form ("acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) : cent >> deltaPsi (10, 0, 50, 100, 0, %f)", pi), "", "colz");
//    c2 -> Print ("deltaPsi.C");
    TCanvas *c3 = new TCanvas ();
    ch -> Draw (Form ("acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) : cent >> deltaPsi2bins (10, 0, 50, 2, 0, %f)", pi), "", "colz");
    c3 -> Print (path + "deltaPsi2bins.C");

    TCanvas *c4 = new TCanvas ();
    ch -> Draw ("cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0])) : cent >> Res_OG(10, 0, 50)", "", "profile");
    Res_OG
    c4 -> Print (path + "Res_OG_" + setNames [set] + ".C");
}
