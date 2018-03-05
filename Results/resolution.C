#include <TMath.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TProfile.h>
#include <iostream>

#define FACTOR 0.797884561

using namespace std;

Double_t SphericalBesselI (Int_t order, Double_t arg) {
  // compute half-integer modified Bessel functions
  // order >0, for order>5/2, interpolation is used
  // (improve this by using recurrence!!!)

  if (order<0) return 0;
  if (arg<1e-7) return 0;

  switch(order) {
    case 0:  // 1/2
      return FACTOR*sqrt(arg)*TMath::SinH(arg)/arg;
      break;

    case 1:  // 3/2
      return FACTOR*sqrt(arg)*( -TMath::SinH(arg)/(arg*arg) + TMath::CosH(arg)/arg );
      break;

    case 2:  // 5/2
      return FACTOR*sqrt(arg)*( (3./(arg*arg)+1.)*TMath::SinH(arg)/arg - 3.*TMath::CosH(arg)/(arg*arg) );
      break;

    default:
      break;

  }

  return 0.5*(TMath::BesselI(order,arg)+TMath::BesselI(order+1,arg));   // use average of integer-order Bessel
}

Double_t ComputeResolutionOllitraut(Double_t ratio, Int_t n) {
 // compute Ollitrault correction factor

    Double_t chisq = -2.*TMath::Log(2.*ratio);
    Double_t chi = sqrt(chisq);
    cout << chi << "\t\t";

  Int_t n1 = (n-1)/2;
  Int_t n2 = (n+1)/2;
  Double_t sumBessel;
  if (n%2==1) sumBessel = TMath::BesselI(n1,0.5*chisq) + TMath::BesselI(n2,0.5*chisq);  // integer order
  // half-integer order approximated via interpolation of nearest integer-order modified Bessel functions
  // (Exact solution is not implemented in Root!), checked ok vs. Ollitrault's plot in arXiv:nucl-ex/9711003v2
  //  else sumBessel = 0.5*( TMath::BesselI(n1,0.5*chisq) + 2*TMath::BesselI(n2,0.5*chisq) + TMath::BesselI(n2+1,0.5*chisq) ) ;
  else sumBessel = SphericalBesselI(n1,0.5*chisq) + SphericalBesselI(n2,0.5*chisq);  // exact solution for order 1/2, 3/2, 5/2
  Double_t fac = sqrt(TMath::Pi()) / 2. * chi * exp(-0.5*chisq) * sumBessel;
  return fac;
}

void calculate (TH2F *deltaPsi2bins, TString flowPath) {
//   deltaPsi2bins->RebinX (2);

    int nBins = deltaPsi2bins -> GetNbinsX ();
    float res, chi;
    float xMax = deltaPsi2bins -> GetXaxis () -> GetXmax ();
    float xMin = deltaPsi2bins -> GetXaxis () -> GetXmin ();
    TString setName = deltaPsi2bins -> GetTitle ();

    TH1F *hResolution = new TH1F ("hR_Oli_" + setName, "hR_Oli_" + setName, nBins, xMin, xMax);
    TH1F *hChi = new TH1F ("hChi_" + setName, "hChi_" + setName, nBins, xMin, xMax);
    cout << "Bin\tRatio\t\tChi\t\tResolution\n";
    for (int bin = 1; bin <= nBins; bin++) {
        Double_t ratio = deltaPsi2bins -> GetBinContent (bin, 2) / (deltaPsi2bins -> GetBinContent (bin, 2) + deltaPsi2bins -> GetBinContent (bin, 1));
        cout << bin << "\t" << ratio << "\t";
        chi = sqrt (- 2. * log (2. * ratio));
        hChi -> SetBinContent (bin, chi);
        res = ComputeResolutionOllitraut (ratio, 1);
        cout << res << endl;
        hResolution -> SetBinContent (bin, res);
    }
    hResolution -> SetMarkerStyle (20);
    TCanvas *c1 = new TCanvas ();
    hResolution -> Draw ();
    c1 -> Print (flowPath + "hChi.C");
    hResolution -> Write ();
    TCanvas *c2 = new TCanvas ();
    hChi -> Draw ();
    c2 -> Print (flowPath + "hR_Oli.C");
    hChi -> Write ();
}

void resolution (int set = 0) {
    Float_t pi = TMath::Pi();
    TString setNames [] = {"protons", "weight", "no_weight"};
    TString corrPath = "../../HADES_corr_" + setNames [set] + "/";
    TString flowPath = "../../HADES_flow_" + setNames [set] + "/";
    TString convPath = "../../HADES_conv_weight/";
    TString filename = "*FWRS_Q_1.root";
    TFile *resFile = new TFile (flowPath + "resolution.root", "recreate");

    TChain *ch = new TChain ("treeQ");
    ch -> Add (corrPath + filename);

    TCanvas *c1 = new TCanvas ();
    ch -> Draw ("acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) : cent >> deltaPhi2bins_" + setNames [set] + Form ("(10, 0, 50, 2, 0, %f)", pi), "", "colz");
    TH2F *deltaPhi2bins = (TH2F*) c1 -> GetPrimitive ("deltaPhi2bins_" + setNames [set]);
    deltaPhi2bins -> SetTitle (setNames [set]);
    c1 -> Print (flowPath + "deltaPsi2bins.C");
    deltaPhi2bins -> Write ();
    calculate (deltaPhi2bins, flowPath);

    TCanvas *c2 = new TCanvas ();
    ch -> Draw ("cos (atan2(Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0])) : cent >> hR_OG_" + setNames [set] + "(10, 0, 50)", "", "profile");
    TH1D *hR_OG = ((TProfile*)c2->GetPrimitive("hR_OG_" + setNames [set])) -> ProjectionX ("hR_OG_" + setNames [set]);
    hR_OG -> SetTitle ("hR_OG_" + setNames [set]);

    float QQ, QQerr;
    int nBins = hR_OG -> GetNbinsX ();
    for (int i = 1; i <= nBins; i++) {
        QQ = hR_OG -> GetBinContent (i);
        QQerr = hR_OG -> GetBinError (i);
        QQerr = QQerr / sqrt (QQ) * 0.5;
        QQ = sqrt (QQ);
        hR_OG -> SetBinContent (i, QQ);
        hR_OG -> SetBinError (i, QQerr);
    }
    hR_OG -> Draw ();
    c2 -> Print (flowPath + "hR_OG.C");
    hR_OG -> Write ();

    TCanvas *c3;
    TH1D *hR_AS;
    if (set == 2)
    {
    c3 = new TCanvas ();
    delete ch;
    ch = new TChain ("Tree");
    ch -> Add (convPath + "*_.root");
    ch -> Draw ("cos (psiA - psiB) : cent >> hR_AS_" + setNames [set] + "(10, 0, 50)", "", "profile");
    hR_AS = ((TProfile*)c3->GetPrimitive("hR_AS_" + setNames [set])) -> ProjectionX ("hR_AS_" + setNames [set]);
    hR_AS -> SetTitle ("hR_AS_" + setNames [set]);

    nBins = hR_AS -> GetNbinsX ();
    for (int i = 1; i <= nBins; i++) {
        QQ = hR_AS -> GetBinContent (i);
        QQerr = hR_AS -> GetBinError (i);
        QQerr = QQerr / sqrt (QQ) * 0.5;
        QQ = sqrt (QQ);
        hR_AS -> SetBinContent (i, QQ);
        hR_AS -> SetBinError (i, QQerr);
    }
    hR_AS -> Draw ();
    c3 -> Print (flowPath + "hR_AS.C");
    hR_AS -> Write ();
    }
    resFile -> Close ();
}
