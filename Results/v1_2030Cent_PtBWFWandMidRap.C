#include <TInterpreter.h>

#include <iostream>
#include "../style.h"

using namespace std;
// B. Kardan 2017 for QM17
void v1_2030Cent_PtBWFWandMidRap(){

  TString title = "v1_2030Cent_PtBWFWandMidRap_SysError";
  style();
  int n=0;
 
  TCanvas *canvas = new TCanvas(title.Data(), title.Data(),10,10,800,800);
  canvas->cd();

     TString PtFitFunction = "[0]*x+[1]*x^3";
 //      TString PtFitFunction = "[0]*x+[1]*x^2+[2]*x^3";
 //    TString PtFitFunction = "[0]*x+[1]*x^2+[2]*x^3+[3]*sqrt(x)";
//  float maxPt =2000;
  float minPt =0;


  float maxPt =2.000;
 // float maxPtFit =.600;
  float minPtFit =0.800;
  float maxPtFit =1.200;
//  float minPtFit =.800; 

//  hist = new TH1F("v_{2} vs Pt","v_{2} vs Pt", 1, 0., 2100);
  TH1F* hist = new TH1F("v_{2} vs Pt","v_{2} vs Pt", 1, 0., 2.100);
  hist->SetLineColor(0);
  
  TAxis *axis = hist->GetXaxis();
  axis->SetTitle("p_{t} [GeV/c]");
  //axis->CenterTitle(kTRUE);
  axis->SetTitleSize(0.05);
  axis->SetLabelSize(0.045);
  axis->SetTitleOffset(1.3);
  axis->SetNdivisions(10);

  axis = hist->GetYaxis();
  axis->SetTitle("v_{1}");
  //axis->CenterTitle(kTRUE);
  axis->SetTitleSize(0.042);
  axis->SetLabelSize(0.045);
  axis->SetTitleOffset(1.4);
  //axis->SetTitleOffset(0.8);
  axis->SetNdivisions(12);
  
  //  hist->SetTitle("Elliptic Flow");
  hist->SetTitle(""); 
  hist->SetStats(0); 
  hist->SetMaximum(0.25); 
  hist->SetMinimum(-0.55); 
  //hist->SetLabelOffset(-0.01,"X");
  gPad->SetTopMargin(.04);
  gPad->SetBottomMargin(.14);
  gPad->SetLeftMargin(.14);
  gPad->SetRightMargin(.03);

//  canvas->SetLogx(); 
  canvas->SetGridx(); 
  canvas->SetGridy(); 
  hist->Draw(); 

//  TLine *tline=new TLine(0.,0, 2100,0.);
  TLine *tline=new TLine(0.,0, 2.100,0.);
  tline->SetLineWidth(1);
  tline->SetLineStyle(7); //wide dash
  tline->Draw("same");

  //------------------------------------------------------------------  
  TLegend *legend = new TLegend(0.15,0.77,0.48,0.95);  // top left

//  TLegend *legend = new TLegend(0.12,0.17,0.48,0.45);  // bottom left
  legend->SetMargin(0.20);
  legend->SetTextFont(42);
//  legend->SetTextSize(0.035);
  legend->SetBorderSize(0);  //no border for legend
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.028); 

  //------------------------------------------------------------------
  
//    TPaveText *textDefault = new TPaveText(0.6,0.17,0.95,0.35,"NDC NB");  //right bottom
    TPaveText *textDefault = new TPaveText(0.6,0.77,0.95,0.95,"NDC NB");  //right top
//  TPaveText *textDefault = new TPaveText(0.55,0.84,0.95,0.95,"NDC NB");  //right top small
  textDefault->SetFillColor(kWhite);

  //textDefault->SetMargin(0.20);
  textDefault->SetTextFont(42);
  textDefault->SetTextAlign(11);
  textDefault->SetTextSize(0.035);

  textDefault->AddText("#bf{HADES}");
  textDefault->AddText("#bf{Au+Au 1.23 AGeV}");
  textDefault->AddText("#bf{#color[2]{Preliminary}}");  
  textDefault->AddText("#bf{centrality 20-30%}");

  textDefault->Draw();

  //-----------------------------------------------------------------------------------------------------
  //cent vs. EPResolution
  float xxd[4]={5.00,15.00,25.00,35.00};
  float yyd[4]={0.608403,0.815295,0.859224,0.850258};  // final values R1
 //  float yyd[4]={1,1,1,1};

  //-----------------------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------------------
  // v_{1} // systematic error  from max from emb protons

//  TString resMCEPv1pt13cent2Legend = "0.55 < y < 0.65";
  Double_t xxSystem[1]  = { 0.05 };
  Double_t yySystem[1]  = { 0.00 };
  Double_t yyerrSystem[1]  = { 0.02 };
  Double_t xxerrSystem[1]  = {0.04};

  TGraphAsymmErrors* resSystem = new TGraphAsymmErrors(1,xxSystem,yySystem,xxerrSystem,xxerrSystem,yyerrSystem,yyerrSystem);
  resSystem->SetMarkerStyle(MarkerStyle[0]);
  resSystem->SetFillColor(2000+5);
  resSystem->Draw("E2");
  legend->AddEntry(resSystem,"Common Systematic Error","F");
  //-----------------------------------------------------------------------------------------------------


  // from file:
  //   /lustre/nyx/hades/user/bkardan/apr12/flow0816/10cent/day_108_gen8_PID14_v1_MCEP_FEP_FEP2_SPsimpleNoScal_QC_MHsBin_allSectors_10cent_ptCut300_EffData0.98-0.7_PhiRPThetaCentralityPID0

  // v_{1} MCEPcent 20-30 -0.05 < y < 0.05

  TString resMCEPv1pt7cent2Legend = "-0.05 < y < 0.05";
  Double_t xxMCEPv1pt7cent2[35]  = {/* 225.000000 , 275.000000 , 325.000000 ,*/ 375.000000 , 425.000000 , 475.000000 , 525.000000 , 575.000000 , 625.000000 , 675.000000 , 725.000000 , 775.000000 , 825.000000 , 875.000000 , 925.000000 , 975.000000 , 1025.000000 , 1075.000000 , 1125.000000 , 1175.000000 , 1225.000000 , 1275.000000 , 1325.000000 , 1375.000000 , 1425.000000 , 1475.000000 , 1525.000000 , 1575.000000 , 1625.000000 , 1675.000000 , 1725.000000 , 1775.000000 , 1825.000000 , 1875.000000 , 1925.000000 };
  Double_t yyMCEPv1pt7cent2[35]  = {/* -0.020000, -0.015764 , -0.013616 ,*/ -0.013230 , -0.010703 , -0.008203 , -0.007025 , -0.004966 , -0.003843 , -0.002546 , -0.003191 , -0.001178 , -0.000602 , -0.001086 , 0.001124 , 0.000232 , -0.002428 , 0.002801 , -0.000780 , -0.006353 , 0.001338 , -0.002416 , -0.003382 , -0.010856 , -0.007574 , -0.001434 , -0.007460 , -0.006766 , -0.011668 , -0.001566 , -0.025651 , -0.011281 , -0.020988 , -0.035345 , -0.008492 };
  Double_t yyerrMCEPv1pt7cent2[35]  = {/* 0.001440, 0.000658 , 0.000584 ,*/ 0.000559 , 0.000554 , 0.000559 , 0.000571 , 0.000591 , 0.000624 , 0.000666 , 0.000717 , 0.000783 , 0.000866 , 0.000961 , 0.001078 , 0.001220 , 0.001404 , 0.001618 , 0.001869 , 0.002167 , 0.002520 , 0.002940 , 0.003428 , 0.004009 , 0.004629 , 0.005403 , 0.006209 , 0.007246 , 0.008618 , 0.010003 , 0.011693 , 0.013590 , 0.015339 , 0.018538 , 0.027672 };

   n=32;   // skip first 3 points
   for (int i=0;i<n;i++) xxMCEPv1pt7cent2[i] /= 1000.;  // MeV -> GeV
  //  for (int i=0;i<n;i++) yyMCEPv1pt7cent2[i] *= -1;  // Flip Values
   for (int i=0;i<n;i++) yyMCEPv1pt7cent2[i] /= yyd[2];  // correction EPresolution
 //  for (int i=0;i<n;i++) yyerrMCEPv1pt7cent2[i]  += 0.02;   // systematic error  from max from emb protons
   Double_t xxerrMCEPv1pt7cent2[35]  = {0.}; for (int i=0;i<n;i++) xxerrMCEPv1pt7cent2[i]  = 0.;
   TGraphAsymmErrors* resMCEPv1pt7cent2 = new TGraphAsymmErrors(n,xxMCEPv1pt7cent2,yyMCEPv1pt7cent2,xxerrMCEPv1pt7cent2,xxerrMCEPv1pt7cent2,yyerrMCEPv1pt7cent2,yyerrMCEPv1pt7cent2);



   TF1* fit_resMCEPv1pt7cent2 = new TF1("fit_resMCEPv1pt7cent2", PtFitFunction, minPtFit, maxPtFit);
   resMCEPv1pt7cent2->Fit(fit_resMCEPv1pt7cent2,"0","",minPtFit, maxPtFit);
   fit_resMCEPv1pt7cent2 = (TF1*)resMCEPv1pt7cent2->GetListOfFunctions()->At(0)->Clone();
   fit_resMCEPv1pt7cent2->SetRange(minPt, maxPt),
   fit_resMCEPv1pt7cent2->SetLineWidth(2);
   fit_resMCEPv1pt7cent2->SetLineColor(1);
//   fit_resMCEPv1pt7cent2->Draw("same");


   resMCEPv1pt7cent2->SetMarkerStyle(MarkerStyle[0]);
   resMCEPv1pt7cent2->SetMarkerColor(2000+1);
   resMCEPv1pt7cent2->SetLineColor(2000+1);
   resMCEPv1pt7cent2->SetLineWidth(2);

   resMCEPv1pt7cent2->Draw("p");
   legend->AddEntry(resMCEPv1pt7cent2,resMCEPv1pt7cent2Legend.Data(),"Pl");

  //-----------------------------------------------------------------------------------------------------



  //-----------------------------------------------------------------------------------------------------
  // v_{1} MCEPcent 20-30 0.15 < y < 0.25

  TString resMCEPv1pt9cent2Legend = " 0.15 < y < 0.25 (#times -1)";
  Double_t xxMCEPv1pt9cent2[33]  = {/* 275.000000 , 325.000000 ,*/ 375.000000 , 425.000000 , 475.000000 , 525.000000 , 575.000000 , 625.000000 , 675.000000 , 725.000000 , 775.000000 , 825.000000 , 875.000000 , 925.000000 , 975.000000 , 1025.000000 , 1075.000000 , 1125.000000 , 1175.000000 , 1225.000000 , 1275.000000 , 1325.000000 , 1375.000000 , 1425.000000 , 1475.000000 , 1525.000000 , 1575.000000 , 1625.000000 , 1675.000000 , 1725.000000 , 1775.000000 , 1825.000000 , 1875.000000 };
  Double_t yyMCEPv1pt9cent2[33]  = {/* 0.039106, 0.052196 ,*/ 0.067803 , 0.079268 , 0.088850 , 0.099894 , 0.108586 , 0.117950 , 0.123667 , 0.129500 , 0.133732 , 0.138042 , 0.143074 , 0.145475 , 0.151530 , 0.148794 , 0.155685 , 0.154464 , 0.155708 , 0.160968 , 0.163680 , 0.167714 , 0.162766 , 0.172272 , 0.151689 , 0.183461 , 0.151355 , 0.171818 , 0.170472 , 0.166611 , 0.145962 , 0.153220 , 0.144113 };
  Double_t yyerrMCEPv1pt9cent2[33]  = {/* 0.004949, 0.000793 ,*/ 0.000597 , 0.000574 , 0.000574 , 0.000588 , 0.000612 , 0.000647 , 0.000694 , 0.000750 , 0.000817 , 0.000896 , 0.000997 , 0.001121 , 0.001270 , 0.001452 , 0.001662 , 0.001937 , 0.002247 , 0.002609 , 0.003042 , 0.003568 , 0.004209 , 0.004897 , 0.005901 , 0.006838 , 0.008092 , 0.009417 , 0.011004 , 0.013112 , 0.015299 , 0.017232 , 0.033170 };

   n=31;   // skip first 2 points
   for (int i=0;i<n;i++) xxMCEPv1pt9cent2[i] /= 1000.;  // MeV -> GeV
   for (int i=0;i<n;i++) yyMCEPv1pt9cent2[i] *= -1;  // Flip Values
   for (int i=0;i<n;i++) yyMCEPv1pt9cent2[i] /= yyd[2];  // correction EPresolution
 //  for (int i=0;i<n;i++) yyerrMCEPv1pt9cent2[i]  += 0.02;   // systematic error  from max from emb protons
   Double_t xxerrMCEPv1pt9cent2[33]  = {0.}; for (int i=0;i<n;i++) xxerrMCEPv1pt9cent2[i]  = 0.;
   TGraphAsymmErrors* resMCEPv1pt9cent2 = new TGraphAsymmErrors(n,xxMCEPv1pt9cent2,yyMCEPv1pt9cent2,xxerrMCEPv1pt9cent2,xxerrMCEPv1pt9cent2,yyerrMCEPv1pt9cent2,yyerrMCEPv1pt9cent2);


   TF1* fit_resMCEPv1pt9cent2 = new TF1("fit_resMCEPv1pt9cent2", PtFitFunction, minPtFit, maxPtFit);
   resMCEPv1pt9cent2->Fit(fit_resMCEPv1pt9cent2,"0","",minPtFit, maxPtFit);
   fit_resMCEPv1pt9cent2 = (TF1*)resMCEPv1pt9cent2->GetListOfFunctions()->At(0)->Clone();
   fit_resMCEPv1pt9cent2->SetRange(minPt, maxPt),
   fit_resMCEPv1pt9cent2->SetLineWidth(2);
   fit_resMCEPv1pt9cent2->SetLineColor(1);
 //  fit_resMCEPv1pt9cent2->Draw("same");

   resMCEPv1pt9cent2->SetMarkerStyle(MarkerStyle[4]);
   resMCEPv1pt9cent2->SetMarkerColor(2000+2);
   resMCEPv1pt9cent2->SetLineColor(2000+2);
   resMCEPv1pt9cent2->SetLineWidth(2);
   resMCEPv1pt9cent2->Draw("p");


  //-----------------------------------------------------------------------------------------------------

   //-----------------------------------------------------------------------------------------------------
   // v_{1} MCEPcent 20-30 -0.25 < y < -0.15

   TString resMCEPv1pt5cent2Legend = "-0.25 < y < -0.15";
   Double_t xxMCEPv1pt5cent2[37]  = {/* 175.000000 , 225.000000 ,*/ 275.000000 , 325.000000 , 375.000000 , 425.000000 , 475.000000 , 525.000000 , 575.000000 , 625.000000 , 675.000000 , 725.000000 , 775.000000 , 825.000000 , 875.000000 , 925.000000 , 975.000000 , 1025.000000 , 1075.000000 , 1125.000000 , 1175.000000 , 1225.000000 , 1275.000000 , 1325.000000 , 1375.000000 , 1425.000000 , 1475.000000 , 1525.000000 , 1575.000000 , 1625.000000 , 1675.000000 , 1725.000000 , 1775.000000 , 1825.000000 , 1875.000000 , 1925.000000 , 1975.000000 };
   Double_t yyMCEPv1pt5cent2[37]  = {/* -0.059950, -0.065955 ,*/ -0.072304 , -0.079803 , -0.087845 , -0.093261 , -0.101308 , -0.109012 , -0.114637 , -0.120210 , -0.125461 , -0.130932 , -0.132920 , -0.135880 , -0.139478 , -0.141474 , -0.147513 , -0.148622 , -0.150889 , -0.152074 , -0.155630 , -0.161061 , -0.158077 , -0.167572 , -0.169543 , -0.165504 , -0.172074 , -0.167301 , -0.169460 , -0.176439 , -0.179716 , -0.199085 , -0.181666 , -0.198111 , -0.202400 , -0.235728 , -0.199348 };
   Double_t yyerrMCEPv1pt5cent2[37]  = {/* 0.001135, 0.000649 ,*/ 0.000575 , 0.000545 , 0.000530 , 0.000527 , 0.000534 , 0.000550 , 0.000582 , 0.000626 , 0.000675 , 0.000729 , 0.000784 , 0.000856 , 0.000948 , 0.001065 , 0.001199 , 0.001369 , 0.001560 , 0.001809 , 0.002094 , 0.002450 , 0.002856 , 0.003348 , 0.003903 , 0.004606 , 0.005415 , 0.006459 , 0.007631 , 0.008880 , 0.010555 , 0.012219 , 0.014791 , 0.017117 , 0.019458 , 0.024343 , 0.031108 };

    n=35;   // skip first 2 points
    for (int i=0;i<n;i++) xxMCEPv1pt5cent2[i] /= 1000.;  // MeV -> GeV
   //  for (int i=0;i<n;i++) yyMCEPv1pt5cent2[i] *= -1;  // Flip Values
    for (int i=0;i<n;i++) yyMCEPv1pt5cent2[i] /= yyd[2];  // correction EPresolution
 //   for (int i=0;i<n;i++) yyerrMCEPv1pt5cent2[i]  += 0.02;   // systematic error  from max from emb protons
    Double_t xxerrMCEPv1pt5cent2[37]  = {0.}; for (int i=0;i<n;i++) xxerrMCEPv1pt5cent2[i]  = 0.;
    TGraphAsymmErrors* resMCEPv1pt5cent2 = new TGraphAsymmErrors(n,xxMCEPv1pt5cent2,yyMCEPv1pt5cent2,xxerrMCEPv1pt5cent2,xxerrMCEPv1pt5cent2,yyerrMCEPv1pt5cent2,yyerrMCEPv1pt5cent2);
 
    TF1* fit_resMCEPv1pt5cent2 = new TF1("fit_resMCEPv1pt5cent2", PtFitFunction, minPtFit, maxPtFit);
    resMCEPv1pt5cent2->Fit(fit_resMCEPv1pt5cent2,"0","",minPtFit, maxPtFit);
    fit_resMCEPv1pt5cent2 = (TF1*)resMCEPv1pt5cent2->GetListOfFunctions()->At(0)->Clone();
    fit_resMCEPv1pt5cent2->SetRange(minPt, maxPt),
    fit_resMCEPv1pt5cent2->SetLineWidth(2);
    fit_resMCEPv1pt5cent2->SetLineColor(1);
//    fit_resMCEPv1pt5cent2->Draw("same");

    resMCEPv1pt5cent2->SetMarkerStyle(MarkerStyle[0]);
    resMCEPv1pt5cent2->SetMarkerColor(2000+2);
    resMCEPv1pt5cent2->SetLineColor(2000+2);
    resMCEPv1pt5cent2->SetLineWidth(2);
    resMCEPv1pt5cent2->Draw("p");
    legend->AddEntry(resMCEPv1pt5cent2,resMCEPv1pt5cent2Legend.Data(),"Pl");
   legend->AddEntry(resMCEPv1pt9cent2,resMCEPv1pt9cent2Legend.Data(),"Pl");
   //-----------------------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------------------
  // v_{1} MCEPcent 20-30 0.55 < y < 0.65

  TString resMCEPv1pt13cent2Legend = " 0.55 < y < 0.65 (#times -1)";
  Double_t xxMCEPv1pt13cent2[21]  = { /*475.000000 ,*/ 525.000000 , 575.000000 , 625.000000 , 675.000000 , 725.000000 , 775.000000 , 825.000000 , 875.000000 , 925.000000 , 975.000000 , 1025.000000 , 1075.000000 , 1125.000000 , 1175.000000 , 1225.000000 , 1275.000000 , 1325.000000 , 1375.000000 , 1425.000000 , 1475.000000 };
  Double_t yyMCEPv1pt13cent2[21]  = { /*0.237992,*/ 0.278802 , 0.304088 , 0.326275 , 0.342267 , 0.355100 , 0.365803 , 0.369051 , 0.380991 , 0.384617 , 0.390919 , 0.390878 , 0.395467 , 0.397664 , 0.403858 , 0.412303 , 0.401257 , 0.396742 , 0.397958 , 0.401196 , 0.281007 };
  Double_t yyerrMCEPv1pt13cent2[21]  = { /*0.024354,*/ 0.001828 , 0.000924 , 0.000813 , 0.000865 , 0.000958 , 0.001080 , 0.001240 , 0.001431 , 0.001678 , 0.001955 , 0.002306 , 0.002713 , 0.003202 , 0.003795 , 0.004522 , 0.005432 , 0.007048 , 0.009557 , 0.015592 , 0.053041 };

   n=19;   // first and last point removed
   for (int i=0;i<n;i++) xxMCEPv1pt13cent2[i] /= 1000.;  // MeV -> GeV
   for (int i=0;i<n;i++) yyMCEPv1pt13cent2[i] *= -1;  // Flip Values
   for (int i=0;i<n;i++) yyMCEPv1pt13cent2[i] /= yyd[2];  // correction EPresolution
 //  for (int i=0;i<n;i++) yyerrMCEPv1pt13cent2[i]  += 0.02;   // systematic error  from max from emb protons
   Double_t xxerrMCEPv1pt13cent2[21]  = {0.}; for (int i=0;i<n;i++) xxerrMCEPv1pt13cent2[i]  = 0.;
   TGraphAsymmErrors* resMCEPv1pt13cent2 = new TGraphAsymmErrors(n,xxMCEPv1pt13cent2,yyMCEPv1pt13cent2,xxerrMCEPv1pt13cent2,xxerrMCEPv1pt13cent2,yyerrMCEPv1pt13cent2,yyerrMCEPv1pt13cent2);


   TF1* fit_resMCEPv1pt13cent2 = new TF1("fit_resMCEPv1pt13cent2", PtFitFunction, minPtFit, maxPtFit);
   resMCEPv1pt13cent2->Fit(fit_resMCEPv1pt13cent2,"0","",minPtFit, maxPtFit);
   fit_resMCEPv1pt13cent2 = (TF1*)resMCEPv1pt13cent2->GetListOfFunctions()->At(0)->Clone();
   fit_resMCEPv1pt13cent2->SetRange(minPt, maxPt),
   fit_resMCEPv1pt13cent2->SetLineWidth(2);
   fit_resMCEPv1pt13cent2->SetLineColor(1);
//   fit_resMCEPv1pt13cent2->Draw("same");


   resMCEPv1pt13cent2->SetMarkerStyle(MarkerStyle[4]);
   resMCEPv1pt13cent2->SetMarkerColor(2000+3);
   resMCEPv1pt13cent2->SetLineColor(2000+3);
   resMCEPv1pt13cent2->SetLineWidth(2);
   resMCEPv1pt13cent2->Draw("p");


  //-----------------------------------------------------------------------------------------------------



   //-----------------------------------------------------------------------------------------------------
   // v_{1} MCEPcent 20-30 -0.65 < y < -0.55

   TString resMCEPv1pt1cent2Legend = "-0.65 < y < -0.55";
   Double_t xxMCEPv1pt1cent2[36]  = { 225.000000 , 275.000000 , 325.000000 , 375.000000 , 425.000000 , 475.000000 , 525.000000 , 575.000000 , 625.000000 , 675.000000 , 725.000000 , 775.000000 , 825.000000 , 875.000000 , 925.000000 , 975.000000 , 1025.000000 , 1075.000000 , 1125.000000 , 1175.000000 , 1225.000000 , 1275.000000 , 1325.000000 , 1375.000000 , 1425.000000 , 1475.000000 , 1525.000000 , 1575.000000 , 1625.000000 , 1675.000000 , 1725.000000 , 1775.000000 , 1825.000000 , 1875.000000 , 1925.000000 , 1975.000000 };
   Double_t yyMCEPv1pt1cent2[36]  = { -0.149305, -0.174699 , -0.202893 , -0.228460 , -0.252716 , -0.277769 , -0.300009 , -0.319143 , -0.336390 , -0.350275 , -0.362407 , -0.371876 , -0.380465 , -0.386678 , -0.395771 , -0.395560 , -0.405415 , -0.405501 , -0.409030 , -0.407438 , -0.417107 , -0.416087 , -0.414373 , -0.409335 , -0.411932 , -0.429315 , -0.418313 , -0.442867 , -0.437715 , -0.441380 , -0.397447 , -0.418517 , -0.428283 , -0.501960 , -0.360471 , -0.445914 };
   Double_t yyerrMCEPv1pt1cent2[36]  = { 0.003014, 0.001170 , 0.000786 , 0.000634 , 0.000572 , 0.000548 , 0.000549 , 0.000569 , 0.000605 , 0.000659 , 0.000730 , 0.000823 , 0.000936 , 0.001075 , 0.001243 , 0.001462 , 0.001712 , 0.002029 , 0.002432 , 0.002921 , 0.003504 , 0.004157 , 0.005017 , 0.006126 , 0.007406 , 0.008711 , 0.010529 , 0.012295 , 0.015066 , 0.017792 , 0.021481 , 0.028331 , 0.029358 , 0.033630 , 0.040499 , 0.052955 };

    n=36;
    for (int i=0;i<n;i++) xxMCEPv1pt1cent2[i] /= 1000.;  // MeV -> GeV
   //  for (int i=0;i<n;i++) yyMCEPv1pt1cent2[i] *= -1;  // Flip Values
    for (int i=0;i<n;i++) yyMCEPv1pt1cent2[i] /= yyd[2];  // correction EPresolution
 //   for (int i=0;i<n;i++) yyerrMCEPv1pt1cent2[i]  += 0.02;   // systematic error  from max from emb protons
    Double_t xxerrMCEPv1pt1cent2[36]  = {0.}; for (int i=0;i<n;i++) xxerrMCEPv1pt1cent2[i]  = 0.;
    TGraphAsymmErrors* resMCEPv1pt1cent2 = new TGraphAsymmErrors(n,xxMCEPv1pt1cent2,yyMCEPv1pt1cent2,xxerrMCEPv1pt1cent2,xxerrMCEPv1pt1cent2,yyerrMCEPv1pt1cent2,yyerrMCEPv1pt1cent2);


    TF1* fit_resMCEPv1pt1cent2 = new TF1("fit_resMCEPv1pt1cent2", PtFitFunction, minPtFit, maxPtFit);
    resMCEPv1pt1cent2->Fit(fit_resMCEPv1pt1cent2,"0","",minPtFit, maxPtFit);
    fit_resMCEPv1pt1cent2 = (TF1*)resMCEPv1pt1cent2->GetListOfFunctions()->At(0)->Clone();
    fit_resMCEPv1pt1cent2->SetRange(minPt, maxPt),
    fit_resMCEPv1pt1cent2->SetLineWidth(2);
    fit_resMCEPv1pt1cent2->SetLineColor(1);
//    fit_resMCEPv1pt1cent2->Draw("same");


    resMCEPv1pt1cent2->SetMarkerStyle(MarkerStyle[0]);
    resMCEPv1pt1cent2->SetMarkerColor(2000+3);
    resMCEPv1pt1cent2->SetLineColor(2000+3);
    resMCEPv1pt1cent2->SetLineWidth(2);
    resMCEPv1pt1cent2->Draw("p");
    legend->AddEntry(resMCEPv1pt1cent2,resMCEPv1pt1cent2Legend.Data(),"Pl");
   legend->AddEntry(resMCEPv1pt13cent2,resMCEPv1pt13cent2Legend.Data(),"Pl");
   //-----------------------------------------------------------------------------------------------------


  
  legend->SetTextFont(62); // 22 = Times New Roman (bold)
  legend->Draw();



  //-------------------------- print -----------
  canvas->Print(Form("./plots/%s.png",title.Data()));
  canvas->Print(Form("./plots/%s.pdf",title.Data()));
//  canvas->Print(Form("%s.eps",title.Data()));
    
} 

void rest(){
    /*


    */
}
