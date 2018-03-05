/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

#define AliFlowAnalysisWithFittingEventPlane_CXX

// root includes
#include "TFile.h"      
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector2.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"


// aliroot includes
#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowVector.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithFittingEventPlane.h"

// Description: FIXME Template maker to serve as a starting point for flow analysis
// Author:      Redmer Alexander Bertens, Utrecht University, 2013
//              rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 


ClassImp(AliFlowAnalysisWithFittingEventPlane)

AliFlowAnalysisWithFittingEventPlane::AliFlowAnalysisWithFittingEventPlane() :
    fBookOnlyBasicCommonHist    (kFALSE),
    fDebug                      (kFALSE),
    fUsePhiWeights              (kFALSE),
    fApplyCorrectionForNUA      (kFALSE),
    fEPresolutionMethod         (0),
    fHarmonic                   (2),
    fHarmonicRP                 (0),
    fPOItype                    (1),
    fWeightsList                (0x0),
    fHistList                   (0x0),
    fCommonHists                (0x0),
    fCommonHistsRes             (0x0)
{ /* constructor */ }
//_____________________________________________________________________________
AliFlowAnalysisWithFittingEventPlane::~AliFlowAnalysisWithFittingEventPlane() 
{
  // destructor
//    delete fHistList;
  //  delete fHistPhiVsPtEta;


}
//_____________________________________________________________________________
void AliFlowAnalysisWithFittingEventPlane::Init() 
{
    //Define all histograms
    if(fHarmonicRP==0) fHarmonicRP = fHarmonic;
    cout<<"---Analysis with FEP of v"<< fHarmonic << " relativ Psi" << fHarmonicRP <<" by fitting of the Angular Distribution---"<<endl;

    //save old value and prevent histograms from being added to directory
    //to avoid name clashes in case multiple analaysis objects are used
    //in an analysis
    Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    TString sName;
    /*
    Int_t iNbinsCent = AliFlowCommonConstants::GetMaster()->GetNbinsCent();
    Double_t  dCentMin = AliFlowCommonConstants::GetMaster()->GetCentMin();            
    Double_t  dCentMax = AliFlowCommonConstants::GetMaster()->GetCentMax();

    Int_t iNbinsMult = AliFlowCommonConstants::GetMaster()->GetNbinsMult();
    Double_t  dMultMin = AliFlowCommonConstants::GetMaster()->GetMultMin();            
    Double_t  dMultMax = AliFlowCommonConstants::GetMaster()->GetMultMax();
    */
    Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
    Double_t dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
    Double_t dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  
    Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
    Double_t dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
    Double_t dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();  
    
    Int_t iNbinsPhi = AliFlowCommonConstants::GetMaster()->GetNbinsPhi();
    Double_t  dPhiMin = AliFlowCommonConstants::GetMaster()->GetPhiMin();	     
    Double_t  dPhiMax = AliFlowCommonConstants::GetMaster()->GetPhiMax();

    // initialize the histograms
    fHistList = new TList();
    fHistList->SetName("cobjFEP");
    fHistList->SetOwner();

    // common histogram container
    fCommonHists = new AliFlowCommonHist("AliFlowCommonHistFEP","AliFlowCommonHist");
    (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 
    fHistList->Add(fCommonHists);

    // common results container
    fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResultsFEP","",fHarmonic);
    fHistList->Add(fCommonHistsRes);
    
    fHistConfig = new TProfile("Flow_FlagsFEP","Flow_FlagsFEP",6,0.5,6.5,"s");
    fHistConfig->GetXaxis()->SetBinLabel(1,"fApplyCorrectionForNUA");
    fHistConfig->GetXaxis()->SetBinLabel(2,"fNormalizationType");
    fHistConfig->GetXaxis()->SetBinLabel(3,"fUsePhiWeights");
    fHistConfig->GetXaxis()->SetBinLabel(4,"fHarmonic");
    fHistConfig->GetXaxis()->SetBinLabel(5,"fHarmonicRP");
    fHistConfig->GetXaxis()->SetBinLabel(6,"fEPresolutionMethod");
    fHistConfig->Fill(1,fApplyCorrectionForNUA);
    //fHistConfig->Fill(2,fNormalizationType);
    fHistConfig->Fill(3,fUsePhiWeights);
    fHistConfig->Fill(4,fHarmonic);
    fHistConfig->Fill(5,fHarmonicRP);
    fHistConfig->Fill(6,fEPresolutionMethod);
    fHistList->Add(fHistConfig);
    
    
    sName = "PtRapidityVSdeltaPhi";
    fHistPhiVsPtEta = new TH3D(sName.Data(), sName.Data(), iNbinsEta, dEtaMin, dEtaMax,
                                                           iNbinsPt,  dPtMin,  dPtMax,
                                                           iNbinsPhi, dPhiMin, dPhiMax);
    fHistPhiVsPtEta->SetXTitle("Rapidity");
    fHistPhiVsPtEta->SetYTitle("P_{t} [MeV/c]");
    fHistPhiVsPtEta->SetZTitle("#phi-#Psi_{RP} [rad]");
    fHistPhiVsPtEta->Sumw2();
    fHistList->Add(fHistPhiVsPtEta);
    
    fHistQaQb = new TH1D("Flow_QaQb","Flow_QaQb",20000,-100.,100.);
    fHistQaQb->SetYTitle("dN/dQaQb");
    fHistQaQb->SetXTitle("dQaQb");
    fHistQaQb->StatOverflows(kTRUE);
    fHistList->Add(fHistQaQb);
    
    fHistDeltaPhiQaQb = new TH1D("Flow_deltaPhiQaQb","Flow_deltaPhiQaQb",180,0,TMath::Pi());
    fHistDeltaPhiQaQb->SetYTitle("dN/d(#phi_a - #phi_b)");
    fHistDeltaPhiQaQb->SetXTitle("(#phi_a - #phi_b)");
    fHistList->Add(fHistDeltaPhiQaQb);
    
     TH1::AddDirectory(oldHistAddStatus);
}
//_____________________________________________________________________________
void AliFlowAnalysisWithFittingEventPlane::Make(AliFlowEventSimple* anEvent) 
{
  // core method, called for each event
  if (!anEvent) return;
  // test statement
  //printf("Numer of POIs %i", anEvent->NumberOfTracks());
     //            cout << " AliFlowAnalysisWithFittingEventPlane::Make   "<< endl;
  
  
  Double_t dPhi   = 0.;
  Double_t dPt    = 0.;
  Double_t dEta   = 0.;
  Double_t wTrack = 1.;
  Double_t wPhi   = 1.; // phi weight
  Double_t weight = 1.;  //total weight
  Int_t    nBinsPhi = 1.;
  Double_t aRP    =0.;
  
  TH1F *phiWeights     = NULL;

  if(fWeightsList){
      if(fUsePhiWeights){
          phiWeights = dynamic_cast<TH1F *>(fWeightsList->FindObject("phi_weights"));
          if(phiWeights) nBinsPhi = phiWeights->GetNbinsX();
      }
  } // end of if(weightsList)
  
  // Get Q vectors for the subevents
  AliFlowVector* vQarray = new AliFlowVector[2];
  if (fUsePhiWeights)
    anEvent->Get2Qsub(vQarray,fHarmonicRP,fWeightsList,kTRUE);
  else
    anEvent->Get2Qsub(vQarray,fHarmonicRP);
  // Subevent a
  AliFlowVector vQa = vQarray[0];
  // Subevent b
  AliFlowVector vQb = vQarray[1];
  delete [] vQarray;

  // multiplicity here corresponds to the V0 equalized multiplicity
  Double_t dMa = vQa.GetMult();
  if( dMa < 2 ) return;
  Double_t dMb = vQb.GetMult();
  if( dMb < 2 ) return;
  
  if(anEvent->IsSetMCReactionPlaneAngle()){
      // get the MC reaction plane angle
      aRP = anEvent->GetMCReactionPlaneAngle();  
  }
  else{
      aRP = TVector2::Phi_0_2pi((anEvent->GetQ(fHarmonicRP)).Phi());
  }

   fHistQaQb->Fill(vQa*vQb/dMa/dMb);
//   fHistDeltaPhiQaQb->Fill(TMath::Abs(vQa.Phi() - vQb.Phi())); //FIXME
   fHistDeltaPhiQaQb->Fill(TMath::Abs(vQa.DeltaPhi(vQb))); //FIXME


  //fill control histograms     
  if (fUsePhiWeights) {
    fCommonHists->FillControlHistograms(anEvent,fWeightsList,kTRUE);
  } else {
    fCommonHists->FillControlHistograms(anEvent);
  }
  
  //calculate flow
  //loop over the tracks of the event
  Int_t iNumberOfTracks = anEvent->NumberOfTracks(); 
//  Int_t iNumberOfRPs = anEvent->GetEventNSelTracksRP(); 
  for (Int_t i=0;i<iNumberOfTracks;i++) {
      AliFlowTrackSimple* pTrack = anEvent->GetTrack(i); 
      if (pTrack){
        if (pTrack->InPOISelection(fPOItype)) {
               dPhi = TVector2::Phi_0_2pi(pTrack->Phi()-aRP);
               dPt  = pTrack->Pt();
               dEta = pTrack->Eta();
               wTrack   = pTrack->Weight();
               if(fUsePhiWeights && phiWeights && nBinsPhi) // determine phi weight for this particle:
               {
                wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi()))); // for case of different bining
               }
               //the total weight is the product
               weight = wTrack*wPhi; //dWPhi*dWPt*dWEta; 
               fHistPhiVsPtEta->Fill(dEta, dPt, dPhi, weight);
           }	      
       }//track selected
    }//loop over tracks
  
  
}
//_____________________________________________________________________________
void AliFlowAnalysisWithFittingEventPlane::GetOutputHistograms(TList *outputListHistos)
{
    //get pointers to all output histograms (called before Finish())
    fHistList = outputListHistos;
    fCommonHists = (AliFlowCommonHist*) fHistList->FindObject("AliFlowCommonHistFEP");
    fCommonHistsRes = (AliFlowCommonHistResults*) fHistList->FindObject("AliFlowCommonHistResultsFEP");
    fHistConfig = (TProfile*)  fHistList->FindObject("Flow_FlagsFEP");
    fHistPhiVsPtEta = (TH3D*)  fHistList->FindObject("PtRapidityVSdeltaPhi");
    fHistQaQb = (TH1D*)        fHistList->FindObject("Flow_QaQb");
    fHistDeltaPhiQaQb = (TH1D*)fHistList->FindObject("Flow_deltaPhiQaQb");
    
    for(Int_t n=0;n<=8;n++){
             fHistParameterPtEta[n] = (TH2D*)gDirectory->Get(Form("%s_%d",fHistPhiVsPtEta->GetName(),n) );
             if (fHistParameterPtEta[n]) delete fHistParameterPtEta[n];
             fHistParameterEta[n] = (TH1D*)gDirectory->Get(Form("%s_xz_%d",fHistPhiVsPtEta->GetName(),n) );
             if (fHistParameterEta[n]) delete fHistParameterEta[n];
    }

    fHistChi2PtEta       = (TH2D*)gDirectory->Get(Form("%s_chi2",fHistPhiVsPtEta->GetName()) );
    if (fHistChi2PtEta) delete fHistChi2PtEta;
    fHistChi2Eta       = (TH1D*)gDirectory->Get(Form("%s_xz_chi2",fHistPhiVsPtEta->GetName()) );
    if (fHistChi2Eta) delete fHistChi2Eta;
    
}   
//_____________________________________________________________________________
void AliFlowAnalysisWithFittingEventPlane::Finish() {
    // calculate flow and fill the AliFlowCommonHistResults
    
    //*************make histograms etc. 
    //if (fDebug) 
    cout<<"*****************************************************"<<endl;
    cout<<"AliFlowAnalysisWithFittingEventPlane::Finish()"<<endl;
    
 //   Int_t iNbinsPhi = AliFlowCommonConstants::GetMaster()->GetNbinsPhi();
    Int_t iNbinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();  
    Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta(); 
 //   Int_t iNbinsCent = AliFlowCommonConstants::GetMaster()->GetNbinsCent();
 //   Int_t iNbinsMult = AliFlowCommonConstants::GetMaster()->GetNbinsMult();
   
    // access harmonic:
    if(fCommonHists && fCommonHists->GetHarmonic())
    {
     fHarmonic = (Int_t)(fCommonHists->GetHarmonic())->GetBinContent(1); // to be improved (moved somewhere else?)
    } 
    
    if(fHistConfig)
    {
     fHarmonicRP = (Int_t)fHistConfig->GetBinContent(5); // to be improved (moved somewhere else?)
     fEPresolutionMethod = (Int_t)fHistConfig->GetBinContent(6); // to be improved (moved somewhere else?)
    }
    if(fHarmonicRP==0) fHarmonicRP = fHarmonic;
    cout<<"---with FEP of v"<< fHarmonic << " relativ Psi" << fHarmonicRP << "---"<<endl;
    if(fEPresolutionMethod==0)      cout<<"---used Method: Ollitrault, no correction only saving  ---" <<endl;
    else if(fEPresolutionMethod==1) cout<<"---used Method: Ollitrault    ---" <<endl;
    else if(fEPresolutionMethod==2) cout<<"---used Method: reference Flow---" <<endl;
        
    Double_t dV =0;
    Double_t dChi = 0.;
    Double_t dVerr = 0.;
    
    if(fEPresolutionMethod==0 || fEPresolutionMethod==1){   
            //Calculate alternative Method(Ollitrault):
            if(!fHistDeltaPhiQaQb) {
                printf(" PANIC: run with full booking !");
                return;
            }
            Int_t i1 = fHistDeltaPhiQaQb->FindBin(TMath::Pi()/2.);
            Int_t i2 = fHistDeltaPhiQaQb->FindBin(TMath::Pi());

            Double_t ratio =   fHistDeltaPhiQaQb->Integral(i1,i2) 
                             / fHistDeltaPhiQaQb->Integral(1 ,i2);
            dChi = sqrt(-2.*TMath::Log(2.*ratio));


            cout << "Ratio = " << ratio << "   chi = " <<  dChi  << endl;
            dV = ComputeResolutionOllitraut(ratio, fHarmonic);
            printf("An estimate of the event plane resolution with Ollitrault Method is: %f\n", dV );
            printf("\nSummary:\n");
            printf("\nRatio = %f   chi = %f   Res{%d} = %f\n",ratio, dChi ,fHarmonicRP,dV);
    }else if(fEPresolutionMethod==2){
            //Calculate reference flow
            //----------------------------------
            //weighted average over (QaQb/NaNb) with weight (NaNb)
            if(!fHistQaQb) {
                printf(" PANIC: run with full booking !");
                return;
            }
    
            Double_t dEntriesQaQb = fHistQaQb->GetEntries();
            if( dEntriesQaQb < 1 )
              return;
            fHistQaQb->GetXaxis()->SetRangeUser(-10.,10.);
            Double_t dQaQb  = fHistQaQb->GetMean();
            if( dQaQb < 0 )
              return;
            Double_t dSpreadQaQb = fHistQaQb->GetRMS();

            // this is the `resolution'
            if(dQaQb <= .0 ) {
              printf(" Panic! the average of QaQb <= 0! Probably you need to run on more events !\n");
              printf("  \tusing dummy value 1 to avoid segfault \n");
              dQaQb = 1.;
              dSpreadQaQb = 1.;
            }  
            printf("QaQb = %f +- %f\n", dQaQb, (dSpreadQaQb/TMath::Sqrt(dEntriesQaQb)) );
            Double_t dV = TMath::Sqrt(dQaQb);

            Int_t fTotalQvector = 3; //FIXME 
            printf("ResSub (sqrt of scalar product of sub-event qvectors) = %f\n", dV );
            printf("fTotalQvector %d \n",fTotalQvector);
            // we just take the spread as the uncertainty

              if(fTotalQvector>2) {
                dChi = FindXi(dV,1e-6);
                printf("An estimate chi of the event plane is: %f\n", dChi );
        
                dV = ComputeResolution( TMath::Sqrt2()*dChi );
                printf("An estimate of the event plane resolution is: %f\n", dV );
              }
  
            Double_t dStatErrorQaQb = dSpreadQaQb;
            if(dQaQb > 0.) dVerr = (1./(2.*pow(dQaQb,0.5)))*dStatErrorQaQb;
            printf("v{RP:%d}(subevents) = %f +- %f\n",fHarmonicRP,dV,dVerr);
    }
    fCommonHistsRes->FillChi(dChi);
    fCommonHistsRes->FillIntegratedFlow(dV,dVerr);

        //--------------------------------------------------------------------
        //fit-function & clean up:

    TF1  *fit_v1v2;
 	fit_v1v2 = new TF1("fit_v1v2","[0]*(1+2*[1]*TMath::Cos(x)"
                                         "+2*[2]*TMath::Cos(2*x)"
                                         "+2*[3]*TMath::Cos(3*x)"
                                         "+2*[4]*TMath::Cos(4*x)"
                                         "+2*[5]*TMath::Cos(5*x)"
                                         "+2*[6]*TMath::Cos(6*x)"
                                         "+2*[7]*TMath::Sin(x)"
                                         "+2*[8]*TMath::Sin(2*x))" ,0,TMath::TwoPi());

    for(Int_t n=0;n<=8;n++){
          fit_v1v2->SetParameter(n,0);
     //     fHistParameterPtEta[n] = (TH2D*)gDirectory->Get(Form("%s_%d",fHistPhiVsPtEta->GetName(),n) );
     //     if (fHistParameterPtEta[n]) delete fHistParameterPtEta[n];
     //     fHistParameterEta[n] = (TH1D*)gDirectory->Get(Form("%s_%d",fHistPhiVsEta->GetName(),n) );
     //     if (fHistParameterEta[n]) delete fHistParameterEta[n];
    }

    fit_v1v2->SetParameter(0,1000);   // FIXME why??
 //   fHistChi2PtEta       = (TH2D*)gDirectory->Get(Form("%s_chi2",fHistPhiVsPtEta->GetName()) );
 //   if (fHistChi2PtEta) delete fHistChi2PtEta;
 //   fHistChi2Eta       = (TH1D*)gDirectory->Get(Form("%s_chi2",fHistPhiVsEta->GetName()) );
 //   if (fHistChi2Eta) delete fHistChi2Eta;


        //--------------------------------------------------------------------
        //PtEta:
        // Looping over all (pt,eta) bins and calculating differential flow with fit-function: 
        cout<<">   Fitting...."<<  flush;
        fHistList->SetOwner();
        fHistPhiVsPtEta->FitSlicesZ(fit_v1v2, 1, 0, 1, 0, 10);
        cout<<"\r>   Fitting finisched.\n"<<  endl;

        for(Int_t n=0;n<=8;n++){
              fHistParameterPtEta[n] = (TH2D*)gDirectory->Get(Form("%s_%d",fHistPhiVsPtEta->GetName(),n) );
              //if( TMath::Abs(dV!=0.) && (n < 7) && (n > 0 ) ) { fHistParameterPtEta[n]->Scale(1/dV,"width"); }  // global correction for RP resolution  v = v/dV;
              fHistList->Add(fHistParameterPtEta[n]);

        }
        fHistChi2PtEta       = (TH2D*)gDirectory->Get(Form("%s_chi2",fHistPhiVsPtEta->GetName()) );
        fHistList->Add(fHistChi2PtEta);

        // Looping over all (pt,eta) bins and calculating differential flow: 
        if(fHistParameterPtEta[fHarmonic])
        for(Int_t p=1;p<=iNbinsPt;p++)
        {
         cout<<"\r> "<<  (float)p/iNbinsPt*100   <<flush;
         for(Int_t e=1;e<=iNbinsEta;e++)
         {
             //    if(h3_1->GetBinEntries(h3_1->GetBin(p,e))>0)
                 Double_t v = fHistParameterPtEta[fHarmonic]->GetBinContent(e,p);
                 if(fEPresolutionMethod>0 && TMath::Abs(dV!=0.) ) v = v/dV;
                 Double_t vError = fHistParameterPtEta[fHarmonic]->GetBinError(e,p);
                 fCommonHistsRes->FillDifferentialFlowPtEtaPOI(e, p, v, vError);
         }
        }
        // Looping over all (pt,eta) bins and calculating sin-term: 
        if(fHistParameterPtEta[7])
        {
            for(Int_t p=1;p<=iNbinsPt;p++)
            {
             for(Int_t e=1;e<=iNbinsEta;e++)
             {
                 Double_t v = fHistParameterPtEta[7]->GetBinContent(e,p);
                 if(fEPresolutionMethod>0 && TMath::Abs(dV!=0.) ) v = v/dV;
                 Double_t vError = fHistParameterPtEta[7]->GetBinError(e,p);
                 fCommonHistsRes->FillDifferentialSinPtEtaPOI(e, p, v, vError);
             } // end of for(Int_t e=1;e<=fnBinsEta;e++)
            } // end of for(Int_t p=1;p<=fnBinsPt;p++)
        }

        //--------------------------------------------------------------------
        // Looping over all eta-bins and calculating differential flow with fit-function: 
        cout<<">   Fitting Eta Bins...."<<  flush;
        fHistList->SetOwner();
        fHistPhiVsEta = (TH2D*)fHistPhiVsPtEta->Project3D("xz");
        fHistList->Add(fHistPhiVsEta);  //for testing
        fHistPhiVsEta->FitSlicesX(fit_v1v2, 0, -1, 10);
        cout<<"\r>   Fitting Eta Bins finisched.\n"<<  endl;

        for(Int_t n=0;n<=8;n++){
              fHistParameterEta[n] = (TH1D*)gDirectory->Get(Form("%s_%d",fHistPhiVsEta->GetName(),n) );
              //if( TMath::Abs(dV!=0.) && (n < 7) && (n > 0 ) ) { fHistParameterEta[n]->Scale(1/dV,"width"); }  // global correction for RP resolution  v = v/dV;
              fHistList->Add(fHistParameterEta[n]);

        }
        fHistChi2Eta       = (TH1D*)gDirectory->Get(Form("%s_chi2",fHistPhiVsEta->GetName()) );
        fHistList->Add(fHistChi2Eta);

        // Looping over all eta bins and calculating differential flow: 
        if(fHistParameterEta[fHarmonic])
        for(Int_t e=1;e<=iNbinsEta;e++)
        {
             //    if(h3_1->GetBinEntries(h3_1->GetBin(p,e))>0)
                 Double_t v = fHistParameterEta[fHarmonic]->GetBinContent(e);
                 if(fEPresolutionMethod>0 && TMath::Abs(dV!=0.) ) v = v/dV;
                 Double_t vError = fHistParameterEta[fHarmonic]->GetBinError(e);
                 fCommonHistsRes->FillDifferentialFlowEtaPOI(e, v, vError);
        }
    
cout<<"\n*****************************************************"<<endl;
    
}
//_____________________________________________________________________________
void AliFlowAnalysisWithFittingEventPlane::WriteHistograms(TDirectoryFile *outputFileName) const {
    // store the final results in output .root file
    outputFileName->Add(fHistList);
    outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}
//_____________________________________________________________________________

//--------------------------------------------------------------------
Double_t AliFlowAnalysisWithFittingEventPlane::ComputeResolution( Double_t x ) const {
  // Computes resolution for Event Plane method
  if(x > 51.3) {
    printf("Warning: Estimation of total resolution might be WRONG. Please check!");
    return 0.99981;
  }
  Double_t a = x*x/4;
  Double_t b = TMath::Exp(-a)*TMath::BesselI0(a)+TMath::Exp(-a)*TMath::BesselI1(a);
  return TMath::Sqrt(TMath::PiOver2())/2*x*b;
}
//--------------------------------------------------------------------
Double_t AliFlowAnalysisWithFittingEventPlane::FindXi( Double_t res, Double_t prec ) const {
  // Computes x(res) for Event Plane method
  if(res > 0.99981) {
    printf("Warning: Resolution for subEvent is high. You reached the precision limit.");
    return 51.3;
  }
  int nSteps =0;
  Double_t xtmp=0, xmin=0, xmax=51.3, rtmp=0, delta=2*prec;
  while( delta > prec ) {
    xtmp = 0.5*(xmin+xmax);
    rtmp = ComputeResolution(xtmp);
    delta = TMath::Abs( res-rtmp );
    if(rtmp>res) xmax = xtmp;
    if(rtmp<res) xmin = xtmp;
    nSteps++;
  }
  return xtmp;
}
//--------------------------------------------------------------------
#define FACTOR 0.797884561

Double_t AliFlowAnalysisWithFittingEventPlane::SphericalBesselI(Int_t order, Double_t arg) const {  
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
//--------------------------------------------------------------------
Double_t AliFlowAnalysisWithFittingEventPlane::ComputeResolutionOllitraut(Double_t ratio, Int_t n) const {  
 // compute Ollitrault correction factor

    Double_t chisq = -2.*TMath::Log(2.*ratio);
    Double_t chi = sqrt(chisq);

  Int_t n1 = (n-1)/2;
  Int_t n2 = (n+1)/2;
  Double_t sumBessel;
  if (n%2==1) sumBessel = TMath::BesselI(n1,0.5*chisq) + TMath::BesselI(n2,0.5*chisq);  // integer order
  // half-integer order approximated via interpolation of nearest integer-order modified Bessel functions
  // (Exact solution is not implemented in Root!), checked ok vs. Ollitrault's plot in arXiv:nucl-ex/9711003v2 
  //  else sumBessel = 0.5*( TMath::BesselI(n1,0.5*chisq) + 2*TMath::BesselI(n2,0.5*chisq) + TMath::BesselI(n2+1,0.5*chisq) ) ;
  else sumBessel = SphericalBesselI(n1,0.5*chisq) + SphericalBesselI(n2,0.5*chisq);  // exact solution for order 1/2, 3/2, 5/2
  Double_t fac = sqrt(TMath::Pi())/2. * chi * exp(-0.5*chisq) * sumBessel;
  return fac;
}


