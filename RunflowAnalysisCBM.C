#include <iostream>
#include <TROOT.h>
#include <TString.h>
#include <TStopwatch.h>
#include "FlowReconstructor/CFlowReconstructor.h"
#include "ManualFunctions.h"

using namespace std;

int RunCBMflowAnalysis (TString option = "all", TString histFileName = "", TString nonUniformInputFileName = "", TString uniformInputFileName = "") {
// options: "correlations", "flow", "all"
// histFileName: route to files with histograms without ".root"
// nonUniformInputFileName: route to file with data without ".root"
// nonUniformInputFileName: route to file with simulations without ".root"

//    option = "all";
//    option = "flow";
//    option = "correlations";

    nonUniformInputFileName = "/mnt/pool/rhic/2/ovgol/CBM_conv/CBM";
    histFileName = "/mnt/pool/rhic/2/ovgol/CBM_flow/CBM";
//    nonUniformInputFileName = "Data/1";
//    histFileName = "Data/1";
//    histFileName = "C:/Users/Administrator/Desktop/Flow_new/RS_E_SS/NEW/NA49_full_E";



	CFlowReconstructor flowReconstructor;
	flowReconstructor.SetNonUniformInputFileName (nonUniformInputFileName);
//	flowReconstructor.SetUniformInputFileName (uniformInputFileName);
	flowReconstructor.SetHistFileName (histFileName);
	flowReconstructor.UseZeroSubevents (0); // do not touch
//	flowReconstructor.SetSamplingMethod (kBootStrapping);
//	flowReconstructor.SetSamplingMethod (kSubsampling);
	flowReconstructor.SetSamplingMethod (kNoSampling);
	flowReconstructor.PropagateResolutionSign ();
	flowReconstructor.CalculateRP ();
//    flowReconstructor.CalculateEP ();
    flowReconstructor.SetFirstStep (kNoCorrections);
    flowReconstructor.SetLastStep (kTwistAndRescale);
//	flowReconstructor.SetComment (comment); // not implemented

	flowReconstructor.AddHarmonic (1);
//	flowReconstructor.SetVariable("eta");
	flowReconstructor.SetVariable("y");
	flowReconstructor.SetNbinsBS (10); // number of samples

//	flowReconstructor.SetReferenceOption (1, "x");
//	flowReconstructor.SetReferenceOption (2, "y");
//    flowReconstructor.SetHarmonicFunction (1, v1); // gen_const
//    flowReconstructor.SetHarmonicFunction (2, v2); // gen_const

	flowReconstructor.SetNrunRange (1508954659, 1508954659);

	flowReconstructor.SetMhRange (10, 350);
	flowReconstructor.SetNbinsMh (10);
	flowReconstructor.SetMhRangeForFlow (10, 350);

	flowReconstructor.SetPtRange (0.0, 2.5);
	flowReconstructor.SetNbinsPt (10);
	flowReconstructor.SetEtaRange (-1.0, 3.0);
	flowReconstructor.SetNbinsEta (8);
	flowReconstructor.SetNbinsEtaRefl (0);

// THREE SUBEVENT WITH PSD
	flowReconstructor.SetResolutionMethod (kThreeSubevents);
	flowReconstructor.SetCentRange (0, 100);
	flowReconstructor.SetNbinsCent (20);
	flowReconstructor.SetCentRangeForFlow (10, 40);
	flowReconstructor.SetPtSubeventsLimits (1, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0);
	flowReconstructor.SetEtaSubeventsLimits (1, 0.0, 1.5, 1.5, 2.5, 2.5, 3.5);
//    flowReconstructor.AddResolutionParticle (1, 1, kPionMinus);
//    flowReconstructor.AddResolutionParticle (1, 1, kPionPlus);
//    flowReconstructor.AddResolutionParticle (1, 1, kProton);
//    flowReconstructor.AddResolutionParticle (1, 1, kAntiProton);
//    flowReconstructor.AddResolutionParticle (1, 2, kPionMinus);
//    flowReconstructor.AddResolutionParticle (1, 2, kPionPlus);
//    flowReconstructor.AddResolutionParticle (1, 2, kProton);
//    flowReconstructor.AddResolutionParticle (1, 2, kAntiProton);
    flowReconstructor.AddResolutionParticle (1, 1, kPSD);
    flowReconstructor.AddResolutionParticle (1, 2, kPSD);
    flowReconstructor.AddResolutionParticle (1, 3, kPSD);

	flowReconstructor.AddFlowParticle (kPionMinus);
	flowReconstructor.SetPtAveragingRange (1, 0.1, 3.0);
	flowReconstructor.SetEtaAveragingRange (1, 0.0, 1.0);
    flowReconstructor.SetResolutionSigns(1, -1, -1, -1);
// END THREE SUBEVENT WITH PSD

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    if (option == "analysis") {
        flowReconstructor.SetNbinsPt (10);
        flowReconstructor.SetNbinsEta (10);
        flowReconstructor.AnalyzeTree ();
    }
	if (option == "all" || option == "correlations") flowReconstructor.GetCorrelations ();
	if (option == "all" || option == "flow") flowReconstructor.GetFlow ();
//	if (option == "reference" || option == "all" || option == "flow") flowReconstructor.Reference (0.0, 2.5, 1.4, 5.0);

    timer.Stop();
    printf("Real time: %f\n",timer.RealTime());
    printf("CPU time: %f\n",timer.CpuTime());

    cout << "Tada!\n";
    return 0;
}

# ifndef __CINT__
int main (int argc, char *argv []) {
    RunCBMflowAnalysis ();
    return 0;
}
# endif
