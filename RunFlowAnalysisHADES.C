#include <iostream>
#include <TROOT.h>
#include <TString.h>
#include <TStopwatch.h>
#include "FlowReconstructor/CFlowReconstructor.h"
#include "ManualFunctions.h"

using namespace std;

int RunFlowAnalysisHADES (TString option = "all", TString histFileName = "", TString nonUniformInputFileName = "", TString uniformInputFileName = "") {
// options: "correlations", "flow", "all"
// histFileName: route to files with histograms without ".root"
// nonUniformInputFileName: route to file with data without ".root"
// nonUniformInputFileName: route to file with simulations without ".root"

//	nonUniformInputFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_conv/12108163524_1_";
//	histFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_corr/12108172955_1_FWRS";
//	histFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_flow_protons/HADES_MDC3S";
	CFlowReconstructor flowReconstructor;
	flowReconstructor.SetNonUniformInputFileName (nonUniformInputFileName);
//	flowReconstructor.SetUniformInputFileName (uniformInputFileName);
	flowReconstructor.SetHistFileName (histFileName);
	flowReconstructor.UseZeroSubevents (0); // do not touch
	flowReconstructor.SetSamplingMethod (kNoSampling); //kSubSampling, kBootStrapping
	flowReconstructor.PropagateResolutionSign ();
    	flowReconstructor.SetFirstStep (kRecentering); // kNoCorrections, kTwistAndRescale
    	flowReconstructor.SetLastStep (kRecentering);

	flowReconstructor.AddHarmonic (1);
//	flowReconstructor.AddHarmonic (2);
//	flowReconstructor.SetVariable("eta"); // peudorapidity
	flowReconstructor.SetVariable ("y"); // rapidity
	flowReconstructor.SetNbinsBS (5); // number of samples
	flowReconstructor.CalculateEP (); // activate Event Plane method

//    flowReconstructor.SetHarmonicFunction (1, v1); // for simulations
//    flowReconstructor.SetHarmonicFunction (2, v2); // for simulations

//	flowReconstructor.SetNrunRange (4200, 4800); // Run number range
//	flowReconstructor.ExcludeRun (3141); // exclude bad runs

	flowReconstructor.SetMhRange (10, 120); // multiplicity range
	flowReconstructor.SetNbinsMh (10); // number of multiplicity bins, currently not used
	flowReconstructor.SetMhRangeForFlow (30, 80); // multiplicity range for flow analysis, applied on "flow" stage, currently not used

	flowReconstructor.SetCentRange (0, 50); // centrality range
	flowReconstructor.SetNbinsCent (10); // number of centrality bins
//	flowReconstructor.SetCentRangeForFlow (0, 20);
	flowReconstructor.SetCentRangeForFlow (20, 30); // centrality range for flow analysis, applied on "flow" stage
//	flowReconstructor.SetCentRangeForFlow (30, 50);

	flowReconstructor.SetPtRange (0.1, 2.0); // pt range for flow analysis
	flowReconstructor.SetNbinsPt (19); // number of pt bins
	flowReconstructor.SetEtaRange (-0.7, 0.9); // (pseudo)rapidity range for flow analysis
	flowReconstructor.SetNbinsEta (16); // number of (pseudo)rapidity bins
	flowReconstructor.SetNbinsEtaRefl (0); // number of reflected (pseudo)rapidity bins on the final plot

//THREE SUBEVENT FW
	flowReconstructor.SetResolutionMethod (kThreeSubevents);
	flowReconstructor.SetPtSubeventsLimits (1, 20.0, 30.0, 20.0, 30.0, 20.0, 30.0); // (harmonic number, ptMinA, ptMaxA, ptMinB, ptMaxB, ptMinC, ptMaxC)
	flowReconstructor.SetEtaSubeventsLimits (1, -0.5, 4.5, 4.5, 6.5, 6.5, 8.5); // (harmonic number, etaMinA, etaMaxA, etaMinB, etaMaxB, etaMinC, etaMaxC)
    flowReconstructor.AddResolutionParticle (1, 1, kFW); // (harmonic number, subevent number, particle to construct subevent - enumerator located in config.h)
    flowReconstructor.AddResolutionParticle (1, 2, kFW);
    flowReconstructor.AddResolutionParticle (1, 3, kFW);
//END THREE SUBEVENT FW


// RANDOM SUBEVENT FW
//    flowReconstructor.SetResolutionMethod (kRandomSubevent);
//    flowReconstructor.SetPtSubeventsLimits (1, 20.0, 30.0, 20.0, 30.0);
//    flowReconstructor.SetEtaSubeventsLimits (1, -0.5, 9.5, -0.5, 9.5);
//    flowReconstructor.AddResolutionParticle (1, 1, kFW);
//    flowReconstructor.AddResolutionParticle (1, 2, kFW);
// END RANDOM SUBEVENT FW

// THREE SUBEVENT MDC
//    flowReconstructor.SetResolutionMethod (kThreeSubevents);
//    flowReconstructor.SetPtSubeventsLimits (1, 0.0, 2.0, 0.0, 2.0, 20, 30);
//    flowReconstructor.SetEtaSubeventsLimits (1, -0.5, -0.3, 0.3, 0.5, -0.5, 9.5);
//    flowReconstructor.AddResolutionParticle (1, 1, kProton);
//    flowReconstructor.AddResolutionParticle (1, 2, kProton);
//    flowReconstructor.AddResolutionParticle (1, 3, kFW);
// END THREE SUBEVENT MDC

// RANDOM SUBEVENT MDC
//    flowReconstructor.SetResolutionMethod (kRandomSubevent);
//    flowReconstructor.SetPtSubeventsLimits (1, 0.0, 2.0, 0.0, 2.0);
//    flowReconstructor.SetEtaSubeventsLimits (1, -0.8, 1.1, -0.8, 1.1);
//    flowReconstructor.AddResolutionParticle (1, 1, kPionMinus);
//    flowReconstructor.AddResolutionParticle (1, 1, kPionPlus);
//    flowReconstructor.AddResolutionParticle (1, 2, kPionMinus);
//    flowReconstructor.AddResolutionParticle (1, 2, kPionPlus);
// END RANDOM SUBEVENT MDC

	flowReconstructor.AddFlowParticle (kProton);
	flowReconstructor.SetPtAveragingRange (1, 0.8, 0.85); // pt integration range for v(cent) and v(rapidity)
	flowReconstructor.SetEtaAveragingRange (1, -0.25, -0.15); // (pseudo)rapidity integration range for v(cent) and v(pt)
    flowReconstructor.SetResolutionSigns (1, 1, 1, 1); // (harmonic number, signA, signB, signC) - use -1 to change flow sign

    TStopwatch timer;
    timer.Reset();
    timer.Start();

	if (option == "all" || option == "correlations") flowReconstructor.GetCorrelations ();
	if (option == "all" || option == "flow") flowReconstructor.GetFlow ();

    timer.Stop();
    printf("Real time: %f\n",timer.RealTime());
    printf("CPU time: %f\n",timer.CpuTime());

    cout << "Tada!\n";
    return 0;
}

# ifndef __CINT__
int main (int argc, char *argv []) {
    RunFlowAnalysisHADES ();
    return 0;
}
# endif
