#include <Riostream.h>
#include <TMath.h>
#include <TString.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "ManualFunctions.h"
#include "Generator/CEventGenerator.h"
#include "Generator/CTreeBuilder.h"
#include "Generator/CAcceptance.h"

using namespace std;

void RunGeneration (TString outputFileName = "../new") {

    outputFileName = "../Generated/gen_const"; // test

    TStopwatch timer;
    timer.Reset();
    timer.Start();

	const Int_t nHarmonics = 2;
	const Int_t nEvents = 10000;
	CEvent* event_;

	CEventGenerator generator (nHarmonics);
//	generator.SetHarmonicFunction (1, v1_new);
//	generator.SetHarmonicFunction (2, v2_new);
	generator.SetHarmonicFunction (1, v1);
	generator.SetHarmonicFunction (2, v2);
	generator.SetPsiFunction (GetPsi);
	generator.SetEventFunction (GetEventVariables);
	generator.SetTrackFunction (GetTrackVariables);
	generator.Init ();

	CTreeBuilder treeBuilder;
	treeBuilder.SetOutputFileName (outputFileName + ".root");
	//treeBuilder.SetOutputFileName ("../Generated_test.root"); // test
	treeBuilder.Init ();

	cout << "\nGenerating:\n";

	for (Int_t i = 0; i < nEvents; i++) {
		cout << "\rEvent " << i + 1 << " from " << nEvents;
		event_ = generator.GenerateEvent ();
		treeBuilder.AddEvent (event_);
	}
	treeBuilder.Finish ();

	cout << "\nApplying acceptance:\n";

    acc -> SetParameters (fitParameters);
    norm =  2.0 * fitParameters [0];
	CAcceptance acceptance (AcceptanceFunction4);
	acceptance.SetInputFile (outputFileName + ".root");
	acceptance.SetOutputFile (outputFileName + "_acc.root");
	acceptance.ProcessTree ();

	cout << "\nFinish" << endl;

    timer.Stop();
    printf("Real time: %f\n",timer.RealTime());
    printf("CPU time: %f\n",timer.CpuTime());
}


//# ifndef __CINT__
//int main (int argc, char *argv []) {
//    //TString s(argv[1]);
//    Generate (argv[1]);
//    return 0;
//}
//# endif
