#include <iostream>
#include <TStopwatch.h>
#include "TreeConverter/CTreeConverter.h"

void RunConvertation (TString inputFileName = "../Source", TString outputFileName = "../Converted") {
	clock_t begin = clock();
//	inputFileName = "/lustre/nyx/cbm/users/ogolosov/NA49_data/3154";
//	outputFileName = "/lustre/nyx/cbm/users/ogolosov/NA49_conv/3154";
//	inputFileName = "Data/3154";
//	outputFileName = "Converted/3154";

    Int_t dEdxSource = 3;
    Int_t centMethod = 1;
	TString source, method;
	if (dEdxSource == 2) source = "MTPC";
	if (dEdxSource == 3) source = "full";
	if (centMethod == 1) method = "mh";
	if (centMethod == 2) method = "E";

	CTreeConverter treeConverter;
	treeConverter.SetInputFileName (inputFileName);
	treeConverter.SetOutputFileName (outputFileName + "_" + source + "_" + method);
	treeConverter.SetInputTreeName ("data");
	treeConverter.SetdEdxSource (dEdxSource);
	treeConverter.SetCentralityMethod (centMethod);  // 1 - multiplicity, 2 - Eveto

	treeConverter.SetSNN (8.8); // GeV
	treeConverter.SetNbinsMh (510);
	treeConverter.SetMhRange (10, 520);
	treeConverter.SetNbinsCent (20);
	treeConverter.SetCentRange (0.0, 100.0);


    TStopwatch timer;
    timer.Reset();
    timer.Start();

	treeConverter.ConvertTree ();

    timer.Stop();
    printf("/nReal time: %f\n",timer.RealTime());
    printf("CPU time: %f\n",timer.CpuTime());
    printf("\nConverted!\n");
}
