#include <iostream>
#include <TStopwatch.h>
#include "TreeConverter/CTreeConverterHADES.h"

void RunConvertationHADES (TString inputFileName = "../Source", TString outputFileName = "../Converted") {
//	inputFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_data/treeMaker/output/Feb_13_13_34/AuAu_1_23AGev_gen9_108.list/tree_12108160806_1";
//    outputFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_conv/12108160806_1";
//	inputFileName = "/lustre/nyx/cbm/users/ogolosov/HADES_data/treeMaker/output/Feb_19_23_17/AuAu_1_23AGev_gen9_108.list/tree_12108160806";
//    outputFileName = "/lustre/nyx/cbm/users/ogolosov/test/12108160806";

	Int_t dEdxSource = 3;
	Int_t centMethod = 1;
	TString source, method;
	if (dEdxSource == 2) source = "MTPC";
	if (dEdxSource == 3) source = "full";
	if (centMethod == 1) method = "mh";
	if (centMethod == 2) method = "E";

	CTreeConverterHADES treeConverter;
	treeConverter.SetInputFileName (inputFileName);
	treeConverter.SetOutputFileName (outputFileName);// + "_" + source + "_" + method);
	treeConverter.SetInputTreeName ("tree");
	treeConverter.SetdEdxSource (dEdxSource);
	treeConverter.SetCentralityMethod (centMethod);  // 1 - multiplicity, 2 - Eveto

	treeConverter.SetSNN (1.23); // GeV
	treeConverter.SetNbinsMh (10);
	treeConverter.SetMhRange (0, 280);
	treeConverter.SetNbinsCent (10);
	treeConverter.SetCentRange (0, 50);


  	TStopwatch timer;
   	timer.Reset();
   	timer.Start();

	treeConverter.ConvertTree ();

	timer.Stop();
	printf("\nReal time: %f\n",timer.RealTime());
	printf("CPU time: %f\n",timer.CpuTime());
	printf("\nConverted!\n");
}

# ifndef __CINT__
int main (int argc, char *argv []) {
    RunConvertationHADES ();
    return 0;
}
# endif
