#include <TSystem.h>
#include <TROOT.h>
#include <Riostream.h>

using std::cout;
using std::endl;

void rootlogon () {
    TString location = gSystem -> WorkingDirectory ();
	location = "/home/basov/ovgol/FlowAnalysis"; // explicit
    location += "/";
    std::cout << location << std::endl;

	TString debugString="+g";

	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsLog.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsEventClassVariable.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsEventClassVariablesSet.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutsBase.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutAbove.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutBelow.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutOutside.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutSetBit.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutsSet.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutValue.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCutWithin.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsQnVector.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsHistogramBase.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsHistogram.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsHistogramChannelized.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsHistogramChannelizedSparse.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsHistogramSparse.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfile.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfile3DCorrelations.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfileChannelized.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfileChannelizedIngress.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfileComponents.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfileCorrelationComponents.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsProfileCorrelationComponentsHarmonics.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDataVector.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDataVectorChannelized.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsQnVectorBuild.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCorrectionStepBase.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCorrectionsSetOnInputData.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCorrectionsSetOnQvector.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCorrectionOnInputData.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsCorrectionOnQvector.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDetectorConfigurationBase.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDetectorConfigurationsSet.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDetectorConfigurationChannels.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDetectorConfigurationTracks.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsDetector.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsManager.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsInputGainEqualization.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsQnVectorRecentering.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsQnVectorAlignment.cxx" + debugString);
	gROOT -> LoadMacro (location + "QnCorrections/QnCorrectionsQnVectorTwistAndRescale.cxx" + debugString);

	gROOT -> LoadMacro (location + "Accessory/accessory.cxx" + debugString);
	gROOT -> LoadMacro (location + "MyDataTree/CTrack.cxx" + debugString);
	gROOT -> LoadMacro (location + "MyDataTree/CEvent.cxx" + debugString);
	gROOT -> LoadMacro (location + "Generator/CEventGenerator.cxx" + debugString);
	gROOT -> LoadMacro (location + "Generator/CTreeBuilder.cxx" + debugString);
	gROOT -> LoadMacro (location + "Generator/CAcceptance.cxx" + debugString);
	gROOT -> LoadMacro (location + "TreeConverter/CTreeConverter.cxx" + debugString);
	gROOT -> LoadMacro (location + "FlowReconstructor/CFlowReconstructor.cxx" + debugString);

	cout << "Compiled???" << endl;
}
