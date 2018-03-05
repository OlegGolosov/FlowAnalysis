#include <TFile.h>
#include <TString.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void check (TString ext = "", TString path = "") {
    TString command;
    command = "ls *" + ext + "_corr.root>|filelist.txt";
    path = "../HADES_corr_protons/";
    gSystem -> cd (path);
    gSystem -> Exec (command);
    gSystem -> Exec ("rm -f corrupted.txt");
    gSystem -> Exec ("touch corrupted.txt");
    gSystem -> Exec ("rm -r -f Corrupted");
    gSystem -> Exec ("mkdir Corrupted");
    gSystem -> Exec ("mkdir Corrupted/log");
    ifstream fileList ("filelist.txt");
    ofstream corruptedList;
    corruptedList.open ("corrupted.txt");
    TFile *f;
    TProfile *p;
    string fileName;
    TString folder = "Recentered/Source Histograms/";
    TString objectNames [36] = {"px1X1aCent_SP", "px1X1bCent_SP", "px1X1cCent_SP", "py1Y1aCent_SP", "py1Y1bCent_SP", "py1Y1cCent_SP",
                                "p2x1X1aPtCent_SP", "p2x1X1bPtCent_SP", "p2x1X1cPtCent_SP", "p2y1Y1aPtCent_SP", "p2y1Y1bPtCent_SP", "p2y1Y1cPtCent_SP",
                                "p2x1X1aEtaCent_SP", "p2x1X1bEtaCent_SP", "p2x1X1cEtaCent_SP", "p2y1Y1aEtaCent_SP", "p2y1Y1bEtaCent_SP", "p2y1Y1cEtaCent_SP",
                                "px1X1aCent_EP", "px1X1bCent_EP", "px1X1cCent_EP", "py1Y1aCent_EP", "py1Y1bCent_EP", "py1Y1cCent_EP",
                                "p2x1X1aPtCent_EP", "p2x1X1bPtCent_EP", "p2x1X1cPtCent_EP", "p2y1Y1aPtCent_EP", "p2y1Y1bPtCent_EP", "p2y1Y1cPtCent_EP",
                                "p2x1X1aEtaCent_EP", "p2x1X1bEtaCent_EP", "p2x1X1cEtaCent_EP", "p2y1Y1aEtaCent_EP", "p2y1Y1bEtaCent_EP", "p2y1Y1cEtaCent_EP"};

    cout << endl;
    if (!fileList) cout << "\nFilelist not found!!!\n";
    while (getline(fileList, fileName)) {
        f = new TFile ((Char_t*)fileName.c_str(), "read");
        if (f) {
            for (Int_t i = 0; i < 36; i++) {
                p = (TProfile*) f -> Get (folder + objectNames [i]);
                if (p -> GetMaximum () > 1. || p -> GetMinimum () < -1.) {
                    cout << fileName << " is corrupted!\n";
                    corruptedList << fileName << endl;
                    command = "mv " + fileName + " Corrupted/" + fileName;
                    gSystem -> Exec (command);
                    fileName=fileName.erase (fileName.find ("_corr.root"), 10);
                    command = "mv log/" + fileName + "* Corrupted/log/" + fileName + "*";
                    break;
                }
            }
        }
        delete f;
    }
    corruptedList.close ();
    if (path != "") gSystem -> cd ("-");
}
