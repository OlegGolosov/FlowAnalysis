#ifndef DEFINES_H
#define DEFINES_H

#include <map>
#include <TString.h>
#include "MyDataTree/CTrack.h"

typedef void (*voidFunction)(std::map <TString, Float_t> &variables);
typedef Float_t (*floatFunction)(std::map <TString, Float_t> &variables);
typedef Bool_t (*boolFunction)(CTrack track);
#define PI TMath::Pi()
#define MAXNHARMONICS 5

enum Particles {
    kNID = 0,
	kPionPlus,
	kPionMinus,
	kProton,
	kAntiProton,
	kElectron,
	kPositron,
	kVeto1,
	kVeto2,
	kVeto3,
	kVeto4,
	kFW,
	kPSD,
	kNPartTypes
};

const Int_t pdgID [kNPartTypes] = {0, 211, -211, 2212, -2212, 11, -11, 81, 82, 83, 84, 85, 86};

const Float_t particleMass [kNPartTypes] = {-999.0, 139.57018, 139.57018, 938.2720814, 938.2720814, 0.5109989461, 0.5109989461};

#endif // DEFINES_H
