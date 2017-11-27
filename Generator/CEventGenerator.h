#ifndef CEVENTGENERATOR_H
#define CEVENTGENERATOR_H

#include <map>
#include <TString.h>
#include "../MyDataTree/CEvent.h"
#include "../config.h"

using namespace std;

class CEventGenerator {

public:
	CEventGenerator (Int_t nHarmonics);
	virtual ~CEventGenerator ();
	void SetHarmonicFunction (Int_t n, floatFunction func);
	void SetPsiFunction (voidFunction func);
	void SetEventFunction (voidFunction func);
	void SetTrackFunction (voidFunction func);
	void Init ();
	CEvent* GenerateEvent ();

private:
	Int_t nHarmonics_;
	Bool_t initFlag_;
	map <TString, Float_t> variables_;
	Bool_t* activeHarmonics_;
	Float_t* v_;
	Float_t* psi_;
	CEvent event_;

	floatFunction* HarmonicFunctions;
	voidFunction GetPsi;
	voidFunction GetEventVariables;
	voidFunction GetTrackVariables;

    Float_t rho (Float_t phi);
};


#endif // CEVENTGENERATOR_H
