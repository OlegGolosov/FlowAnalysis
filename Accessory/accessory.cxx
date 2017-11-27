#ifndef ACCESSORY_CXX
#define ACCESSORY_CXX

#include "accessory.h"

Bool_t Acc::IsInfNaN (Double_t x) {
    if (x < 0.0 || x > 0.0) return 0;
    else return 1;
}


Bool_t Acc::CheckHistogram (TH1 *h) {
    Bool_t flag = 0;
    Int_t nBins = h -> GetNbinsX ();
    for (Int_t i = 1; i <= nBins; i++) {
        if (IsInfNaN (h -> GetBinContent (i)) || IsInfNaN (h -> GetBinError (i))) {
            flag = 1;
            h -> SetBinContent (i, 1.0);
            h -> SetBinError (i, 0.0);
        }
    }
    return flag;
}

#endif // ACCESSORY_CXX
