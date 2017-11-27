#ifndef ACCESSORY_H
#define ACCESSORY_H

#include <TH1.h>

namespace Acc {
    Bool_t IsInfNaN (Double_t x);
    Bool_t CheckHistogram (TH1 *h);
}

#endif // ACCESSORY_H
