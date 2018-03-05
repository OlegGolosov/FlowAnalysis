#include <TRandom3.h>
#include <TH1F.h>
#include <vector>

void randomizer (int nTracks = 200, int nEvents = 10000) {
    TRandom2 r (0);

    TH2F *h = new TH2F ("h", "h", nTracks, 0, nTracks, nTracks, 0, nTracks);
    int a, b;
    vector <int> index;
    for (int i = 0; i < nEvents; i++) {
        random_shuffle (myvector.begin (), myvector.end ());
        a = 0;
        b = 0;
        for (int j = 0; j < nTracks; j++) {
            index.push_back (j);
        }
        random_shuffle (myvector.begin (), myvector.end ());
        for (int j = 0; j < nTracks; j++) {
            index.push_back (j);
            if (r.Rndm () <= 0.5) a++;
            else b++;
        }
    h -> Fill (a, b);
    }
    h -> Draw ("colz");
}


