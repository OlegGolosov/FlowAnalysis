void resolution2 (int nEvents = 500000) {
    float QQ, QQerr;
    int nBins;
    TH1D *resolution;
    TChain *ch = new TChain ("treeQ");
    TString filename = "../../HADES_corr_no_weight/*FWRS_Q_1.root";
    ch->Add (filename);
    TCanvas *c1 = new TCanvas ();
    ch -> Draw ("cos (atan2(Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0])) : cent >> R_OG(10, 0, 50)", "", "profile", nEvents);
    resolution = ((TProfile*)c1->GetPrimitive("R_OG")) -> ProjectionX ();
    nBins = resolution -> GetNbinsX ();
    for (int i = 1; i <= nBins; i++) {
        QQ = resolution -> GetBinContent (i);
        QQerr = resolution -> GetBinError (i);
        QQerr = QQerr / sqrt (QQ) * 0.5;
        QQ = sqrt (QQ);
        resolution -> SetBinContent (i, QQ);
        resolution -> SetBinError (i, QQerr);
    }
    resolution -> Draw ();
}
