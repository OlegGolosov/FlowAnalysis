{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Tue Feb 20 19:07:27 2018) by ROOT version5.34/14
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",10,10,700,500);
   c1_n2->Range(0,0,1,1);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetFrameBorderMode(0);
   
   TH2F *deltaPsi2bins = new TH2F("deltaPsi2bins","acos (cos (atan2 (Ya [0], Xa[0]) - atan2 (Yb[0], Xb[0]))) : cent ",10,0,50,2,0,3.141593);
   deltaPsi2bins->SetBinContent(13,2770148);
   deltaPsi2bins->SetBinContent(14,3405154);
   deltaPsi2bins->SetBinContent(15,3764644);
   deltaPsi2bins->SetBinContent(16,4120336);
   deltaPsi2bins->SetBinContent(17,4136848);
   deltaPsi2bins->SetBinContent(18,4309073);
   deltaPsi2bins->SetBinContent(19,4066883);
   deltaPsi2bins->SetBinContent(20,3990297);
   deltaPsi2bins->SetBinContent(21,3489538);
   deltaPsi2bins->SetBinContent(22,301141);
   deltaPsi2bins->SetBinContent(25,1962364);
   deltaPsi2bins->SetBinContent(26,1612975);
   deltaPsi2bins->SetBinContent(27,1186324);
   deltaPsi2bins->SetBinContent(28,972623);
   deltaPsi2bins->SetBinContent(29,835922);
   deltaPsi2bins->SetBinContent(30,830154);
   deltaPsi2bins->SetBinContent(31,812183);
   deltaPsi2bins->SetBinContent(32,885917);
   deltaPsi2bins->SetBinContent(33,943864);
   deltaPsi2bins->SetBinContent(34,108134);
   deltaPsi2bins->SetEntries(4.450452e+07);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   deltaPsi2bins->SetLineColor(ci);
   deltaPsi2bins->GetXaxis()->SetLabelFont(42);
   deltaPsi2bins->GetXaxis()->SetLabelSize(0.035);
   deltaPsi2bins->GetXaxis()->SetTitleSize(0.035);
   deltaPsi2bins->GetXaxis()->SetTitleFont(42);
   deltaPsi2bins->GetYaxis()->SetLabelFont(42);
   deltaPsi2bins->GetYaxis()->SetLabelSize(0.035);
   deltaPsi2bins->GetYaxis()->SetTitleSize(0.035);
   deltaPsi2bins->GetYaxis()->SetTitleFont(42);
   deltaPsi2bins->GetZaxis()->SetLabelFont(42);
   deltaPsi2bins->GetZaxis()->SetLabelSize(0.035);
   deltaPsi2bins->GetZaxis()->SetTitleSize(0.035);
   deltaPsi2bins->GetZaxis()->SetTitleFont(42);
   deltaPsi2bins->Draw("colz");
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
}
