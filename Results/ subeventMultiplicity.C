{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Feb 19 14:21:04 2018) by ROOT version5.34/14
   TCanvas *c1 = new TCanvas("c1", "c1",10,64,700,500);
   c1->Range(-3.75,-3.75,33.75,33.75);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH2F *subeventMult = new TH2F("subeventMult","mha[0]:mhb[0]",30,0,30,30,0,30);
   subeventMult->SetBinContent(33,130423);
   subeventMult->SetBinContent(34,154196);
   subeventMult->SetBinContent(66,328792);
   subeventMult->SetBinContent(67,638424);
   subeventMult->SetBinContent(99,1057689);
   subeventMult->SetBinContent(100,1564179);
   subeventMult->SetBinContent(132,2146197);
   subeventMult->SetBinContent(133,2781765);
   subeventMult->SetBinContent(165,3415082);
   subeventMult->SetBinContent(166,3963460);
   subeventMult->SetBinContent(198,4317154);
   subeventMult->SetBinContent(199,4398510);
   subeventMult->SetBinContent(231,4189508);
   subeventMult->SetBinContent(232,3728898);
   subeventMult->SetBinContent(264,3112156);
   subeventMult->SetBinContent(265,2431519);
   subeventMult->SetBinContent(297,1788128);
   subeventMult->SetBinContent(298,1238304);
   subeventMult->SetBinContent(330,812532);
   subeventMult->SetBinContent(331,505240);
   subeventMult->SetBinContent(363,298083);
   subeventMult->SetBinContent(364,167562);
   subeventMult->SetBinContent(396,90429);
   subeventMult->SetBinContent(397,46513);
   subeventMult->SetBinContent(429,23305);
   subeventMult->SetBinContent(430,11190);
   subeventMult->SetBinContent(462,5160);
   subeventMult->SetBinContent(463,2423);
   subeventMult->SetBinContent(495,1124);
   subeventMult->SetBinContent(496,529);
   subeventMult->SetBinContent(528,244);
   subeventMult->SetBinContent(529,105);
   subeventMult->SetBinContent(561,53);
   subeventMult->SetBinContent(562,26);
   subeventMult->SetBinContent(594,14);
   subeventMult->SetBinContent(595,12);
   subeventMult->SetBinContent(627,2);
   subeventMult->SetBinContent(628,1);
   subeventMult->SetBinContent(660,2);
   subeventMult->SetBinContent(661,1);
   subeventMult->SetBinContent(727,1);
   subeventMult->SetEntries(4.334894e+07);
   subeventMult->SetContour(20);
   subeventMult->SetContourLevel(0,0);
   subeventMult->SetContourLevel(1,219925.5);
   subeventMult->SetContourLevel(2,439851);
   subeventMult->SetContourLevel(3,659776.5);
   subeventMult->SetContourLevel(4,879702);
   subeventMult->SetContourLevel(5,1099628);
   subeventMult->SetContourLevel(6,1319553);
   subeventMult->SetContourLevel(7,1539478);
   subeventMult->SetContourLevel(8,1759404);
   subeventMult->SetContourLevel(9,1979330);
   subeventMult->SetContourLevel(10,2199255);
   subeventMult->SetContourLevel(11,2419180);
   subeventMult->SetContourLevel(12,2639106);
   subeventMult->SetContourLevel(13,2859032);
   subeventMult->SetContourLevel(14,3078957);
   subeventMult->SetContourLevel(15,3298882);
   subeventMult->SetContourLevel(16,3518808);
   subeventMult->SetContourLevel(17,3738734);
   subeventMult->SetContourLevel(18,3958659);
   subeventMult->SetContourLevel(19,4178584);
   
   TPaletteAxis *palette = new TPaletteAxis(30.1875,0,31.875,30,subeventMult);
palette->SetLabelColor(1);
palette->SetLabelFont(42);
palette->SetLabelOffset(0.005);
palette->SetLabelSize(0.035);
palette->SetTitleOffset(1);
palette->SetTitleSize(0.035);
   palette->SetFillColor(100);
   palette->SetFillStyle(1001);
   subeventMult->GetListOfFunctions()->Add(palette,"br");
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.695,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("subeventMult");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries =   4.334894e+07");
   text = ptstats->AddText("Mean x =  5.648");
   text = ptstats->AddText("Mean y =  5.148");
   text = ptstats->AddText("RMS x =  1.981");
   text = ptstats->AddText("RMS y =  1.978");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   subeventMult->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(subeventMult);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   subeventMult->SetLineColor(ci);
   subeventMult->GetXaxis()->SetLabelFont(42);
   subeventMult->GetXaxis()->SetLabelSize(0.035);
   subeventMult->GetXaxis()->SetTitleSize(0.035);
   subeventMult->GetXaxis()->SetTitleFont(42);
   subeventMult->GetYaxis()->SetLabelFont(42);
   subeventMult->GetYaxis()->SetLabelSize(0.035);
   subeventMult->GetYaxis()->SetTitleSize(0.035);
   subeventMult->GetYaxis()->SetTitleFont(42);
   subeventMult->GetZaxis()->SetLabelFont(42);
   subeventMult->GetZaxis()->SetLabelSize(0.035);
   subeventMult->GetZaxis()->SetTitleSize(0.035);
   subeventMult->GetZaxis()->SetTitleFont(42);
   subeventMult->Draw("colz");
   
   TPaveText *pt = new TPaveText(0.3714655,0.9365254,0.6285345,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("mha[0]:mhb[0]");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
