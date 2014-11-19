void f_DrawProfile(Char_t *pda_p, Char_t *pda_s)
{
  // Analytical projection
  TVpPointDetectorArrayCylindrical *pdapPtr = 
    new TVpPointDetectorArrayCylindrical(pda_p);
  TH1F *hp = pdapPtr->GetYProfile(pdapPtr->GetNx()/2);
  TCanvas *c1 = new TCanvas("c1", "primary");
  c1->cd();
  c1->SetLogy();
  hp->SetTitle("Primary");
  hp->Draw("hist");

  // MC projection
  TVpPointDetectorArrayCylindrical *pdasPtr = 
    new TVpPointDetectorArrayCylindrical(pda_s);
  TH1F *hs = pdasPtr->GetYProfile(pdasPtr->GetNx()/2);
  TCanvas *c2 = new TCanvas("c2", "scatter");
  c2->cd();
  c2->SetLogy();
  hs->SetTitle("Scatter");
  hs->Draw("E");

  // Scatter to primary ratio
  TCanvas *c3 = new TCanvas("c3", "scatter to primary");
  c3->cd();
  c3->SetLogy(1);
  TH1F *hsp = new TH1F((*hs)/(*hp));
  hsp->GetYaxis()->SetTitle("R");
  hsp->SetTitle("Scatter to primary");
  hsp->Draw("hist");
}
