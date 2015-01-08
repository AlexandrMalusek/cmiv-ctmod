void x_DrawProfile()
{
  gROOT->SetStyle("Plain");
  gSystem->Load("libRCTmod.so");
  gStyle->SetOptStat(0);
  gROOT->LoadMacro("f_DrawProfile.C");
  f_DrawProfile("ctmod_p.pda", "ctmod_s.pda");
}
