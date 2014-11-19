{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gSystem->Load("libRCTmod.so");
  
  gROOT->LoadMacro("TInputFile.cc");
  gROOT->LoadMacro("getMaterialManager.cc");
  gROOT->LoadMacro("getGeometry.cc");
  gROOT->LoadMacro("getTomograph.cc");

  TInputFile inputFile;
}
