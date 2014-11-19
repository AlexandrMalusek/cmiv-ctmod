TInputFile inputFile;

void x_DrawGeometry()
{ 
  Double_t angle = 0.0;

  // Read parameters from an input file
  std::ifstream run_inp("input.txt");

  inputFile.Read(run_inp);
  inputFile.Write();

  TVpMaterialManager *materialManagerPtr = getMaterialManager();
  TVpSetupTomograph *setupTomographPtr = getTomograph(materialManagerPtr);
  setupTomographPtr->SetPosition(angle);
  setupTomographPtr->fGeometryPtr->SetViewRange(TVpVector3(70, 70, 70));
  setupTomographPtr->fGeometryPtr->Draw();
  setupTomographPtr->fSource1Ptr->Draw(TVpSource::kDrawAxes);
    // (TVpSource::kDrawAxes | TVpSource::kDrawRandomParticles);
  setupTomographPtr->fDetectorPtr->Draw
    (TVpPointDetectorArray::kDrawElementAxes |
     TVpPointDetectorArray::kDrawElementPoints |
     TVpPointDetectorArray::kDrawPdaAxes);
}
