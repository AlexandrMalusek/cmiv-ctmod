TInputFile inputFile;

void x_DrawGeometry()
{ 
  Double_t angle = 0.0;

  // Read parameters from an input file
  std::ifstream run_inp("input.txt");

  inputFile.Read(run_inp);
  inputFile.Write();

  // Define the rectangle cut
  TVpMatrix3x3 rot = rotationMatrix(TVpVector3(1, 0, 0), angle);
  TVpRectangleCut *rcut = new  TVpRectangleCut
    (rot*TVpVector3(0, 0, 40),   // center
     rot*TVpVector3(1, 0, 0),    // baseX
     rot*TVpVector3(0, 1, 0),    // baseY
     20,                         // sizeX [cm]
     40,                         // sizeY [cm]
     255,                        // dimX
     255);                       // dimY

  TVpSetupTomograph *setupTomographPtr = getTomograph();
  setupTomographPtr->SetPosition(angle);
  setupTomographPtr->fGeometryPtr->SetViewRange(TVpVector3(60, 60, 60));
  setupTomographPtr->fGeometryPtr->Draw();
  setupTomographPtr->fSource1Ptr->Draw();
  setupTomographPtr->fDetectorPtr->Draw
    (TVpPointDetectorArray::kDrawElementAxes |
     TVpPointDetectorArray::kDrawElementPoints |
     TVpPointDetectorArray::kDrawPdaAxes);
  rcut->Draw();
  rcut->DrawSourceParticle(setupTomographPtr->fSource1Ptr, 300);
}
