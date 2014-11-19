#include <string>
#include <fstream>
#include "TVpGeometry.h"
#include "TVpSetupTomograph.h"
#include "TVpRunManagerTomograph.h"
#include "ctmod.h"
#include "TInputFile.h"

TInputFile inputFile;

int main()
{
  // Generated files
  std::ofstream logFile("ctmod_status.xml");
  
  inputFile.Read();
  inputFile.Write();

  // Materials
  TVpMaterialManager *materialManagerPtr = getMaterialManager();
  materialManagerPtr->InitializeMaterialData();
  
  TVpSetupTomograph *setupTomographPtr = getTomograph(materialManagerPtr);
  setupTomographPtr->PrintStatus(logFile);

  // Calculate the analytical primary projection
  setupTomographPtr->AnalyticProjection();
  setupTomographPtr->GetPointDetectorArrayPtr()->WritePdaFile("ctmod_p.pda");

  // Calculate the analytical primary projection without the phantom
  setupTomographPtr->GetPointDetectorArrayPtr()->Zero();
  setupTomographPtr->AnalyticProjection(0, 1, 0);
  setupTomographPtr->GetPointDetectorArrayPtr()->WritePdaFile("ctmod_e.pda");

  // MC simulation of the projection
  // Use the alternative source.
  setupTomographPtr->ActivateSource(1);
  setupTomographPtr->GetPointDetectorArrayPtr()->Zero();
  TVpRunManagerTomograph *trmPtr = new TVpRunManagerTomograph(setupTomographPtr);

  trmPtr->Run(inputFile.GetNumOfPhotons());
  trmPtr->fSetupTomographPtr->GetPointDetectorArrayPtr()->
    WritePdaFile("ctmod_s.pda", 1);
  
  return 0;
}
