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

  TVpSetupTomograph *setupTomographPtr = getTomograph();
  setupTomographPtr->PrintStatus(logFile);

  // Calculate the analytical primary projection
  setupTomographPtr->AnalyticProjection();
  setupTomographPtr->GetPointDetectorArrayPtr()->WritePdaFile("ctmod_p.pda");

  // MC simulation of the projection
  setupTomographPtr->ActivateSource(1);  // Use the alternative source.
  setupTomographPtr->GetPointDetectorArrayPtr()->Zero();
  TVpRunManagerTomograph *trmPtr = new TVpRunManagerTomograph(setupTomographPtr);
  trmPtr->Run(inputFile.GetNumOfPhotons());
  trmPtr->fSetupTomographPtr->GetPointDetectorArrayPtr()->
    WritePdaFile("ctmod_s.pda", 1);

  return 0;
}
