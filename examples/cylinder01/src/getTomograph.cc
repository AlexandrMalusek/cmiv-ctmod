#include <cmath>
#include "TInputFile.h"

#ifndef __CINT__
#include <iostream>
#include <fstream>
#include "TVpPitchSequential.h"
#include "TVpSetupTomograph.h"
#include "TVpSourceCylinderFanBeam.h"
#include "TVpPointDetectorArrayCylindrical.h"
#include "TVpDimensionsTomograph.h"
#include "ctmod.h"
#endif

extern TInputFile inputFile;

TVpSetupTomograph *getTomograph(TVpMaterialManager *mmPtr)
{
  ////////////////////////////////////////////////////////////////////////
  // The tomograph
  ////////////////////////////////////////////////////////////////////////
  
  TVpDimensionsTomograph *dimensionsPtr = new TVpDimensionsTomograph();
  ifstream in("dimTom.txt");
  dimensionsPtr->Read(in);
  in.close();

  // Point detector array
  TVpDetectorResponse *aerPtr = new TVpDetectorResponse();
  aerPtr->ReadAerFile("detector/YGdOEu_02.aer");

  TVpPointDetectorArrayCylindrical *pdaPtr = new TVpPointDetectorArrayCylindrical
    (dimensionsPtr->GetDaw(),  // PDA size in x-direction (axial) in cm
     dimensionsPtr->GetDal(),  // PDA size in y-direction in cm
     dimensionsPtr->GetSdd(),  // PDA-source distance = radius
     dimensionsPtr->GetNdw(),  // Num. of PD in x-direction
     dimensionsPtr->GetNdl(),  // Num. of PD in y-direction
     1,                        // Histogram's number of channels
     0.0,                      // Histogram's min energy in keV
     1e+30,                    // Histogram's max energy in keV
     aerPtr);                  // Pointer to the detector response class
  pdaPtr->SetScoredQuantity(TVpPointDetector::kDetectorResponse);
  
  // Spectrum and bowtie filter
  TVpSpectrum *spectrumPtr =
    new TVpSpectrum("spectra/w120_15_Al35.spe");
  AVpBowTieFilter *bfPtr = 0;

  // The geometry
  TVpGeometry *geometryPtr = getGeometry(mmPtr);
  
  // Scanning mode
  TVpPitchSequential *pitchPtr = new TVpPitchSequential
    (1.0);   // distance between slices
  
  // Isotropic source emitting particles into a cone - primary projection
  TVpSourceCylinderFanBeam *source0Ptr = new TVpSourceCylinderFanBeam
    (spectrumPtr,               // X-ray spectrum
     dimensionsPtr->GetSdd(),   // source-detector distance = radius
     inputFile.GetBeamWidth(),  // beam width
     dimensionsPtr->GetDal());  // beam length
  source0Ptr->SetSolid(geometryPtr->GetSolid(0));
  source0Ptr->SetBowTieFilter(bfPtr);
  
  // Isotropic source emitting particles into a cone - MC simulation
  Double_t optSlit =  2 * dimensionsPtr->GetSdd() *
    asin(inputFile.GetCylinderRadius() / dimensionsPtr->GetSad());
  TVpSourceCylinderFanBeam *source1Ptr = new TVpSourceCylinderFanBeam
    (spectrumPtr,              // X-ray spectrum
     dimensionsPtr->GetSdd(),  // source-detector distance = radius
     inputFile.GetBeamWidth(), // beam width
     optSlit);                 // beam length
  source1Ptr->SetSolid(geometryPtr->GetSolid(0));
  source1Ptr->SetBowTieFilter(bfPtr);
  
  TVpSetupTomograph *setupTomographPtr =
    new TVpSetupTomograph(source0Ptr, geometryPtr, pdaPtr, pitchPtr, dimensionsPtr);
  setupTomographPtr->SetSource1(source1Ptr);
  setupTomographPtr->SetPosition(0, 0);

  return setupTomographPtr;
}
