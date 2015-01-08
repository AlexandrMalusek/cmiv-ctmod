#include <cmath>
#include "TInputFile.h"

#ifndef __CINT__
#include <iostream>
#include <fstream>
#include "TVpSetupTomograph.h"
#include "TVpSourceCylinderFanBeam.h"
#include "TVpPointDetectorArrayCylindrical.h"
#include "TVpPitchHelical.h"
#include "ctmod.h"
#endif

extern TInputFile inputFile;

TVpSetupTomograph *getTomograph()
{
  ////////////////////////////////////////////////////////////////////////
  // The tomograph
  ////////////////////////////////////////////////////////////////////////
  
  TVpDimensionsTomograph *dimensionsPtr = new TVpDimensionsTomograph();
  ifstream in("dimTom.txt");
  dimensionsPtr->Read(in);
  in.close();
  
  // Detector response function
  TVpDetectorResponse *aerPtr = new TVpDetectorResponse();
  aerPtr->ReadAerFile("detector/YGdOEu_02.aer");
  aerPtr->SetAntiscatterGrid(0);

  TVpPointDetectorArrayCylindrical *pdaPtr = new TVpPointDetectorArrayCylindrical
    (dimensionsPtr->GetDaw(),  // PDA size in x-direction (axial) in cm
     dimensionsPtr->GetDal(),  // PDA size in y-direction in cm
     dimensionsPtr->GetSdd(),  // PDA-source distance = radius
     dimensionsPtr->GetNdw(),  // Num. of PD in x-direction
     dimensionsPtr->GetNdl(),  // Num. of PD in y-direction
     1,              // Number of channels
     0.0,            // Min registered-photon energy in keV
     1e+30,          // Max registered-photon energy in keV
     aerPtr);        // Pointer to the detector response class
  pdaPtr->SetScoredQuantity(TVpPointDetector::kDetectorResponse);

  // Spectrum
  TVpSpectrum *spectrumPtr;
  Double_t e_eps = 1e-6;  // we deal with floating point numbers here
  if (std::abs(inputFile.GetEnergy() - 120) < e_eps)
    spectrumPtr = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  else if (std::abs(inputFile.GetEnergy() - 140) < e_eps)
    spectrumPtr = new TVpSpectrum("spectra/w140_16_Cu05.spe");
  else
    {
      std::cerr << "Error: getTomograph: only 120 and 140 kV allowed, using 140 kV\n";
      spectrumPtr = new TVpSpectrum("spectra/w140_16_Cu05.spe");
    }
  AVpBowTieFilter *bfPtr = 0;
  
  // The geometry
  TVpGeometry *geometryPtr = getGeometry();

  // The tomograph pitch helical
  TVpPitchHelical *pitchPtr = new TVpPitchHelical(1.0);

  // Isotropic source emitting particles in a cone - primary projection
  TVpSourceCylinderFanBeam *source0Ptr = new TVpSourceCylinderFanBeam
    (spectrumPtr,               // X-ray spectrum
     dimensionsPtr->GetSdd(),   // source-detector distance = radius
     inputFile.GetBeamWidth(),  // beam width
     dimensionsPtr->GetDal());  // beam length
  source0Ptr->SetSolid(geometryPtr->GetSolid(0));
  source0Ptr->SetBowTieFilter(bfPtr);

  // Isotropic source emitting particles in a cone - MC simulation
  TVpSourceCylinderFanBeam *source1Ptr = new TVpSourceCylinderFanBeam
    (spectrumPtr,              // X-ray spectrum
     dimensionsPtr->GetSdd(),  // source-detector distance = radius
     inputFile.GetBeamWidth(), // beam width
     30.0);                    // beam length
  source1Ptr->SetSolid(geometryPtr->GetSolid(0));
  source1Ptr->SetBowTieFilter(bfPtr);
  
  TVpSetupTomograph *setupTomographPtr =
    new TVpSetupTomograph(source0Ptr, geometryPtr, pdaPtr, pitchPtr, dimensionsPtr);
  setupTomographPtr->SetSource1(source1Ptr);
  setupTomographPtr->SetPosition(0.0);

  return setupTomographPtr;
}
