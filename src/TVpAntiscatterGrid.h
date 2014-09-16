#ifndef TVpAntiscatterGrid_h
#define TVpAntiscatterGrid_h

#include "TObject.h"
#include "TH2.h"
#include <iosfwd>
#include "TVpMaterial.h"
#include "TVpVector3.h"

class TVpAntiscatterGrid
{
 public:
  Double_t fStripHeight;             // Strip height in cm
  Double_t fStripWidth;              // Strip width in cm
  Double_t fInterspaceWidth;         // Interspace width in cm
  Double_t fUpperCoverWidth;         // Upper cover width in cm
  Double_t fLowerCoverWidth;         // Lower cover width in cm
  TVpMaterial *fMaterialCover;       //! Cover material
  TVpMaterial *fMaterialStrip;       //! Strip material
  TVpMaterial *fMaterialInterspace;  //! Interspace material
  
  TVpAntiscatterGrid();
  TVpAntiscatterGrid(Double_t gridHeight, Double_t stripWidth, Double_t interspaceWidth,
		     Double_t upperCoverWidth, Double_t lowerCoverWidth,
		     TVpMaterial *materialCover, TVpMaterial *materialStrip,
		     TVpMaterial *materialInterspace);
  virtual ~TVpAntiscatterGrid();
  virtual void       PrintStatus(std::ostream &out = std::cout) const;
  Double_t           Transmission(Double_t energy, const TVpVector3 &directionL) const;
  TH2D *GetHistTransmissionThetaPhi(Double_t energy, Int_t numTheta = 64, Int_t numPhi = 64);

  ClassDef(TVpAntiscatterGrid,1) // Attenuation of an antiscatter grid
};

#endif  // TVpAntiscatterGrid_h
