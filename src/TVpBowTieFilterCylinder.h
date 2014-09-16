#ifndef TVpBowTieFilterCylinder_h
#define TVpBowTieFilterCylinder_h

#include "TObject.h"
#include "AVpBowTieFilter.h"
#include "TVpMaterial.h"

class TVpBowTieFilterCylinder : public AVpBowTieFilter
{
 private:
  Double_t      fRadius2;       // fRadius2 = fRadius * fRadius
  Double_t      fSad2;          // fSad2 = fSad * fSad
  Double_t      fMaxThickness;  // Max thickness

 public:
  Double_t      fRadius;        // Cylinder radius in cm
  Double_t      fSad;           // source-axis distance in cm
  Double_t      fBeamAngle;     // Polar angle in the x-direction
  TVpMaterial  *fMaterialPtr;   //! filter material, water by default

  TVpBowTieFilterCylinder();
  TVpBowTieFilterCylinder(Double_t radius, Double_t sad, Double_t beamAngle,
			  TVpMaterial *materialPtr);
  virtual ~TVpBowTieFilterCylinder();

  Double_t         GetThickness(Double_t energy, const TVpVector3 &directionL) const;
  virtual Double_t GetTransmission(Double_t energy, const TVpVector3 &directionL) const;
  virtual void     PrintStatus(std::ostream &out = std::cout) const;

  ClassDef(TVpBowTieFilterCylinder,1) // Bow-tie filter, cylinder phantom
};


#endif  // TVpBowTieFilterCylinder_h
