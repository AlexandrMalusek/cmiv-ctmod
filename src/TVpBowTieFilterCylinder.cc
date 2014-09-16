//______________________________________________________________________________
//
// TVpBowTieFilterCylinder represents an ideal (i.e, non-scattering) bowtie
// filter which compensates for a cylindrical phantom.  The radius and material
// of the cylinder can be specified.
//
// Method: The bowtie filter decreses the weight of photons emitted from the
// source.  Note that the decrease in weight may be quite large.  For a 32 cm
// water phantom and the linear attenuation coefficient of 0.2 1/cm, the
// attenuation factor is exp(-32*0.2) = exp(-6.4) = 1.7e-3.
//______________________________________________________________________________

#include <math.h>
#include "TVpBowTieFilterCylinder.h"

ClassImp(TVpBowTieFilterCylinder)

//______________________________________________________________________________
TVpBowTieFilterCylinder::TVpBowTieFilterCylinder()
{
  // Default constructor.

  fMaterialPtr = 0;
  fRadius = fRadius2 = 0.0;
  fSad = fSad2 = 0.0;
  fMaxThickness = 0.0;
  fBeamAngle = 0.0;
}

//______________________________________________________________________________
TVpBowTieFilterCylinder::TVpBowTieFilterCylinder(Double_t radius, Double_t sad,
						 Double_t beamAngle,
						 TVpMaterial *materialPtr)
{
  // Create the bow-tie filter object.

  fRadius = radius;
  fSad = sad;
  fBeamAngle = beamAngle;
  fMaterialPtr = materialPtr;
  fRadius2 = fRadius * fRadius;
  fSad2 = fSad * fSad;
  fMaxThickness = 2.0 * fRadius / cos(beamAngle);
}

//______________________________________________________________________________
TVpBowTieFilterCylinder::~TVpBowTieFilterCylinder()
{
  // Destructor.  The material is not deleted.
}

//______________________________________________________________________________
Double_t TVpBowTieFilterCylinder::GetThickness(Double_t energy,
					       const TVpVector3 &directionL) const
{
  // Get thickness of the water (by default) filter in a given direction and
  // for a given photon energy.  Currently, the result does not depend on the
  // energy.
  //
  // Input parameters:
  // - energy - photon energy in keV
  // - directionL - direction from source
  //
  // Method:
  // The thickness is defined as the difference between fMaxThickness and the
  // the pathlength through the water cylinder.

  Double_t thickness = fMaxThickness;
  Double_t u2 = directionL.GetY() * directionL.GetY() 
    + directionL.GetZ() * directionL.GetZ(); 
  Double_t D = u2 * fRadius2 - directionL.GetY() * directionL.GetY()* fSad2;
  if (D > 0)
    thickness -= 2.0 * sqrt(D) / u2;
  return thickness;
}

//______________________________________________________________________________
Double_t TVpBowTieFilterCylinder::GetTransmission(Double_t energy,
						  const TVpVector3 &directionL) const
{
  // Calculate the transmission in a given direction and for a given photon
  // energy.
  //
  // Input parameters:
  // - energy - photon energy in keV
  // - directionL - direction from source
  //
  // Method:
  // The transmission is defined as the probability of passing the bow-tie filter
  // without any interaction, i.e. exp(-radiological_path).
  
  return exp(-fMaterialPtr->GetLac(energy) * GetThickness(energy, directionL));
}

//______________________________________________________________________________
void TVpBowTieFilterCylinder::PrintStatus(std::ostream &out) const
{
  // Print the object's status
  
  out << "<TVpBowTieFilterCylinder>\n"
      << "Radius: " << fRadius << '\n'
      << "SAD: " << fSad << '\n'
      << "Beam angle: " << fBeamAngle << '\n'
      <<  "<TVpBowTieFilterCylinder>\n";
}
