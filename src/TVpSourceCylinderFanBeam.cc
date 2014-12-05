//______________________________________________________________________________
//
// TVpSourceCylinderFanBeam defines a source which emits particles
// isotropically into a fan beam which projects into a rectangle on a
// cylindrical surface.  The source may have a bow-tie filter attached.
//
// See TVpSource for information about the coordinate system.

#include <cmath>
#include <stdio.h>
#include "TVpSourceCylinderFanBeam.h"
#include "misc.h"

ClassImp(TVpSourceCylinderFanBeam)

//______________________________________________________________________________
TVpSourceCylinderFanBeam::TVpSourceCylinderFanBeam
(TVpSpectrum *spectrum, Double_t radius, Double_t slitX, Double_t slitY)
  : TVpSource(spectrum)
{
  // Construct the source and set it to the default position. An example:
  //
  //   TVpSourceCylinderFanBeam = new TVpSourceCylinderFanBeam
  //     (spectrum,
  //      100.0,    // Radius of the cylindrical surface [cm]
  //      2.0,      // Slit size in X (axial) direction on the cylindrical surface
  //      90.0)     // Slit size in Y direction on the cylindrical surface
  
  f4PiR = 4 * M_PI * radius;
  fRadius = radius;
  fSlitWidth = slitX / 2.0;
  fAngle = slitY / (2 * radius);
  fRadius2 = radius * radius;
  fCosAngle = cos(fAngle);
  fBias = GetSolidAngle() / (4 * M_PI);
}

//______________________________________________________________________________
Double_t TVpSourceCylinderFanBeam::GetSolidAngle() const
{
  // Return solid angle
  
  return 4 * fAngle * fSlitWidth / sqrt(fSlitWidth * fSlitWidth + fRadius 
					* fRadius);
}

//______________________________________________________________________________
void TVpSourceCylinderFanBeam::GetParticle(TVpParticle *particle) const
{
  // Set particle's kinematic parameters.
  //
  // Output parameters:
  // - particle - particle whose parameters are to be set
  //
  // Method:
  // The position of the point source in LOC is (0,0,0) and its beam axis is
  // (0,0,-1).  Sample proper direction and energy.  Transform the direction
  // and position into UNI and set the particle's kinematic parameters.

  // Generate random direction. The symetry axis is (0,0,1).
  Double_t alpha = fAngle * (2*getRand()-1.0);
  Double_t x, rnd;
  do
    {
      x = (2.0*getRand()-1.0) * fSlitWidth;
      rnd = getRand();
    }
  while (rnd > j(x));
  
  Double_t nu = 1 / sqrt(x*x + fRadius2);  // = 1/norm
  TVpVector3 dirL(nu*x, nu*fRadius * sin(alpha), -nu*fRadius * cos(alpha));
  particle->fDirection = DirLocToUni(dirL);
  particle->fPosition = fTraVecL2u; // = PosLocToUni(TVpVector3(0,0,0));
  particle->fEnergy = fSpectrumPtr->GetRandEnergy();
  particle->SetSolid(fSolidPtr);
  particle->fWeight = fBias;
  if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
void TVpSourceCylinderFanBeam::GetParticleHeadedToPoint
(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel)
{
  // Set a particle's kinematic parameters so that it originates in the X-ray
  // source and heads towards a point.
  //
  // Input parameters:
  // - point - the end point
  // - channel - energy channel in the X-ray spectrum
  //
  // Output parameters:
  // - particle - particle whose parameters are to be set

  particle->fPosition = fTraVecL2u;  // = PosLocToUni(TVpVector3(0,0,0));
  particle->fDirection = normalize(*pointUni - particle->fPosition);  // in Uni
  fSpectrumPtr->GetEnergyAndWeight(channel, particle->fEnergy, particle->fWeight);
  particle->SetSolid(fSolidPtr);

  // Account for the beam shape and the bowtie filter
  TVpVector3 dirL = DirUniToLoc(particle->fDirection);
  if (! IsInsideBeam(dirL))
    particle->fWeight = 0.0;
  else if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
void TVpSourceCylinderFanBeam::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream

  out << "<TVpSourceCylinderFanBeam>\n"
      << "$Id: TVpSourceCylinderFanBeam.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Slit size X: " << fSlitWidth << " cm\n"
      << "Cone angle: " << fAngle << " rad\n"
      << "Solid angle: " << GetSolidAngle() << " sr\n"
      << "Bias: " << fBias << '\n';
  fSpectrumPtr->PrintStatus(out);
  out << "</TVpSourceCylinderFanBeam>\n";
}

//______________________________________________________________________________
Double_t TVpSourceCylinderFanBeam::j(Double_t x) const
{
  // Return the value of the "j" function (see theory).
  //
  // Input parameters:
  // - x - the x-value (see theory)
  
  Double_t v = fRadius / sqrt(x*x + fRadius2);
  return v*v*v;
}


//______________________________________________________________________________
void TVpSourceCylinderFanBeam::Draw(Int_t option, Int_t nParticles) const
{
  // Draw a filled circle at the centre of the source and a line segment in
  // the direction of the beam axis.

  TVpSource::Draw(option, nParticles);
  if ((option & TVpSource::kDrawShape) != 0)
    {
      // FIX
    }
}
