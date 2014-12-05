//______________________________________________________________________________
//
// TVpSourcePlaneFanBeam2 defines a source which emits particles
// isotropically into a fan beam which projects into a rectangle on a
// planar surface.

#include <math.h>
#include <stdio.h>
#include "TVpSourcePlaneFanBeam2.h"
#include "misc.h"

ClassImp(TVpSourcePlaneFanBeam2)

//______________________________________________________________________________
TVpSourcePlaneFanBeam2::TVpSourcePlaneFanBeam2()
{
  // Default constructor.

  f4Pi = 4 * M_PI;
  fBias = GetSolidAngle() / (4 * M_PI);
}


//______________________________________________________________________________
TVpSourcePlaneFanBeam2::TVpSourcePlaneFanBeam2
(TVpSpectrum *spectrumPtr, Double_t distToSlit, Double_t slitWidthX,
 Double_t slitLengthY)
{
  // Constructor

  f4Pi = 4 * M_PI;
  fSpectrumPtr = spectrumPtr;
  fDistToSlit = distToSlit;
  fSlitLength = slitLengthY / 2.0;
  fSlitWidth = slitWidthX / 2.0;
  fDistToSlit2 = fDistToSlit * fDistToSlit;
  fBias = GetSolidAngle() / (4 * M_PI);
}

//______________________________________________________________________________
TVpSourcePlaneFanBeam2::~TVpSourcePlaneFanBeam2()
{
  // Destructor

}

//______________________________________________________________________________
void TVpSourcePlaneFanBeam2::GetParticle(TVpParticle *particle) const
{
  // Generete new particle in the source.


  // Generate random direction. The symetry axis is (0,0,1). 
  Double_t y = InvKx(fSlitWidth, (2*getRand()-1) * Kx(fSlitWidth, fSlitLength));
  Double_t x = InvJy((2*getRand()-1) * Jy(fSlitWidth, y), y);
  Double_t nu = 1 / sqrt(x*x + y*y + fDistToSlit2);
  //cout << '\t' << nu*x << '\t' << nu*y << '\t' << -nu*fDistToSlit << endl;
  TVpVector3 dirL(nu*x, nu*y, -nu*fDistToSlit);
  particle->fDirection = DirLocToUni(dirL);
  particle->fPosition = fTraVecL2u; // = PosLocToUni(TVpVector3(0,0,0));
  particle->fEnergy = fSpectrumPtr->GetRandEnergy();
  particle->SetSolid(fSolidPtr);
  particle->fWeight = fBias;
  if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
void TVpSourcePlaneFanBeam2::GetParticleHeadedToPoint
(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel)
{
  // Generate new particle headed towards the point "point".

  particle->fPosition = fTraVecL2u;  // = PosLocToUni(TVpVector3(0,0,0));
  particle->fDirection = normalize(*pointUni - particle->fPosition);
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
Double_t TVpSourcePlaneFanBeam2::GetSolidAngle() const
{
  // Return solid angle

  return 4.0 * atan(fSlitLength * fSlitWidth /
		    (fDistToSlit * sqrt(fDistToSlit2 + fSlitLength*fSlitLength
					+ fSlitWidth*fSlitWidth)));
}


//______________________________________________________________________________
void TVpSourcePlaneFanBeam2::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpSourcePlaneFanBeam2>\n"
      << "$Id: TVpSourcePlaneFanBeam2.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Slit size X: " << fSlitWidth << '\n'
      << "Solid angle: " << GetSolidAngle() << " sr" << '\n'
      << "Bias: " << fBias << '\n';
  fSpectrumPtr->PrintStatus(out);
  out << "</TVpSourcePlaneFanBeam2>\n";
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam2::Jy(Double_t x, Double_t y) const
{
  // The Jy function

  Double_t v1 = fDistToSlit2 + y*y;
  return fDistToSlit * x /(f4Pi * v1 * sqrt(v1 + x*x));
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam2::InvJy(Double_t x, Double_t y) const
{
  // The inverse Jy function

  Double_t v1 = fDistToSlit2 + y*y;
  Double_t v2 = f4Pi * x * v1;
  Double_t v3 = fDistToSlit / v2;
  Double_t v4 = sqrt(v1 / ( v3*v3 -1));
  return (x < 0) ? -v4 : v4;
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam2::Kx(Double_t x, Double_t y) const
{
  // The Kx function

  return atan(x*y / (fDistToSlit*sqrt(fDistToSlit2 + x*x + y*y))) / f4Pi;
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam2::InvKx(Double_t x, Double_t y) const
{
  // The inverse Kx function

  Double_t v1 = x / (fDistToSlit * tan(f4Pi*y));
  Double_t v2 = sqrt((x*x + fDistToSlit2) / (v1*v1 -1));
  return (y < 0) ? -v2 : v2;
}


 //______________________________________________________________________________
void TVpSourcePlaneFanBeam2::Draw(Int_t option, Int_t nParticles) const
{
  // Draw the source symbol: a magenta full circle with a short line segment
  // in the direction of the beam's axis.
  
  TVpSource::Draw(option, nParticles);
  if ((option & TVpSource::kDrawShape) != 0)
    {
      // FIX
    }
}
