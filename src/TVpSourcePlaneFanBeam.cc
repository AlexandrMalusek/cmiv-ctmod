//______________________________________________________________________________
//
// TVpSourcePlaneFanBeam defines a source which emits particles isotropically
// into a fan beam which cuts a rectangle on a planar surface that is
// perpendicular to the beam's central axis.
//
// Example:
// TVpSourcePlaneFanBeam *source0 = new TVpSourcePlaneFanBeam
//   (spectrumPtr,   // X-ray spectrum
//    100,           // distance to slit in cm
//    30,            // slit length in cm (x-axis)
//    40);           // slit width in cm (y-axis)
// source0->SetSolid(geometry->GetSolid(0));

#include <cmath>
#include <cstdio>
#include "TVpMatrix3x3.h"
#include "TVpSourcePlaneFanBeam.h"
#include "misc.h"

ClassImp(TVpSourcePlaneFanBeam)

//______________________________________________________________________________
TVpSourcePlaneFanBeam::TVpSourcePlaneFanBeam()
{
  // Default constructor.  Deprecated.

  f4Pi = 4 * M_PI;
  fBias = GetSolidAngle() / (4 * M_PI);
}

//______________________________________________________________________________
TVpSourcePlaneFanBeam::TVpSourcePlaneFanBeam
(TVpSpectrum *spectrumPtr, Double_t distToSlit, Double_t slitWidthX,
 Double_t slitLengthY)
{
  // Constructor.  The beam is defined by a rectangular opening (a slit) at a
  // given distance.  The slit is centered at the beam axis and is
  // perpendicular to it.
  //
  // Input:
  // - spectrumPtr - energy spectrum of the source
  // - distToSlit - distance to the slit in cm
  // - slitWidthX - full width of the slit in cm (x-axis)
  // - slitLengthY - full length of the slit in cm (y-axis)

  f4Pi = 4 * M_PI;
  fSpectrumPtr = spectrumPtr;
  fDistToSlit = distToSlit;
  fSlitLength = slitLengthY / 2.0;
  fSlitWidth = slitWidthX / 2.0;
  fDistToSlit2 = fDistToSlit * fDistToSlit;
  fBias = GetSolidAngle() / (4 * M_PI);
}

//______________________________________________________________________________
TVpSourcePlaneFanBeam::~TVpSourcePlaneFanBeam()
{
  // Destructor

}

//______________________________________________________________________________
void TVpSourcePlaneFanBeam::GetParticle(TVpParticle *particle) const
{
  // Setup particle's parameters as if it was generated in the source.
  //
  // Input:
  // particle - pointer to an already existing particle

  // Generate random direction. The symetry axis is (0,0,1). 
  Double_t x, y, rnd;
  do
    {
      x = (2.0*getRand()-1.0) * fSlitWidth;
      y = (2.0*getRand()-1.0) * fSlitLength;
      rnd = getRand();
    }
  while (rnd > j(x,y));

  Double_t nu = 1 / sqrt(x*x + y*y + fDistToSlit2);
  TVpVector3 dirL(nu*x, nu*y, -nu*fDistToSlit);
  particle->fDirection =  DirLocToUni(dirL);
  particle->fPosition = fTraVecL2u; // = PosLocToUni(TVpVector3(0,0,0));
  particle->fEnergy = fSpectrumPtr->GetRandEnergy();
  particle->SetSolid(fSolidPtr);
  particle->fWeight = fBias;
  if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
void TVpSourcePlaneFanBeam::GetParticleHeadedToPoint
(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel)
{
  // Set particle's parameters as if it was generated in the source, had
  // energy given by the energy channel of the spectrum and was directed
  // towards a given point.
  //
  // Input:
  // - particle - pointer to an already existing particle
  // - point - point towards which the particle is heading
  // - channel - energy channel of the spectrum

  particle->fPosition = fTraVecL2u;  // = PosLocToUni(TVpVector3(0,0,0));
  particle->fDirection = normalize(*pointUni - particle->fPosition);
  fSpectrumPtr->GetEnergyAndWeight(channel, particle->fEnergy, particle->fWeight);
  particle->SetSolid(fSolidPtr);

  // Don't waste time with zero weight particles
  if (particle->fWeight == 0.0)
    return;
  
  // Account for the beam shape and the bowtie filter
  TVpVector3 dirL = DirUniToLoc(particle->fDirection);
  if (! IsInsideBeam(dirL))
    particle->fWeight = 0.0;
  else if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam::GetSolidAngle() const
{
  // Return solid angle.

  return 4.0 * atan(fSlitLength * fSlitWidth /
		    (fDistToSlit * sqrt(fDistToSlit2 + fSlitLength*fSlitLength
					+ fSlitWidth*fSlitWidth)));
}


//______________________________________________________________________________
void TVpSourcePlaneFanBeam::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input:
  // - out - output stream

  out << "<TVpSourcePlaneFanBeam>\n"
      << "$Id: TVpSourcePlaneFanBeam.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Slit size X: " << fSlitWidth << '\n'
      << "Solid angle: " << GetSolidAngle() << " sr" << '\n'
      << "Bias: " << fBias << '\n';
  fSpectrumPtr->PrintStatus(out);
  out << "</TVpSourcePlaneFanBeam>\n";
}

//______________________________________________________________________________
Double_t TVpSourcePlaneFanBeam::j(Double_t x, Double_t y) const
{
  // The j(x,y) function

  Double_t v1 = fDistToSlit / sqrt(x*x + y*y + fDistToSlit2);
  return v1*v1*v1;
}


//______________________________________________________________________________
void TVpSourcePlaneFanBeam::Draw(Int_t option, Int_t nParticles) const
{
  // Draw the source symbol: a magenta full circle with a short line segment
  // in the direction of the beam's axis.
  
  TVpSource::Draw(option, nParticles);
  if (option & TVpSource::kDrawShape != 0)
    {
      // FIX
    }
}
