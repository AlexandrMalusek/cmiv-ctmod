//______________________________________________________________________________
//
// TVpSourceIso defines a point isotropic source which emits particles into a
// cone. It produces a circular beam.
//
// Example:
// TVpSpectrum *spectrumPtr = new TVpSpectrum("spectra/w120_16_Cu05.spe");
// TVpSourceIso *source0Ptr = new TVpSourceIso
//   (spectrumPtr,      // X-ray spectrum
//    6 * d2r);         // 6 deg cone angle
// source0->SetSolid(geometry->GetSolid(0));  // source location

#include <cmath>
#include <cstdio>
#include "TVpSourceIso.h"
#include "misc.h"

ClassImp(TVpSourceIso)

//______________________________________________________________________________
  TVpSourceIso::TVpSourceIso(TVpSpectrum *spectrumPtr, Double_t coneAngle)
    : TVpSource(spectrumPtr)
{
  // Constructor.  The source position is undefined.
  //
  // Input parameters:
  // - spectrumPtr - energy distribution, or spectrum, of the source
  // - coneAngle - the cone angle in rad measured from (0,0,-1)

  fCosConeAngle = cos(coneAngle);
  fBias = GetSolidAngle() / (4 * M_PI);
}

//______________________________________________________________________________
void TVpSourceIso::GetParticle(TVpParticle *particle) const
{
  // Generete new particle in the source.

  // Position
  // GenerateRandomPosition(particle);

  // Generate random direction. The symetry axis is (0,0,1). 
  Double_t cosTheta = -fCosConeAngle - getRand() * (1.0 - fCosConeAngle); 
  Double_t sinTheta = sqrt(1.0 - cosTheta*cosTheta);
  Double_t phi = 2 * M_PI * getRand();
  Double_t sinPhi = sin(phi);
  Double_t cosPhi = cos(phi);
  TVpVector3 dirL(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
  particle->fDirection = DirLocToUni(dirL);
  particle->fPosition = fTraVecL2u; // = PosLocToUni(TVpVector3(0,0,0));
  particle->fEnergy = fSpectrumPtr->GetRandEnergy();
  particle->SetSolid(fSolidPtr);
  particle->fWeight = fBias;
  if (fBowTieFilterPtr != 0)
    particle->fWeight *= fBowTieFilterPtr->GetTransmission(particle->fEnergy, dirL);
}

//______________________________________________________________________________
void TVpSourceIso::GetParticleHeadedToPoint
(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel)
{
  // Generate new particle headed towards the point "point".
  // Return the directional weight (i.e. 0 for a particle outside the cone beam)
  // If the return weight is 0 then the particle is undefined.

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
Double_t TVpSourceIso::GetSolidAngle() const
{
  // Return solid angle

  return 2 * M_PI * (1.0 - fCosConeAngle);
}


//______________________________________________________________________________
void TVpSourceIso::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpSourceIso>\n"
      << "$Id: TVpSourceIso.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Cone angle: " << acos(fCosConeAngle)*180.0/M_PI << '\n'
      << "Bias: " << fBias << '\n';
  fSpectrumPtr->PrintStatus(out);
  out << "</TVpSourceIso>\n";
}


#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
  
//______________________________________________________________________________
void TVpSourceIso::Draw(Int_t option, Int_t nParticles) const
{
  // Draw the source.

  TVpSource::Draw(option, nParticles);
  if ((option & TVpSource::kDrawShape) != 0)
    {
      // FIX
    }
}
