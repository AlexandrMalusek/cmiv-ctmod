//______________________________________________________________________________
//
// TVpParticle stores kinematic parameters of a particle in a geometry.
//
// Interaction types which may appear in DED or PED files:
// 0 - kKilled                  4 - kCoherentScattering
// 1 - kTakenFromStack          5 - kPairProduction
// 2 - kPhotoefect              6 - kTrackEnd
// 3 - kComptonScattering

#include <cmath>
#include <iomanip>
#include <stdio.h>
#include "misc.h"
#include "TVpParticle.h"
#include "TVpGeometry.h"

ClassImp(TVpParticle)

//______________________________________________________________________________
void TVpParticle::NewDirection()
{
  // Calculate new direction of a particle after scattering.  Scattering
  // angles are given by data members fSinTheta, fCosTheta, fSinPhi, fCosPhi

#define DOUBLE_LIMIT 1.0e-10

  Double_t p, q;
  
  TVpVector3 oldDirection = fDirection;
  Double_t *oDir = oldDirection.fR;
  
  q = 1.0 - oDir[2]*oDir[2];
  if (q>DOUBLE_LIMIT)
    {
      p = fSinTheta / sqrt(q);
      
      fDirection.Set(oDir[0]*fCosTheta+(oDir[1]*fSinPhi-oDir[0]*oDir[2]*fCosPhi)*p,
		     oDir[1]*fCosTheta-(oDir[0]*fSinPhi+oDir[1]*oDir[2]*fCosPhi)*p,
		     oDir[2]*fCosTheta+q*fCosPhi*p);
    }
  else
    fDirection.Set(fSinTheta*fCosPhi, 
		   fSinTheta*fSinPhi,
		   (oDir[2]>0.0) ? fCosTheta: -fCosTheta);
  
#undef DOUBLE_LIMIT
}

//______________________________________________________________________________
void TVpParticle::SetIsotropicSinPhiCosPhi()
{
  // Sample azimuthal scattering angle from the uniform distribution and
  // calculate fSinPhi and fCosPhi accordingly.

  Double_t phi = 2 * M_PI * getRand();
  fSinPhi = sin(phi);
  fCosPhi = cos(phi);
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& s, TVpParticle& p)
{
  // Print particle attributes: position, direction, energy, weight, and solid
  // index in the universe coordinate system.

  s << std::scientific 
    << std::setw(14) << p.GetUniversalPosition() << ' '
    << std::setw(14) << p.GetUniversalDirection() << ' '
    << std::setw(14) << p.fEnergy << ' '
    << std::setw(14) << p.fWeight << ' '
    << ' ' << p.fSolid->GetIndex() << ' ' << p.fSolid->GetSubIndex(p.GetPos());
  return s;
}

//______________________________________________________________________________
TVpVector3 TVpParticle::GetUniversalPosition() const
{
  // Return the particle position in Universe coordinates

  TVpVector3  pos = fSolid->fRotMatL2u * fPosition + fSolid->fTraVecL2u;
  return pos;
}

//______________________________________________________________________________
TVpVector3 TVpParticle::GetUniversalDirection() const
{
  // Return the particle direction in Universe coordinates

  TVpVector3 dir = fSolid->fRotMatL2u * fDirection;
  return dir;
}

#ifdef USE_HISTORY
//______________________________________________________________________________
void TVpParticle::AddToHistory()
{
  // Add the particle to its history

  fHistory.push_back(*this);
}

//______________________________________________________________________________
void TVpParticle::ClearHistory()
{
  // Clear the particle's history
  
  fHistory.clear();
}
#endif // USE_HISTORY

#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

#ifdef USE_HISTORY
//______________________________________________________________________________
void TVpParticle::DrawTrajectory()
{
  // Draw the particle's trajectory

  std::cout << "Making polyline: " << fHistory.size() << std::endl;
  TVpVector3 pos;
  TPolyLine3D *tra = new TPolyLine3D(fHistory.size());
  for (UInt_t i = 0; i < fHistory.size(); i++)
    {
      // Use universal coordinates
      pos = fHistory[i].GetUniversalPosition();
      tra->SetPoint(i, pos.fR[0], pos.fR[1], pos.fR[2]);
    }
  tra->Draw();
}
#endif // USE_HISTORY
