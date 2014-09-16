//______________________________________________________________________________
//
// TVpSource is an abstract class implementing a source of photons.
//
// Source intensity:
// 1.  All results correspond to an isotropic source which emits 1 photon into
// the solid angle 4Pi.
// 2.  If a derived class implements a source which emits photons into a solid
// angle Omega which is smaller than 4Pi (e.g. to a cone) than it MUST set
// fBias = Omega/(4Pi).  The weight w' of generated photons is then decreased
// as w' = fBias * w, where w is the original photon weight which may be set
// by the TVpSpectrum class.  In an analog simulation, photons emited outside
// Omega would be killed.  Here, we don't simulate them but we adjust the
// weight of the ones emitted into Omega instead.
//
// Coordinate system:
// 1. Particle positions and directions are given in the universe coordinate
// system and the source must be located in the universe, i.e. fSolidPtr must
// point to the universe.  Other setup will break the coordinate
// transformation routines.
// 2. When a new particle is generated from the source, it is first assumed
// that the source is in its default position, i.e. its center is (0,0,SAD)
// and the beam axis points at (0,0,-1).  Then the particle's position is set
// to fCenter and its direction is transformed using the rotation matrix
// fRotMatrix.  Note:  This concept may change in the future.
//
// Update: Two concepts of transformation.
//
// 1. Originally, the position of the source was defined in derived classes.
// These might have defined e.g. the center and the beam axis.  This concept
// is not good when (in principle) anisotropic elements are added to sources.
// Thus a new concept has been developed.  Old routines:
// - SetPosition(const TVpVector3& center, const TVpVector3& direction) = 0;
// - Move(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec) = 0;
//
// 2. Each object derives from the TVpObjectLocation class.  Transformation
// are defined via active translations and rotations.  New routines:
// - SetActiveTranslation(TVpVector3 *traVec);
// - SetActiveRotation(TVpMatrix3x3 *rotMat);

#include <stdio.h>
#include "TVpSource.h"

ClassImp(TVpSource)

//______________________________________________________________________________
TVpSource::TVpSource()
{
  // Default constructor
  
  fSolidPtr = 0;
  fSpectrumPtr = 0;
  fBias = 1.0;
  fBowTieFilterPtr = 0;
}

//______________________________________________________________________________
TVpSource::TVpSource(TVpSpectrum *spectrumPtr)
{
  // Constructor with partial initialization. 
  //
  // Input parameters:
  // - spectrumPtr - pointer to the X-ray spectrum
  
  fSolidPtr = 0;
  fSpectrumPtr = spectrumPtr;
  fBias = 1.0;
  fBowTieFilterPtr = 0;
}

//______________________________________________________________________________
TVpSource::~TVpSource()
{
  // Destructor.  The X-ray spectrum is not deleted.
}

//______________________________________________________________________________
void TVpSource::GetParticlePosition
(TVpVector3& position, Int_t pointIndex, Int_t numOfPoints) const
{
  // Get the particle position.  Limitation: Only the center is considered
  // now.
  //
  // Input parameters:
  // - pointIndex - index of a multipoint source.  Currently ignored.
  // - numOfPoints - number of points in the multipont source.  Curently
  //   ignored.
  //
  // Output parameters:
  // - position - particle's position

  position = fTraVecL2u;
}

#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
//______________________________________________________________________________
void TVpSource::Draw(Int_t option, Int_t nParticles) const
{
  // Draw a filled circle at the centre of the source and a line segment in
  // the direction of the beam axis.
  
  TVpVector3 posSrcU = PosLocToUni(TVpVector3(0,0,0));  // position
  TVpVector3 dirSrcU = DirLocToUni(TVpVector3(0,0,-1));  // beam axis direction 
  
  if ((option & TVpSource::kDrawAxes) != 0)
    DrawLcsAxes(5, 5, 5);
  if ((option & TVpSource::kDrawBeamAxis) != 0)
    {
      TVpVector3 posNewU = posSrcU + 10*dirSrcU;
      TPolyLine3D *sdir = new TPolyLine3D(2);
      sdir->SetPoint(0, posSrcU.GetX(), posSrcU.GetY(), posSrcU.GetZ());
      sdir->SetPoint(1, posNewU.GetX(), posNewU.GetY(), posNewU.GetZ());
      sdir->Draw();
    }
  if ((option & TVpSource::kDrawRandomParticles) != 0)
    {
      TPolyLine3D *line;
      TVpVector3 src;
      TVpVector3 end;
      TVpParticle particle;
      for (Int_t i = 0; i < nParticles; i++)
	{
	  // Set particle's parameters
	  GetParticle(&particle);
	  src = particle.fPosition;
	  // The length is proportional to the photon's weight
	  end = particle.fPosition + 100 * (particle.fWeight / fBias)
	    * particle.fDirection;
	  line = new TPolyLine3D(2);
	  line->SetPoint(0, src.GetX(), src.GetY(), src.GetZ());
	  line->SetPoint(1, end.GetX(), end.GetY(), end.GetZ());
	  line->SetLineColor(kBlack);
	  line->Draw();
	}
    }
}
