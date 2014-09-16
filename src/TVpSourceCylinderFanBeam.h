#ifndef TVpSourceCylinderFanBeam_h
#define TVpSourceCylinderFanBeam_h

#include "TObject.h"
#include "TVpSource.h"
#include "TVpSpectrum.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"

class TVpSourceCylinderFanBeam : public TVpSource
{
 protected:
  Double_t      f4PiR;           // f4Pi = 4 * M_PI * fRadius
  Double_t      fRadius2;        // fRadius2 = fRadius * fRadius
  Double_t      fCosAngle;       // fCosAngle = cos(fAngle)

  Double_t j(Double_t x) const;

 public:
  Double_t      fRadius;         // = sdd, source-detector distance in cm
  Double_t      fAngle;          // -fAngle <= alpha <= fAngle
  Double_t      fSlitWidth;      // -fSlitWidth <= x <= fSlitWidth

  inline TVpSourceCylinderFanBeam();
  TVpSourceCylinderFanBeam
    (TVpSpectrum *spectrum, Double_t radius, Double_t slitX, Double_t slitY);
  virtual inline ~TVpSourceCylinderFanBeam();
  
  void     GetParticle(TVpParticle *particle) const;
  void     GetParticleHeadedToPoint(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel);
  Double_t GetSolidAngle() const;
  void     PrintStatus(std::ostream &out = std::cout) const;
  virtual Int_t GetType() const {return 2;}
  inline virtual Int_t IsInsideBeam(const TVpVector3& dirL) const;

  void     Draw
    (Int_t option = (Int_t) (TVpSource::kDrawAxes | TVpSource::kDrawShape),
     Int_t nParticles = 100) const;
  ClassDef(TVpSourceCylinderFanBeam,1) // Cylindrical fan beam source, rejection sampling
};

//______________________________________________________________________________
inline TVpSourceCylinderFanBeam::TVpSourceCylinderFanBeam()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpSourceCylinderFanBeam::~TVpSourceCylinderFanBeam()
{
  // Destructor
}

//______________________________________________________________________________
inline Int_t TVpSourceCylinderFanBeam::IsInsideBeam(const TVpVector3& dirL) const
{
  // Return 1 if the direction is inside the beam.  Return 0 otherwise.
  //
  // Input parameters:
  // - dirL - direction in the local coordinate system

  Double_t u1 = dirL.GetY();
  Double_t u2 = dirL.GetZ();
  Double_t spu = sqrt(u1*u1 + u2*u2);  // size of u projected to x=0
  Double_t cosAlpha = -u2 / spu;
  if (cosAlpha < fCosAngle)
    return 0;
  Double_t t = fRadius / spu;
  Double_t x = t * dirL.GetX();
  return -fSlitWidth <= x && x <= fSlitWidth;
}

#endif  // TVpSourceCylinderFanBeam_h
