#ifndef TVpSourcePlaneFanBeam_h
#define TVpSourcePlaneFanBeam_h

#include "TObject.h"
#include "TVpSource.h"
#include "TVpSpectrum.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"

class TVpSourcePlaneFanBeam : public TVpSource
{
 protected:
  Double_t      f4Pi;            // f4Pi = 4 * M_PI
  Double_t      fDistToSlit2;    // fDistToSlit2 = fDistToSlit * fDistToSlit

  Double_t j(Double_t x, Double_t y) const;

 public:
  Double_t      fDistToSlit;     //
  Double_t      fSlitLength;     // fSlitLength = l / 2, where l is the length
  Double_t      fSlitWidth;      // fSlitWidth = w / 2, where w is the width

  TVpSourcePlaneFanBeam();
  TVpSourcePlaneFanBeam
    (TVpSpectrum *spectrumPtr, Double_t distToSlit, Double_t slitWidthX, Double_t slitLengthY);
  virtual ~TVpSourcePlaneFanBeam();
  
  void     GetParticle(TVpParticle *particle) const;
  void     GetParticleHeadedToPoint(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel);
  Double_t GetSolidAngle() const;
  void     PrintStatus(std::ostream &out = std::cout) const;
  virtual Int_t GetType() const {return 1;}
  inline virtual Int_t IsInsideBeam(const TVpVector3& dirL) const;
  
  void     Draw
    (Int_t option = (Int_t) (TVpSource::kDrawAxes | TVpSource::kDrawShape),
     Int_t nParticles = 100) const;
  ClassDef(TVpSourcePlaneFanBeam,1) // Planar fan beam source, rejection sampling
};

//______________________________________________________________________________
inline Int_t TVpSourcePlaneFanBeam::IsInsideBeam(const TVpVector3& dirL) const
{
  // Return 1 if the direction is inside the beam.  Return 0 otherwise.
  //
  // Input parameters:
  // - dirL - direction in the local coordinate system

  Double_t t = -fDistToSlit / dirL.GetZ();
  if (t < 0)
    return 0;
  Double_t x = t * dirL.GetX();
  if (x < -fSlitWidth || fSlitWidth < x)
    return 0;
  Double_t y = t * dirL.GetY();
  return (-fSlitLength <= y && y <= fSlitLength);
}

#endif  // TVpSourcePlaneFanBeam_h
