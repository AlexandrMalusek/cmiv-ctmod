#ifndef TVpSourceIso_h
#define TVpSourceIso_h

#include "TObject.h"
#include "TVpSource.h"
#include "TVpSpectrum.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"

class TVpSourceIso : public TVpSource
{
 public:
  Double_t      fCosConeAngle;    // cosine of the cone angle

  inline TVpSourceIso();
  TVpSourceIso(TVpSpectrum *spectrumPtr, Double_t coneAngle);
  virtual inline ~TVpSourceIso();
  
  void     GetParticle(TVpParticle *particle) const;
  void     GetParticleHeadedToPoint(TVpParticle *particle, TVpVector3 *pointUni, Int_t channel);
  Double_t GetSolidAngle() const;
  void     PrintStatus(std::ostream &out = std::cout) const;
  virtual Int_t GetType() const {return 3;}
  inline virtual Int_t IsInsideBeam(const TVpVector3& dirL) const;

  void     Draw
    (Int_t option = (Int_t) (TVpSource::kDrawAxes | TVpSource::kDrawShape),
     Int_t nParticles = 100) const;

  ClassDef(TVpSourceIso,1) // Isotropic source
};

//______________________________________________________________________________
inline TVpSourceIso::TVpSourceIso()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpSourceIso::~TVpSourceIso()
{
  // Destructor
}

//______________________________________________________________________________
inline Int_t TVpSourceIso::IsInsideBeam(const TVpVector3& dirL) const
{
  // Return 1 if the direction is inside the beam.  Return 0 otherwise.
  //
  // Input parameters:
  // - dirL - direction in the local coordinate system

  return (-dirL.GetZ() <= fCosConeAngle) ? 0 : 1;
}


#endif  // TVpSourceIso_h
