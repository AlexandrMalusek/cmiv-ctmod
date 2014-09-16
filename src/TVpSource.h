#ifndef TVpSource_h
#define TVpSource_h

#include "TObject.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpObjectLocation.h"
#include "TVpParticle.h"
#include "TVpSolid.h"
#include "AVpBowTieFilter.h"

class TVpSource : public TVpObjectLocation
{
 public:
  TVpSolid        *fSolidPtr;          //! Pointer to the solid where the source is
  TVpSpectrum     *fSpectrumPtr;       //! Pointer to the energy spectrum
  Double_t         fBias;              //  Bias factor
  AVpBowTieFilter *fBowTieFilterPtr;   //! Pointer to the bow-tie filter

 public:
  TVpSource();
  TVpSource(TVpSpectrum *spectrumPtr);
  virtual ~TVpSource();

  inline void      SetSolid(TVpSolid *solid);
  inline void      SetBowTieFilter(AVpBowTieFilter *bowTieFilterPtr);
  virtual void     GetParticle(TVpParticle *particle) const = 0;
  virtual void     GetParticleHeadedToPoint
    (TVpParticle *particle, TVpVector3 *point, Int_t channel) = 0;
  virtual Double_t GetBias() const;
  virtual Double_t GetSolidAngle() const = 0;
  virtual void     PrintStatus(std::ostream &out = std::cout) const = 0;
  virtual Int_t    GetType() const {return 0;}
  virtual Int_t    IsInsideBeam(const TVpVector3& dirL) const = 0;
  void             GetParticlePosition(TVpVector3& position, Int_t pointIndex,
				       Int_t numOfPoints) const;

  enum EDrawOption {kDrawBeamAxis = 1, kDrawAxes = 2, kDrawShape = 4, kDrawRandomParticles = 8};
  virtual void     Draw
    (Int_t option = (Int_t)(kDrawAxes | kDrawShape), Int_t nParticles = 100) const;
  ClassDef(TVpSource,1) // Source base class
};


//______________________________________________________________________________
inline void TVpSource::SetSolid(TVpSolid *solidPtr)
{
  // Set the solid where the source is located.
  //
  // Input parameters:
  // - solidPtr - pointer to the solid containg the source
  
  fSolidPtr = solidPtr;
}

//______________________________________________________________________________
inline void TVpSource::SetBowTieFilter(AVpBowTieFilter *bowTieFilterPtr)
{
  // Set the pointer to a bow-tie filter
  //
  // Input parameters:
  // - bowTieFilterPtr - pointer to a bow-tie filter

  fBowTieFilterPtr = bowTieFilterPtr;
}

//______________________________________________________________________________
inline Double_t TVpSource::GetBias() const
{
  // Return the bias of the source.

  return fBias;
}

#endif   // TVpSource_h
