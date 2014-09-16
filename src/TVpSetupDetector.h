#ifndef TVpSetupDetector_h
#define TVpSetupDetector_h

#include "TObject.h"
#include <cmath>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpSource.h"
#include "TVpGeometry.h"
#include "AVpSetup.h"

class TVpSetupDetector : public AVpSetup
{
 public:
  TVpSource             *fSource;         //! Source of photons
  TVpGeometry           *fGeometry;       //! Geometry

  TVpVector3             fSrcTraVec;      //  Source translation vector
  TVpMatrix3x3           fSrcRotMat;      //  Source rotation matrix
  
  TVpSetupDetector();
  TVpSetupDetector(TVpSource *source, TVpGeometry *geometry);
  virtual ~TVpSetupDetector();
  inline void         SetSrcMatVec(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec);
  void                SetPosition(TVpVector3& helTraVec, TVpMatrix3x3& helRotMat);
  inline TVpSource   *GetSourcePtr();
  inline TVpGeometry *GetGeometryPtr() const;

  inline void         SetSourceParticle(TVpParticle& particle);

  ClassDef(TVpSetupDetector,1) // Detector simulation setup
};

//______________________________________________________________________________
inline void TVpSetupDetector::SetSrcMatVec(TVpMatrix3x3 const& rotMat,
					   TVpVector3 const& traVec)
{
  fSrcRotMat = rotMat;
  fSrcTraVec = traVec;
}

//______________________________________________________________________________
inline TVpSource *TVpSetupDetector::GetSourcePtr()
{
  // Return pointer to the tomograph's source
 
  return fSource;
}

//______________________________________________________________________________
inline TVpGeometry *TVpSetupDetector::GetGeometryPtr() const
{
  // Return pointer to the tomograph's geometry
 
  return fGeometry;
}

//______________________________________________________________________________
inline void TVpSetupDetector::SetSourceParticle(TVpParticle& particle)
{
  // Set particle parameters

  fSource->GetParticle(&particle);
}

#endif  // TVpSetupDetector_h
