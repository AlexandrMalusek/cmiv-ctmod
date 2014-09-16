#ifndef TVpSetupTomograph_h
#define TVpSetupTomograph_h

#include "TObject.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpSource.h"
#include "TVpPointDetectorArray.h"
#include "TVpPitchHelical.h"
#include "AVpSetup.h"
#include "TVpVolumeIntegrator.h"
#include "TVpDimensionsTomograph.h"

class TVpSetupTomograph : public AVpSetup
{
 protected:
  TVpSource              *fSourcePtr;      //! Active source (fSource0Ptr or fSource1Ptr)

 public:
  enum EPdaType {kPdaUnknown = 0, kPdaPlanar = 1, kPdaCylindrical = 2};
  enum ESourceType{kSourceUnknown = 0, kSourcePlanar = 1, kSourceCylindrical = 2, kSourceIso = 3};

  TVpSource              *fSource0Ptr;     //! Primary source of photons
  TVpSource              *fSource1Ptr;     //! Secondary source of photons
  TVpGeometry            *fGeometryPtr;    //! Geometry
  TVpPointDetectorArray  *fDetectorPtr;    //! Detector
  TVpPitch               *fPitchPtr;       //! The trajectory of the source
  TVpDimensionsTomograph *fDimensionsPtr;  //! Tomograph dimensions

  EPdaType                fPdaType;        //  PDA type
  ESourceType             fSourceType;     //  Source type
  TVpVector3              fSrcTraVec;      //  Source active translation vector
  TVpMatrix3x3            fSrcRotMat;      //  Source active rotation matrix
  TVpVector3              fPdaTraVec;      //  PDA active translation vector
  TVpMatrix3x3            fPdaRotMat;      //  PDA active rotation matrix
  
  TVpSetupTomograph();
  TVpSetupTomograph
    (TVpSource *sourcePtr, TVpGeometry *geometryPtr, TVpPointDetectorArray *detectorPtr,
     TVpPitch *pitchPtr, TVpDimensionsTomograph *dimensionsPtr);
  virtual ~TVpSetupTomograph();
  inline void SetSrcMatVec(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec);
  inline void SetPdaMatVec(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec);
  void        SetPosition(Double_t rotationAngle, Int_t sliceIndex = 0);
  void        AnalyticProjection(TVpPointDetectorArray *pda = 0, Int_t numOfPoints = 1,
				 Int_t withPhantom = 1);
  void        AnalyticSingleScatterProjection(TVpVolumeIntegrator *viPtr,
					      TVpPointDetectorArray *pdaPtr = 0,
					      Int_t numOfPoints = 1);
  //void PedFileProjection(std::istream &pedIn, TVpPointDetectorArray *pda = 0);
  inline TVpSource             *GetSourcePtr() const;
  inline TVpGeometry           *GetGeometryPtr() const;
  inline TVpPointDetectorArray *GetPointDetectorArrayPtr() const;
  inline TVpPitch              *GetPitchPtr() const;
  inline void                   SetSourceParticle(TVpParticle& particle);
  void        PrintStatus(std::ostream &out = std::cout) const;
  void        ActivateSource(Int_t index);
  void        SetSource1(TVpSource *source);

  ClassDef(TVpSetupTomograph,1) // Tomograph simulation setup
};

//______________________________________________________________________________
inline void TVpSetupTomograph::SetSrcMatVec(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec)
{
  // Set the source position in space. The default object position is that the
  // object's LCS is the same as the UCS.  Active rotation via rotMat is
  // applied first, then active translation via traVec.
  //
  // Input parameters:
  // - rotMat - active rotation
  // - traVec - active translation

  fSrcRotMat = rotMat;
  fSrcTraVec = traVec;
}

//______________________________________________________________________________
inline void TVpSetupTomograph::SetPdaMatVec(TVpMatrix3x3 const& rotMat, TVpVector3 const& traVec)
{
  // Set the PDA position in space. The default object position is that the
  // object's LCS is the same as the UCS.  Active rotation via rotMat is
  // applied first, then active translation via traVec.
  //
  // Input parameters:
  // - rotMat - active rotation
  // - traVec - active translation
 
  fPdaRotMat = rotMat;
  fPdaTraVec = traVec;
}

//______________________________________________________________________________
inline TVpSource *TVpSetupTomograph::GetSourcePtr() const
{
  // Return pointer to the tomograph's source
 
  return fSourcePtr;
}

//______________________________________________________________________________
inline TVpGeometry *TVpSetupTomograph::GetGeometryPtr() const
{
  // Return pointer to the tomograph's geometry
 
  return fGeometryPtr;
}

//______________________________________________________________________________
inline TVpPointDetectorArray *TVpSetupTomograph::GetPointDetectorArrayPtr() const
{
  // Return pointer to the tomograph's PDA
 
  return fDetectorPtr;
}

//______________________________________________________________________________
inline TVpPitch *TVpSetupTomograph::GetPitchPtr() const
{
  // Return pointer to the tomograph's pitch (source trajectory)
 
  return fPitchPtr;
}

//______________________________________________________________________________
inline void TVpSetupTomograph::SetSourceParticle(TVpParticle& particle)
{
  // Set particle parameters

  fSourcePtr->GetParticle(&particle);
}

#endif  // TVpSetupTomograph_h
