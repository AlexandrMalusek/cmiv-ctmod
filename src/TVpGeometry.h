#ifndef TVpGeometry_h
#define TVpGeometry_h

#include "TObject.h"
#include "TView3D.h"
#include <iosfwd>
#include "TVpVoxelArray.h"
#include "TVpVector3.h"

class TVpGeometry
{
 private:
  Double_t      fInfinity;
  Double_t      fShift;

  void          InsertSolidPtr(TVpSolid *solid);

 public:
  enum EStepResult {kOverlapSolidEntered, kCurrentSolidLeft, kNoChange,
		    kTerminated};

  TVpSolid     *fUniverse;    //! The solid containing all other solids
  Int_t         fSolidNumber; //  Number of solids
  TVpSolidPtr  *fSolidPtr;    //! Array of solid pointers

  TVpGeometry();
  TVpGeometry(TVpSolid *universe);
  virtual ~TVpGeometry();
  void           Initialize();
  void           PrintStatus(std::ostream &out = std::cout) const;
  void           PrintStatistics(std::ostream &out = std::cout,
				 Long_t numOfHist = 1) const;
  EStepResult    Step(TVpParticle *particle, const Double_t l, Double_t& t) const;
  TVpSolid      *GetSolid(Int_t indexOfSolid) const;
  TVpSolid      *GetSolid(TVpVector3& posLoc, TVpSolid *base = 0) const;
  Double_t       GetOpticalPath(TVpParticle particle, Double_t distance) const;
  void           Voxelize(TVpVoxelArray *voxelArrayPtr, TVpVector3 traVec,
			  TVpMatrix3x3 rotMat) const;
  Double_t       GetNativeVolume(Int_t indexOfSolid, Int_t numInteg = 0) const;
  Double_t       GetNativeMass(Int_t indexOfSolid, Int_t numInteg = 0) const;
  Double_t       GetAbsorbedDose(Int_t indexOfSolid, Long_t numOfHist,
				 Int_t numInteg = 0) const;
  Int_t          TraceRay(TVpVector3& startPosUni, TVpVector3& endPosUni) const;

  // TVpVector3     fViewRange;                     // obsolete. the view range
  Double_t       fViewRangeMin[3];
  Double_t       fViewRangeMax[3];
  TView3D*       Draw(TView3D *view = 0) const;  // Draw the geometry in a canvas
  void           SetViewRange(const TVpVector3& viewRange);
  void           SetViewRange(const Double_t* rmin, const Double_t* rmax);

  ClassDef(TVpGeometry,1) // Simulation geometry
};

#endif  // TVpGeometry_h
