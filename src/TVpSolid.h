#ifndef TVpSolid_h
#define TVpSolid_h

#include "TObject.h"
#include <iosfwd>
#include <string>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpMaterial.h"
#include "TVpObjectLocation.h"
#include "TVpSolidNeighbor.h"
#include "TVpSolidStatistics.h"

class TVpSolidNeighbor;
class TVpSolid;
typedef TVpSolid* TVpSolidPtr;

class TVpSolid : public TVpObjectLocation, public TVpSolidStatistics
{
 public:
  Char_t           *fName;            //! Descriptive name of the solid
  Int_t             fIndex;           //  Index of the solid, must be unique in the geometry
  TVpMaterial      *fMaterial;        //! Pointer to material (homogenous solids)

  Int_t             fOverlapNumber;   //  Number of overlaping solids
  TVpSolidNeighbor *fOverlap;         //! Array of overlaping solids
  Int_t             fBaseNumber;      //  Number of base solids
  TVpSolidNeighbor *fBase;            //! Array of base solids

  TVpSolid(const Char_t *name = "", Int_t index = 0);
  virtual ~TVpSolid();

  inline void          SetMaterial(TVpMaterial *materialPtr);
  void                 SetOverlap(TVpSolidPtr *overlapList);
  void                 SetBase(TVpSolidPtr *baseList);
  void                 InitializeNeighbors();
  virtual void         PrintStatus(std::ostream &out = std::cout) const;
  virtual Int_t        RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
					 const Double_t l, Double_t& t) const = 0;
  virtual Int_t        RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
					  const Double_t l, Double_t& t) const = 0;
  virtual Int_t        IsInside(TVpVector3& pos) const = 0;
  inline Int_t         GetIndex() const;
  inline virtual TVpMaterial*  GetMaterial(Double_t rLoc[]) const;
  virtual Int_t        GetSubIndex(Double_t r[]) const;
  virtual Double_t     GetPathLength(Double_t r[], Double_t u[], Double_t energy,
				     Double_t opticalPath) const;
  virtual Double_t     GetOpticalPathInside(Double_t r[], Double_t u[], Double_t energy, 
					    Double_t distance) const;
  virtual Double_t     GetVolume() const = 0;
  Double_t             GetDensity() const;
  virtual void         Draw(Int_t color = -1) const = 0;

  ClassDef(TVpSolid,1) // Geometry: base class of all solids
};

//______________________________________________________________________________
inline void TVpSolid::SetMaterial(TVpMaterial *materialPtr)
{
  // Set material of the solid.
  //
  // Input parameters:
  // - materialPtr - pointer to the material

  fMaterial = materialPtr;
}

//______________________________________________________________________________
inline Int_t TVpSolid::GetIndex() const
{
  // Return the solid index

  return fIndex;
}

//______________________________________________________________________________
inline TVpMaterial* TVpSolid::GetMaterial(Double_t rLoc[]) const
{
  // Return the pointer to the material.  The position "rLoc" is used by voxel
  // arrays and ignored by solids which do not have internal structure.
  //
  // Input parameters:
  // - rLoc[3] - position in local coordinates. 

  return fMaterial;
}

#endif  // TVpSolid_h
