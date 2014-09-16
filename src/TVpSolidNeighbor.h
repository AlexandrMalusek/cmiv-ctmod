#ifndef TVpSolidNeighbor_h
#define TVpSolidNeighbor_h

#include "TObject.h"
#include <ostream>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpSolid.h"

class TVpSolid;

class TVpSolidNeighbor
{
 public:
  TVpSolid     *fSolidPtr;      //! Solid
  TVpVector3   *fTraVecC2nPtr;  //! Translation vector
  TVpMatrix3x3 *fRotMatC2nPtr;  //! Rotation matrix

  TVpSolidNeighbor();
  virtual ~TVpSolidNeighbor();

  void          PrintStatus(std::ostream &out = std::cout) const;

  ClassDef(TVpSolidNeighbor,1) // Geometry: solid neighbours - coordinate transformations
};

#endif  // TVpSolidNeighbor_h
