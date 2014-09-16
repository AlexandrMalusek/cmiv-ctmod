#ifndef TVpObjectLocation_h
#define TVpObjectLocation_h

#include "TObject.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"

class TVpObjectLocation
{
 public:
  TVpVector3        fTraVecL2u;  //  Active translation vector, Loc -> Uni
  TVpMatrix3x3      fRotMatL2u;  //  Active rotation matrix, Loc -> Uni
  TVpMatrix3x3      fRotMatL2uT; //  fRotMatL2uT = transpose(fRotMatL2u)
  
  TVpObjectLocation();
  virtual ~TVpObjectLocation();

  void           SetActiveTranslation(TVpVector3 *traVec);
  void           SetActiveRotation(TVpMatrix3x3 *rotMat);
  void           SetActiveTransformation(TVpMatrix3x3 *rotMat, TVpVector3 *traVec);
  TVpVector3     PosLocToUni(const TVpVector3& posLoc) const;
  TVpVector3     DirLocToUni(const TVpVector3& dirLoc) const;
  TVpVector3     PosUniToLoc(const TVpVector3& posUni) const;
  TVpVector3     DirUniToLoc(const TVpVector3& dirUni) const;
  void           PrintLocation(std::ostream &out = std::cout) const;

  void           DrawLcsAxes(Double_t lengthX, Double_t lengthY, Double_t lengthZ) const;
  ClassDef(TVpObjectLocation,1)  // ObjectLocation base class
};

#endif   // TVpObjectLocation_h
