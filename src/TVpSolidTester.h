#ifndef TVpSolidTester_h
#define TVpSolidTester_h

#include "TObject.h"
#include "TVpVector3.h"
#include "TVpSolid.h"
#include "TVpVoxelArray.h"
#include "TH1.h"

class TVpSolidTester
{
 public:
  Double_t    fXMin;
  Double_t    fXMax;
  Double_t    fYMin;
  Double_t    fYMax;
  Double_t    fZMin;
  Double_t    fZMax;

  TVpSolidTester(Double_t fXMin, Double_t fXMax, Double_t fYMin, Double_t fYMax,
		 Double_t fZMin, Double_t fZMax);
  virtual ~TVpSolidTester();
  void  DrawPointsInside(TVpSolid *solid, Int_t nSamples);
  void  DrawPointsRayIntersectionIn(TVpSolid *solid, Int_t nSamples);
  void  DrawPointsRayIntersectionOut(TVpSolid *solid, Int_t nSamples);
  void  PrintVoxelArrayPathTable(TVpVoxelArray *va, Double_t energy, Int_t nSamples);
  void  PrintVoxelArrayPathDifference(TVpVoxelArray *va, Double_t energy, Int_t nSamples);
  void  PrintVoxelArrayOpticalPathTwoDirections(TVpVoxelArray *va, Double_t energy,
						Int_t nSamples);
  TVpVector3  GetRandomPoint();
  TVpVector3  GetRandomDirection();

  ClassDef(TVpSolidTester,1) // Geometry: test ray tracing routines of solids
};

#endif // TVpSolidTester_h
