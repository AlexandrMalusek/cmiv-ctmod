#ifndef TVpMath_h
#define TVpMath_h

#include "TObject.h"
#include "TVectorF.h"

class TVpMath
{
 public:
  TVpMath();
  virtual ~TVpMath() {};
  static Double_t GetLinLinInterpolation(Double_t x1, Double_t y1,
					 Double_t x2, Double_t y2, Double_t x);
  static Double_t GetLinLinInterpolation(Double_t x, Int_t n, Double_t x0, Double_t dx,
					 Double_t *y);
  static Double_t GetLogLogInterpolation(Double_t x1, Double_t y1,
					 Double_t x2, Double_t y2, Double_t x);
  static void SetPoly1(Double_t& p0, Double_t& p1, Double_t x1, Double_t y1,
		       Double_t x2, Double_t y2);
  static Int_t FindIndexByBinarySearch(Int_t n, Double_t *g, Double_t x);
  static Int_t FindIndexByLinearSearch(Int_t n, Double_t *g, Double_t x);

  static inline Int_t Get3dArrayIndex(Int_t ix, Int_t iy, Int_t iz,
				      Int_t nx, Int_t ny, Int_t nz);
  static void RebinCtArray(TVectorF &cto, Int_t nxo, Int_t nyo, Int_t nzo,
			   Int_t mx, Int_t my, Int_t mz,
			   TVectorF &ctn, Int_t& nxn, Int_t& nyn, Int_t& nzn);
				

  ClassDef(TVpMath,0)   // Math interpolation routines
};

//______________________________________________________________________________
inline Int_t TVpMath::Get3dArrayIndex(Int_t ix, Int_t iy, Int_t iz,
				      Int_t nx, Int_t ny, Int_t nz)
{
  // Return index of an element in an nx * ny * nz array.

  return iz*ny*nx + iy*nx + ix;
}

#endif  // TVpMath_h
