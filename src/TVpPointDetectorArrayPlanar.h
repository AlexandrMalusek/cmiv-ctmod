#ifndef TVpPointDetectorArrayPlanar_h
#define TVpPointDetectorArrayPlanar_h

#include "TVpPointDetectorArray.h"
#include "TObject.h"

class TVpPointDetectorArrayPlanar : public TVpPointDetectorArray
{
 public:
  Double_t   fArraySizeX;    // Size of the array in X-direction
  Double_t   fArraySizeY;    // Size of the array in Y-direction
  
  TVpPointDetectorArrayPlanar();
  TVpPointDetectorArrayPlanar
    (Double_t arraySizeX, Double_t arraySizeY, Int_t nx, Int_t ny, Int_t numOfChannels,
     Double_t minEnergy, Double_t maxEnergy, TVpDetectorResponse *detectorResponse = 0);
  TVpPointDetectorArrayPlanar(const Char_t *fileName);
  virtual ~TVpPointDetectorArrayPlanar() {};
  void     PrintStatus(std::ostream &out = std::cout) const;
  Int_t    GetType() const {return 1;}
  Int_t    WritePdaFile(const Char_t *fileName, Int_t withErrors = 0);
  Int_t    ReadPdaFile(const Char_t *fileName);
  TVpPointDetectorArray* Clone(Int_t numDivX = 1, Int_t numDivY = 1) const;
  void     SetDefaultPosition();

  ClassDef(TVpPointDetectorArrayPlanar,1) // Point detector array, planar shape
};

#endif // TVpPointDetectorArrayPlanar_h
