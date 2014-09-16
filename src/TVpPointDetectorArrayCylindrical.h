#ifndef TVpPointDetectorArrayCylindrical_h
#define TVpPointDetectorArrayCylindrical_h

#include "TVpPointDetectorArray.h"

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

class TVpPointDetectorArrayCylindrical : public TVpPointDetectorArray
{
 public:
  Double_t   fArrayRadius;    // cylinder radius in cm
  Double_t   fArrayWidth;     // PDA width = cylinder height
  Double_t   fArrayLength;    // PDA lenght = length of the circle segment

  TVpPointDetectorArrayCylindrical();
  TVpPointDetectorArrayCylindrical
    (Double_t arrayWidth, Double_t arrayLength, Double_t arrayRadius,
     Int_t na, Int_t nr, Int_t numOfChannels, Double_t minEnergy,
     Double_t maxEnergy, TVpDetectorResponse *detectorResponse = 0);
  TVpPointDetectorArrayCylindrical(Char_t *fileName);
  virtual ~TVpPointDetectorArrayCylindrical() {};
  void     PrintStatus(std::ostream &out = std::cout) const;
  Int_t    GetType() const {return 2;}
  Int_t    WritePdaFile(Char_t *fileName, Int_t withErrors = 0);
  Int_t    ReadPdaFile(Char_t *fileName);
  TVpPointDetectorArray* Clone(Int_t numDivX = 1, Int_t numDivY = 1) const;
  void     SetDefaultPosition();

  TH2F   *GetImage(Int_t isRotated = 0) const;
  TH1F   *GetYProfile(Int_t sliceX, Int_t isCm = 0) const;

  ClassDef(TVpPointDetectorArrayCylindrical,1) // Point detector array, cylindrical shape
};

#endif // TVpPointDetectorArrayCylindrical_h
