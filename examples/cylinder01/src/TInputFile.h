#ifndef TInputFile_h
#define TInputFile_h

#include <iostream>
#include "TObject.h"

class TInputFile
{
 private:
  ULong_t     fNumOfPhotons;
  Double_t    fCylinderRadius;
  Double_t    fCylinderHeight;
  Double_t    fBeamWidth;
 
 public:
  TInputFile() {}
  ~TInputFile() {};
  
  void Read(std::istream &in = std::cin);
  void Write(std::ostream &out = std::cout) const;
  
  ULong_t  GetNumOfPhotons() const {return fNumOfPhotons;}
  Double_t GetCylinderRadius() const {return fCylinderRadius;}
  Double_t GetCylinderHeight() const {return fCylinderHeight;}
  Double_t GetBeamWidth() const {return fBeamWidth;}
};

#endif   // TInputFile
