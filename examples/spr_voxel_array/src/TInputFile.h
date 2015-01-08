#ifndef TInputFile_h
#define TInputFile_h

#include <iostream>
#include "TObject.h"

class TInputFile
{
 private:
  ULong_t     fNumOfPhotons;
  Double_t    fEnergy;
  Double_t    fVoxelArrayCenterX;
  Double_t    fVoxelArrayCenterY;
  Double_t    fVoxelArrayCenterZ;
  Double_t    fBeamWidth;

 public:
  TInputFile() {}
  ~TInputFile() {};
  
  void Read(std::istream &in = std::cin);
  void Write(std::ostream &out = std::cout) const;
  
  ULong_t  GetNumOfPhotons() const {return fNumOfPhotons;}
  Double_t GetEnergy() const {return fEnergy;}
  Double_t GetVoxelArrayCenterX() const {return fVoxelArrayCenterX;}
  Double_t GetVoxelArrayCenterY() const {return fVoxelArrayCenterY;}
  Double_t GetVoxelArrayCenterZ() const {return fVoxelArrayCenterZ;}
  Double_t GetBeamWidth() const {return fBeamWidth;}
};

#endif   // TInputFile
