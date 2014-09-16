#ifndef TVpBox_h
#define TVpBox_h

#include "TObject.h"
#include "TVpSolid.h"

class TVpBox : public TVpSolid
{
 public:
  Double_t fSizeX;
  Double_t fSizeY;
  Double_t fSizeZ;
  
  TVpBox();
  TVpBox(const Char_t *name, Int_t index, Double_t sizeX, Double_t sizeY, Double_t sizeZ);
  virtual ~TVpBox() {};

  Int_t      RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
			       const Double_t l, Double_t& t) const;
  Int_t      RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				const Double_t l, Double_t& t) const;
  Int_t      IsInside(TVpVector3& pos) const;
  Double_t   GetVolume() const;
  void       PrintStatus(std::ostream &out = std::cout) const;
  inline Double_t GetSizeX() const;
  inline Double_t GetSizeY() const;
  inline Double_t GetSizeZ() const;
  void       Draw(Int_t color = -1) const;

  ClassDef(TVpBox,1) // Geometry solid: a box
};

//______________________________________________________________________________
inline Double_t TVpBox::GetSizeX() const
{
  // Return the box size in x-direction in cm.

  return fSizeX;
}

//______________________________________________________________________________
inline Double_t TVpBox::GetSizeY() const
{
  // Return the box size in y-direction in cm.
  
  return fSizeY;
}

//______________________________________________________________________________
inline Double_t TVpBox::GetSizeZ() const
{
  // Return the box size in z-direction in cm.

  return fSizeZ;
}

#endif  // TVpBox_h
