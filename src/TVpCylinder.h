#ifndef TVpCylinder_h
#define TVpCylinder_h

#include "TObject.h"
#include "TVpSolid.h"

class TVpCylinder : public TVpSolid
{
 private:
  Double_t   fRadius2;   // radius squared

 public:
  Double_t   fRadius;    // radius
  Double_t   fHeight;    // height

  TVpCylinder();
  TVpCylinder(const Char_t *name, Int_t index, Double_t radius, Double_t height);
  virtual ~TVpCylinder() {};

  Int_t        RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				 const Double_t l, Double_t& t) const;
  Int_t        RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				  const Double_t l, Double_t& t) const;
  Int_t        IsInside(TVpVector3& pos) const;
  Double_t     GetVolume() const;
  void         PrintStatus(std::ostream &out = std::cout) const;

  void       Draw(Int_t color = -1) const;
  ClassDef(TVpCylinder,1) // Geometry solid: a cylinder
};

#endif  // TVpCylinder_h
