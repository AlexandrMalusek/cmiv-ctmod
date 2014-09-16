#ifndef TVpSphere_h
#define TVpSphere_h

#include "TObject.h"
#include "TVpSolid.h"

class TVpSphere : public TVpSolid
{
 private:
  Double_t   fRadius2;   // radius squared
  
 public:
  Double_t   fRadius;    // radius

  TVpSphere();
  TVpSphere(const Char_t *name, Int_t index, Double_t radius);
  virtual ~TVpSphere() {};

  Int_t        RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				 const Double_t l, Double_t& t) const;
  Int_t        RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				  const Double_t l, Double_t& t) const;
  Int_t        IsInside(TVpVector3& pos) const;
  Double_t     GetVolume() const;
  void         PrintStatus(std::ostream &out = std::cout) const;
  void         Draw(Int_t color = -1) const;

  ClassDef(TVpSphere,1) // Geometry solid: a sphere
};

#endif  // TVpSphere_h
