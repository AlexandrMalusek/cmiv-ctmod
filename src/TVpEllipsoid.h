#ifndef TVpEllipsoid_h
#define TVpEllipsoid_h

#include "TObject.h"
#include "TVpSolid.h"

class TVpEllipsoid : public TVpSolid
{
 private:
  Double_t   fia2;  // 1/(fa*fa)
  Double_t   fib2;  // 1/(fb*fb)
  Double_t   fic2;  // 1/(fc*fc)
  
 public:
  Double_t   fa;    // length
  Double_t   fb;    // length
  Double_t   fc;    // length
  
  TVpEllipsoid();
  TVpEllipsoid(const Char_t *name, Int_t index, Double_t a, Double_t b, Double_t c);
  virtual ~TVpEllipsoid() {};

  Int_t        RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				 const Double_t l, Double_t& t) const;
  Int_t        RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				  const Double_t l, Double_t& t) const;
  Int_t        IsInside(TVpVector3& pos) const;
  Double_t     GetVolume() const;
  void         PrintStatus(std::ostream &out = std::cout) const;
  void         Draw(Int_t color = -1) const;

  ClassDef(TVpEllipsoid,1) // Geometry solid: an ellipsoid
};

#endif  // TVpEllipsoid_h
