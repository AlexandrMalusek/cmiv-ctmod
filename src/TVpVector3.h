#ifndef TVpVector3_h
#define TVpVector3_h

#include <cstdio>
#include <iostream>
#include <ostream>

#include "TObject.h"

class TVpVector3;
typedef TVpVector3* TVpVector3Ptr;

class TVpVector3
{
 public:
  Double_t fR[3];

  inline TVpVector3();
  inline TVpVector3(Double_t x, Double_t y, Double_t z);
  virtual ~TVpVector3() {};

  inline Double_t *Get();
  inline void      Set(Double_t x, Double_t y, Double_t z);
  inline Double_t  GetX() const;
  inline Double_t  GetY() const;
  inline Double_t  GetZ() const;
  void             SetPolar(Double_t size, Double_t theta, Double_t phi);
  void             PrintValue();
  void             PrintValue(FILE *fp);
  void             PrintStatus(std::ostream &out = std::cout) const;
  Int_t            IsEqual(const TVpVector3& a, Double_t eps) const;

  friend TVpVector3 operator - (const TVpVector3& a);
  friend TVpVector3 operator + (const TVpVector3& a, const TVpVector3& b);
  friend TVpVector3 operator - (const TVpVector3& a, const TVpVector3& b);
  friend TVpVector3 operator * (const Double_t d, const TVpVector3& a);
  friend Double_t operator | (const TVpVector3& a, const TVpVector3& b);
  friend TVpVector3 operator * (const TVpVector3& a, const TVpVector3& b);
  friend Int_t operator == (const TVpVector3& a, const TVpVector3& b);
  friend Int_t operator != (const TVpVector3& a, const TVpVector3& b);
  friend Double_t norm(const TVpVector3& a);
  friend TVpVector3 normalize(const TVpVector3& a);
  friend std::ostream& operator<<(std::ostream& s, const TVpVector3& a);
  friend std::istream& operator>>(std::istream& s, TVpVector3& a);

  ClassDef(TVpVector3,1) // Basic 3 vector operations
};

//______________________________________________________________________________
inline TVpVector3::TVpVector3()
{
  // Default constructor.  All vector components are initialized to 0.
  //
  // Example:
  // root [] TVpVector3 v;

  fR[0] = fR[1] = fR[2] = 0.0;
}

//______________________________________________________________________________
inline TVpVector3::TVpVector3(Double_t x, Double_t y, Double_t z)
{
  // Constructor with full initialization.
  //
  // Input parameters:
  // - x,y,z - vector components
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);

  fR[0] = x; fR[1] = y; fR[2] = z;
}

//______________________________________________________________________________
inline void TVpVector3::Set(Double_t x, Double_t y, Double_t z)
{
  // Set x, y, and z components of a vector.
  //
  // Input parameters;
  // - x,y,z - vector components
  //
  // Example:
  // root [] TVpVector3 v;
  // root [] v.Set(1, 2, 3);

  fR[0] = x; fR[1] = y; fR[2] = z;
}

//______________________________________________________________________________
inline Double_t *TVpVector3::Get()
{
  // Return a pointer to the first component of the internal 3-component
  // array.  See the internal representation of TVpVector3.
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << v.Get()[0] << endl;
  // 1.00000000e+00

  return fR;
}

//______________________________________________________________________________
inline Double_t TVpVector3::GetX() const
{
  // Return x-component of the vector
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << v.GetX() << endl;
  // 1.00000000e+00

  return fR[0];
}

//______________________________________________________________________________
inline Double_t TVpVector3::GetY() const
{
  // Return y-component of the vector
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << v.GetY() << endl;
  // 2.00000000e+00
  
  return fR[1];
}

//______________________________________________________________________________
inline Double_t TVpVector3::GetZ() const
{
  // Return z-component of the vector
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << v.GetZ() << endl;
  // 3.00000000e+00 

  return fR[2];
}

#endif // TVpVector3_h
