#ifndef TVpMatrix3x3_h
#define TVpMatrix3x3_h

#include "TObject.h"
#include "TVpVector3.h"

class TVpMatrix3x3;
typedef TVpMatrix3x3* TVpMatrix3x3Ptr;

class TVpMatrix3x3
{
 public:
  Double_t fM[3][3];
  
  inline TVpMatrix3x3();
  inline TVpMatrix3x3(const Double_t m00, const Double_t m01, const Double_t m02,
		      const Double_t m10, const Double_t m11, const Double_t m12,
		      const Double_t m20, const Double_t m21, const Double_t m22);
  TVpMatrix3x3(const Double_t omega, const Double_t theta, const Double_t phi);
  virtual ~TVpMatrix3x3();
  
  inline void Set(const Double_t m00, const Double_t m01, const Double_t m02,
		  const Double_t m10, const Double_t m11, const Double_t m12,
		  const Double_t m20, const Double_t m21, const Double_t m22);
  void PrintStatus(std::ostream &out = std::cout) const;
  Double_t Determinant() const;
  friend TVpMatrix3x3 operator - (const TVpMatrix3x3& a);
  friend TVpMatrix3x3 operator + (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
  friend TVpMatrix3x3 operator - (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
  friend TVpMatrix3x3 operator * (const Double_t d, const TVpMatrix3x3& a);
  friend TVpVector3 operator * (const TVpMatrix3x3& a, const TVpVector3& b);
  friend TVpMatrix3x3 operator * (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
  friend Int_t operator == (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
  friend Int_t operator != (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
  friend TVpMatrix3x3 transpose(const TVpMatrix3x3& a);
  friend std::ostream& operator<<(std::ostream& s, const TVpMatrix3x3& a);

  ClassDef(TVpMatrix3x3,1)   // Basic 3x3 matrix operations
};

TVpMatrix3x3 rotationMatrix(const TVpVector3& n, Double_t angle);
TVpMatrix3x3 rotationMatrix(const TVpVector3& oldDir, const TVpVector3& newDir);

//______________________________________________________________________________
inline TVpMatrix3x3::TVpMatrix3x3()
{
  // Default constructor.  Matrix elements are NOT initialized.
}

//______________________________________________________________________________
inline TVpMatrix3x3::TVpMatrix3x3(const double m00, const double m01, const double m02,
				  const double m10, const double m11, const double m12,
				  const double m20, const double m21, const double m22)
{
  // Constructor with full initialization
  //
  // Input parameters:
  // - m01, ..., m22 - matrix elements
  //
  // Example:
  // root[] TVpMatrix3x3 m(11, 12, 13, 21, 22, 23, 31, 32, 33);
  // root[] m.PrintStatus();
  // <TVpMatrix3>
  //    1.10000000e+01   1.20000000e+01   1.30000000e+01
  //    2.10000000e+01   2.20000000e+01   2.30000000e+01
  //    3.10000000e+01   3.20000000e+01   3.30000000e+01
  // </TVpMatrix3>

  fM[0][0]=m00; fM[0][1]=m01; fM[0][2]=m02;
  fM[1][0]=m10; fM[1][1]=m11; fM[1][2]=m12;
  fM[2][0]=m20; fM[2][1]=m21; fM[2][2]=m22;
}

//______________________________________________________________________________
inline TVpMatrix3x3::~TVpMatrix3x3()
{
  // Destructor
}

//______________________________________________________________________________
inline void TVpMatrix3x3::Set(const Double_t m00, const Double_t m01, const Double_t m02,
			      const Double_t m10, const Double_t m11, const Double_t m12,
			      const Double_t m20, const Double_t m21, const Double_t m22)
{
  // Set the matrix data
  //
  // Input parameters:
  // - m00, ..., m22 - matrix elements
  //
  // Example:
  // root[] TVpMatrix3x3 m;
  // root[] m.Set(11, 12, 13, 21, 22, 23, 31, 32, 33);
  // root[] m.PrintStatus();
  // <TVpMatrix3>
  //    1.10000000e+01   1.20000000e+01   1.30000000e+01
  //    2.10000000e+01   2.20000000e+01   2.30000000e+01
  //    3.10000000e+01   3.20000000e+01   3.30000000e+01
  // </TVpMatrix3>

  fM[0][0]=m00; fM[0][1]=m01; fM[0][2]=m02;
  fM[1][0]=m10; fM[1][1]=m11; fM[1][2]=m12;
  fM[2][0]=m20; fM[2][1]=m21; fM[2][2]=m22;
}

#endif // TVpMatrix3x3
