//______________________________________________________________________________
//
// TVpVector3 implements basic 3-dimensional vector operations.
//
// The following functions and operators are defined via friend functions.
// Terminology follows Wolfram MathWorld, http://mathworld.wolfram.com/).
// - Unary minus
//   friend TVpVector3 operator - (const TVpVector3& a);
// - Vector addition
//   friend TVpVector3 operator + (const TVpVector3& a, const TVpVector3& b);
// - Vector subtraction
//   friend TVpVector3 operator - (const TVpVector3& a, const TVpVector3& b);
// - Scalar multiplication
//   friend TVpVector3 operator * (const Double_t d, const TVpVector3& a);
// - Dot product
//   friend Double_t operator | (const TVpVector3& a, const TVpVector3& b);
// - Cross product
//   friend TVpVector3 operator * (const TVpVector3& a, const TVpVector3& b);
// - Is equal (precise)
//   friend Int_t operator == (const TVpVector3& a, const TVpVector3& b);
// - Is not equal (precise)
//   friend Int_t operator != (const TVpVector3& a, const TVpVector3& b);
// - Vector Euclidean norm
//   friend Double_t norm(const TVpVector3& a);
// - Normalized vector
//   friend TVpVector3 normalize(const TVpVector3& a);
// - Write to a stream
//   friend std::ostream& operator<<(std::ostream& s, const TVpVector3& a);
// - Read from a stream
//   friend std::istream& operator>>(std::istream& s, TVpVector3& a);
//
// Note: TVpVector3 was written before the support of 3-dimensional vectors
// (see ROOT::Math::PositionVector3D and ROOT::Math::DisplacementVector3D) was
// added to ROOT.  In other words, I did not reinvent the wheel.
//
// Examples:
// gSystem->Load("libRCTmod.so");
// TVpVector3 v1(1, 2, 3);
// TVpVector3 v2(3, 2, 1);
// TVpVector3 v;
// v = -v1;
// v = v1 + v2;
// v = v1 - v2;
// v = 3.0 * v1;
// Double_t sp = (v1|v2);
// v = v1 * v2;
// cout << (v1 == TVpVector3(1, 2, 4)) << endl;
// cout << (v1 != TVpVector3(1, 2, 3)) << endl;
// cout << v1.IsEqual(TVpVector(1, 2, 3, 1e-16)) << endl;
// cout << norm(v1) << endl;
// v = normalize(v1);

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include "TVpVector3.h"

ClassImp(TVpVector3)

//______________________________________________________________________________
void TVpVector3::SetPolar(Double_t size, Double_t theta, Double_t phi)
{
  // Set vector from polar components.
  //
  // Input parameters:
  // - size - vector length
  // - theta - polar angle in radians
  // - phi - azimuthal angle in radians
  //
  // Example:
  // root [] TVpVector3 e1, e2, e3;
  // root [] e1.SetPolar(1.0, TMath::Pi()/2, 0.0);
  // root [] e1.PrintStatus();
  // <TVpVector3>    1.00000000e+00   0.00000000e+00   6.12303177e-17 </TVpVector3>
  // root [] e2.SetPolar(1.0, TMath::Pi()/2, TMath::Pi()/2);
  // root [] e2.PrintStatus();
  // <TVpVector3>    6.12303177e-17   1.00000000e+00   6.12303177e-17 </TVpVector3>
  // root [] e3.SetPolar(1.0, 0.0, 0.0);
  // root [] e3.PrintStatus();
  // <TVpVector3>    0.00000000e+00   0.00000000e+00   1.00000000e+00 </TVpVector3>

  fR[0] = size * sin(theta) * cos(phi);
  fR[1] = size * sin(theta) * sin(phi);
  fR[2] = size * cos(theta);
}

//______________________________________________________________________________
TVpVector3 operator * (const Double_t d, const TVpVector3& a)
{
  // Multiply a vector by a number
  
  TVpVector3 v;
  
  v.fR[0] = d * a.fR[0];
  v.fR[1] = d * a.fR[1];
  v.fR[2] = d * a.fR[2];
  return v;
}

//______________________________________________________________________________
TVpVector3 operator - (const TVpVector3& a)
{
  // Unary minus
  
  TVpVector3 v;
  
  v.fR[0] = -a.fR[0];
  v.fR[1] = -a.fR[1];
  v.fR[2] = -a.fR[2];
  return v;
}

//______________________________________________________________________________
TVpVector3 operator + (const TVpVector3& a, const TVpVector3& b)
{
  // Binary plus
  
  TVpVector3 v;
  
  v.fR[0] = a.fR[0] + b.fR[0];
  v.fR[1] = a.fR[1] + b.fR[1];
  v.fR[2] = a.fR[2] + b.fR[2];
  return v;
}

//______________________________________________________________________________
TVpVector3 operator - (const TVpVector3& a, const TVpVector3& b)
{
  // Binary minus

  TVpVector3 v;
  
  v.fR[0] = a.fR[0] - b.fR[0];
  v.fR[1] = a.fR[1] - b.fR[1];
  v.fR[2] = a.fR[2] - b.fR[2];
  return v;
}

//______________________________________________________________________________
Double_t operator | (const TVpVector3& a, const TVpVector3& b)
{
  // Dot product
  
  return a.fR[0] * b.fR[0] + a.fR[1] * b.fR[1] + a.fR[2] * b.fR[2];
}

//______________________________________________________________________________
TVpVector3 operator * (const TVpVector3& a, const TVpVector3& b)
{
  // Vector product

  TVpVector3 v;
  
  v.fR[0] = a.fR[1] * b.fR[2] - a.fR[2] * b.fR[1];
  v.fR[1] = a.fR[2] * b.fR[0] - a.fR[0] * b.fR[2];
  v.fR[2] = a.fR[0] * b.fR[1] - a.fR[1] * b.fR[0];
  return v;
}

//______________________________________________________________________________
Int_t operator == (const TVpVector3& a, const TVpVector3& b)
{
  // Is equal

  return &a == &b ||
    a.fR[0] == b.fR[0] && a.fR[1] == b.fR[1] && a.fR[2] == b.fR[2];
}

//______________________________________________________________________________
Int_t operator != (const TVpVector3& a, const TVpVector3& b)
{
  // Is not equal

  return !(a == b);
}

//______________________________________________________________________________
Double_t norm(const TVpVector3& a)
{
  // Return Euclidean norm of the vector (L2-norm).
  //
  // Input parameters:
  // - a - the vector
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << norm(v) << endl;
  // 3.74165739e+00

  return std::sqrt(a.fR[0] * a.fR[0] + a.fR[1] * a.fR[1] + a.fR[2] * a.fR[2]);
}

//______________________________________________________________________________
TVpVector3 normalize(const TVpVector3& a)
{
  // Return normalized vector.  The original vector is not changed.
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] TVpVector3 n = normalize(v);
  // root [] n.PrintStatus();
  // <TVpVector3>    2.67261242e-01   5.34522484e-01   8.01783726e-01 </TVpVector3>

  TVpVector3 v;
  v = (1.0 / norm(a)) * a;
  return v;
}

//______________________________________________________________________________
void TVpVector3::PrintValue()
{
  // Print vector components into stdout.  Then print the newline character.
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] v.PrintValue();
  //     1.000000     2.000000     3.000000

  printf("%12f %12f %12f\n", fR[0], fR[1], fR[2]);
}

//______________________________________________________________________________
void TVpVector3::PrintValue(FILE *fp)
{
  // Print vector components into the file fp.  Then print the newline
  // character.
  //
  // Input parameters:
  // - fp - file structure
  //
  // Example:
  // root [] FILE *fp = fopen("test.dat", "w");
  // root [] v.PrintValue(fp);
  // root [] fclose(fp);

  fprintf(fp, "%12f %12f %12f\n", fR[0], fR[1], fR[2]);
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& s, const TVpVector3& a)
{
  // Write the vector to a stream.  No newline.  Operator overloading for
  // streams does not work in a ROOT interactive session.
  //
  // Input parameters:
  // - s - ouput stream (default = cout)
  // - a - the vector
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] cout << v << endl;
  //    1.00000000e+00   2.00000000e+00   3.00000000e+00
  
  s << std::scientific << std::setprecision(8) 
    << std::setw(17) << a.fR[0]
    << std::setw(17) << a.fR[1]
    << std::setw(17) << a.fR[2];
  return s;
}

//______________________________________________________________________________
std::istream& operator>>(std::istream& s, TVpVector3& a)
{
  // Read a vector from a stream.  Operator overloading for streams does not
  // work in a ROOT interactive session.
  //
  // Input parameters:
  // - s - input stream (default = cin)
  // - a - the vector

  s >> a.fR[0] >> a.fR[1] >> a.fR[2];
  return s;
}

//______________________________________________________________________________
void TVpVector3::PrintStatus(std::ostream &out) const
{
  // Print the vector status to a stream.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpVector3 v(1, 2, 3);
  // root [] v.PrintStatus();
  // <TVpVector3>    1.00000000e+00   2.00000000e+00   3.00000000e+00 </TVpVector3>

  out << "<TVpVector3>\n"
      << "$Id: TVpVector3.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << *this
      << "\n</TVpVector3>\n"; 
}

//______________________________________________________________________________
Int_t TVpVector3::IsEqual(const TVpVector3& a, Double_t eps) const
{
  // Return 1 if all vector components are within the eps distance and 0
  // otherwise.
  //
  // Input parameters:
  // - a - tested vector
  // - eps - epsilon distance
  //
  // Example:
  // root [] e1.SetPolar(1.0, TMath::Pi()/2, 0.0);
  // root [] e1.PrintStatus();
  // <TVpVector3>    1.00000000e+00   0.00000000e+00   6.12303177e-17 </TVpVector3>
  // root [] cout << e1.IsEqual(TVpVector3(1,0,0), 1e-17) << endl;
  // 0
  // root [] cout << e1.IsEqual(TVpVector3(1,0,0), 1e-16) << endl;
  // 1
  
  return std::abs(fR[0] - a.fR[0]) < eps &&
    std::abs(fR[1] - a.fR[1]) < eps &&
    std::abs(fR[2] - a.fR[2]) < eps; 
}
