//______________________________________________________________________________
//
// TVpMatrix3x3 implements basic 3x3 matrix operations.
//
// The following functions and operators are defined via friend functions.
// - Unary minus
//   friend TVpMatrix3x3 operator - (const TVpMatrix3x3& a);
// - Matrix addition
//   friend TVpMatrix3x3 operator + (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
// - Matrix subtraction
//   friend TVpMatrix3x3 operator - (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
// - Scalar multiplication
//   friend TVpMatrix3x3 operator * (const Double_t d, const TVpMatrix3x3& a);
// - Matrix-vector multiplication
//   friend TVpVector3 operator * (const TVpMatrix3x3& a, const TVpVector3& b);
// - Matrix multiplication
//   friend TVpMatrix3x3 operator * (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
// - Rotation matrix (axis and angle given)
//   friend TVpMatrix3x3 rotationMatrix(const TVpVector3& n, Double_t angle);
// - Rotation matrix (old and new directions given)
//   friend TVpMatrix3x3 rotationMatrix(const TVpVector3& oldDir, const TVpVector3& newDir);
// - Is equal (up to eps absolute difference allowed)
//   friend int operator == (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
// - Is not equal (absolute difference larger than eps)
//   friend int operator != (const TVpMatrix3x3& a, const TVpMatrix3x3& b);
// - Transposed matrix
//   friend TVpMatrix3x3 transpose(const TVpMatrix3x3& a);
// - Write to a stream
//   friend std::ostream& operator<<(std::ostream& s, const TVpMatrix3x3& a);
//
// Note: TVpMatrix3x3 was written before matrix support was added to ROOT.  In
// other words, I did not reinvent the wheel.

#include <cmath>
#include <cstdio>
#include <iomanip>
#include "TVpMatrix3x3.h"

ClassImp(TVpMatrix3x3)

//______________________________________________________________________________
TVpMatrix3x3::TVpMatrix3x3(const Double_t omega, const Double_t theta, const Double_t phi)
{
  // Constructor based on Euler's angles.  Rotations are active and positive.
  //
  // First, a rotation of angle omega about the z-axis, followed by a rotation
  // of theta about the y-axis and, finally, a rotation of angle phi about the
  // z-axis.
  //
  // R(omega, theta, phi) = Rz(phi) Ry(theta) Rz(omega)
  //
  // Method:
  // The matrix was calculated using Mathematica:
  // In[]:= R1 = RotationMatrix3D[-Omega, 0, 0]
  // Out[]:= {{Cos[Omega], -Sin[Omega], 0}, {Sin[Omega], Cos[Omega], 0}, {0, 0, 1}}
  // In[]:= R2 = RotationMatrix3D[Pi/2, -Theta, -Pi/2]
  // Out[]:= {{Cos[Theta], 0, Sin[Theta]}, {0, 1, 0}, {-Sin[Theta], 0, Cos[Theta]}}
  // In[]:= R3 = RotationMatrix3D[-Phi, 0, 0]
  // Out:= {{Cos[Phi], -Sin[Phi], 0}, {Sin[Phi], Cos[Phi], 0}, {0, 0, 1}}
  // In[]:= R3.R2.R1
  // Out[]:= {{Cos[Omega] Cos[Phi] Cos[Theta] - Sin[Omega] Sin[Phi],
  //           -Cos[Phi] Cos[Theta] Sin[Omega] - Cos[Omega] Sin[Phi], 
  //           Cos[Phi] Sin[Theta]},
  //          {Cos[Phi] Sin[Omega] + Cos[Omega] Cos[Theta] Sin[Phi],
  //           Cos[Omega] Cos[Phi] - Cos[Theta] Sin[Omega] Sin[Phi],
  //           Sin[Phi] Sin[Theta]},
  //          {-Cos[Omega] Sin[Theta],
  //           Sin[Omega] Sin[Theta],
  //           Cos[Theta]}}
  //
  // Example:
  // // active positive (right-hand) rotation around x-axis by Pi/2
  // root[] TVpMatrix3x3 m(TMath::Pi()/2, TMath::Pi()/2, -TMath::Pi()/2);
  // root[] m.PrintStatus();
  // <TVpMatrix3>
  //    1.00000000e+00   6.12303177e-17   6.12303177e-17
  //    6.12303177e-17   6.12303177e-17  -1.00000000e+00
  //   -6.12303177e-17   1.00000000e+00   6.12303177e-17
  // </TVpMatrix3>
  // root[] (m * TVpVector3(0,1,0)).PrintStatus();
  // <TVpVector3>    6.12303177e-17   6.12303177e-17   1.00000000e+00 </TVpVector3>

  Double_t cosOmega = cos(omega);
  Double_t sinOmega = sin(omega);
  Double_t cosPhi = cos(phi);
  Double_t sinPhi = sin(phi);
  Double_t cosTheta = cos(theta);
  Double_t sinTheta = sin(theta);
  fM[0][0] =  cosPhi*cosTheta*cosOmega - sinPhi*sinOmega;
  fM[0][1] = -cosPhi*cosTheta*sinOmega - sinPhi*cosOmega;
  fM[0][2] =  cosPhi*sinTheta;
  fM[1][0] =  sinPhi*cosTheta*cosOmega + cosPhi*sinOmega;
  fM[1][1] = -sinPhi*cosTheta*sinOmega + cosPhi*cosOmega;
  fM[1][2] =  sinPhi*sinTheta;
  fM[2][0] = -sinTheta*cosOmega;
  fM[2][1] =  sinTheta*sinOmega;
  fM[2][2] =  cosTheta;
}

//______________________________________________________________________________
Double_t TVpMatrix3x3::Determinant() const
{
  // Return the 3rd order determinant of the 3x3 matrix.
  
  Double_t det = fM[0][0]*fM[1][1]*fM[2][2] + fM[0][1]*fM[1][2]*fM[2][0]
    + fM[0][2]*fM[1][0]*fM[2][1] - fM[0][0]*fM[1][2]*fM[2][1]
    - fM[0][1]*fM[1][0]*fM[2][2] - fM[0][2]*fM[1][1]*fM[2][0];
  return det;
}

//______________________________________________________________________________
TVpMatrix3x3 operator - (const TVpMatrix3x3& a)
{
  // Unary minus
  
  TVpMatrix3x3 b;

  b.fM[0][0] = -a.fM[0][0];
  b.fM[0][1] = -a.fM[0][1];
  b.fM[0][2] = -a.fM[0][2];
  b.fM[1][0] = -a.fM[1][0];
  b.fM[1][1] = -a.fM[1][1];
  b.fM[1][2] = -a.fM[1][2];
  b.fM[2][0] = -a.fM[2][0];
  b.fM[2][1] = -a.fM[2][1];
  b.fM[2][2] = -a.fM[2][2];
  return b;
}

//______________________________________________________________________________
TVpMatrix3x3 operator + (const TVpMatrix3x3& a, const TVpMatrix3x3& b)
{
  // Binary plus

  TVpMatrix3x3 C;

  C.fM[0][0] = a.fM[0][0] + b.fM[0][0];
  C.fM[0][1] = a.fM[0][1] + b.fM[0][1];
  C.fM[0][2] = a.fM[0][2] + b.fM[0][2];
  C.fM[1][0] = a.fM[1][0] + b.fM[1][0];
  C.fM[1][1] = a.fM[1][1] + b.fM[1][1];
  C.fM[1][2] = a.fM[1][2] + b.fM[1][2];
  C.fM[2][0] = a.fM[2][0] + b.fM[2][0];
  C.fM[2][1] = a.fM[2][1] + b.fM[2][1];
  C.fM[2][2] = a.fM[2][2] + b.fM[2][2];
  return C;
}

//______________________________________________________________________________
TVpMatrix3x3 operator - (const TVpMatrix3x3& a, const TVpMatrix3x3& b)
{
  // Binary minus

  TVpMatrix3x3 C;
  
  C.fM[0][0] = a.fM[0][0] - b.fM[0][0];
  C.fM[0][1] = a.fM[0][1] - b.fM[0][1];
  C.fM[0][2] = a.fM[0][2] - b.fM[0][2];
  C.fM[1][0] = a.fM[1][0] - b.fM[1][0];
  C.fM[1][1] = a.fM[1][1] - b.fM[1][1];
  C.fM[1][2] = a.fM[1][2] - b.fM[1][2];
  C.fM[2][0] = a.fM[2][0] - b.fM[2][0];
  C.fM[2][1] = a.fM[2][1] - b.fM[2][1];
  C.fM[2][2] = a.fM[2][2] - b.fM[2][2];
  return C;
}

//______________________________________________________________________________
TVpMatrix3x3 operator *(const Double_t d, const TVpMatrix3x3& a)
{
  // Multiplication by a number

  TVpMatrix3x3 b;

  b.fM[0][0] = d * a.fM[0][0];
  b.fM[0][1] = d * a.fM[0][1];
  b.fM[0][2] = d * a.fM[0][2];
  b.fM[1][0] = d * a.fM[1][0];
  b.fM[1][1] = d * a.fM[1][1];
  b.fM[1][2] = d * a.fM[1][2];
  b.fM[2][0] = d * a.fM[2][0];
  b.fM[2][1] = d * a.fM[2][1];
  b.fM[2][2] = d * a.fM[2][2];
  return b;
}

//______________________________________________________________________________
TVpVector3 operator * (const TVpMatrix3x3& a, const TVpVector3& x)
{
  // Multiply a vector by a matrix

  TVpVector3 y;
  
  y.fR[0] = a.fM[0][0]*x.fR[0] + a.fM[0][1]*x.fR[1] + a.fM[0][2]*x.fR[2];
  y.fR[1] = a.fM[1][0]*x.fR[0] + a.fM[1][1]*x.fR[1] + a.fM[1][2]*x.fR[2];
  y.fR[2] = a.fM[2][0]*x.fR[0] + a.fM[2][1]*x.fR[1] + a.fM[2][2]*x.fR[2];
  return y;
}

//______________________________________________________________________________
TVpMatrix3x3 operator * (const TVpMatrix3x3& a, const TVpMatrix3x3& b)
{
  // Multiply 2 matrices

  Int_t i, j, k;
  const Double_t *pa, *pb;
  Double_t *pC, sum;
  TVpMatrix3x3 C;
  
  for (i = 0; i < 3; i++)
    for (pC = C.fM[i], j = 0; j < 3; j++)
      {
	for (pa = a.fM[i], pb = &b.fM[0][j], sum = 0.0, k = 0; k < 3; k++)
	  sum += pa[k] * pb[3*k];
	pC[j] = sum;
      }
  return C;
}

//______________________________________________________________________________
Int_t operator == (const TVpMatrix3x3& a, const TVpMatrix3x3& b)
{
  // Is equal
  
  Double_t eps = 1e-15;

  return &a == &b ||
    std::abs(b.fM[0][0] - a.fM[0][0]) < eps &&
    std::abs(b.fM[0][1] - a.fM[0][1]) < eps &&
    std::abs(b.fM[0][2] - a.fM[0][2]) < eps &&
    std::abs(b.fM[1][0] - a.fM[1][0]) < eps &&
    std::abs(b.fM[1][1] - a.fM[1][1]) < eps &&
    std::abs(b.fM[1][2] - a.fM[1][2]) < eps &&
    std::abs(b.fM[2][0] - a.fM[2][0]) < eps &&
    std::abs(b.fM[2][1] - a.fM[2][1]) < eps &&
    std::abs(b.fM[2][2] - a.fM[2][2]) < eps;
}

//______________________________________________________________________________
int operator != (const TVpMatrix3x3& a, const TVpMatrix3x3& b)
{
  // Is not equal

  return !(a == b);
}

//______________________________________________________________________________
TVpMatrix3x3 transpose(const TVpMatrix3x3& a)
{
  // Return transposed matrix

  TVpMatrix3x3 b;
  
  b.fM[0][0] = a.fM[0][0];
  b.fM[0][1] = a.fM[1][0];
  b.fM[0][2] = a.fM[2][0];
  b.fM[1][0] = a.fM[0][1];
  b.fM[1][1] = a.fM[1][1];
  b.fM[1][2] = a.fM[2][1];
  b.fM[2][0] = a.fM[0][2];
  b.fM[2][1] = a.fM[1][2];
  b.fM[2][2] = a.fM[2][2];
  return b;
}

//______________________________________________________________________________
TVpMatrix3x3 rotationMatrix(const TVpVector3& n, double angle)
{
  // Calculate the active rotation matrix when axis and rotation angle are
  // given.  For example: Let n points in the z-direction and the angle is
  // Pi/2.  Then the rotation matrix m transforms a vector pointing to the
  // x-direction into a vector pointing to the y-direction.  In other words,
  // e2 = m * e1
  //
  // Input parameters:
  // - n - axis, must be a unit vector, i.e. norm(n) == 1
  // - angle - rotation angle in radians
  //
  // Example:
  // root [] TVpMatrix3x3 m = rotationMatrix(TVpVector3(0,0,1), TMath::Pi()/2);
  // root [] m.PrintStatus();
  // <TVpMatrix3>
  //    1.57019580e-16  -1.00000000e+00   0.00000000e+00
  //    1.00000000e+00   1.57019580e-16   0.00000000e+00
  //    0.00000000e+00   0.00000000e+00   1.00000000e+00
  // </TVpMatrix3>
  // root [] (m * TVpVector3(1,0,0)).PrintStatus()
  // <TVpVector3>    1.57019580e-16   1.00000000e+00   0.00000000e+00 </TVpVector3>

  Double_t q[4];
  Double_t sn;
  TVpMatrix3x3 C;
  
  angle /= 2.0;
  q[0] = cos(angle);
  sn = sin(angle);
  q[1] = sn * n.fR[0];
  q[2] = sn * n.fR[1];
  q[3] = sn * n.fR[2];
  C.fM[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  C.fM[0][1] = 2*(q[1]*q[2]-q[0]*q[3]);
  C.fM[0][2] = 2*(q[1]*q[3]+q[0]*q[2]);
  C.fM[1][0] = 2*(q[1]*q[2]+q[0]*q[3]);
  C.fM[1][1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  C.fM[1][2] = 2*(q[2]*q[3]-q[0]*q[1]);
  C.fM[2][0] = 2*(q[1]*q[3]-q[0]*q[2]);
  C.fM[2][1] = 2*(q[2]*q[3]+q[0]*q[1]);
  C.fM[2][2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  return C;
}

//______________________________________________________________________________
TVpMatrix3x3 rotationMatrix(const TVpVector3& oldDir, const TVpVector3& newDir)
{
  // Calculate the active rotation matrix which rotates the oldDir vector into
  // the newDir vector.
  //
  // Input parameters:
  // - oldDir - old direction vector
  // - newDir - new direction vector
  //
  // Example:
  // root [] (rotationMatrix(TVpVector3(1,0,0), TVpVector3(0,1,0))).PrintStatus()
  // <TVpMatrix3>
  //    1.57019580e-16  -1.00000000e+00   0.00000000e+00
  //    1.00000000e+00   1.57019580e-16   0.00000000e+00
  //    0.00000000e+00   0.00000000e+00   1.00000000e+00
  // </TVpMatrix3>

  
  Double_t eps = 1e-12;
  Double_t rotAngle = acos(oldDir | newDir);
  TVpVector3 axis = oldDir * newDir;
  if (norm(axis) < eps)
    {
      // The vectors are almost parallel or antiparallel and the rotation axis
      // is not well defined.  Use the x-axis if oldDir is not parallel with
      // it or the y-axis otherwise.
      if (oldDir.IsEqual(TVpVector3(1,0,0), 1e-4))
	axis.Set(0,1,0);
      else
	axis.Set(1,0,0);
    }
  else
    axis = normalize(axis);
  return rotationMatrix(axis, rotAngle);
}

//______________________________________________________________________________
void TVpMatrix3x3::PrintStatus(std::ostream &out) const
{
  // Print the matrix status to a stream.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpMatrix3x3 m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
  // root [] m1.PrintStatus();
  // <TVpMatrix3>
  //    1.00000000e+00   2.00000000e+00   3.00000000e+00
  //    4.00000000e+00   5.00000000e+00   6.00000000e+00
  //    7.00000000e+00   8.00000000e+00   9.00000000e+00
  // </TVpMatrix3>

  out << "<TVpMatrix3>\n"
      << "$Id: TVpMatrix3x3.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << *this
      << "</TVpMatrix3>\n"; 
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& s, const TVpMatrix3x3& a)
{
  // Write matrix elements in a tabular form to a stream.  Operator
  // overloading for streams does not work in a ROOT interactive session.
  //
  // Input parameters:
  // - s - ouput stream (default = cout)
  // - a - the matrix
  //
  // Example:
  // root [] TVpMatrix3x3 m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
  // root [] cout << m1;
  //    1.00000000e+00   2.00000000e+00   3.00000000e+00
  //    4.00000000e+00   5.00000000e+00   6.00000000e+00
  //    7.00000000e+00   8.00000000e+00   9.00000000e+00

  s << std::scientific << std::setprecision(8);
  for (Int_t i = 0; i < 3; ++i)
    for (Int_t j = 0; j < 3; ++j)
      {
	s << std::setw(17) << a.fM[i][j];
	if (j == 2)
	  s << '\n';
      }
  return s;
}
