//______________________________________________________________________________
//
// TVpEllipsoid defines a homogenous ellipsoid.
//
// Local coordinate system: The origin is in the center of the ellipsoid,
// principal semiaxes are identical with main axes of the ellipsoid; see the
// figure:
//Begin_Html
/*
<img src="png/solid_2.png">
*/
//End_Html
// 
// Example:
// gSystem->Load("libRCTmod.so");
// TVpEllipsoid *ell2 = new TVpEllisoid(
//   "ell2",    // name
//   2,         // unique solid index
//   10,        // principal semiaxis length in cm (x-axis)
//   7,         // principal semiaxis length in cm (y-axis)
//   5);        // principal semiaxis length in cm (z-axis)
// ell2->SetActiveTranslation(new TVpVector3(originX, originY, originZ));
// ell2->SetActiveRotation(&rotMat);
// ell2->SetMaterial(matWater);
//
// See TVpSolid for the description of the translation and rotation matrices.

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TVpEllipsoid.h"

ClassImp(TVpEllipsoid)

//______________________________________________________________________________
TVpEllipsoid::TVpEllipsoid()
  : TVpSolid("ellipsoid", 0)
{
  // Default constructor.  Initialize name to "ellipsoid" and principal
  // semiaxis lengths to 0.

  fa = fb = fc = fia2 = fib2 = fic2 = 0.0;
}

//______________________________________________________________________________
TVpEllipsoid::TVpEllipsoid(const Char_t *name, Int_t index,
			   Double_t a, Double_t b, Double_t c)
  : TVpSolid(name, index)
{
  // Constructor.
  //
  // Input parameters:
  // - name - user defined descriptive name
  // - index - solid index.  Must be unique in the geometry (range = 0, ...)
  // - a - principal semiaxis length (x-axis) in cm
  // - b - principal semiaxis length (y-axis) in cm
  // - c - principal semiaxis length (z-axis) in cm

  fa = a;
  fb = b;
  fc = c;

  fia2 = 1/(a*a);
  fib2 = 1/(b*b);
  fic2 = 1/(c*c);
}

//______________________________________________________________________________
Int_t TVpEllipsoid::RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				      const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the ellipsoid
  // and a line segment described by the point "pos", direction "dir", and
  // length "l".  Return 1 if the intersection exists and 0 otherwise.
  // Suppose the point "pos" is inside the ellipsoid.
  //
  // Input parameters:
  // - pos - starting point of the line segment in LCS
  // - dir - direction of the line segment in LCS
  // - l - length of the line segment in cm
  //
  // Output parameters:
  // - t - the distance in cm to the first intersection from the point "pos".
  //       If the intersection does not exist then "t" is undefined.
  //
  // Method:
  // Intersections are found by solving a quadratic equation.
  //
  // Example:
  // root [] gSystem->Load("libRCTmod.so");
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << elli->RayIntersectionIn(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 11.5739
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //     7.682213     7.682213     7.682213

  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  
#ifdef ELLIPSOID_DEBUG
  if (IsInside(pos) == kFALSE)
    {
      std::cerr << "TVpEllipsoid::RayIntersectionIn: The point r is outside the ellipsoid "
		<< fName << "\n.";
      exit(1);
    }
#endif

  // at^2 + 2bt + c = 0
  Double_t a = u[0]*u[0]*fia2 + u[1]*u[1]*fib2 + u[2]*u[2]*fic2;
  Double_t b = r[0]*u[0]*fia2 + r[1]*u[1]*fib2 + r[2]*u[2]*fic2;
  Double_t c = r[0]*r[0]*fia2 + r[1]*r[1]*fib2 + r[2]*r[2]*fic2 - 1.0;
  t = (-b + sqrt(b*b - a*c))/a;    // D>=0 because the point r is inside the ellipsoid
  if (t < l)
    return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________
Int_t TVpEllipsoid::RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				       const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the ellipsoid
  // and a line segment described by the point "pos", direction "dir", and
  // length "l".  Return 1 if the intersection exists and 0 otherwise.
  // Suppose the point "pos" is outside the ellipsoid.
  //
  // Input parameters:
  // - pos - starting point of the line segment in LCS
  // - dir - direction of the line segment in LCS
  // - l - length of the line segment in cm
  //
  // Output parameters:
  // - t - the distance in cm to the first intersection from the point "pos".
  //       If the intersection does not exist then "t" is undefined.
  //
  // Method:
  // Intersections are found by solving a quadratic equation.
  //
  // Example:
  // root [] gSystem->Load("libRCTmod.so");
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] TVpVector3 pos(-20, -20, -20);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << elli->RayIntersectionOut(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 21.335
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //    -7.682213    -7.682213    -7.682213

  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  Double_t D;
  
  // at^2 + 2bt + c = 0
  Double_t a = u[0]*u[0]*fia2 + u[1]*u[1]*fib2 + u[2]*u[2]*fic2;
  Double_t b = r[0]*u[0]*fia2 + r[1]*u[1]*fib2 + r[2]*u[2]*fic2;
  Double_t c = r[0]*r[0]*fia2 + r[1]*r[1]*fib2 + r[2]*r[2]*fic2 - 1.0;
  if ((D = b*b - a*c) <= 0.0)
    return kFALSE;
  if ((t = (-b - sqrt(D))/a) > 0.0 && t < l)
    return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________
Int_t TVpEllipsoid::IsInside(TVpVector3& pos) const
{
  // Return 1 if the point "pos" is inside the ellipsoid (the ellipsoid's
  // boundary is included) and 0 otherwise.
  //
  // Input parameters:
  // - pos - tested point
  //
  // Example:
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] TVpVector3 pos1(1, 2, 3);          // this point is inside the ellipsoid
  // root [11] cout << elli->IsInside(pos1) << endl;
  // 1
  // root [12] TVpVector3 pos2(-20, -20, -20);  // this point is outside the ellipsoid
  // root [13] cout << elli->IsInside(pos2) << endl;
  // 0
  // root [14] TVpVector3 pos3(10, 0, 0);       // this point is on the boundary 
  // root [16] cout << elli->IsInside(pos3) << endl;
  // 1

  Double_t *r = pos.Get();
  Double_t d = r[0]*r[0]*fia2 + r[1]*r[1]*fib2 + r[2]*r[2]*fic2;
  return d <= 1.0;
}

//______________________________________________________________________________
Double_t TVpEllipsoid::GetVolume() const
{
  // Return the volume in cm^3.
  //
  // Example:
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] cout << elli->GetVolume() << endl;
  // 12566.4

  return 4.0 / 3.0 * M_PI * fa * fb * fc;
}

//______________________________________________________________________________
void TVpEllipsoid::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  //  Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] elli->PrintStatus()
  // <TVpEllipsoid>
  // Radius [cm]:    10      15      20
  // <TVpSolid>
  // Name:   ellipsoid
  // Index:  1
  // Material:       Undefined or voxel array
  // Position of origin in Universe coordinates [cm]
  //    0.00000000e+00   0.00000000e+00   0.00000000e+00
  // Rotation matrix:
  //    1.00000000e+00   0.00000000e+00   0.00000000e+00
  //    0.00000000e+00   1.00000000e+00   0.00000000e+00
  //    0.00000000e+00   0.00000000e+00   1.00000000e+00
  // Overlap solid indices:  None
  // Base solid indices:     None
  // </TVpSolid>
  // </TVpEllipsoid>

  out << "<TVpEllipsoid>\n"
      << "$Id: TVpEllipsoid.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Radius [cm]:\t" << fa << '\t' << fb << '\t' << fc << '\n';
  TVpSolid::PrintStatus(out);
  out << "</TVpEllipsoid>"  << std::endl;
}

#include "TPolyLine3D.h"

//______________________________________________________________________________
void TVpEllipsoid::Draw(Int_t color) const
{
  // Draw 3 main ellipses.
  //
  // Input parameters:
  // - color - line color of main circles (default = solid_index+20)
  //
  // Example: 
  // root [] TVpEllipsoid *elli = new TVpEllipsoid("ellipsoid", 1, 10, 15, 20);
  // root [] elli->Draw(1);  // the canvas is automatically created

  if (color == -1)
    color = fIndex+20;

  const Int_t nSegments = 100; // number of line segments per one ellipse
  Double_t dAlpha = 2 * M_PI / nSegments;
  TPolyLine3D *xellipse = new TPolyLine3D(nSegments+1);
  TPolyLine3D *yellipse = new TPolyLine3D(nSegments+1);
  TPolyLine3D *zellipse = new TPolyLine3D(nSegments+1);
  TVpVector3 lp;  // point in local coordinates
  TVpVector3 gp;  // point in global coordinates
  
  // x-y, zellipse
  for (Int_t i = 0; i <= nSegments; i++)
    {
      Double_t alpha = i * dAlpha;
      Double_t rca = fa * cos(alpha);
      Double_t rsa = fb * sin(alpha);
      lp.Set(rca, rsa, 0);
      gp = PosLocToUni(lp);
      zellipse->SetPoint(i, gp.GetX(), gp.GetY(), gp.GetZ());
    }

  // x-z, yellipse
  for (Int_t i = 0; i <= nSegments; i++)
    {
      Double_t alpha = i * dAlpha;
      Double_t rca = fa * cos(alpha);
      Double_t rsa = fc * sin(alpha);
      lp.Set(rca, 0, rsa);
      gp = PosLocToUni(lp);
      yellipse->SetPoint(i, gp.GetX(), gp.GetY(), gp.GetZ());
    }

  // y-z, xellipse
  for (Int_t i = 0; i <= nSegments; i++)
    {
      Double_t alpha = i * dAlpha;
      Double_t rca = fb * cos(alpha);
      Double_t rsa = fc * sin(alpha);
      lp.Set(0, rca, rsa);
      gp = PosLocToUni(lp);
      xellipse->SetPoint(i, gp.GetX(), gp.GetY(), gp.GetZ());
    }

  xellipse->SetLineColor(color);
  yellipse->SetLineColor(color);
  zellipse->SetLineColor(color);
  xellipse->Draw();
  yellipse->Draw();
  zellipse->Draw();

  DrawLcsAxes(fa/10, fb/10, fc/10);
}
