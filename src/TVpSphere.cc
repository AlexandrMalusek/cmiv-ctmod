//______________________________________________________________________________
//
// TVpSphere defines a homogenous sphere.
//
// Local coordinate system: The origin is in the center of the sphere, see the
// figure:
//Begin_Html
/*
<img src="png/solid_1.png">
*/
//End_Html
// 
// Example:
// gSystem->Load("libRCTmod.so");
// TVpSphere *sph1 = new TVpSphere(
//   "sph1",    // name
//   1,         // unique solid index
//   10);       // radius in cm
// sph1->SetActiveTranslation(new TVpVector3(originX, originY, originZ));
// sph1->SetActiveRotation(&rotMat);
// sph1->SetMaterial(matWater);
//
// See TVpSolid for the description of the translation and rotation matrices.

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TVpSphere.h"

ClassImp(TVpSphere)

//______________________________________________________________________________
TVpSphere::TVpSphere()
  : TVpSolid("sphere", 0)
{
  // Default constructor.  Initialize name to "sphere" and set radius to 0.

  fRadius = fRadius2 = 0.0;
}

//______________________________________________________________________________
TVpSphere::TVpSphere(const Char_t *name, Int_t index, Double_t radius)
  : TVpSolid(name, index)
{
  // Constructor.
  //
  // Input parameters:
  // - name - user defined descriptive name
  // - index - solid index.  Must be unique in the geometry (range = 0, ...)
  // - radius - radius in cm
  
  fRadius = radius;
  fRadius2 = fRadius*fRadius;
}

//______________________________________________________________________________
Int_t TVpSphere::RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				   const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the sphere and a
  // line segment described by the point "pos", direction "dir", and length
  // "l".  Return 1 if the intersection exists and 0 otherwise.  Suppose the
  // point "pos" is inside the sphere.
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
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << sphere->RayIntersectionIn(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 8.26795
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //     5.773503     5.773503     5.773503
  // root [] cout << norm(new) << endl;
  // 10
  
  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  double b,c;
  
  #ifdef SPHERE_DEBUG
  if (IsInside(pos) == kFALSE)
    {
      std::cerr << "TVpSphere::RayIntersectionIn: The point r is outside the sphere "
		<< fName << "\n.";
      exit(1);
    }
  #endif
  
  b = r[0]*u[0] + r[1]*u[1] + r[2]*u[2];
  c = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  t = -b + sqrt(b*b - c + fRadius2);  // D>=0 because the point r is inside the sphere
  if (t < l)
    return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________
Int_t TVpSphere::RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				    const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the sphere and a
  // line segment described by the point "pos", direction "dir", and length
  // "l".  Return 1 if the intersection exists and 0 otherwise.  Suppose the
  // point "pos" is outside the sphere.
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
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] TVpVector3 pos(-20, -20, -20);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << sphere->RayIntersectionOut(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 24.641
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //    -5.773503    -5.773503    -5.773503
  // root [] cout << norm(new) << endl;
  // 10

  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  Double_t b, c, D;

  b = r[0]*u[0] + r[1]*u[1] + r[2]*u[2];
  c = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  if ((D = b*b - c + fRadius2) <= 0.0)
    return kFALSE;
  if ((t = -b - sqrt(D)) > 0.0 && t < l)
    return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________
Int_t TVpSphere::IsInside(TVpVector3& pos) const
{
  // Return 1 if the point "pos" is inside the sphere (the sphere's boundary
  // is included) and 0 otherwise.
  //
  // Input parameters:
  // - pos - tested point
  //
  // Example:
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] TVpVector3 pos1(1, 2, 3);        // this point is inside the sphere
  // root [] cout << sphere->IsInside(pos1) << endl;
  // 1
  // root [] TVpVector3 pos2(-20, -20, -20);  // this point is outside the sphere
  // root [] cout << sphere->IsInside(pos2) << endl;
  // 0
  // root [] TVpVector3 pos3(10, 0, 0);       // this point is on the boundary
  // root [] cout << sphere->IsInside(pos3) << endl;
  // 1

  Double_t *r = pos.Get();
  return r[0]*r[0] + r[1]*r[1] + r[2]*r[2] <= fRadius2;
}

//______________________________________________________________________________
Double_t TVpSphere::GetVolume() const
{
  // Return the volume in cm^3.
  //
  // Example:
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] cout << sphere->GetVolume() << endl;
  // 4188.79
  
  return 4.0 / 3.0 * M_PI * fRadius2 * fRadius;
}

//______________________________________________________________________________
void TVpSphere::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  //  Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] sphere->PrintStatus()
  // <TVpSphere>
  // Radius: 10 cm
  // <TVpSolid>
  // Name:   sphere
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
  // </TVpSphere>

  out << "<TVpSphere>\n"
      << "$Id: TVpSphere.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Radius: " << fRadius << " cm\n";
  TVpSolid::PrintStatus(out);
  out << "</TVpSphere>" << std::endl;
}

#include "TPolyLine3D.h"

//______________________________________________________________________________
void TVpSphere::Draw(Int_t color) const
{
  // Draw 3 main circles.
  //
  // Input parameters:
  // - color - line color of main circles (default = solid_index+20)
  //
  // Example:
  // root [] TVpSphere *sphere = new TVpSphere("sphere", 1, 10);
  // root [] sphere->Draw(1);  // the canvas is automatically created

  if (color == -1)
    color = fIndex+20;

  const Int_t nSegments = 100; // number of line segments per one circle
  Double_t dAlpha = 2 * M_PI / nSegments;
  TPolyLine3D *xcircle = new TPolyLine3D(nSegments+1);
  TPolyLine3D *ycircle = new TPolyLine3D(nSegments+1);
  TPolyLine3D *zcircle = new TPolyLine3D(nSegments+1);

  for (Int_t i = 0; i <= nSegments; i++)
    {
      Double_t alpha = i * dAlpha;
      Double_t rca = fRadius * cos(alpha);
      Double_t rsa = fRadius * sin(alpha);
      xcircle->SetPoint(i, fTraVecL2u.GetX(), fTraVecL2u.GetY() + rca,
			fTraVecL2u.GetZ() + rsa);
      ycircle->SetPoint(i, fTraVecL2u.GetX() + rca, fTraVecL2u.GetY(),
			fTraVecL2u.GetZ() + rsa);
      zcircle->SetPoint(i, fTraVecL2u.GetX() + rca, fTraVecL2u.GetY() + rsa,
			fTraVecL2u.GetZ());
    }
  xcircle->SetLineColor(color);
  ycircle->SetLineColor(color);
  zcircle->SetLineColor(color);
  xcircle->Draw();
  ycircle->Draw();
  zcircle->Draw();

  DrawLcsAxes(fRadius/10, fRadius/10, fRadius/10);
}
