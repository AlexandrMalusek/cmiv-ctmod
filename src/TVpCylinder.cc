//______________________________________________________________________________
//
// TVpCylinder defines a homogenous cylinder.
//
// Local coordinate system: The origin is in the centre of the lower base, the
// z-axis is identical with the cylinder axis; see the figure:
//Begin_Html
/*
<img src="png/solid_3.png">
*/
//End_Html
// 
// Example:
// gSystem->Load("libRCTmod.so");
// TVpCylinder *cyl3 = new TVpCylinder(
//   "cyl3",    // name
//   3,         // unique solid index
//   8,         // radius in cm
//   16);       // height in cm
// cyl3->SetActiveTranslation(new TVpVector3(originX, originY, originZ));
// cyl3->SetActiveRotation(&rotMat);
// cyl3->SetMaterial(matWater);
//
// See TVpSolid for the description of the translation and rotation matrices.

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "TVpCylinder.h"

ClassImp(TVpCylinder)

//______________________________________________________________________________
TVpCylinder::TVpCylinder()
  : TVpSolid("cylinder", 0)
{
  // Default constructor.  Initialize name to "cylinder" and radius and height
  // to 0.

  fRadius = fRadius2 = fHeight = 0.0;
}

//______________________________________________________________________________
TVpCylinder::TVpCylinder(const Char_t *name, Int_t index, Double_t radius, Double_t height)
  : TVpSolid(name, index)
{
  // Constructor.
  //
  // Input parameters:
  // - name - user defined descriptive name
  // - index - solid index.  Must be unique in the geometry (range = 0, ...)
  // - radius - radius of the base in cm
  // - height - height of the cylinder in cm

  fRadius = radius;
  fHeight = height;
  fRadius2 = fRadius*fRadius;
}

//______________________________________________________________________________
Int_t TVpCylinder::RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				     const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the cylinder and
  // a line segment described by the point "pos", direction "dir", and length
  // "l".  Return 1 if the intersection exists and 0 otherwise.  Suppose the
  // point "pos" is inside the cylinder.
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
  // Test bases first, then the cylindrical surface.
  //
  // Example:
  // root [] gSystem->Load("libRCTmod.so");
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << cyl->RayIntersectionIn(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 10.5154
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //     7.071068     7.071068     7.071068
  
  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  Double_t x0, x1, a, b, c;
  
#ifdef CYLINDER_DEBUG
  if (IsInside(pos) == kFALSE)
    {
      std::cerr << "TVpCylinder::RayIntersectionIn: The point r is outside the cylinder "
		<< fName << "\n.";
      exit(1);
    }
#endif
  
  if (u[2] < 0.0)
    t = -r[2]/u[2];			// lower base
  else
    t=(fHeight - r[2])/u[2];     	// upper base
  x0 = r[0] + t*u[0];
  x1 = r[1] + t*u[1];
  if (x0*x0 + x1*x1 < fRadius2)
    {
      if (t < l)
	return kTRUE;
      return kFALSE;
    }
  
  a = u[0]*u[0] + u[1]*u[1];	        // lateral surface
  b = r[0]*u[0] + r[1]*u[1];
  c = r[0]*r[0] + r[1]*r[1] - fRadius2;
  t = (-b + sqrt(b*b - a*c)) / a;       // D>=0 because r is inside the cylinder
  if (t < l)
    return kTRUE;
  return kFALSE;
}


//______________________________________________________________________________
Int_t TVpCylinder::RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
				      const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the box and a
  // line segment described by the point "pos", direction "dir", and length
  // "l".  Return 1 if the intersection exists and 0 otherwise.  Suppose the
  // point "pos" is outside the box.
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
  // Test bases first, then the cylindrical surface.
  //
  // Example:
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // root [] TVpVector3 pos(1, 1, -2);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << cyl->RayIntersectionOut(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 3.4641
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //     3.000000     3.000000     0.000000

  Double_t *r = pos.Get();
  Double_t *u = dir.Get();
  Double_t x0, x1;
  Double_t a, b, c, D;
  
  if (u[2] > 0.0)
    t = -r[2] / u[2];			// lower base
  else
    t = (fHeight - r[2]) / u[2];     	// upper base
  if (t > l)
    return kFALSE;

  x0 = r[0] + t*u[0];
  x1 = r[1] + t*u[1];
  if (t > 0.0 && x0*x0 + x1*x1 < fRadius2)
    return kTRUE;
  a = u[0]*u[0]+u[1]*u[1];			// lateral surface
  b = r[0]*u[0]+r[1]*u[1];
  c = r[0]*r[0]+r[1]*r[1]-fRadius2;
  if ((D = b*b-a*c) <= 0.0)
    return kFALSE;
  t = (-b - sqrt(D)) / a;
  if (t < 0.0 || t > l)
    return kFALSE;
  x0 = r[2] + t*u[2];
  if (0.0 < x0 && x0 < fHeight)
    return kTRUE;
  return kFALSE;
}


//______________________________________________________________________________
Int_t TVpCylinder::IsInside(TVpVector3& pos) const
{
  // Return 1 if the point "pos" is inside the cylinder (the cylinder's
  // boundary is included) and 0 otherwise.
  //
  // Input parameters:
  // - pos - tested point 
  //
  // Example:
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // root [] TVpVector3 pos1(1, 2, 3);     // this point is inside the cylinder
  // root [] cout << cyl->IsInside(pos1) << endl;
  // 1
  // root [] TVpVector3 pos2(-2, -2, -2);  // this point is outside the cylinder
  // root [] cout << cyl->IsInside(pos2) << endl;
  // 0
  // root [] TVpVector3 pos3(1, 2, 20);    // this point is on the boundary
  // root [] cout << cyl->IsInside(pos3) << endl;
  // 1

  Double_t *r = pos.Get();
  return 0.0 <= r[2] && r[2] <= fHeight && r[0]*r[0] + r[1]*r[1] <= fRadius2;
}

//______________________________________________________________________________
Double_t TVpCylinder::GetVolume() const
{
  // Return the box volume in cm^3.
  //
  // Example:
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // root [] cout << cyl->GetVolume() << endl;
  // 6283.19

  return M_PI * fRadius2 * fHeight;
}

//______________________________________________________________________________
void TVpCylinder::PrintStatus(std::ostream &out) const
{
  // Print the object status
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // root [18] cyl->PrintStatus()
  // <TVpCylinder>
  // Radius: 10 cm
  // Height: 20 cm
  // <TVpSolid>
  // Name:   cylinder
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
  // </TVpCylinder>
  
  out << "<TVpCylinder>\n"
      << "$Id: TVpCylinder.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Radius: " << fRadius << " cm\n"
      << "Height: " << fHeight << " cm\n";
  TVpSolid::PrintStatus(out);
  out << "</TVpCylinder>" << std::endl;
}

#include "TPolyLine3D.h"

//______________________________________________________________________________
void TVpCylinder::Draw(Int_t color) const
{
  // Draw 2 circles and the axis.
  // Input parameters:
  // - color - line color of main circles (default = solid_index+12)
  //
  // Example:
  // root [] TVpCylinder *cyl = new TVpCylinder("cylinder", 1, 10, 20);
  // cyl->Draw(1);  // the canvas is automatically created

  if (color == -1)
    color = fIndex+12;

  const Int_t nSegments = 100; // number of line segments per one circle
  Double_t dAlpha = 2 * M_PI / nSegments;
  TPolyLine3D *baseU = new TPolyLine3D(nSegments+1);
  TPolyLine3D *baseL = new TPolyLine3D(nSegments+1);
  TVpVector3 pointLocal;     // local coordinates
  TVpVector3 pU;

  for (Int_t i = 0; i <= nSegments; i++)
    {
      Double_t alpha = i * dAlpha;
      Double_t rca = fRadius * cos(alpha);
      Double_t rsa = fRadius * sin(alpha);
      // Lower base
      pointLocal.Set(rca, rsa, 0);
      pU =  PosLocToUni(pointLocal);
      baseL->SetPoint(i, pU.GetX(), pU.GetY(), pU.GetZ());
      // Upper base
      pointLocal.Set(rca, rsa, fHeight);
      pU =  PosLocToUni(pointLocal);
      baseU->SetPoint(i, pU.GetX(), pU.GetY(), pU.GetZ());
    }
  baseL->SetLineColor(color);
  baseU->SetLineColor(color);
  baseL->Draw();
  baseU->Draw();

  TPolyLine3D *axis = new TPolyLine3D(2);
  pU =  PosLocToUni(TVpVector3(0, 0, 0));
  axis->SetPoint(0, pU.GetX(), pU.GetY(), pU.GetZ());
  pU =  PosLocToUni(TVpVector3(0,0,fHeight));
  axis->SetPoint(1, pU.GetX(), pU.GetY(), pU.GetZ());
  axis->SetLineColor(color);
  axis->Draw();

  DrawLcsAxes(fRadius/10, fRadius/10, fHeight/10);
}
