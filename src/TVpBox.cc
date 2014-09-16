//______________________________________________________________________________
//
// TVpBox defines a homogenous box.
//
// Local coordinate system (LCS): The origin is in the corner, edges coincide
// with the coordinate system axes; see the figure:
//Begin_Html
/*
<img src="png/solid_4.png">
*/
//End_Html
// 
// Example:
// gSystem->Load("libRCTmod.so");
// TVpBox *box4 = new TVpBox(
//   "box4",    // name
//   4,         // unique solid index
//   3,         // edge length in cm (x-axis)
//   4,         // edge length in cm (y-axis)
//   5);        // edge length in cm (z-axis)
// box4->SetActiveTranslation(new TVpVector3(originX, originY, originZ));
// box4->SetActiveRotation(&rotMat);
// box4->SetMaterial(matWater);
//
// See TVpSolid for the description of the translation and rotation matrices.

#include "TVpBox.h"

ClassImp(TVpBox)

//______________________________________________________________________________
TVpBox::TVpBox()
  : TVpSolid("box", 0)
{
  // Default constructor.  Initialize name to "box" and sizes to 0.
  
  fSizeX = fSizeY = fSizeZ = 0.0;
}

//______________________________________________________________________________
TVpBox::TVpBox(const Char_t *name, Int_t index, Double_t sizeX, Double_t sizeY, Double_t sizeZ)
  : TVpSolid(name, index)
{
  // Constructor.
  //
  // Input parameters:
  // - name - user defined descriptive name
  // - index - solid index.  Must be unique in the geometry (range = 0, ...)
  // - sizeX - box edge length (x-axis) in cm
  // - sizeY - box edge length (y-axis) in cm
  // - sizeZ - box edge length (z-axis) in cm
  
  fSizeX = sizeX;
  fSizeY = sizeY;
  fSizeZ = sizeZ;
}

//______________________________________________________________________________
Int_t TVpBox::RayIntersectionIn(TVpVector3& pos, TVpVector3& dir,
				const Double_t l, Double_t& t) const
{
  // Calculate the distance to the first intersection between the box and a
  // line segment described by the point "pos", direction "dir", and length
  // "l".  Return 1 if the intersection exists and 0 otherwise.  Suppose the
  // point "pos" is inside the box.
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
  // Directions are divided into 8 octants.  An itersection with 3 planes
  // corresponding to the proper octant is calculated.
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 100.0;
  // root [] Double_t t;
  // root [] cout << box->RayIntersectionIn(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 15.5885
  // root [] TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //    10.000000    10.000000    10.000000

  Double_t x, y, z;
  Double_t *u = dir.fR;
  Double_t *r = pos.fR;

  if (u[0] > 0.0)
    if (u[1] > 0.0)
      if (u[2] > 0.0)
	{
	  // 1. octant
	  // test S1
	  t = (fSizeX - r[0]) / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S2
	  t = (fSizeY-r[1])/u[1];
	  x = r[0] + t * u[0];						
	  z = r[2] + t * u[2];						
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S5
	  t = (fSizeZ - r[2]) / u[2];
	  return (t < l) ? 1 : 0;
	}
      else
	{
	  // 5. octant
	  // test S1
	  t = (fSizeX - r[0]) / u[0];
	  y = r[1] + t * u[1];						
	  z = r[2] + t * u[2];						
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)			
	    return (t<l) ? 1 : 0;

	  // test S2
	  t =(fSizeY - r[1]) / u[1];
	  x = r[0] + t * u[0];					       
	  z = r[2] + t * u[2];						
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)			
	    return (t<l) ? 1 : 0;                                 	
	  
	  // test S6
	  t = -r[2] / u[2];
	  return (t<l) ? 1 : 0;
	}
    else
      if (u[2] > 0.0)
	{
	  // 4. octant
	  // test S1
	  t = (fSizeX - r[0]) / u[0];
	  y = r[1] + t * u[1];						
	  z = r[2] + t * u[2];						
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)			
	    return (t < l) ? 1 : 0;                                 	
	  
	  // testS4
	  t = -r[1] / u[1];
	  x = r[0] + t * u[0];
	  z = r[2] + t * u[2];
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S5
	  t = (fSizeZ - r[2]) / u[2];
	  return (t < l) ? 1 : 0;
	}
      else
	{
	  // 8. octant
	  // test S1
	  t = (fSizeX - r[0]) / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0 < y && y < fSizeY && 0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S4
	  t = -r[1] / u[1];
	  x = r[0] + t * u[0];
	  z = r[2] + t * u[2];
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;
	  
	  // test S6
	  t = -r[2] / u[2];
	  return (t < l) ? 1 : 0;
	}
  else
    if (u[1] > 0.0)
      if (u[2]>0.0)
	{
	  // 2. octant
	  // test S2
	  t = (fSizeY - r[1]) / u[1];
	  x = r[0] + t * u[0];						
	  z = r[2] + t * u[2];						
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)			
	    return (t < l) ? 1 : 0;                                 	
	  
	  // test S3
	  t = -r[0] / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S5
	  t = (fSizeZ - r[2]) / u[2];
	  return (t < l) ? 1 : 0;
	}
      else
	{
	  // 6. octant
	  // test S2
	  t = (fSizeY - r[1]) / u[1];
	  x = r[0] + t * u[0];						
	  z = r[2] + t * u[2];						
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)			
	    return (t < l) ? 1 : 0;                                  	

	  // test S3
	  t = -r[0] / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S6
	  t = -r[2] / u[2];
	  return (t < l) ? 1 : 0;
	}
    else
      if (u[2] > 0.0)
	{
	  // 3. octant
	  // test S3
	  t = -r[0] / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;                                  	

	  // test S4
	  t  = -r[1] / u[1];
	  x = r[0] + t * u[0];
	  z = r[2] + t * u[2];
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S5
	  t = (fSizeZ - r[2]) / u[2];
	  return (t < l) ? 1 : 0;
	}
      else
	{
	  // 8. octant
	  // test S3
	  t = -r[0] / u[0];
	  y = r[1] + t * u[1];
	  z = r[2] + t * u[2];
	  if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;

	  // test S4
	  t = -r[1] / u[1];
	  x = r[0] + t * u[0];
	  z = r[2] + t * u[2];
	  if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)
	    return (t < l) ? 1 : 0;
	  
	  // test S6
	  t = -r[2] / u[2];
	  return (t < l) ? 1 : 0;
	}
}


#define test_S1                                                          \
  if ((t = fSizeX - r[0]) < 0.0)                                         \
    {                                                                    \
      t /= u[0];                                                         \
      if (t > l)                                                         \
        return 0;                                                        \
      y = r[1] + t * u[1];                                               \
      z = r[2] + t * u[2];                                               \
      if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)                \
        return 1;                                                        \
    }

#define test_S2                                                          \
  if ((t = fSizeY - r[1]) < 0.0)                                         \
    {                                                                    \
      t /= u[1];                                                         \
      if (t > l)                                                         \
        return 0;                                                        \
      x = r[0] + t * u[0];                                               \
      z = r[2] + t * u[2];                                               \
      if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)                \
        return 1;                                                        \
    }

#define test_S3                                                          \
  if (r[0] < 0.0)                                                        \
    {                                                                    \
      t = -r[0] / u[0];                                                  \
      if (t > l)                                                         \
        return 0;                                                        \
      y = r[1] + t * u[1];                                               \
      z = r[2] + t * u[2];                                               \
      if (0.0 < y && y < fSizeY && 0.0 < z && z < fSizeZ)                \
        return 1;                                                        \
    }

#define test_S4                                                          \
  if (r[1] < 0.0)                                                        \
    {                                                                    \
      t = -r[1] / u[1];                                                  \
      if (t > l)                                                         \
        return 0;                                                        \
      x = r[0] + t * u[0];                                               \
      z = r[2] + t * u[2];                                               \
      if (0.0 < x && x < fSizeX && 0.0 < z && z < fSizeZ)                \
         return 1;                                                       \
    }
                                                                         \
#define test_S5                                                          \
  if ((t = fSizeZ - r[2]) < 0.0)                                         \
    {                                                                    \
      t /= u[2];                                                         \
      if (t > l)                                                         \
        return 0;                                                        \
      x = r[0] + t * u[0];                                               \
      y = r[1] + t * u[1];                                               \
      if (0.0 < x && x < fSizeX && 0.0 < y && y < fSizeY)                \
        return 1;                                                        \
    }

#define test_S6                                                          \
  if (r[2] < 0.0)                                                        \
    {                                                                    \
      t = -r[2] / u[2];                                                  \
      if (t > l)                                                         \
        return 0;                                                        \
      x = r[0] + t * u[0];                                               \
      y = r[1] + t * u[1];                                               \
      if (0.0 < x && x < fSizeX && 0.0 < y && y < fSizeY)                \
        return 1;                                                        \
    }

//______________________________________________________________________________
Int_t TVpBox::RayIntersectionOut(TVpVector3& pos, TVpVector3& dir,
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
  // Directions are divided into 8 octants.  An itersection with 3 planes
  // corresponding to the proper octant is calculated.
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] TVpVector3 pos(2, -1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] Double_t l = 10.0;
  // root [] Double_t t;
  // root [] cout << box->RayIntersectionOut(pos, dir, l, t) << endl;
  // 1
  // root [] cout << t << endl;
  // 1.73205081e+00
  // TVpVector3 new = pos + t * dir;  // calculate the intersection
  // root [] new.PrintValue();
  //     3.000000     0.000000     2.000000

  Double_t x, y, z;
  Double_t *u = dir.fR;
  Double_t *r = pos.fR;

  if (u[0] > 0.0)
    if (u[1] > 0.0)
      if (u[2] > 0.0)
	{
	  // 1. octant
	  test_S3
	  test_S4
	  test_S6
	  return 0;
	}
      else
	{
	  // 5. octant
	  test_S3
	  test_S4
	  test_S5
	  return 0;
	}
    else
      if (u[2] > 0)
	{
	  // 4. octant
	  test_S2
	  test_S3
	  test_S6
	  return 0;
	}
      else
	{
	  // 8. octant
	  test_S2
	  test_S3
	  test_S5
	  return 0;
	}
  else
    if (u[1] > 0.0)
      if (u[2] > 0.0)
	{
	  // 2. octant
	  test_S1
	  test_S4
	  test_S6
	  return 0;
	}
      else
	{
	  // 6. octant
	  test_S1
	  test_S4
	  test_S5
          return 0;
	}
    else
      if (u[2] > 0.0)
	{
	  // 3. octant
	  test_S1
	  test_S2
	  test_S6
	  return 0;
	}
      else
	{
	  // 8. octant 
	  test_S1
	  test_S2
	  test_S5
	  return 0;
	}
  
}

#undef test_S1
#undef test_S2
#undef test_S3
#undef test_S4
#undef test_S5
#undef test_S6

//______________________________________________________________________________
Int_t TVpBox::IsInside(TVpVector3& pos) const
{
  // Return 1 if the point "pos" is inside the box (the box's boundary is
  // included) and 0 otherwise.
  //
  // Input parameters:
  // - pos - tested point 
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] TVpVector3 pos1(1, 2, 3);     // this point is inside the box
  // root [] cout << box->IsInside(pos1) << endl;
  // 1
  // root [] TVpVector3 pos2(-2, -2, -2);  // this point is outside the box
  // root [] cout << box->IsInside(pos2) << endl;
  // 0
  // root [] TVpVector3 pos3(10, 15, 20);  // this point is on the boundary
  // root [] cout << box->IsInside(pos3) << endl;
  // 1

  Double_t *r = pos.Get();

  return 0 <= r[0] && r[0] <= fSizeX 
    && 0 <= r[1] && r[1] <= fSizeY
    && 0 <= r[2] && r[2] <= fSizeZ;
}

//______________________________________________________________________________
Double_t TVpBox::GetVolume() const
{
  // Return the box volume in cm^3.
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] cout << box->GetVolume() << endl;
  // 3.00000000e+03

  return fSizeX * fSizeY * fSizeZ;
}

//______________________________________________________________________________
void TVpBox::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->PrintStatus()
  // <TVpBox>
  // Size [cm]:      1.00000000e+01  1.50000000e+01  2.00000000e+01
  // <TVpSolid>
  // Name:   box
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
  // </TVpBox>
  
  out << "<TVpBox>\n"
      << "$Id: TVpBox.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Size [cm]:\t" << fSizeX 
      << '\t' << fSizeY
      << '\t' << fSizeZ << '\n';
  TVpSolid::PrintStatus(out);
  out << "</TVpBox>" << std::endl;
}


#include "TPolyLine3D.h"

//______________________________________________________________________________
void TVpBox::Draw(Int_t color) const
{
  // Draw box edges.
  //
  // Input parameters:
  // - color - line color of main circles (default = solid_index+12)
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->Draw(1);  // the canvas is automatically created

  if (color == -1)
    color = fIndex+12;

  TVpVector3 go = PosLocToUni(TVpVector3(0, 0, 0));
  TVpVector3 gx = DirLocToUni(TVpVector3(fSizeX, 0, 0));
  TVpVector3 gy = DirLocToUni(TVpVector3(0, fSizeY, 0));
  TVpVector3 gz = DirLocToUni(TVpVector3(0, 0, fSizeZ));

  TVpVector3 gp;
  // Lower base
  TPolyLine3D *baseL = new TPolyLine3D(5);
  baseL->SetPoint(0, go.GetX(), go.GetY(), go.GetZ());
  gp = go+gx;
  baseL->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gx+gy;
  baseL->SetPoint(2, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gy;
  baseL->SetPoint(3, gp.GetX(), gp.GetY(), gp.GetZ());
  baseL->SetPoint(4, go.GetX(), go.GetY(), go.GetZ());

  // Upper base
  TPolyLine3D *baseU = new TPolyLine3D(5);
  gp = go+gz;
  baseU->SetPoint(0, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gx+gz;
  baseU->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gx+gy+gz;
  baseU->SetPoint(2, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gy+gz;
  baseU->SetPoint(3, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gz;
  baseU->SetPoint(4, gp.GetX(), gp.GetY(), gp.GetZ());

  // side 1
  TPolyLine3D *side1 = new TPolyLine3D(2);
  gp = go;
  side1->SetPoint(0, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gz;
  side1->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());

  // side 2
  TPolyLine3D *side2 = new TPolyLine3D(2);
  gp = go+gx;
  side2->SetPoint(0, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gx+gz;
  side2->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());

  // side 3
  TPolyLine3D *side3 = new TPolyLine3D(2);
  gp = go+gx+gy;
  side3->SetPoint(0, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gx+gy+gz;
  side3->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());

  // side 4
  TPolyLine3D *side4 = new TPolyLine3D(2);
  gp = go+gy;
  side4->SetPoint(0, gp.GetX(), gp.GetY(), gp.GetZ());
  gp = go+gy+gz;
  side4->SetPoint(1, gp.GetX(), gp.GetY(), gp.GetZ());

  baseL->SetLineColor(color);
  baseU->SetLineColor(color);
  side1->SetLineColor(color);
  side2->SetLineColor(color);
  side3->SetLineColor(color);
  side4->SetLineColor(color);
  baseL->Draw();
  baseU->Draw();
  side1->Draw();
  side2->Draw();
  side3->Draw();
  side4->Draw();

  DrawLcsAxes(fSizeX/10, fSizeY/10, fSizeZ/10);
}
