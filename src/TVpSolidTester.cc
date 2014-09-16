//______________________________________________________________________________
//
// TVpSolidTester tests ray intersection routines

//#include <limits>
#include <math.h>
#include "misc.h"
#include "TVpSolidTester.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

ClassImp(TVpSolidTester)

//______________________________________________________________________________
TVpSolidTester::TVpSolidTester(Double_t xMin, Double_t xMax,
			       Double_t yMin, Double_t yMax,
			       Double_t zMin, Double_t zMax)
{
  // Constructor

  fXMin = xMin;
  fXMax = xMax;
  fYMin = yMin;
  fYMax = yMax;
  fZMin = zMin;
  fZMax = zMax;
}

//______________________________________________________________________________
TVpSolidTester::~TVpSolidTester()
{
  // Destructor
}

//______________________________________________________________________________
TVpVector3 TVpSolidTester::GetRandomPoint()
{
  // Return point located randomly inside the box

  return TVpVector3 (getRand(fXMin, fXMax),
		     getRand(fYMin, fYMax),
		     getRand(fZMin, fZMax));
}

//______________________________________________________________________________
TVpVector3 TVpSolidTester::GetRandomDirection()
{
  // Return random direction

  // The symetry axis is (0,0,1). 
  Double_t cosTheta = 2.0 * getRand() - 1.0;
  Double_t sinTheta = sqrt(1.0 - cosTheta*cosTheta);
  Double_t phi = 2 * M_PI * getRand();
  Double_t sinPhi = sin(phi);
  Double_t cosPhi = cos(phi);
  return TVpVector3(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
}

//______________________________________________________________________________
void TVpSolidTester::DrawPointsInside(TVpSolid *solid, Int_t nSamples)
{
  // Draw points that are inside
  
  TPolyMarker3D *pm;
  for (Int_t i = 0; i < nSamples; i++)
    {
      TVpVector3 pos = GetRandomPoint();
      if (solid->IsInside(pos))
	{
	  pm = new TPolyMarker3D (1, 8,"");
	  pm->SetMarkerColor(6);
	  pm->SetMarkerStyle(8);
	  pm->SetMarkerSize(0.1);
	  pm->SetPoint(0, pos.fR[0], pos.fR[1], pos.fR[2]);
	  pm->Draw();
	}
    }
}

//______________________________________________________________________________
void TVpSolidTester::DrawPointsRayIntersectionIn(TVpSolid *solid, Int_t nSamples)
{
  // Draw intersection points, ray and surface, origin inside

  
  TPolyMarker3D *pm;
  TVpVector3 dir;
  //Double_t l = numeric_limits<Double_t>::max();
  Double_t l = 1e+30;
  Double_t t;

  for (Int_t i = 0; i < nSamples; i++)
    {
      TVpVector3 pos = GetRandomPoint();
      if (solid->IsInside(pos))
	{
	  dir = GetRandomDirection();
	  if (solid->RayIntersectionIn(pos, dir, l, t))
	    {
	      pos = pos + t * dir;
	      pm = new TPolyMarker3D (1, 8,"");
	      pm->SetMarkerColor(6);
	      pm->SetMarkerStyle(8);
	      pm->SetMarkerSize(0.1);
	      pm->SetPoint(0, pos.fR[0], pos.fR[1], pos.fR[2]);
	      pm->Draw();
	    }
	  else
	    std::cerr << "Error: Point is inside but no intersection occurs\n"; 
	}
    }
}

//______________________________________________________________________________
void TVpSolidTester::DrawPointsRayIntersectionOut(TVpSolid *solid, Int_t nSamples)
{
  // Draw intersection points, ray and surface, origin outside

  
  TPolyMarker3D *pm;
  TVpVector3 dir;
  //Double_t l = numeric_limits<Double_t>::max();
  Double_t l = 1e+30;
  Double_t t;

  for (Int_t i = 0; i < nSamples; i++)
    {
      TVpVector3 pos = GetRandomPoint();
      if (! solid->IsInside(pos))
	{
	  dir = GetRandomDirection();
	  if (solid->RayIntersectionOut(pos, dir, l, t))
	    {
	      pos = pos + t * dir;
	      pm = new TPolyMarker3D (1, 8,"");
	      pm->SetMarkerColor(6);
	      pm->SetMarkerStyle(8);
	      pm->SetMarkerSize(0.1);
	      pm->SetPoint(0, pos.fR[0], pos.fR[1], pos.fR[2]);
	      pm->Draw();
	    }
	}
    }
}

//______________________________________________________________________________
void TVpSolidTester::PrintVoxelArrayPathTable(TVpVoxelArray *va, Double_t energy,
					      Int_t nSamples)
{
  // Create a resizeable histogram and fill it with differences ...

  TVpVector3 pos1, pos2;  // positions of two random points
  TVpVector3 dir;         // direction
  Double_t rpt;           // real path true
  Double_t rpc;           // real path calculated
  Double_t opt;           // optical path true
  Double_t opc;           // optical path calculated
  TVpMaterial *matPtr;    // material
  Int_t ti1, ti2;         // tissue indices

  std::cout << std::scientific
	    << "                x                y                z "
	    << "                x                y                z"
	    << "            rpt            rpc            opc            opt" << std::endl;
  for (Int_t i = 0; i < nSamples; i++)
    {
      pos1 = GetRandomPoint();
      pos2 = GetRandomPoint();
      dir = pos2 - pos1;
      rpt = norm(dir);
      dir = (1/rpt)*dir;
      opc = va->GetOpticalPathInside(pos1.Get(), dir.Get(), energy, rpt);
      rpc = va->GetPathLength(pos1.Get(), dir.Get(), energy, opc);
      std::cout << pos1 << ' ' << pos2 << ' '
		<< rpt << ' ' << rpc << ' ' << opc;

      // If inside the solid then calculate the opt
      ti1 =  va->GetTissueIndex(pos1.Get());
      ti2 =  va->GetTissueIndex(pos2.Get());
      std::cout << ' ' << ti1 << ' ' << ti2;
      if (ti1 == 1 && ti2 == 1)
	{
	  matPtr = va->GetMaterial(pos1.Get());
	  opt = rpt * matPtr->GetLac(energy);
	  std::cout << ' ' << opt;
	}
      std::cout << std::endl;
    }
}

//______________________________________________________________________________
void TVpSolidTester::PrintVoxelArrayPathDifference(TVpVoxelArray *va, Double_t energy,
						   Int_t nSamples)
{
  // Create a resizeable histogram and fill it with differences ...

  TVpVector3 pos1, pos2;  // positions of two random points
  TVpVector3 dir;         // direction
  Double_t rpt;           // real path true
  Double_t rpc;           // real path calculated
  Double_t opc;           // optical path calculated

  std::cout << std::scientific;
  for (Int_t i = 0; i < nSamples; i++)
    {
      pos1 = GetRandomPoint();
      pos2 = GetRandomPoint();
      dir = pos2 - pos1;
      rpt = norm(dir);
      dir = (1/rpt)*dir;
      opc = va->GetOpticalPathInside(pos1.Get(), dir.Get(), energy, rpt);
      rpc = va->GetPathLength(pos1.Get(), dir.Get(), energy, opc);
      std::cout << rpc - rpt << std::endl;
    }
}

//______________________________________________________________________________
void TVpSolidTester::PrintVoxelArrayOpticalPathTwoDirections(TVpVoxelArray *va,
							     Double_t energy,
							     Int_t nSamples)
{
  // Create a resizeable histogram and fill it with differences ...

  TVpVector3 pos1, pos2;  // positions of two random points
  TVpVector3 dir1, dir2;  // directions
  Double_t op1, op2;      // optical path
  Double_t rp;            // real path
  
  std::cout << std::scientific;
  for (Int_t i = 0; i < nSamples; i++)
    {
      pos1 = GetRandomPoint();
      pos2 = GetRandomPoint();
      dir1 = pos2 - pos1;
      rp = norm(dir1);
      dir1 = (1/rp)*dir1;
      dir2 = -1.0 * dir1;
      op1 = va->GetOpticalPathInside(pos1.Get(), dir1.Get(), energy, rp);
      op2 = va->GetOpticalPathInside(pos2.Get(), dir2.Get(), energy, rp);
      std::cout << op2 - op1 << ' ' << std::endl;
    }
}
