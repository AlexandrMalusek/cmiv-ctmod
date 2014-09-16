//______________________________________________________________________________
//
// TVpSolid defines attributes common to all solids (spheres, boxes,
// cylinders, ellipsoids, voxel arrays, ...).
//
// In CTmod, a solid is a basic building block of the geometry.  Each solid
// has an index, name, and position in space.  The TVpSolid class also stores
// information about neighbours: a list of overlaping solids and a list of
// base solids.  These are used during particle tracking.  Homogenous solids
// store a pointer to one material, voxel arrays store pointers to all
// materials used in these arrays.
//
// In the following, UCS stands for the Universe Coordinate System and LCS
// stands for the local coordinate system.  The user specifies positions of
// solids in UCS.  Particle tracking routines use the LCS of each solid.
//
//
// To position a solid in the geometry:
//
// - Define the universe.  It may be a sphere, box, or other solid.  Use
// (0,0,0) as the translation vector and the unit matrix as the rotation
// matrix.
//
// - Find the position of the solid in LCS (see class description of TVpBox,
// TVpCylinder, TVpSphere, TVpEllipsoid, or TVpVoxelArray).
//
// - Define the translation vector as an active translation of the LCS of the
// solid into the desired position.
//
// - Define the rotation matrix as an active rotation of the solid in its LCS.
//
// TVpSolid is an abstract class and thus it cannot be instantiated.  In the
// following examples, TVpBox is used instead.  The box in the default
// position is drawn with black lines, the transformed one is drawn with
// magenta lines.
//
// Example 1:
// // Right-hand (positive) Pi/2 rotation about the y-axis.
// gSystem->Load("libRCTmod.so");
// TVpBox *boxA = new TVpBox("boxA", 1, 10, 15, 20);
// TVpBox *boxB = new TVpBox("boxB", 1, 10, 15, 20);
// boxB->SetActiveRotation(new TVpMatrix3x3(0, TMath::Pi()/2, 0));
//Begin_Html
/*
<img src="png/TVpSolid_SetActiveRotation1.png">
*/
//End_Html
//
// Example 2:
// Right-hand (positive) Pi/2 rotation about the z-axis.
// gSystem->Load("libRCTmod.so");
// TVpBox *boxA = new TVpBox("boxA", 1, 10, 15, 20);
// TVpBox *boxB = new TVpBox("boxB", 1, 10, 15, 20);
// boxB->SetActiveRotation(new TVpMatrix3x3(TMath::Pi()/2, 0, 0));
//Begin_Html
/*
<img src="png/TVpSolid_SetActiveRotation2.png">
*/
//End_Html
//
// Example 3:
// // 1. Right-hand (positive) Pi/2 rotation about the z-axis.
// // 2. Right-hand (positive) Pi/2 rotation about the y-axis.
// gSystem->Load("libRCTmod.so");
// TVpBox *boxA = new TVpBox("boxA", 1, 10, 15, 20);
// TVpBox *boxB = new TVpBox("boxB", 1, 10, 15, 20);
// boxB->SetActiveRotation(
//   new TVpMatrix3x3(TMath::Pi()/2, TMath::Pi()/2, 0));
//Begin_Html
/*
<img src="png/TVpSolid_SetActiveRotation3.png">
*/
//End_Html
//
// Example 4:
// // 1. Right-hand (positive) Pi/2 rotation about the z-axis.
// // 2. Right-hand (positive) Pi/2 rotation about the y-axis.
// // 3. Move the box by (2,4,6).
// gSystem->Load("libRCTmod.so");
// TVpBox *boxA = new TVpBox("boxA", 1, 10, 15, 20);
// TVpBox *boxB = new TVpBox("boxB", 1, 10, 15, 20);
// boxB->SetActiveRotation(
//   new TVpMatrix3x3(TMath::Pi()/2, TMath::Pi()/2, 0));
// boxB->SetActiveTranslation(new TVpVector3(2, 4, 6));
//Begin_Html
/*
<img src="png/TVpSolid_SetTransformation.png">
*/
//End_Html

#include <stdio.h>
#include "TVpSolid.h"

ClassImp(TVpSolid)

//______________________________________________________________________________
TVpSolid::TVpSolid(const Char_t *name, Int_t index)
{
  // Constructor also used as the default constructor.  Set undefined members
  // to 0, use unit rotation matrix.
  //
  // Input parameters:
  // - name - user defined descriptive name.
  // - index - solid index.  Must be unique in the geometry (range = 0, ...)

  fOverlapNumber = fBaseNumber = 0;
  fOverlap = fBase = 0;
  fMaterial = 0;
  fIndex = index;

  // Create a copy of the name.  Use a string in the future.
  Int_t strLen = strlen(name) + 1;
  fName = new Char_t[strLen];
  strncpy(fName, name, strLen);

  fTraVecL2u = TVpVector3(0, 0, 0);
  fRotMatL2u = TVpMatrix3x3(1,0,0, 0,1,0, 0,0,1);
  fRotMatL2uT = TVpMatrix3x3(1,0,0, 0,1,0, 0,0,1);  // = transpose(fRotMatL2u)
}

//______________________________________________________________________________
TVpSolid::~TVpSolid()
{
  // Destructor.  Destruct all data related to the solid, do not destruct the
  // material.

  delete[] fName;
  delete[] fOverlap;
  delete[] fBase;
}

//______________________________________________________________________________
void TVpSolid::SetOverlap(TVpSolidPtr *overlapSolidList)
{
  // Set the list of overlapping solids.
  //
  // Input parameters:
  // - overlapSolidList - a null terminated list of overlapping solids.  Use 0
  // to specify that the solid has no overlapping solids.
  //
  // Example:
  // root [] TVpSolidPtr sphere0Overlap[] = {cylinder1, 0};
  // root [] sphere0->SetOverlap(sphere0Overlap);
  // root [] box->SetOverlap(0);  // empty list

  if (overlapSolidList == 0)
    {
      fOverlapNumber = 0;
      fOverlap = 0;
      return;
    }

  // Find the number of solids in the list
  fOverlapNumber = 0;
  for (TVpSolidPtr *oPtr = overlapSolidList; *oPtr != 0; oPtr++)
    fOverlapNumber++;

  // Allocate array of neighbors
  fOverlap = new TVpSolidNeighbor[fOverlapNumber];

  // Initialize the the SolidNeighbor array
  for (Int_t i = 0; i < fOverlapNumber; i++)
    fOverlap[i].fSolidPtr = overlapSolidList[i];
}

//______________________________________________________________________________
void TVpSolid::SetBase(TVpSolidPtr *baseSolidList)
{
  // Set the list of base solids.
  //
  // Input parameters:
  // - baseSolidList - a null terminated list of base solids.  Use 0 to
  // specify that the solid has no base solids.
  //
  // Example:
  // root [] TVpSolidPtr cylinder1Base[] = {sphere0, 0};
  // root [] cylinder1->SetBase(cylinder1Base);
  // root [] box->SetBase(0);  // Empty list

  if (baseSolidList == 0)
    {
      fBaseNumber = 0;
      fBase = 0;
      return;
    }

  // Find the number of solids in the list
  fBaseNumber = 0;
  for (TVpSolidPtr *oPtr = baseSolidList; *oPtr != 0; oPtr++)
    fBaseNumber++;

  // Allocate array of neighbors
  fBase = new TVpSolidNeighbor[fBaseNumber];

  // Initialize the the SolidNeighbor array
  for (Int_t i = 0; i < fBaseNumber; i++)
    fBase[i].fSolidPtr = baseSolidList[i];
}

//______________________________________________________________________________
void TVpSolid::InitializeNeighbors()
{
  // Initialize lists of neighbours.  Solid-to-solid rotation matrices and
  // translation vectors are calculated for all neighbours.

  TVpVector3 traVec;
  TVpMatrix3x3 rotMat;
  TVpVector3 vec0(0,0,0);
  TVpMatrix3x3 matI(1,0,0, 0,1,0, 0,0,1);
  TVpSolidPtr neighborPtr;

  for (Int_t i = 0; i < fOverlapNumber; i++)
    {
      neighborPtr = fOverlap[i].fSolidPtr;
      // Calculate the relative translation vector
      traVec = neighborPtr->fRotMatL2uT * (fTraVecL2u - neighborPtr->fTraVecL2u);
      fOverlap[i].fTraVecC2nPtr = (traVec == vec0) ? 0 : new TVpVector3(traVec);

      // Calculate the relative rotation matrix
      rotMat = neighborPtr->fRotMatL2uT * fRotMatL2u;
      fOverlap[i].fRotMatC2nPtr = (rotMat == matI) ? 0 : new TVpMatrix3x3(rotMat);
    }

  for (Int_t i = 0; i < fBaseNumber; i++)
    {
      neighborPtr = fBase[i].fSolidPtr;
      // Calculate the relative translation vector
      traVec = neighborPtr->fRotMatL2uT * (fTraVecL2u - neighborPtr->fTraVecL2u);
      fBase[i].fTraVecC2nPtr = (traVec == vec0) ? 0 : new TVpVector3(traVec);

      // Calculate the relative rotation matrix
      rotMat = neighborPtr->fRotMatL2uT * fRotMatL2u;
      fBase[i].fRotMatC2nPtr = (rotMat == matI) ? 0 : new TVpMatrix3x3(rotMat);
    }
}

//______________________________________________________________________________
void TVpSolid::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->PrintStatus();
  // ...
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
  // ... 

  out << "<TVpSolid>\n"
      << "$Id: TVpSolid.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Name:\t" << fName << '\n'
      << "Index:\t" << fIndex << '\n';

  // For a voxel array, fMaterial == 0
  if (fMaterial == 0)
    out << "Material:\t" << "Undefined or voxel array" << '\n';
  else
    out << "Material:\t" << fMaterial->GetName() << '\n';

  out << "Position of the solid is defined via T and R so that xUni = R*xLoc + T\n"
      << "Translation vector T:\n"
      << fTraVecL2u << '\n'
      << "Rotation matrix R:\n"
      << fRotMatL2uT;

  // Overlap solids
  if (fOverlapNumber == 0)
    out << "Overlap solid indices:\tNone\n";
  else
    {
      out << "Overlap solid indices:";
      for (Int_t i = 0; i < fOverlapNumber; i++)
	out << ' ' << fOverlap[i].fSolidPtr->fIndex;
      out << '\n';
      for (Int_t i = 0; i < fOverlapNumber; i++)
	fOverlap[i].PrintStatus(out);
    }

  // Base solids
  if (fBaseNumber == 0)
    out << "Base solid indices:\tNone\n";
  else
    {
      out << "Base solid indices:";
      for (Int_t i = 0; i < fBaseNumber; i++)
	out << ' ' << fBase[i].fSolidPtr->fIndex;
      out << '\n';
      for (Int_t i = 0; i < fBaseNumber; i++)
	fBase[i].PrintStatus(out);
    }
  out << "</TVpSolid>\n";
}

//______________________________________________________________________________
Int_t TVpSolid::GetSubIndex(Double_t r[]) const
{
  // Return the solid subindex, e.g the tissue index of a voxel array,
  // corresponding to point "r".  Return 0 for solids with no internal
  // structure.
  //
  // Input parameters:
  // - r[3] - position in LCS
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] Double_t r[3] = {1, 2, 3};
  // root [] cout << box->GetSubIndex(r) << endl;
  // 0

  return 0;
}

//______________________________________________________________________________
Double_t TVpSolid::GetDensity() const
{
  // Return the density in g/cm^3 for homogenous solids and -1 for voxel
  // arrays (when fMaterial == 0).  This is just a TEMPORARY WORKAROUND.
  //
  // Example:
  // root [] TVpMaterial *matWater = new TVpMaterial("material/water.mat");
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->SetMaterial(matWater);
  // root [] cout << box->GetDensity() << endl;
  // 1.00000000e+00

  if (fMaterial != 0)
    return fMaterial->GetDensity();

  return -1.0;
}

//______________________________________________________________________________
Double_t TVpSolid::GetPathLength(Double_t *r, Double_t *u, Double_t energy,
				 Double_t opticalPath) const
{
  // Return path length of a line segment specified by position r[3],
  // direction u[3], and optical path.  The position and direction are in
  // local coordinate system of the solid (LCS). The point r must be inside
  // the solid.
  //
  // For homogenous solids, return optical path divided by the linear
  // attenuation coefficient. For voxel arrays, return the corresponding path
  // length if inside the array or a very big number if the particle leaves
  // the array.
  //
  // Input parameters:
  // - r[3] - ray's position in LCS in cm
  // - u[3] - ray's direction in LCS.  Must be a unit vector.
  // - energy - photon energy in keV
  // - opticalPath - optical path (dimensionless)
  //
  // Example:
  // root [] TVpMaterial *matWater = new TVpMaterial("material/water.mat",
  // "material/water_m80901.cff", "material/water.isf");
  // root [] matWater->Initialize();
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->SetMaterial(matWater);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] cout << box->GetPathLength(pos.Get(), dir.Get(), 60.0, 1.02856) << endl;
  // 4.99999

  Double_t lac = fMaterial->GetToCsg(energy);

  return opticalPath / lac;
}

//______________________________________________________________________________
Double_t TVpSolid::GetOpticalPathInside(Double_t r[], Double_t u[], Double_t energy,
					Double_t distance) const
{
  // Return optical path of a line segment specified by position r[3],
  // direction u[3], and length.  The position and direction are in local
  // coordinate system of the solid (LCS). The line segment must be inside the
  // solid.
  //
  // Input:
  // - r[3] - ray's position in LCS in cm
  // - u[3] - ray's direction in LCS.  Must be a unit vector.
  // - energy - photon energy in keV
  // - distance - path length in cm
  //
  // Example:
  // root [] TVpMaterial *matWater = new TVpMaterial("material/water.mat",
  // "material/water_m80901.cff", "material/water.isf");
  // root [] matWater->Initialize();
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->SetMaterial(matWater);
  // root [] TVpVector3 pos(1, 1, 1);
  // root [] TVpVector3 dir = normalize(TVpVector3(1, 1, 1));
  // root [] cout << box->GetOpticalPathInside(pos.Get(), dir.Get(), 60.0, 5.0) << endl;
  // 1.02856

  Double_t lac = fMaterial->GetToCsg(energy);
  Double_t opticalPath = lac * distance;
  return opticalPath;
}

