//______________________________________________________________________________
//
// TVpGeometry defines the phantom geometry.
//
//
// Example:
// gSystem->Load("libRCTmod.so");
//
// // Definition of materials
// TVpMaterial *matVacuum = new TVpMaterial("material/vacuum.mat");  
// matVacuum->Initialize();
//
// TVpMaterial *matWater = new TVpMaterial
//   ("material/water.mat",
//    "material/water_m80901.cff",  //_m80901
//    "material/water.isf");
// matWater->Initialize();
//
// // Rotation matrix
// // 1. Right-hand (positive) Pi/2 rotation about the z-axis.
// // 2. Right-hand (positive) Pi/2 rotation about the y-axis.
// TVpMatrix3x3 rotMat(TMath::Pi()/2, TMath::Pi()/2, 0);
//
// // sphere, the universe
// TVpSphere *sphere = new TVpSphere("sphere", 0, 80);  // The Universe
// sphere->SetMaterial(matVacuum);
//
// // box inside the sphere
// TVpBox *box = new TVpBox("box", 1, 20, 40, 60);
// box->SetActiveTranslation(new TVpVector3(-60/2.0, -20/2.0, -40/2.0));
// box->SetActiveRotation(&rotMat);
// box->SetMaterial(matWater);
//
// // Base and overlap lists
// TVpSolidPtr sphereOverlap[] = {box, 0};
// sphere->SetOverlap(sphereOverlap);
// sphere->SetBase(0);
// TVpSolidPtr boxBase[] = {sphere, 0};
// box->SetOverlap(0);
// box->SetBase(boxBase);
//
// TVpGeometry *geometry = new TVpGeometry(sphere);
//
//Begin_Html
/*
<img src="png/geometryBox.png">
*/
//End_Html

#include <iomanip>
#include "TVpGeometry.h"
#include "TVpSolid.h"

ClassImp(TVpGeometry)

//______________________________________________________________________________
TVpGeometry::TVpGeometry()
{
  // Default constructor

  fUniverse = 0;
  fInfinity = 1e10;   // cm
  fShift = 1.0e-20;   // cm
  fSolidNumber = 0;

  SetViewRange(TVpVector3(100,100,100));
}

//______________________________________________________________________________
TVpGeometry::TVpGeometry(TVpSolid *universe)
{
  // Constructor with full initialization.

  fUniverse = universe;
  fInfinity = 1e10;   // cm
  fShift = 1.0e-20;   // cm
  fSolidNumber = 0;
  Initialize();

  SetViewRange(TVpVector3(100,100,100));
}

//______________________________________________________________________________
TVpGeometry::~TVpGeometry()
{
  // Destructor
}

//______________________________________________________________________________
void TVpGeometry::Initialize()
{
  // Initialize the geometry structures This implementation uses max 128
  // solids. Will be improved later via the STL vector.
  
  const Int_t maxSolid = 128;
  fSolidPtr = new TVpSolidPtr[maxSolid];

  // Insert all solids
  InsertSolidPtr(fUniverse);

  // Initialize them
  for (Int_t i = 0; i < fSolidNumber; i++)
    fSolidPtr[i]->InitializeNeighbors();
}

//______________________________________________________________________________
void TVpGeometry::InsertSolidPtr(TVpSolid *solid)
{
  // Insert a solid and all its overlapping solids into the fSolidPtr array.
  // This routine is used during geometry initialization to fill up the
  // fSolidPtr array.
  //
  // Input parameters:
  // - solid - the solid pointer
  //
  // Method:
  // If the pointer is not in the fSolidPtr array then insert it there.  Use
  // recursion and check all the overlapping solids.
  
  // Check if it is not there
  Int_t isThere = 0;
  for (Int_t i = 0; i < fSolidNumber; i++)
    if (fSolidPtr[i] == solid)
      isThere = 1;

  if (! isThere)
    fSolidPtr[fSolidNumber++] = solid;  // Put it there

  // Use recursion for overlaping solids
  for (Int_t i = 0; i < solid->fOverlapNumber; i++)
    InsertSolidPtr(solid->fOverlap[i].fSolidPtr);
}

//______________________________________________________________________________
TVpSolid *TVpGeometry::GetSolid(Int_t indexOfSolid) const
{
  // Return a pointer to the solid with a given index.
  //
  // Input parameters:
  // - indexOfSolid - the solid index, range: 0, ...

  for (Int_t i = 0; i < fSolidNumber; i++)
    if (fSolidPtr[i]->fIndex == indexOfSolid)
      return fSolidPtr[i];

  return 0;
}

//______________________________________________________________________________
Double_t TVpGeometry::GetOpticalPath(TVpParticle particle, Double_t distance) const
{
  // Calculate optical path corresponding to real path length "distance".
  // Note that particle is passed by value because it is transported
  // inside this routine.
  //
  // Input parameters:
  // - particle - the particle for which the optical path is calculated
  // - distance - the particle's real path
  
  Double_t t;
  Double_t opticalPath = 0.0;  // Optical path = exp(-lac * t)
  EStepResult stepResult;
  TVpParticle particleOld;     // temporary storage of particle data

  do
    {
      // Remember the particle data
      particleOld = particle;

      // Try to transport the particle over the distance "distance".
      // The actual transport length is in the variable "t".
      stepResult = Step(&particle, distance, t);
      
      // Use the old particle to calculate the optical path
      Double_t dop = particleOld.fSolid->GetOpticalPathInside(particleOld.GetPos(),
							      particleOld.GetDir(),
							      particleOld.GetEnergy(),
							      t);
      opticalPath += dop;
      distance -= t;
    }
  while (stepResult == kOverlapSolidEntered || stepResult == kCurrentSolidLeft);
  return opticalPath;
}

//______________________________________________________________________________
Int_t TVpGeometry::TraceRay(TVpVector3& startPosUni, TVpVector3& endPosUni) const
{
  // Print crossing points between a virtual photon travelling from
  // startPosUni to endPosUni and the geometry.
  //
  // Input parameters:
  // - startPosUni - start possition in universe coordinates
  // - endPosUni - end position in universe coordinates
  
  EStepResult stepResult;
  Char_t *namePtr;
  Double_t infinity = 1e38;
  Double_t t;

  // Get local coordinates
  TVpSolid *solidPtr = GetSolid(startPosUni);
  TVpVector3 dirU = normalize(endPosUni - startPosUni);
  TVpVector3 dirL = solidPtr->DirUniToLoc(dirU);
  TVpVector3 posL = solidPtr->PosUniToLoc(startPosUni);
  TVpVector3 posU;

  // Create virtual particle
  TVpParticle particle;
  particle.SetPosition(posL);
  particle.SetDirection(dirL);
  particle.SetSolid(solidPtr);
  particle.SetEnergy(1.0);
  particle.SetWeight(1.0);
  std::cout << "Start:\n" << particle << std::endl;
  
  do
    {
      // Try to transport the particle over the distance "distance".
      // The actual transport length is in the variable "t".
      stepResult = Step(&particle, infinity, t);
      
      // Use the old particle to calculate the optical path
      switch (stepResult)
	{
	  case kOverlapSolidEntered:
	    std::cout << "o ";
	    break;
	  case kCurrentSolidLeft:
	    std::cout << "b ";
	    break;
	  case kTerminated:
	    std::cout << "t ";
	    break;
	  default:
	    std::cout << "internal error ";
	    break;
	}
      solidPtr = particle.GetSolidPtr();
      posL = particle.GetLocPos();
      posU = solidPtr->PosLocToUni(posL);
      namePtr = solidPtr->GetMaterial(posL.Get())->GetName();

      std::cout << std::scientific
		<< std::setw(3) << solidPtr->GetIndex() << ' '
		<< t << ' '
		<< posU << ' '
		<< posL << ' '
		<< namePtr << std::endl;
    }
  while (stepResult == kOverlapSolidEntered || stepResult == kCurrentSolidLeft);

  return 0;
}

//______________________________________________________________________________
TVpSolid *TVpGeometry::GetSolid(TVpVector3& posLoc, TVpSolid *solid) const
{
  // Return pointer to the solid which contains the point.  Position 'posLoc'
  // is specified in the coordinate system of the solid 'solid'.  If solid is
  // not specified then start from the universe. In this case, if the point
  // isn't in universe then return 0.
  //
  // Input parameters:
  // - posLoc - position of the point in local coordinates
  // - solid - pointer to the solid from which the search starts
  //
  // Method:
  // Overlapping solids are tested first.  If the point is contained in one of
  // them then recursion is used to find if the point is in one of its
  // overlapping solids.  Otherwise the pointer to the current solid is
  // returned.  No test is performed that the point is really contained in the
  // starting solid.

  TVpVector3 locPos;          // Position in local coordinates of the probed solid
  TVpVector3 *traVecPtr;      // translation vector of the probed solid
  TVpMatrix3x3 *rotMatPtr;    // rotation matrix of the probed solid
  TVpSolidNeighbor *neighbor;

  // If no solid specified, use the universe
  if (solid == 0)
    {
      solid = fUniverse;
      if (solid->IsInside(posLoc) == 0)
	return 0;       // We are not in the universe. Where are we?
    }
  
  // Traverse the overlap solid list
  for (Int_t i = 0; i < solid->fOverlapNumber; i++)
    {
      neighbor = &solid->fOverlap[i];
      // Get coordinates in the new solid
      if ((rotMatPtr = neighbor->fRotMatC2nPtr) != 0)
	locPos = *rotMatPtr * posLoc;
      else
	locPos = posLoc;
      if ((traVecPtr = neighbor->fTraVecC2nPtr) != 0)
	locPos = locPos + *traVecPtr;
      
      if (neighbor->fSolidPtr->IsInside(locPos))
	// The point is inside an overlaping solid. Use recursion
	// to investigate the overping solid.
	return GetSolid(locPos, neighbor->fSolidPtr);
    }

  // The point wasn't located in any of the overlaping solids. It must
  // be in the original solid.
  return solid;
}
  
//______________________________________________________________________________
TVpGeometry::EStepResult TVpGeometry::Step(TVpParticle *particle, const Double_t l,
					   Double_t& t) const
{
  // Transport the particle and return the result.
  //
  // Input parameters:
  // - particle - pointer to the particle to be transported
  // - l - requested step length in cm
  // - t - distance travelled in cm (may be smaller than l)
  //
  // Method:
  // First, check if the line segment intersects an overlapping solid.  If yes
  // then select the closest solid and transport the particle to its boundary.
  // Otherwise check the current solid.  If there is an intersection then
  // transport the particle there and find the corresponding base solid.
  // Otherwise transport the particle in the current solid.  Coordinates are
  // transformed into the coordinate system of the entered solid.  Counters of
  // particles which left the current solid and entered a new solid are
  // updated.

  Double_t ts;
  TVpVector3 newSolidPos;   // Position in the new solid coordinates
  TVpVector3 newSolidDir;   // Direction in the new solid coordinates
  TVpSolid *newSolid = 0;

  TVpVector3 locPos;          // Position in local coordinates of the probed solid
  TVpVector3 locDir;          // Direction in local coordinated of the tested solid
  TVpVector3 *traVecPtr;         // translation vector of the probed solid
  TVpMatrix3x3 *rotMatPtr;       // rotation matrix of the probed solid
  TVpSolidNeighbor *neighbor;

  TVpSolid *solid = particle->fSolid;  // Current solid where the particle is
  ts = t = l;
  
  // Check the overlap solid list
  for (Int_t i = 0; i < solid->fOverlapNumber; i++)
    {
      neighbor = &solid->fOverlap[i];
      
      // Get coordinates in the new solid
      if ((rotMatPtr = neighbor->fRotMatC2nPtr) != 0)
	{
	  locPos = *rotMatPtr * particle->fPosition;
	  locDir = *rotMatPtr * particle->fDirection;
	}
      else
	{
	  locPos = particle->fPosition;
	  locDir = particle->fDirection;
	}
      if ((traVecPtr = neighbor->fTraVecC2nPtr) != 0)
	locPos = locPos + *traVecPtr;
      
      // Test intersection
      if (neighbor->fSolidPtr->RayIntersectionOut(locPos, locDir, l, ts))
	if (ts < t)
	  {
	    // Remember the solid
	    t = ts;
	    newSolid = neighbor->fSolidPtr;
	    newSolidPos = locPos;
	    newSolidDir = locDir;
	  }
    }
  
  // If we crossed, return. We asume convex solids here.
  if (newSolid != 0)
    {
      // Move the particle and set up its position in the new solid coordinates
      particle->fPosition = newSolidPos + t * newSolidDir;
      particle->fDirection = newSolidDir;
      particle->fSolid = newSolid;
      return kOverlapSolidEntered;
    }
  
  // Check the current solid
  if (solid->RayIntersectionIn(particle->fPosition, particle->fDirection, l, t))
    {
      // The particle leaves the current solid.
      particle->fPosition = particle->fPosition + 
	(t + fShift) * particle->fDirection;      // needs some tuning

      // Test if the particle leaves the universe
      if (solid == fUniverse)
	return kTerminated;
    

      // Find the base solid
      for (Int_t i = 0; i < solid->fBaseNumber; i++)
	{
	  neighbor = &solid->fBase[i];
	  
	  // Move to the local coordinate system.
	  if ((rotMatPtr = neighbor->fRotMatC2nPtr) != 0)
	    {
	      locPos = *rotMatPtr * particle->fPosition;
	      locDir = *rotMatPtr * particle->fDirection;
	    }
	  else
	    {
	      locPos = particle->fPosition;
	      locDir = particle->fDirection;
	    }
	  if ((traVecPtr = neighbor->fTraVecC2nPtr) != 0)
	    locPos = locPos + *traVecPtr;
	  
	  if (neighbor->fSolidPtr->IsInside(locPos))
	    {
	      // The point is inside the base solid.
	      particle->fPosition = locPos;
	      particle->fDirection = locDir;
	      particle->fSolid = neighbor->fSolidPtr;
	      return kCurrentSolidLeft;
	    } 
	}
    }
  
  // We stay in the current solid
  t = l;
  particle->fPosition = particle->fPosition + t * particle->fDirection;
  return kNoChange;
}

//______________________________________________________________________________
Double_t TVpGeometry::GetNativeVolume(Int_t indexOfSolid, Int_t numInteg) const
{
  // Return native (non-overlapped) volume of the solid with index
  // indexOfSolid or
  // -1 ... invalid indexOfSolid
  // -2 ... cannot calculate analytically and numeric integration not requested
  // -3 ... numerical integration not yet implemented (= internal error)
  //
  // Input parameters:
  // - indexOfSolid - index of the solid
  // - numInteg - perform numerical integration (range=0,1, default = 0)
  
  TVpSolid *solidPtr = GetSolid(indexOfSolid);
  if (solidPtr == 0)
    return -1.0;
  Double_t volume = solidPtr->GetVolume();
  
  // No overlapping solids?
  if (solidPtr->fOverlapNumber == 0)
    return volume;
  
  // Are all overlapping solids inside and separated?
  Int_t additive = 1;
  Double_t volumeOfOverlap = 0.0;
  for (Int_t i = 0; i < solidPtr->fOverlapNumber; i++)
    {
      TVpSolid *sPtr = solidPtr->fOverlap[i].fSolidPtr;
      if (sPtr->fBaseNumber != 1)  // the original solid only
	{
	  additive = 0;
	  continue;
	}
      volumeOfOverlap += sPtr->GetVolume();
    }
  if (additive == 1)
    return volume - volumeOfOverlap;
  
  // Do not perform numerical integration
  if (numInteg == 0)
    return -2;
  
  // Numerical integration not yet implemented!
  return -3;
}
 
//______________________________________________________________________________
Double_t TVpGeometry::GetNativeMass(Int_t indexOfSolid, Int_t numInteg) const
{
  // Return native (non-overlapped) mass of the solid with index
  // indexOfSolid or
  // -1 ... invalid indexOfSolid
  // -2 ... cannot calculate analytically and numeric integration not requested
  // -3 ... numerical integration not yet implemented (= internal error)
  //
  // Input parameters:
  // - indexOfSolid - index of the solid
  // - numInteg - perform numerical integration (range=0,1, default = 0)

  Double_t volume = GetNativeVolume(indexOfSolid, numInteg);
  if (volume < 0)
    return volume;
  TVpSolid *solidPtr = GetSolid(indexOfSolid);

  // HOT FIX
  if (solidPtr->fMaterial == 0)
    return -4;

  Double_t density = solidPtr->fMaterial->GetDensity();
  return volume * density;
}

//______________________________________________________________________________
Double_t TVpGeometry::GetAbsorbedDose(Int_t indexOfSolid, Long_t numOfHist,
				      Int_t numInteg) const
{
  // Return absorbed dose

  Double_t mass = GetNativeMass(indexOfSolid, numInteg);
  if (mass < 0)
    return mass;
  TVpSolid *solidPtr = GetSolid(indexOfSolid);
  Double_t energyImparted = solidPtr->GetEnergyImparted(numOfHist);
  return energyImparted / mass;
}

//______________________________________________________________________________
void TVpGeometry::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpGeometry>\n"
      << "$Id: TVpGeometry.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Solids:\t" << fSolidNumber << '\n';

  for (Int_t i = 0; i < fSolidNumber; i++)
    fSolidPtr[i]->PrintStatus(out);
  out << "</TVpGeometry>\n";
}

//______________________________________________________________________________
void TVpGeometry::PrintStatistics(std::ostream &out, Long_t numOfHist) const
{
  // Print statistics

  out << "<statistics_geometry>\n";
  for (Int_t i = 0; i < fSolidNumber; i++)
    {
      out << "<statistics_solid>\n"
	  << "  solid number: " << i << '\n';
      fSolidPtr[i]->PrintStatistics(out, numOfHist);
      out << "</statistics_solid>\n";
    }
  out << "</statistics_geometry>\n";

  // Volume, mass, average absorbed dose
  out << std::scientific;
  out << "<statistics_geometry_table>\n";
  out << "index\tvolume       mass    D [keV/g] per hist. \n";
  for (Int_t i = 0; i < fSolidNumber; i++)
    out << i << '\t'
	<< GetNativeVolume(i) << ' '
	<< GetNativeMass(i) << ' '
	<< GetAbsorbedDose(i, numOfHist) << '\n';
  out << "</statistics_geometry_table>\n";
}

//______________________________________________________________________________
void TVpGeometry::Voxelize(TVpVoxelArray *voxelArrayPtr, TVpVector3 traVec,
			   TVpMatrix3x3 rotMat) const
{
  // Voxelize the geometry.  The dimension and size of the voxel array
  // must be defined.  The function sets up the voxel array data which
  // can be saved to a file by the voxel array member functions.
  //
  // Now, only one sample point in the voxel center is used, in the
  // future, a grid of points may be used.
  //
  // At present, only the tissue index is set. The solid index is
  // stored there.
  //
  // voxelArrayPtr->WriteVADFile ... can be used
  // voxelArrayPtr->WriteVAMFile ... cannot be used
  // voxelArrayPtr->WriteVANFile ... cannot be used

  TVpVector3 posLoc;   // position in the voxel array local coordinates
  TVpVector3 posUni;   // position in the universe coordinates
  TVpSolid *solidPtr;

  // Set the tissue index
  for (Int_t ix = 0; ix < voxelArrayPtr->GetNX(); ix++)
      for (Int_t iy = 0; iy < voxelArrayPtr->GetNY(); iy++)
	for (Int_t iz = 0; iz < voxelArrayPtr->GetNZ(); iz++)
	  {
	    posLoc.Set((ix+0.5) * voxelArrayPtr->fDpX,
		       (iy+0.5) * voxelArrayPtr->fDpY,
		       (iz+0.5) * voxelArrayPtr->fDpZ);
	    posUni = rotMat * posLoc + traVec;
#if 0
	    if (ix == 0 && iy == 0 && iz == 0)
	      std::cout << "posLoc = " << posLoc.fR[0] << ' '
			<< posLoc.fR[1] << ' '
			<< posLoc.fR[2] << '\n'
			<< "posUni = " << posUni.fR[0] << ' '
			<< posUni.fR[1] << ' '
			<< posUni.fR[2] << '\n';
#endif
	    solidPtr = GetSolid(posUni);
	    voxelArrayPtr->fTissue[voxelArrayPtr->GetIndex(iz, iy, ix)] =
	      solidPtr->fIndex;
	  }
}

#include "TCanvas.h"
#include "TPolyLine3D.h"

//______________________________________________________________________________
TView3D* TVpGeometry::Draw(TView3D *viewPtr) const
{
  // Draw the geometry.  Return a pointer to the "view" used.

  TCanvas *c1;
  if (viewPtr == 0)
    {
      c1 = new TCanvas("geometry", "geometry", 800, 800);
      viewPtr = new TView3D(1, fViewRangeMin, fViewRangeMax);
      c1->Draw();
    }

  // Draw x-axis
  TPolyLine3D *xaxis = new TPolyLine3D(2);
  xaxis->SetPoint(0, 0, 0, 0);
  xaxis->SetPoint(1, fViewRangeMax[0], 0, 0);
  xaxis->SetLineColor(2);
  xaxis->Draw();

  // Draw y-axis
  TPolyLine3D *yaxis = new TPolyLine3D(2);
  yaxis->SetPoint(0, 0, 0, 0);
  yaxis->SetPoint(1, 0,  fViewRangeMax[1], 0);
  yaxis->SetLineColor(3);
  yaxis->Draw();

  // Draw z-axis
  TPolyLine3D *zaxis = new TPolyLine3D(2);
  zaxis->SetPoint(0, 0, 0, 0);
  zaxis->SetPoint(1, 0, 0, fViewRangeMax[2]);
  zaxis->SetLineColor(4);
  zaxis->Draw();

  // Draw all solids
  for (Int_t i = 0; i < fSolidNumber; i++)
    fSolidPtr[i]->Draw();

  return viewPtr;
}

//______________________________________________________________________________
void TVpGeometry::SetViewRange(const TVpVector3& viewRange)
{
  // Set the data members only.  Deprecated function since ROOT v. 5.16/00.

  for (Int_t i = 0; i < 3; i++)
    {
      fViewRangeMin[i] = -viewRange.fR[i];
      fViewRangeMax[i] = viewRange.fR[i];
    }
}

//______________________________________________________________________________
void TVpGeometry::SetViewRange(const Double_t* rmin, const Double_t* rmax)
{
  // Set the data members only.

  for (Int_t i = 0; i < 3; i++)
    {
      fViewRangeMin[i] = rmin[i];
      fViewRangeMax[i] = rmax[i];
    }
}
