#include <cmath>
#include "TInputFile.h"

#ifndef __CINT__
#include <iostream>
#include "TVpMaterial.h"
#include "TVpMatrix3x3.h"
#include "TVpSolid.h"
#include "TVpSphere.h"
#include "TVpVoxelArray.h"
#include "TVpGeometry.h"
#endif

extern TInputFile inputFile;

TVpGeometry *getGeometry()
{
  ////////////////////////////////////////////////////////////////////////
  // Materials

  TVpMaterial *matVacuum = new TVpMaterial("material/vacuum.mat");  
  matVacuum->Initialize();

  ////////////////////////////////////////////////////////////////////////
  // Solids

  TVpMatrix3x3 matI(1,0,0, 0,1,0, 0,0,1);
  // e3 -> e1
  TVpMatrix3x3 rotMat = transpose(rotationMatrix(TVpVector3(0,0,1), TVpVector3(1,0,0)));

  TVpSphere *sphere0 = new TVpSphere("sphere0", 0, 250);  // The Universe
  TVpVector3 traVec0 = TVpVector3(0.0, 0.0, 0.0);
  sphere0->SetActiveTranslation(&traVec0);
  sphere0->SetMaterial(matVacuum);

  TVpVoxelArray *voxelArray1 = new TVpVoxelArray("phantom/phantom2/phantom2.vam",
  						 "voxelArray", 1);
  voxelArray1->PrintMatTable();
  
  voxelArray1->SetActiveTranslation(new TVpVector3(-inputFile.GetVoxelArrayCenterX(),
						   -inputFile.GetVoxelArrayCenterY(),
						   -inputFile.GetVoxelArrayCenterZ()));

  ////////////////////////////////////////////////////////////////////////
  // Set base and overlap lists

  // sphere0
  TVpSolidPtr sphere0Overlap[] = {voxelArray1, 0};
  sphere0->SetOverlap(sphere0Overlap);
  sphere0->SetBase(0);

  // voxel array
  TVpSolidPtr voxelArray1Base[] = {sphere0, 0};
  voxelArray1->SetOverlap(0);
  voxelArray1->SetBase(voxelArray1Base);

  TVpGeometry *geometry = new TVpGeometry(sphere0);
  
  return geometry;
}
