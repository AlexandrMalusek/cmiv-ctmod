// Cylindrical phantom

#include <cmath>
#include "TInputFile.h"

#ifndef __CINT__
#include <iostream>
#include "TVpMaterial.h"
#include "TVpMatrix3x3.h"
#include "TVpSolid.h"
#include "TVpSphere.h"
#include "TVpCylinder.h"
#include "TVpGeometry.h"
#include "TVpMaterialManager.h"
#endif

extern TInputFile inputFile;

TVpGeometry *getGeometry(TVpMaterialManager *mmPtr)
{
  ////////////////////////////////////////////////////////////////////////
  // Materials


  ////////////////////////////////////////////////////////////////////////
  // Solids

  TVpSphere *sphere0Ptr = new TVpSphere("sphere0", 0, 110);  // The Universe
  sphere0Ptr->SetMaterial(mmPtr->GetMaterial("vacuum"));

  Double_t cylHeight = inputFile.GetCylinderHeight();
  Double_t cylRadius = inputFile.GetCylinderRadius();
  TVpMatrix3x3 rotMat = rotationMatrix(TVpVector3(0,0,1), TVpVector3(1,0,0));
  TVpCylinder *cylinder1Ptr = new TVpCylinder
    ("cylinder1", // descriptive name
     1,           // index
     cylRadius,   // radius in cm
     cylHeight);  // height in cm
  TVpVector3 traVec1 = TVpVector3(-cylHeight/2, 0.0, 0.0);
  cylinder1Ptr->SetActiveTransformation(&rotMat,&traVec1);
  cylinder1Ptr->SetMaterial(mmPtr->GetMaterial("water"));

  ////////////////////////////////////////////////////////////////////////
  // Set base and overlap lists

  // sphere0
  TVpSolidPtr sphere0Overlap[] = {cylinder1Ptr, 0};
  sphere0Ptr->SetOverlap(sphere0Overlap);
  sphere0Ptr->SetBase(0);

  // cylinder1
  TVpSolidPtr cylinder1Base[] = {sphere0Ptr, 0};
  cylinder1Ptr->SetOverlap(0);
  cylinder1Ptr->SetBase(cylinder1Base);

  TVpGeometry *geometryPtr = new TVpGeometry(sphere0Ptr);

  return geometryPtr;
}
