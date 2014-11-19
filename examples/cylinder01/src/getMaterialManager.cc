// Cylindrical phantom

#ifndef __CINT__
#include "TVpMaterial.h"
#include "TVpMaterialManager.h"
#endif

TVpMaterialManager *getMaterialManager()
{
  TVpMaterialManager *materialManagerPtr = new TVpMaterialManager;

  TVpMaterial *mat_vacuum_Ptr = new TVpMaterial("material/vacuum.mat");
  materialManagerPtr->RegisterMaterial("vacuum", mat_vacuum_Ptr);
  
  TVpMaterial *mat_water_Ptr = new TVpMaterial
    ("material/water.mat",
     "material/water_m80901.cff",
     "material/water.isf");
  materialManagerPtr->RegisterMaterial("water", mat_water_Ptr);

  return materialManagerPtr;
}
