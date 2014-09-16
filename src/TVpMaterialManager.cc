//______________________________________________________________________________
//
// TVpMaterialManager manages materials.

#include <fstream>
#include <iostream>
#include "TVpMaterialManager.h"

ClassImp(TVpMaterialManager)

//______________________________________________________________________________
TVpMaterialManager::TVpMaterialManager()
{
  // Constructor

  TVpMaterialDefaults matDef;
  
  fMinEnergy = matDef.GetEnergyMinIUB();
  fMaxEnergy = matDef.GetEnergyMaxIUB();
  fEnergyStep = (fMaxEnergy - fMinEnergy) / matDef.GetDimIub();
  fUseCff = 1;
  fUseIsf = 1;
}

//______________________________________________________________________________
TVpMaterialManager::TVpMaterialManager(const Char_t *fileName)
{
  // Constructor

  ReadMMGFile(fileName);
}

//______________________________________________________________________________
Int_t TVpMaterialManager::ReadMMGFile(const Char_t *fileName)
{
  // Read MMG file

  std::string recName;                 // The record name e.g. "# Size"
  std::string recRest;                 // The rest of the record

  std::ifstream mmg(fileName);
  if (!mmg)
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Cannot open file: " << fileName << '\n';
      return 1;
    }

  // Read and process the file header
  // # Format:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Format"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Format:\" is missing.\n";
      return 2;
    }
  std::getline(mmg, recRest);
  if (recRest != std::string(" MMG 2.0"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Format version should be:\n"
	   << "\" MMG 2.0\"" << " but is:\n"
	   << '"' << recRest << "\"\n";
      return 3;
    }

  // # Min energy:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Min energy"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Min energy:\" is missing.\n";
      return 4;
    }
  mmg >> fMinEnergy;
  std::getline(mmg, recRest);

  // # Max energy:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Max energy"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Max energy:\" is missing.\n";
      return 5;
    }
  mmg >> fMaxEnergy;
  std::getline(mmg, recRest);

  // # Energy step:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Energy step"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Energy step:\" is missing.\n";
      return 6;
    }
  mmg >> fEnergyStep;
  std::getline(mmg, recRest);

  // # Use CFF:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Use CFF"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Use CFF:\" is missing.\n";
      return 7;
    }
  mmg >> fUseCff;
  std::getline(mmg, recRest);

  // # Use ISF:
  std::getline(mmg, recName, ':');
  if (recName != std::string("# Use ISF"))
    {
      std::cerr << "Error: VoxelArray::ReadMMGFile: Record \"# Use ISF:\" is missing.\n";
      return 8;
    }
  mmg >> fUseIsf;
  std::getline(mmg, recRest);

  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterialManager::WriteMMGFile(const Char_t *fileName)
{
  // Write the MMG file.
  
  std::ofstream mmg(fileName);
  if (!mmg)
    {
      std::cerr << "TVpVoxelArray::WriteMMGFile: Cannot open file: " << fileName << '\n';
      return 1;
    }
  
  mmg << "# Format: MMG 2.0\n";
  mmg << "# Min energy: " << fMinEnergy << " keV\n";
  mmg << "# Max energy: " << fMaxEnergy << " keV\n";
  mmg << "# Energy step: " << fEnergyStep << " keV\n";
  mmg << "# Use CFF: " << fUseCff << '\n';
  mmg << "# Use ISF: " << fUseIsf << '\n';

  if (!mmg.good())
    {
      std::cerr << "Error: TVpVoxelArray::WriteMMGFile: An output error occured.\n";
      return 2;
    }

  return 0;
}

//______________________________________________________________________________
TVpMaterial *TVpMaterialManager::RegisterMaterial(const Char_t *fileNameMAT)
{
  // Register material file

  const Char_t *fileNameCff, *fileNameIsf;
  std::string strFileNameMat = std::string(fileNameMAT);
  std::string::size_type idx = strFileNameMat.find('.');
  std::string strBaseName = strFileNameMat.substr(0, idx);
  std::string strFileNameCff = strBaseName + std::string(".cff");
  std::string strFileNameIsf = strBaseName + std::string(".isf");
  
  fileNameCff = (fUseCff) ? strFileNameCff.c_str() : 0;
  fileNameIsf = (fUseIsf) ? strFileNameIsf.c_str() : 0;

  TVpMaterial *mptr;
  if ((mptr = fMaterialPtr[fileNameMAT]) == 0)
    {
      // Register the new material
      mptr = fMaterialPtr[fileNameMAT] = 
	new TVpMaterial(fileNameMAT, fileNameCff, fileNameIsf);

      // Set material defaults
      TVpMaterialDefaults matDef;

      matDef.fDimIub = (int) ((fMaxEnergy - fMinEnergy) / fEnergyStep);
      matDef.fEnergyMinIUB = fMinEnergy;
      matDef.fEnergyMaxIUB = fMinEnergy +  matDef.fDimIub * fEnergyStep;

      mptr->Initialize(matDef);
      mptr->SetUseFf(fUseCff);
      mptr->SetUseSf(fUseIsf);
    }
  
  return mptr;
}

//______________________________________________________________________________
void TVpMaterialManager::RegisterMaterial
(const Char_t *materialName, TVpMaterial *materialPtr)
{
  // Register material.
  //
  // Input parameters:
  // - materialName - unique name of the material.
  // - materialPtr - pointer to an existing material

  if (fMaterialPtr[materialName] == 0)
    fMaterialPtr[materialName] = materialPtr;  // Register the new material
  else
    {
      std::cerr << "Error: TVpMaterialManager::RegisterMaterial: "
		<< "materialName = " << materialName
		<< " has already been registered.\n";
      return;
    }
}

//______________________________________________________________________________
TVpMaterial *TVpMaterialManager::GetMaterial(const Char_t *fileNameMAT)
{
  // Return a pointer to the material defined by the file fileNameMAT

  return fMaterialPtr[fileNameMAT];
}

//______________________________________________________________________________
void TVpMaterialManager::PrintStatus(std::ostream &out) const
{
  // Print status.

  typedef TVpStringMaterialPtrMap::const_iterator CI;
  TVpMaterial *matPtr;

  out << "<TVpMaterialManager>\n";
  for (CI p = fMaterialPtr.begin(); p!= fMaterialPtr.end(); ++p)
    {
      matPtr = p->second;
      out << "Material name: \"" << p->first << "\"\n"
	  << "  MAT file: \"" << matPtr->fFileNameMAT << "\"\n"
	  << "  CFF file: \"" << matPtr->fFileNameCFF << "\"\n"
	  << "  ISF file: \"" << matPtr->fFileNameISF << "\"\n";
    }
  out << "</TVpMaterialManager>\n";
}

//______________________________________________________________________________
void TVpMaterialManager::InitializeMaterialData()
{
  // Initialize material data
  
  typedef TVpStringMaterialPtrMap::const_iterator CI;
  std::cerr << "<Info: TVpMaterialManager::InitializeMaterialData: started.>\n";
  for (CI p = fMaterialPtr.begin(); p!= fMaterialPtr.end(); ++p)
    {
      std::cerr << p->first << "...\n";
      p->second->Initialize();
    }
  std::cerr << "<Info: TVpMaterialManager::InitializeMaterialData: finished.>\n";
}
