#ifndef TVpMaterialManager_h
#define TVpMaterialManager_h

#include "TObject.h"
#include <map>
#include <string>
#include "TVpMaterial.h"

typedef std::map<std::string, TVpMaterialPtr> TVpStringMaterialPtrMap;

class TVpMaterialManager
{
 public:
  Double_t                fMinEnergy;    // min energy in material data energy grid
  Double_t                fMaxEnergy;    // max energy in material data energy grid
  Double_t                fEnergyStep;   // energy step in material data energy grid
  Int_t                   fUseCff;       // use coherent scattering form factors
  Int_t                   fUseIsf;       // use incoherent scattering scattering functions
  TVpStringMaterialPtrMap fMaterialPtr;  //!

  TVpMaterialManager();
  TVpMaterialManager(const Char_t *fileName);
  virtual ~TVpMaterialManager() {};

  Int_t        ReadMMGFile(const Char_t *fileName);
  Int_t        WriteMMGFile(const Char_t *fileName);
  TVpMaterial *RegisterMaterial(const Char_t *fileNameMAT);
  void         RegisterMaterial(const Char_t *name, TVpMaterial *materialPtr);
  TVpMaterial *GetMaterial(const Char_t *fileNameMAT);
  void         PrintStatus(std::ostream &out = std::cout) const;
  void         InitializeMaterialData();

  ClassDef(TVpMaterialManager,1) // Material manager which handles voxel array materials
};

#endif  // TVpMaterialManager_h
