#ifndef TVpBowTieFilterSomatomDefinition_h
#define TVpBowTieFilterSomatomDefinition_h

#include "TObject.h"
#include "AVpBowTieFilter.h"
#include "TVpBowTieFilter1dTable.h"

class TVpBowTieFilterSomatomDefinition : public TVpBowTieFilter1dTable
{
 public:
  Double_t      fTiSlabThickness; // Ti slab thickness in cm
  TVpMaterial  *fMatTiPtr;        //! filter material, titanium

  TVpBowTieFilterSomatomDefinition();
  TVpBowTieFilterSomatomDefinition(Char_t *fileNameBFT, TVpMaterial *matTiPtr = 0);
  virtual ~TVpBowTieFilterSomatomDefinition();

  virtual Double_t GetTransmission(Double_t energy, const TVpVector3 &directionL) const;
  virtual void     PrintStatus(std::ostream &out = std::cout) const;
  
  ClassDef(TVpBowTieFilterSomatomDefinition,1) // Bowtie filter, Siemens Somatom Definition
};


#endif  // TVpBowTieFilterSomatomDefinition_h
