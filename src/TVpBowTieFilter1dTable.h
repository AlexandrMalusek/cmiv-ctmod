#ifndef TVpBowTieFilter1dTable_h
#define TVpBowTieFilter1dTable_h

#include "TObject.h"
#include <string>
#include "AVpBowTieFilter.h"
#include "TVpMaterial.h"

class TVpBowTieFilter1dTable : public AVpBowTieFilter
{
 protected:
  Double_t      fFanAngleCosStep;// Step length in equidistant grid.  Undefined otherwise.

 public:
  std::string   fName;          // Descriptive name
  Int_t         fFormatVersion; // Format version
  Int_t         fDim;           // Number of points
  Double_t     *fFanAngleCos;   //! Cosine of the fan angle
  Double_t     *fThickness;     //! filter thickness in cm
  TVpMaterial  *fMaterialPtr;   //! filter material, water by default

  TVpBowTieFilter1dTable();
  TVpBowTieFilter1dTable(const Char_t *fileNameBFT);
  virtual ~TVpBowTieFilter1dTable();

  Double_t         GetThickness(Double_t energy, const TVpVector3 &directionL) const;
  virtual Double_t GetTransmission(Double_t energy, const TVpVector3 &directionL) const;
  virtual void     PrintStatus(std::ostream &out = std::cout) const;
  Int_t            ReadBftFile(const Char_t *fileNameBFT);

  TGraph          *GetBftThickness() const;

  ClassDef(TVpBowTieFilter1dTable,1) // Bowtie filter, transmission data in 1D table
};


#endif  // TVpBowTieFilter1dTable_h
