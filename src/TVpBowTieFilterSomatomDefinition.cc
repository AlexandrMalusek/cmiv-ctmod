//______________________________________________________________________________
// The bowtie filter of the Siemens Somatom Definition CT scanner consist of a
// 0.6 mm thick Ti plate and an aluminum formfilter.  The latter is defined via
// a table in a BFT file.
//
// data/bowtie_filters/Siemens_SomatomDefinition_fixed.bft
// data/bowtie_filters/Siemens_SomatomDefinition_fixedAndMovable.bft
//
//______________________________________________________________________________

#include <math.h>
#include "TVpBowTieFilterSomatomDefinition.h"

ClassImp(TVpBowTieFilterSomatomDefinition)

//______________________________________________________________________________
TVpBowTieFilterSomatomDefinition::TVpBowTieFilterSomatomDefinition()
{
  // Default constructor.  Initialize default dimensions, set pointers to
  // material classes to 0.

  fMatTiPtr = 0;
  fTiSlabThickness = 0.06;  // cm
  fTiSlabThickness = 0.0;
}

//______________________________________________________________________________
TVpBowTieFilterSomatomDefinition::TVpBowTieFilterSomatomDefinition(Char_t *fileNameBFT,
								   TVpMaterial *matTiPtr)
{
  // Create the bowtie filter object.
  //
  // Input parameters:
  // - fileNameBFT - BFT file name
  // - matTiPtr - pointer to titanium

  fMatTiPtr = (matTiPtr != 0) ? matTiPtr : new TVpMaterial("material/Ti.mat");
  fMatTiPtr->Initialize();
  fTiSlabThickness = 0.06;  // cm
  ReadBftFile(fileNameBFT);
}

//______________________________________________________________________________
TVpBowTieFilterSomatomDefinition::~TVpBowTieFilterSomatomDefinition()
{
  // Destructor.  Instances of material classes are not deleted.
}

//______________________________________________________________________________
Double_t TVpBowTieFilterSomatomDefinition::GetTransmission(Double_t energy,
						  const TVpVector3 &directionL) const
{
  // Calculate the transmission in a given direction and for a given photon
  // energy.
  //
  // Input parameters:
  // - energy - photon energy in keV
  // - directionL - direction from source
  //
  // Method:
  // The transmission is defined as the probability of passing the bow-tie filter
  // without any interaction, i.e. exp(-radiological_path).

  Double_t thicknessAl = GetThickness(energy, directionL);
  Double_t thicknessTi = fTiSlabThickness / directionL.GetZ();
  Double_t transmission = exp(-fMaterialPtr->GetLac(energy)*thicknessAl -
			      fMatTiPtr->GetLac(energy)*thicknessTi);
  return transmission;
}

//______________________________________________________________________________
void TVpBowTieFilterSomatomDefinition::PrintStatus(std::ostream &out) const
{
  // Print the object's status
  
  out << "<TVpBowTieFilterSomatomDefinition>\n"
      <<  "<TVpBowTieFilterSomatomDefinition>\n";
}
