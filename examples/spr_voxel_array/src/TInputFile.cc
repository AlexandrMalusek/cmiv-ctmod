#include "TInputFile.h"
#include <fstream>
#include "TString.h"

//______________________________________________________________________________
void TInputFile::Read(std::istream &in)
{
  TString restOfLine;   // STL strings still don't work in CINT 

  in >> fNumOfPhotons >> restOfLine
     >> fEnergy >> restOfLine
     >> fVoxelArrayCenterX >> restOfLine
     >> fVoxelArrayCenterY >> restOfLine
     >> fVoxelArrayCenterZ >> restOfLine
     >> fBeamWidth >> restOfLine;
}

//______________________________________________________________________________
void TInputFile::Write(std::ostream &out) const
{
  out << "<InputFile>\n"  
      << fNumOfPhotons << '\t' << "# Number of simulated photons\n"
      << fEnergy << '\t' << "# Energy in keV\n"
      << fVoxelArrayCenterX << '\t' << "# Voxel array center X in cm\n"
      << fVoxelArrayCenterY << '\t' << "# Voxel array center Y in cm\n"
      << fVoxelArrayCenterZ << '\t' << "# Voxel array center Z in cm\n"
      << fBeamWidth <<  '\t'<< "#_Beam_width_in_cm\n"
      << "</InputFile>\n";
}
