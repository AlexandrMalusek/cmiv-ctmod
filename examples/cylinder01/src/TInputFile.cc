#include "TInputFile.h"
#include <fstream>
#include "TString.h"

//______________________________________________________________________________
void TInputFile::Read(std::istream &in)
{
  TString restOfLine;   // STL strings still don't work in CINT 

  in >> fNumOfPhotons >> restOfLine
     >> fCylinderRadius >> restOfLine
     >> fCylinderHeight >> restOfLine
     >> fBeamWidth >> restOfLine;
}

//______________________________________________________________________________
void TInputFile::Write(std::ostream &out) const
{
  out << "<InputFile>\n"  
      << fNumOfPhotons << '\t' << "#_Number_of_simulated_photons\n"
      << fCylinderRadius << '\t' << "#_Cylinder_radius_in_cm\n"
      << fCylinderHeight << '\t' << "#_Cylinder_height_in_cm\n"
      << fBeamWidth <<  '\t' << "#_Beam_width_in_cm\n"
      << "</InputFile>\n";
}
