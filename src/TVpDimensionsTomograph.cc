//______________________________________________________________________________
//
// TVpDimensionsTomograph defines dimensions of a CT scanner.  The purpose of
// this class is to provide a centralized storage of technical data.  This is
// especially useful for interactive scripts.
//
// Names are converted according to CTmod's conventions, for example: SDD, the
// source-axis distance, is named fSdd in the case of the data memeber.  The
// corresponding getter and setter are GetSdd() and SetSdd(), respectively.
//
// Front view:
//Begin_Html
/*
<img src="png/geometry_1.png">
*/
//End_Html
// Side view:
//Begin_Html
/*
<img src="png/geometry_2.png">
*/
//End_Html
// Top view:
//Begin_Html
/*
<img src="png/geometry_3.png">
*/
//End_Html
//
// An example:
// TVpDimensionsTomograph *dtPtr = new TVpDimensionsTomograph();
// dtPtr->SetNdl(32);
// dtPtr->SetNdw(32);
// dtPtr->SetSad(70.0);
// dtPtr->SetSdd(100.0);
// dtPtr->SetDal(90.0);
// dtPtr->SetDaw(40.0);

#include <string>
#include "TVpDimensionsTomograph.h"

ClassImp(TVpDimensionsTomograph)

//______________________________________________________________________________
TVpDimensionsTomograph::TVpDimensionsTomograph()
{
  // Default constructor.  Initialize values to zero.
  
  fNdl = fNdw = 0;
  fSad = fSdd = fDal = fDaw = 0.0;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetNdl(Int_t ndl)
{
  // Set the number of detector elements in the length direction

  fNdl = ndl;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetNdw(Int_t ndw)
{
  // Set the number of detector elements in the width direction

  fNdw = ndw;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetSad(Double_t sad)
{
  // Set the source - axis distance in cm.

  fSad = sad;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetSdd(Double_t sdd)
{
  // Set the source - detector distance in cm.

  fSdd = sdd;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetDal(Double_t dal)
{
  // Set the detector array length in cm.

  fDal = dal;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::SetDaw(Double_t daw)
{
  // Set the detector array width in cm.

  fDaw = daw;
}

//______________________________________________________________________________
void TVpDimensionsTomograph::Read(std::istream &in)
{
  std::string restOfLine;

  std::getline(in, restOfLine);
  in >> fNdl;  std::getline(in, restOfLine);
  in >> fNdw;  std::getline(in, restOfLine);
  in >> fSad;  std::getline(in, restOfLine);
  in >> fSdd;  std::getline(in, restOfLine);
  in >> fDal;  std::getline(in, restOfLine);
  in >> fDaw;  std::getline(in, restOfLine);
}

//______________________________________________________________________________
void TVpDimensionsTomograph::Write(std::ostream &out) const
{
  out << "<TVpDimensionsTomograph>\n"  
      << fNdl << '\t' << "# NDL - num. of det. elements in the perp. dir.\n"
      << fNdw << '\t' << "# NDW - num. of det. elements in the axial dir.\n"
      << fSad << '\t' << "# SAD - source-axis distance in cm\n"
      << fSdd << '\t' << "# SDD - source-detector distance in cm\n"
      << fDal << '\t' << "# DAL - det. array length in the  perp. dir. in cm\n"
      << fDaw << '\t' << "# DAW - det. array width (in the axial dir.) in cm\n"
      << "</TVpDimensionsTomograph>\n";
}
