//______________________________________________________________________________
//
// TVpSolidNeighbor references a solid and its spatial orientation.  An array
// of TVpSolidNeighbor objects is used to describe base and overlapping solids
// of a given solid.
//
// If fTraVecC2nPtr = 0 then a null vector, (0,0,0), is used.
// If fRotMatC2nPtr = 0 then an identity matrix, {1,0,0; 0,1,0; 0,0,1}, is used.

#include <stdio.h>
#include "TVpSolidNeighbor.h"

ClassImp(TVpSolidNeighbor)

//______________________________________________________________________________
TVpSolidNeighbor::TVpSolidNeighbor()
{
  // Default constructor.  Initialize pointers to 0.
  
  fSolidPtr = 0;
  fTraVecC2nPtr = 0;
  fRotMatC2nPtr = 0;
}

//______________________________________________________________________________
TVpSolidNeighbor::~TVpSolidNeighbor()
{
  // Destructor.  Delete transformation matrices.
  
  delete fTraVecC2nPtr;
  delete fRotMatC2nPtr;
}

//______________________________________________________________________________
void TVpSolidNeighbor::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  
  out << "<TVpSolidNeighbor>\n"
      << "$Id: TVpSolidNeighbor.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Name:\t" << fSolidPtr->fName << '\n'
      << "Index:\t" << fSolidPtr->fIndex << '\n'
      << "Transformation matrices are defined so that xNew = R*xCur + T\n"
      << "Translation vector T:\n";
  
  if (fTraVecC2nPtr == 0)
    out << TVpVector3(0, 0, 0) << " , fTraVec == 0\n";
  else
    out << *fTraVecC2nPtr << '\n';

  out << "Rotation matrix R:\n";
  if (fRotMatC2nPtr == 0)
    out << TVpMatrix3x3(1,0,0, 0,1,0, 0,0,1) << "fRotMat == 0\n";
  else
    out << *fRotMatC2nPtr;
  out << "</TVpSolidNeighbor>\n";
}
