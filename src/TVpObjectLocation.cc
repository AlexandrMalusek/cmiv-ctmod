//______________________________________________________________________________
//
// TVpObjectLocation defines a position of an object in space.

#include <iostream>
#include "TVpObjectLocation.h"

ClassImp(TVpObjectLocation)

//______________________________________________________________________________
TVpObjectLocation::TVpObjectLocation()
{
  // Constructor.  Set the active translation vector equal to the null vector
  // and the active rotation matrix equal to the identity matrix.

  fTraVecL2u = TVpVector3(0, 0, 0);
  fRotMatL2u = TVpMatrix3x3(1,0,0, 0,1,0, 0,0,1);
  fRotMatL2uT = TVpMatrix3x3(1,0,0, 0,1,0, 0,0,1); // = transpose(fRotMatL2u)
}

//______________________________________________________________________________
TVpObjectLocation::~TVpObjectLocation()
{
  // Destructor.

}

//______________________________________________________________________________
void TVpObjectLocation::SetActiveTranslation(TVpVector3 *traVec)
{
  // Set the translation vector, i.e. the position of the origin of the
  // solid's LCS in UCS.
  //
  // Input parameters:
  // - traVec - position of the solid's LCS in UCS
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // root [] box->SetActiveTranslation(new TVpVector3(2.0, 4.0, 8.0));

  fTraVecL2u = *traVec;
}

//______________________________________________________________________________
void TVpObjectLocation::SetActiveRotation(TVpMatrix3x3 *rotMat)
{
  // Set the rotation matrix, i.e. the orientation of the solid in UCS.  The
  // matrix must be unitary, no check is performed.
  //
  // Input parameters:
  // - rotMat - active rotation which transforms coordinate axes of the LCS of
  // the solid.
  //
  // Example:
  // root [] TVpBox *box = new TVpBox("box", 1, 10, 15, 20);
  // // 1. Right-hand (positive) Pi rotation about the z-axis.
  // // 2. Right-hand (positive) Pi/2 rotation about the y-axis.
  // root [] box->SetActiveRotation(new TVpMatrix3x3(TMath::Pi(), TMath::Pi()/2, 0));

  fRotMatL2u = *rotMat;
  fRotMatL2uT = transpose(fRotMatL2u);
}

//______________________________________________________________________________
void TVpObjectLocation::SetActiveTransformation(TVpMatrix3x3 *rotMat, TVpVector3 *traVec)
{
  // Set the object position.  See SetActiveRotation() and
  // SetActiveTranslation().
  //
  // Input parameters:
  // - rotMat - active rotation
  // - traVec - active translation

  SetActiveRotation(rotMat);
  SetActiveTranslation(traVec);
}

//______________________________________________________________________________
TVpVector3 TVpObjectLocation::PosLocToUni(const TVpVector3& posLoc) const
{
  // Return position of a point in the UCS.  The point is specified in the
  // LCS of the current solid.
  //
  // Input parameters:
  // - posUni - position of the point in the UCS

  TVpVector3 posUni = fRotMatL2u * posLoc + fTraVecL2u;
  return posUni;
}

//______________________________________________________________________________
TVpVector3 TVpObjectLocation::DirLocToUni(const TVpVector3& dirLoc) const
{
  // Return direction of a free vector in UCS.  The direction is specified in
  // LCS of the current solid.
  //
  // Input parameters:
  // - dirLoc - direction of the free vector in the LCS

  TVpVector3 dirUni = fRotMatL2u * dirLoc;
  return dirUni;
}

//______________________________________________________________________________
TVpVector3 TVpObjectLocation::PosUniToLoc(const TVpVector3& posUni) const
{
  // Return position of a point in the LCS of the current solid.  The point is
  // specified in the UCS.
  //
  // Input parameters:
  // - posUni - position of the point in the UCS

  TVpVector3 posLoc = fRotMatL2uT * (posUni - fTraVecL2u);
  return posLoc;
}

//______________________________________________________________________________
TVpVector3 TVpObjectLocation::DirUniToLoc(const TVpVector3& dirUni) const
{
  // Return direction of a free vector in the LCS of the current solid.  The
  // direction is specified in the UCS.
  //
  // Input parameters:
  // - dirUni - direction of the free vector in the UCS

  TVpVector3 dirLoc = fRotMatL2uT * dirUni;
  return dirLoc;
}

//______________________________________________________________________________
void TVpObjectLocation::PrintLocation(std::ostream &out) const
{
  // Print translation and rotation matrices
  
  out << "<TVpObjectLocation>\n"
      << "TraVecL2u: " << fTraVecL2u << '\n'
      << "RotMatL2u:\n" << fRotMatL2u
      << "</TVpObjectLocation>\n";
}

#include "TPolyLine3D.h"

//______________________________________________________________________________
void TVpObjectLocation::DrawLcsAxes(Double_t lengthX, Double_t lengthY, Double_t lengthZ) const
{
  // Draw LCS axes.
  //
  // Input parameters:
  // lengthX, lengthY, lengthZ - lengths of line segments

  TVpVector3 up;

  TPolyLine3D *ex = new TPolyLine3D(2);
  up = PosLocToUni(TVpVector3(0, 0, 0));
  ex->SetPoint(0, up.fR[0], up.fR[1], up.fR[2]);
  up = PosLocToUni(TVpVector3(lengthX, 0, 0));
  ex->SetPoint(1, up.fR[0], up.fR[1], up.fR[2]);
  ex->SetLineColor(2);
  ex->Draw();

  TPolyLine3D *ey = new TPolyLine3D(2);
  up = PosLocToUni(TVpVector3(0, 0, 0));
  ey->SetPoint(0, up.fR[0], up.fR[1], up.fR[2]);
  up = PosLocToUni(TVpVector3(0, lengthY, 0));
  ey->SetPoint(1, up.fR[0], up.fR[1], up.fR[2]);
  ey->SetLineColor(3);
  ey->Draw();

  TPolyLine3D *ez = new TPolyLine3D(2);
  up = PosLocToUni(TVpVector3(0, 0, 0));
  ez->SetPoint(0, up.fR[0], up.fR[1], up.fR[2]);
  up = PosLocToUni(TVpVector3(0, 0, lengthZ));
  ez->SetPoint(1, up.fR[0], up.fR[1], up.fR[2]);
  ez->SetLineColor(4);
  ez->Draw();
}
