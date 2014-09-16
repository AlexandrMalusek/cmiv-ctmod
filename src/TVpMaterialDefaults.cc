//______________________________________________________________________________
//
// TVpMaterialDefaults defines default values of internals used in TVpMaterial
// and related classes.  Their arrays must be re-allocated and re-calculated
// in a consistent way and thus all changes should be done via the
// TVpMaterialDefaults class.
//
// The defaults are:
//   fDimIub = 1024;
//   fDimA = 1024;
//   fDimLnSfg = 1024;
//   fDimLnFfg = 1024;
//   fEnergyMinIUB = 10;     // [keV]
//   fEnergyMaxIUB = 150.0;  // [keV]
//
// Due to speed reasons, lin-lin interpolation is used to approximate data
// which are better approximated by a log-log interpolation.  The higher the
// array dimensions given by fDimIub, fDimA, fDimLnSfg, and fDimLnFfg are, the
// more precise the interpolation is.  But, on the other hand, large arrays
// consume more memory and their initialization time is longer.  The material
// data energy range is given by fEnergyMinIUB and fEnergyMaxIUB.  Cross
// sections down to 1 keV are often available but, since characteristic
// photons are not simulated in CTmod, the results in the range from 1 keV to
// 10 keV are not very reliable.  The practical upper energy limit is about
// 200 keV.  The code may work for higher energies but it has not been tested
// in this range.
//
// Data members fDimA, fDimLnSfg, fDimLnFfg, fEnergyMinIUB, and fEnergyMaxIUB
// are explained in TVpMaterialGridData.h
//
// Example:
//
// TVpMaterialDefaults *mdPtr = new TVpMaterialDefaults();
// mdPtr->SetEnergyMinIUB(1.0);  // The new value is 1 keV

#include "TVpMaterialDefaults.h"

ClassImp(TVpMaterialDefaults)

//______________________________________________________________________________
TVpMaterialDefaults::TVpMaterialDefaults()
{
  // Default constructor - set the defaults

  fDimIub = 1024;
  fDimA = 1024;
  fDimLnSfg = 1024;
  fDimLnFfg = 1024;
  fEnergyMinIUB = 10;     // [keV]
  fEnergyMaxIUB = 150.0;  // [keV]
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetDimIub(Int_t dimIub)
{
  // Set the dimension of IUB arrays.  They store cross sections and
  // interaction probability fractions.

  if (dimIub < 2)
    {
      std::cerr << "Error: TVpMaterialDefaults::SetDimIub: dimIub < 2, dimIub = "
		<< dimIub << ". Setting dimIub = 1024\n";
      dimIub = 1024;
    }
  fDimIub = dimIub;
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetDimA(Int_t dimA)
{
  // Set the dimension of the A array which used by the coherent scattering
  // sampling routine.  A is derived from form factors.

  fDimA = dimA;
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetDimLnSfg(Int_t dimLnSfg)
{
  // Set the dimension of the LnSfg arrays which approximate the scattering
  // function of the incoherent scattering.

  fDimLnSfg = dimLnSfg;
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetDimLnFfg(Int_t dimLnFfg)
{
  // Set the dimension of the LnFfg arrays which approximate form factors of
  // the coherent scattering.

  fDimLnFfg = dimLnFfg;
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetEnergyMinIUB(Double_t energyMinIUB)
{
  // Set the lower cutoff value of material data cross sections.

  fEnergyMinIUB = energyMinIUB;
}

//______________________________________________________________________________
void TVpMaterialDefaults::SetEnergyMaxIUB(Double_t energyMaxIUB)
{
  // Set the upper cutoff value of material data cross sections.

  fEnergyMaxIUB = energyMaxIUB;
}

//______________________________________________________________________________
void TVpMaterialDefaults::PrintStatus(std::ostream &out) const
{
  // Print the object status.

  out << "<TVpMaterialDefaults>\n"
      << "$Id: TVpMaterialDefaults.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "DimIub: " << fDimIub << '\n'
      << "DimA: " << fDimA << '\n'
      << "DimLnSfg:" << fDimLnSfg << '\n'
      << "DimLnFfg:" << fDimLnFfg << '\n'
      << "EnergyMinIUB" << fEnergyMinIUB << '\n'
      << "EnergyMaxIUB" << fEnergyMaxIUB << '\n'
      << "</TVpMaterialDefaults>\n";
}
