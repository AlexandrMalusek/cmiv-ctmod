//______________________________________________________________________________
//
// TVpMaterial class defines methods which return material density, cross
// section data, and histograms of various quantities.  It provides material
// related data to other classes via functions which were inherited from
// TVpMaterialGridData and TVpMaterialFileData classes.
//
// Simulations use the following functions:
// - TVpMaterialFileData::GetDensity()
// - TVpMaterialGridData::GetToCsg()
// - TVpMaterialGridData::GetSfg()
// - TVpMaterialGridData::GetInIub()
// - TVpMaterialGridData::GetPhIub()
// - TVpMaterial::GetA()
// - TVpMaterial::GetInvA()
// - TVpMaterial::GetInPdfOmega()
// - TVpMaterial::GetCoPdfOmega()
//
// ******
// Plots:
// ******
//
// Probability density functions of coherent and incoherent scattering:
// - GetHistCoPdfTheta(Double_t energy, Int_t numTheta)
// - GetHistInPdfTheta(Double_t energy, Int_t numTheta)
// - GetHistCoPdfOmega(Double_t energy, Int_t numTheta)
// - GetHistInPdfOmega(Double_t energy, Int_t numTheta)
//
// Cross sections (cross sections in barns):
// - GetHistKnCsTheta(Double_t energy, Int_t numTheta);  // Klein-Nishina
// - GetHistInCsTheta(Double_t energy, Int_t numTheta);
// - GetHistCoCsTheta(Double_t energy, Int_t numTheta);
// Notes:
// 1. Absolute values of the differential cross section are correct only if
// absolute values of the form factors or incoherent scattering functions are
// correct.  For compounds and mixtures only relative values are often
// available.  CTmod does not know which normalization procedure was used and
// therefore it always writes out a warning message.
// 2. Absolute values of differential cross sections are not used in MC
// simulations.
//
//
// CINT may have problems in an interactive session with commands which extend
// over several lines, use the block { } in this case.
//
// Example:
// {
//   gSystem->Load("libRCTmod.so");           // load the library 
//   TVpMaterial *matWater = new TVpMaterial
//     ("material/water.mat",                 // cross sections
//      "material/water_m80901.cff",          // form factors
//      "material/water.isf");                // incoherent scattering function
//   matWater->Initialize();                  // creation of internal data tables
// }

#include "TVpMaterial.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "misc.h"
#include "TVpMath.h"

ClassImp(TVpMaterial)

//______________________________________________________________________________
TVpMaterial::TVpMaterial()
  : TVpMaterialFileData(), TVpMaterialGridData()
{
  // Default constructor. Initialize data members to 0

  fValU = fValA = fPolyA0 = 0;
}

//______________________________________________________________________________
TVpMaterial::TVpMaterial(const Char_t *fileNameMat, const Char_t *fileNameCff,
			 const Char_t *fileNameIsf)
  : TVpMaterialFileData(), TVpMaterialGridData()
{
  // Constructor which loads cross sections, form factors and the incoherent
  // scattering function from data files.  See
  // TVpMaterialFileData::ReadMatFile(), TVpMaterialFileData::ReadCffFile(),
  // and TVpMaterialFileData::ReadIsfFile() for format description.  Info
  // messages are printed when reading starts.
  //
  // Input:
  // - fileNameMat - file name of a file in MAT format containing cross sections
  // - fileNameCff - file name of a file in CFF format containing form factors
  // - fileNameIsf - file name of a file in ISF format containing scattering function
  //
  // Example:
  // {
  //   TVpMaterial *matWater = new TVpMaterial
  //     ("material/water.mat",
  //      "material/water_m80901.cff",
  //      "material/water.isf");
  // }
  // <Info TVpMaterialFileData::ReadMatFile(material/water.mat)>
  // <Info TVpMaterialFileData::ReadCffFile(material/water_m80901.cff)>
  // <Info TVpMaterialFileData::ReadIsfFile(material/water.isf)>

  fValU = fValA = fPolyA0 = 0;  // Form factor related data

  ReadMatFile(fileNameMat);
  if (fileNameCff != 0)
    ReadCffFile(fileNameCff);
  if (fileNameIsf != 0)
    ReadIsfFile(fileNameIsf);
}

//______________________________________________________________________________
TVpMaterial::~TVpMaterial()
{
  // Destructor. Delete all allocated arrays.

  delete [] fValU;
  delete [] fValA;
  delete [] fPolyA0;
}

//______________________________________________________________________________
Int_t TVpMaterial::Initialize()
{
  // Initialize TVpMaterial internals with default values.  Default values are
  // defined by the class TVpMaterialDefaults.

  TVpMaterialDefaults materialDefaults;  // filled with defaults
  return Initialize(materialDefaults);
}

//______________________________________________________________________________
Int_t TVpMaterial::Initialize(TVpMaterialDefaults& materialDefaults)
{
  // Initialize TVpMaterial internals with user-supplied values and return 0
  // if no errors were detected.  Retun a non-zero value otherwise.  See the
  // source code and the class TVpMaterialDefaults for more details.
  //
  // Input:
  // - materialDefaults - material defaults
  //
  // Method:
  // This function initializes internal arrays which are used during
  // simulations.  It is CPU intensive and therefore it is not performed
  // automatically when material data are read from files.

  fDimA         = materialDefaults.fDimA;
  fDimIub       = materialDefaults.fDimIub;
  fDimA         = materialDefaults.fDimA;
  fDimLnSfg     = materialDefaults.fDimLnSfg;
  fDimLnFfg     = materialDefaults.fDimLnFfg;
  fEnergyMinIUB = materialDefaults.fEnergyMinIUB;
  fEnergyMaxIUB = materialDefaults.fEnergyMaxIUB;
  
  // Csg and Iub initialization (grid cross section and upper-boundary arrays)
  // ************************************************************************
  if (fEnergyMinIUB < fEnergyMinCsf)
    {
      std::cerr << "Error: TVpMaterial::Initialize: "
		<< "fEnergyMinIUB < fEnergyMinCsf for material "
		<< fName << '\n';
      return 1;
    }
  
  if (fEnergyMaxIUB > fEnergyMaxCsf)
    {
      std::cerr << "Error: TVpMaterial::Initialize: "
		<< "fEnergyMaxIUB > fEnergyMaxCsf for material "
		<< fName << '\n';
      return 2;
    }
  InitializeCsgIub();

  // Sf initialization (scattering function arrays)
  if (fUseSf != 0)
    {
      if (fDimLnSfg <= 2)
	{
	  std::cerr << "Error: TVpMaterial::Initialize: "
		    << "fDimLnSfg = " << fDimLnSfg << ", it must be > 2";
	  return 3;
	}
      InitializeSfg();
    }
  
  // Ff initialization (form factor arrays)
  if (fUseFf != 0)
    {
      if (fDimLnFfg <= 2)
	{
	  std::cerr << "Error: TVpMaterial::InitializeInSFG: fDimLnFfg = " << fDimLnFfg
		    << ", it must be > 2";
	  return 4;
	}

      InitializeFfg();
      InitializeATable();
    }

  // Normalization factors initialization (Used by pdfs of In and Co)
  InitializeNog();
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterial::InitializeCsgIub()
{
  // Set Csg and Iub arrays.  Old arays are deleted.

  fChannelToEnergy = (fEnergyMaxIUB - fEnergyMinIUB) / fDimIub;
  fEnergyToChannel = 1.0 / fChannelToEnergy;

  // De-allocate old arrays to avoid memory leaks.  It's safe because these
  // pointers are initialized to 0 by the constructor.
  delete [] fEnergyIUB;
  delete [] fToCsg0; delete [] fToCsg1;
  delete [] fInIub0; delete [] fInIub1;
  delete [] fPhIub0; delete [] fPhIub1;
  delete [] fCoIub0; delete [] fCoIub1;
  delete [] fInCsg0; delete [] fInCsg1;
  delete [] fPhCsg0; delete [] fPhCsg1;
  delete [] fCoCsg0; delete [] fCoCsg1;

  // Allocate data arrays
  fEnergyIUB = new Double_t[fDimIub];
  fToCsg0 = new Double_t[fDimIub]; fToCsg1 = new Double_t[fDimIub];
  fInIub0 = new Double_t[fDimIub]; fInIub1 = new Double_t[fDimIub];
  fPhIub0 = new Double_t[fDimIub]; fPhIub1 = new Double_t[fDimIub];
  fCoIub0 = new Double_t[fDimIub]; fCoIub1 = new Double_t[fDimIub];
  fInCsg0 = new Double_t[fDimIub]; fInCsg1 = new Double_t[fDimIub];
  fPhCsg0 = new Double_t[fDimIub]; fPhCsg1 = new Double_t[fDimIub];
  fCoCsg0 = new Double_t[fDimIub]; fCoCsg1 = new Double_t[fDimIub];

  // Check if the allocation is OK
  if ( fEnergyIUB == 0 || fToCsg0 == 0 || fToCsg1 == 0 ||
       fInIub0 == 0 || fInIub1 == 0 || fPhIub0 == 0 || fPhIub1 == 0 ||
       fCoIub0 == 0 || fCoIub1 == 0 || fInCsg0 == 0 || fInCsg1 == 0 ||
       fPhCsg0 == 0 || fPhCsg1 == 0 || fCoCsg0 == 0 || fCoCsg1 == 0)
    {
      std::cerr << "Error: TVpMaterial::InitializeCsgIub: "
		<< "Allocation of Csg and Iub arrays failed.\n";
      return 1;
    }

  // Calculate probability intervals
  Double_t energy1, energy2, toCsf1, toCsf2;
  Double_t inCsf1, inCsf2, coCsf1, coCsf2, phCsf1, phCsf2;  // Csf values
  Double_t inIub1, inIub2, coIub1, coIub2, phIub1, phIub2;  // Iub values
  for (Int_t i = 0; i < fDimIub; i++)
    {
      energy1 = GetEnergy(i);
      energy2 = GetEnergy(i+1);
      fEnergyIUB[i] = energy1;

      // Get file data values in 2 points
      inCsf1 = GetInCsf(energy1);  inCsf2 = GetInCsf(energy2);
      coCsf1 = GetCoCsf(energy1);  coCsf2 = GetCoCsf(energy2);
      phCsf1 = GetPhCsf(energy1);  phCsf2 = GetPhCsf(energy2);
      toCsf1 = inCsf1 + coCsf1 + phCsf1;  toCsf2 = inCsf2 + coCsf2 + phCsf2;
      // Calculate 1-st degree polynoms of Csg
      TVpMath::SetPoly1(fInCsg0[i], fInCsg1[i],	energy1, inCsf1, energy2, inCsf2);
      TVpMath::SetPoly1(fCoCsg0[i], fCoCsg1[i], energy1, coCsf1, energy2, coCsf2);
      TVpMath::SetPoly1(fPhCsg0[i], fPhCsg1[i], energy1, phCsf1, energy2, phCsf2);
      // Calculate 1-st degree polynoms of Iub.
      inIub1 = inCsf1 / toCsf1;
      inIub2 = inCsf2 / toCsf2;  
      phIub1 = (inCsf1 + phCsf1) / toCsf1;
      phIub2 = (inCsf2 + phCsf2) / toCsf2;
      coIub1 = (inCsf1 + phCsf1 + coCsf1) / toCsf1;
      coIub2 = (inCsf2 + phCsf2 + coCsf2) / toCsf2;
      TVpMath::SetPoly1(fInIub0[i], fInIub1[i], energy1, inIub1, energy2, inIub2);
      TVpMath::SetPoly1(fPhIub0[i], fPhIub1[i], energy1, phIub1, energy2, phIub2);
      TVpMath::SetPoly1(fCoIub0[i], fCoIub1[i], energy1, coIub1, energy2, coIub2);
      
      // Convert ToCsg from [cm^2/g] to [cm^-1]
      toCsf1 *= fDensity;  toCsf2 *= fDensity;
      TVpMath::SetPoly1(fToCsg0[i], fToCsg1[i], energy1, toCsf1, energy2, toCsf2);
    }
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterial::InitializeSfg()
{
  // Set the fSfg0 and fSfg1 arrays used by Sfg.  The first bin uses lin-lin
  // interpolation, the rest uses log-log.  See TVpMaterialGridData.cc for a
  // detailed description.

  fLnSfgXMin = log(fInXf[1]);         // the first non-zero value
  fLnSfgXMax = log(fInXf[fDimSff-1]); // the last value
  fDLnSfgX = (fLnSfgXMax - fLnSfgXMin) / (fDimLnSfg - 1);

  // De-allocate old arrays to avoid memory leaks.  It's safe because these
  // pointers are initialized to 0 by the constructor.
  delete [] fLnSfg0; delete [] fLnSfg1;
  
  // Alocate arrays
  fLnSfg0 = new Double_t[fDimLnSfg]; fLnSfg1 = new Double_t[fDimLnSfg];
  if (fLnSfg0 == 0 || fLnSfg1 == 0)
    {
      std::cerr << "Error:TVpMaterial::InitializeSfg: "
		<< "Cannot allocate fLnSfg0 or fLnSfg1 array\n";
      return 1;
    }
  
  // Bin 0 uses lin-lin interpolation
  fLnSfg1[0] = (fInSff[1]-fInSff[0]) / (fInXf[1]-fInXf[0]);
  fLnSfg0[0] = fInSff[0] - fInXf[0] * fLnSfg1[0];

  Double_t lnX0, lnX1, lnY0, lnY1;
  // Bins 1, 2, ... use log-log interpolation
  for (Int_t i = 1; i < fDimLnSfg; i++)
    {
      lnX0 = fLnSfgXMin + (i-1) * fDLnSfgX;
      lnX1 = lnX0 + fDLnSfgX;
      lnY0 = log(GetSff(exp(lnX0)));
      lnY1 = log(GetSff(exp(lnX1)));
      fLnSfg1[i] = (lnY1 - lnY0) / fDLnSfgX;
      fLnSfg0[i] = lnY1 - lnX1 * fLnSfg1[i];
    }
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterial::InitializeFfg()
{
  // Set the fFfg0 and fFfg1 arrays used by Ffg.  The first bin uses lin-lin
  // interpolation, the rest uses log-log.  See TVpMaterialGridData.cc for a
  // detailed description.

  fLnFfgXMin = log(fCoXf[1]);         // the first non-zero value
  fLnFfgXMax = log(fCoXf[fDimFff-1]); // the last value
  fDLnFfgX = (fLnFfgXMax - fLnFfgXMin) / (fDimLnFfg - 1);

  // De-allocate old arrays to avoid memory leaks.  It's safe because these
  // pointers are initialized to 0 by the constructor.
  delete [] fLnFfg0; delete [] fLnFfg1;
  
  // Alocate arrays
  fLnFfg0 = new Double_t[fDimLnFfg];
  fLnFfg1 = new Double_t[fDimLnFfg];
  if (fLnFfg0 == 0 || fLnFfg1 == 0)
    {
      std::cerr << "Error:TVpMaterial::InitializeFfg: "
		<< "Cannot allocate fLnFfg0 or fLnFfg1 array\n";
      return 1;
    }
  
  // Bin 0 uses lin-lin interpolation
  fLnFfg1[0] = (fCoFff[1]-fCoFff[0]) / (fCoXf[1]-fCoXf[0]);
  fLnFfg0[0] = fCoFff[0] - fCoXf[0] * fLnFfg1[0];

  Double_t lnX0, lnX1, lnY0, lnY1;
  // Bins 1, 2, ... use log-log interpolation
  for (Int_t i = 1; i < fDimLnFfg; i++)
    {
      lnX0 = fLnFfgXMin + (i-1) * fDLnFfgX;
      lnX1 = lnX0 + fDLnFfgX;
      lnY0 = log(GetFff(exp(lnX0)));
      lnY1 = log(GetFff(exp(lnX1)));
      fLnFfg1[i] = (lnY1 - lnY0) / fDLnFfgX;
      fLnFfg0[i] = lnY1 - lnX1 * fLnFfg1[i];
    }
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterial::InitializeNog()
{
  // Initialize Ntg and Nog arrays.  Lin-lin interpolation is used.

  // De-allocate old arrays to avoid memory leaks.  It's safe because these
  // pointers are initialized to 0 by the constructor.
  delete [] fInNtg0; delete [] fInNtg1;
  delete [] fCoNtg0; delete [] fCoNtg1;  
  delete [] fInNog0; delete [] fInNog1;
  delete [] fCoNog0; delete [] fCoNog1;  

  fInNtg0 = new Double_t[fDimIub]; fInNtg1 = new Double_t[fDimIub];
  fCoNtg0 = new Double_t[fDimIub]; fCoNtg1 = new Double_t[fDimIub];
  fInNog0 = new Double_t[fDimIub]; fInNog1 = new Double_t[fDimIub];
  fCoNog0 = new Double_t[fDimIub]; fCoNog1 = new Double_t[fDimIub];

  // To avoid errors, store calculated values in temporary arrays and evaluate
  // polynomial coefficients from these arrays.
  Double_t *inNtg = new Double_t[fDimIub+1];  // point data - thus + 1
  Double_t *coNtg = new Double_t[fDimIub+1];
  Double_t *inNog = new Double_t[fDimIub+1];
  Double_t *coNog = new Double_t[fDimIub+1];

  if (inNtg == 0 || coNtg == 0 || inNog == 0 || coNog == 0 || 
      fInNtg0 == 0 || fInNtg1 == 0 || fCoNtg0 == 0 || fCoNtg1 == 0 ||
      fInNog0 == 0 || fInNog1 == 0 || fCoNog0 == 0 || fCoNog1 == 0)
    {
      std::cerr << "Error: TVpMaterial::InitializeNog: "
		<< "Allocation of Ntg and Nog arrays failed\n";
      return 1;
    }

  Double_t energy1, energy2;
  // First, set the normalization factors to 1.0, i.e. inNtg[i] = coNtg[i] = 1.0
  for (Int_t i = 0; i < fDimIub; i++)
    {
      fInNtg0[i] =  fCoNtg0[i] = 1.0;
      fInNtg1[i] =  fCoNtg1[i] = 0.0;
      fInNog0[i] =  fCoNog0[i] = 1.0;
      fInNog1[i] =  fCoNog1[i] = 0.0;
    }
  
  // Second, calculate new normalization factors so that the integral returns 1.0
  for (Int_t i = 0; i <= fDimIub; ++i)
    {
      Double_t energy = fEnergyMinIUB + i * fChannelToEnergy;
      inNtg[i] = 1.0 / GetInPdfIntegral(energy);
      coNtg[i] = 1.0 / GetCoPdfIntegral(energy);
      inNog[i] = 1.0 / GetInPdfOmegaIntegral(energy);
      coNog[i] = 1.0 / GetCoPdfOmegaIntegral(energy);
    }

  // Third, calculate polynomial coefficients
  for (Int_t i = 0; i < fDimIub; i++)
    {
      energy1 = GetEnergy(i);
      energy2 = GetEnergy(i+1);
      TVpMath::SetPoly1(fInNtg0[i], fInNtg1[i], energy1, inNtg[i], energy2, inNtg[i+1]);
      TVpMath::SetPoly1(fCoNtg0[i], fCoNtg1[i], energy1, coNtg[i], energy2, coNtg[i+1]);
      TVpMath::SetPoly1(fInNog0[i], fInNog1[i], energy1, inNog[i], energy2, inNog[i+1]);
      TVpMath::SetPoly1(fCoNog0[i], fCoNog1[i], energy1, coNog[i], energy2, coNog[i+1]);
    }
  delete [] inNtg;  delete [] coNtg;
  delete [] inNog;  delete [] coNog;
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterial::GetBin(Double_t x, Int_t edgeDim, Double_t *lowEdge)
{
  // Return the bin number. Start from 0, -1 underflow, -2 overflow
  // Sequential search, a faster algorithm will be used later.

  if (x < lowEdge[0])
    return -1;        // underflow

  for (Int_t i = 1; i < edgeDim; i++)
    {
      if (x < lowEdge[i])
	return i-1;
    }

  return -2;          // overflow
}

//______________________________________________________________________________
Double_t TVpMaterial::GetA(Double_t u)
{
  // Return the A function value.

  Int_t i = GetBin(u, fDimA, fValU);
  Double_t val = fValA[i] + fPolyA0[i] * (u - fValU[i]);
  return val;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInvA(Double_t au)
{
  // Return the invA function value.

  Int_t i = GetBin(au, fDimA, fValA);
  Double_t val = (au - fValA[i]) / fPolyA0[i] + fValU[i];
  return val;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoXSample(Double_t energy)
{
  // Return x.  See coherent scattering sampling theory.

  Double_t umax = energy / TVpConstant::hc;
  umax *= umax;
  Double_t u = GetInvA(getRand() * GetA(umax));
  return sqrt(u);
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoPdfTheta(Double_t energy, Double_t cosTheta)
{
  // Return the value of the "per-theta" coherent scattering PDF for a given
  // angle and energy.
  //
  // Input:
  // energy ... energy of the incident photon in keV
  // cosTheta ....... cosine of the scattering angle, cosTheta = cos(Theta) 
  //
  // Example:
  // cout << matWater->GetCoPdfTheta(30, 0.5) << endl;
  // 0.287408

  Double_t F;  // Form factor
  Double_t cosTheta2 = cosTheta * cosTheta;
  if (fUseFf == 0)
    F = 1.0;
  else    
    {
      Double_t x = GetMomentumTransferX(energy, cosTheta);
      F = GetFff(x);
    }
  Double_t sinTheta = sqrt(1.0 - cosTheta2);
  Double_t value = (1 + cosTheta2) * sinTheta * F * F;
  Double_t norm = GetCoNtg(energy);  // Normalization factor
  return norm * value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoPdfOmega(Double_t energy, Double_t cosTheta)
{
  // Return the value of the "per-omega" coherent scattering PDF for a given
  // angle and energy.
  //
  // Input:
  // energy ... energy of the incident photon in keV
  // cosTheta ....... cosine of the scattering angle, cosTheta = cos(Theta) 
  //
  // Example:
  // cout << matWater->GetCoPdfOmega(30, 0.5) << endl;
  // 0.0528188

  Double_t F;  // Form factor
  Double_t cosTheta2 = cosTheta * cosTheta;
  if (fUseFf == 0)
    F = 1.0;
  else    
    {
      Double_t x = GetMomentumTransferX(energy, cosTheta);
      F = GetFff(x);
    }
  Double_t value = (1 + cosTheta2) * F * F;
  Double_t norm = GetCoNog(energy);  // Normalization factor
  return norm * value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInPdfTheta(Double_t energy, Double_t cosTheta)
{
  // Return the value of the "per-theta" incoherent scattering PDF for a given
  // angle and energy.
  //
  // Input:
  // energy ..... energy of the incident photon in keV
  // cosTheta ... cosine of the scattering angle
  //
  // Example:
  // cout << matWater->GetInPdfTheta(30, 0.5) << endl;
  // 0.426402

  Double_t S;  // Scattering function
  if (fUseSf == 0)
    S = 1.0;
  else
    {
      Double_t x = GetMomentumTransferX(energy,  cosTheta);
      S = GetSfg(x);
    }
  Double_t value = GetKnCsTheta(energy, cosTheta) * S;
  Double_t norm = GetInNtg(energy);  // Normalization factor
  return norm * value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInPdfOmega(Double_t energy, Double_t cosTheta)
{
  // Return the value of the "per-omega" incoherent scattering PDF for a given
  // angle and energy.
  //
  // Input:
  // energy ..... energy of the incident photon in keV
  // cosTheta ... cosine of the scattering angle
  //
  // Example:
  // cout << matWater->GetInPdfOmega(30, 0.5) << endl;
  // 0.0783625

  Double_t S;  // Scattering function
  if (fUseSf == 0)
    S = 1.0;
  else
    {
      Double_t x = GetMomentumTransferX(energy,  cosTheta);
      S = GetSfg(x);
    }
  Double_t value = GetKnCsOmega(energy, cosTheta) * S;
  Double_t norm = GetInNog(energy);  // Normalization factor
  return norm * value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInPdfIntegral(Double_t energy)
{
  // Return the integral of GetInPdfTheta(energy, cosTheta).  This function is
  // used to calculate the normalization arrays fInNtg0[] and fInNtg1[].  If
  // already normalized, this function returns 1.0.
  //
  // Input:
  // - energy - incident photon energy in keV
  //
  // Example:
  // cout << matWater->GetInPdfIntegral(30) << endl;
  // 1

  return Integral(0, energy, 0, M_PI);
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoPdfIntegral(Double_t energy)
{
  // Return the integral of GetCoPdfTheta(energy, cosTheta).  This function is
  // used to calculate the normalization arrays fCoNtg0[] and fCoNtg1[].  If
  // already normalized, this function returns 1.0.
  //
  // Input:
  // - energy - incident photon energy in keV
  //
  // Example:
  // cout << matWater->GetCoPdfIntegral(30) << endl;
  // 0.999986

  return Integral(1, energy, 0, M_PI);
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInPdfOmegaIntegral(Double_t energy)
{
  // Return the integral of GetInPdfOmega(energy, cosTheta).  This function is
  // used to calculate the normalization arrays fInNog0[] and fInNog1[].  If
  // already normalized, this function returns 1.0.
  //
  // Input:
  // - energy - incident photon energy in keV
  //
  // Example:
  // cout << matWater->GetInPdfOmegaIntegral(30) << endl;
  // 1

  return Integral(2, energy, 0, M_PI);
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoPdfOmegaIntegral(Double_t energy)
{
  // Return the integral of GetCoPdfOmega(energy, cosTheta).  This function is
  // used to calculate the normalization arrays fCoNog0[] and fCoNog1[].  If
  // already normalized, this function returns 1.0.
  //
  // Input:
  // - energy - incident photon energy in keV
  //
  // Example:
  // cout << matWater->GetCoPdfOmegaIntegral(30) << endl;
  // 0.999986

  return Integral(3, energy, 0, M_PI);
}

//______________________________________________________________________________
Double_t TVpMaterial::GetKnCsTheta(Double_t energy, Double_t cosTheta)
{
  // Return the "per-theta" Klein-Nishina differential cross section in cm^2.
  // It corresponds to the scattering of a photon on a free electron.
  //
  // Latex formulas:
  // ds / d\theta = (r_e^2/2) (E'/E)^2 (E/E' + E'/E - \sin^2 \theta)
  // E' = E / (1 + E / (m_e c^2) (1 - \cos \theta)
  //
  // Parameters:
  // energy ...... energy of incident photon [keV]
  // cosTheta .... cosine of the scattering angle

  Double_t cosTheta2 = cosTheta * cosTheta;
  Double_t sinTheta = sqrt(1.0 - cosTheta2);
  Double_t rat = 1.0 / (1.0 + energy / TVpConstant::fElectronRestEnergy
			* (1.0 - cosTheta));  // = E' / E
  Double_t knd = TVpConstant::fClassicalElectronRadius22
    * rat * rat * (rat + 1.0/rat + cosTheta2 - 1.0) * 2 * M_PI * sinTheta;
  return knd;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetKnCsOmega(Double_t energy, Double_t cosTheta)
{
  // Return the "per-omega" Klein-Nishina differential cross section in cm^2.
  // It corresponds to the scattering of a photon on a free electron.
  //
  // Latex formulas:
  // ds / d\theta = (r_e^2/2) (E'/E)^2 (E/E' + E'/E - \sin^2 \theta)
  // E' = E / (1 + E / (m_e c^2) (1 - \cos \theta)
  //
  // Parameters:
  // energy ...... energy of incident photon [keV]
  // cosTheta .... cosine of the scattering angle

  Double_t cosTheta2 = cosTheta * cosTheta;
  Double_t rat = 1.0 / (1.0 + energy / TVpConstant::fElectronRestEnergy
			* (1.0 - cosTheta));  // = E' / E
  Double_t knd = TVpConstant::fClassicalElectronRadius22
    * rat * rat * (rat + 1.0/rat + cosTheta2 - 1.0);
  return knd;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetKnCsI(Double_t energy)
{
  // Return the Klein-Nishina cross section.
  // Status: CHECK ME

  Double_t pre2 = 2* M_PI * 7.940775e-26;
  Double_t k = energy / TVpConstant::fElectronRestEnergy;

  Double_t h1 = 1.0 + 2 * k;
  Double_t kni = pre2 * ( (1.0 + k) / (k*k) * (2*(1.0+k)/h1 - log(h1)/k) +
			log(h1)/(2*k) - (1+3*k)/(h1*h1));
  return kni;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetInCsTheta(Double_t energy, Double_t cosTheta)
{
  // Return either the "per-theta" differential cross section for incoherent
  // scattering of a photon on an atom in cm^2/radian if the incoherent
  // scattering function is used or the per-theta Klein-Nishina differential
  // cross section in cm^2/radian otherwise.  It is used to (1) calculate the
  // total coherent scattering cross section, (2) by GetInPdfTheta() which is
  // used in MC simulation.
  //
  // Parameters:
  // energy ..... energy of the incident photon in [keV]
  // cosTheta ... cosine of the scattering angle
  
  // value is calculated according to GetKnCsTheta()
  Double_t cosTheta2 = cosTheta * cosTheta;
  Double_t sinTheta = sqrt(1.0 - cosTheta2);
  Double_t rat = 1.0 / (1.0 + energy / TVpConstant::fElectronRestEnergy
			* (1.0 - cosTheta));  // = E' / E
  Double_t value = TVpConstant::fClassicalElectronRadius22
    * rat * rat * (rat + 1.0/rat + cosTheta2 - 1.0) * 2 * M_PI * sinTheta;
  
  // Form factor correction
  if (fUseSf != 0)
    {
      Double_t x = GetMomentumTransferX(energy, cosTheta);
      Double_t isf = GetSfg(x);
      value *= isf;
    }
  return value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCoCsTheta(Double_t energy, Double_t cosTheta)
{
  // Return either the "per-theta" differential cross section for coherent
  // scattering of a photon on an atom in cm^2/radian if the form factor is
  // used or the per-theta Thompson differential cross section in cm^2/radian
  // otherwise.  It is used to (1) calculate the total coherent scattering
  // cross section, (2) by GetCoPdfTheta() which is used in MC simulation.
  //
  // Parameters:
  // energy ..... energy of the incident photon in [keV]
  // cosTheta ... cosine of the scattering angle

  Double_t cosTheta2 = cosTheta * cosTheta;
  Double_t sinTheta = sqrt(1.0 - cosTheta2);
  // 
  Double_t value = TVpConstant::fClassicalElectronRadius22
    * (1.0 + cosTheta2) * 2 * M_PI * sinTheta;
  
  // Form factor correction
  if (fUseFf != 0)
    {
      Double_t sinTheta_2 = sqrt((1.0 - cosTheta)/2.0);
      Double_t x = energy * sinTheta_2 / TVpConstant::hc;
      Double_t cff = GetFff(x);
      value *= cff * cff;
    }

  return value;
}

//______________________________________________________________________________
Double_t TVpMaterial::GetCSSum(Double_t energy, Int_t csSum)
{
  // Return the sum of cross sections
  
  Double_t sum = 0.0;

  if (csSum & kInCS)
    sum = GetInCsf(energy);
  if (csSum & kPhCS)
    sum += GetPhCsf(energy);
  if (csSum & kCoCS)
    sum += GetCoCsf(energy);

  return sum;
}


//______________________________________________________________________________
void TVpMaterial::SetDensity(Double_t density)
{
  // Set the material density [g/cm^3]
  
  if (fFileNameMAT == 0)
    {
      fprintf(stderr, 
	      "Warning: TVpMaterial::SetDensity: No material data, Setting the density has no efect.\n");
      return;
    }
  Double_t denRat = density / fDensity;
  fDensity = density;
  for (Int_t i = 0; i < fDimIub; i++)
    {
      fToCsg0[i] *= denRat;
      fToCsg1[i] *= denRat;
    }
}

//______________________________________________________________________________
void TVpMaterial::InitializeATable()
{
  // Initialize the A function arrays

  // The function A is very steep around 0. Both A and its inversion
  // function invA are used which makes the problem very numerically
  // sensitive. Therefore the original grid points of FF are preserved
  // and linear interpolation is used. As a result, the identity
  // InvA(A(x)) == x is preserved but the algorithm is slow.

  fDimA = fDimIub;                  // Preserve the original grid
  fValU = new Double_t[fDimA];
  fValA = new Double_t[fDimA];     // These values are used for checking only

  fPolyA0 = new Double_t[fDimA-1]; // fDimA points and fDimA-1 segments

  // Initialize fPolyA0 array. It contains average values.
  for (Int_t i = 0; i < fDimA-1; i++)
    {
      Double_t F1 = GetFfgValue(i);
      Double_t F2 = GetFfgValue(i+1);
      fPolyA0[i] = 0.5 * (F1*F1 + F2*F2);
    }

  // Initialize fValU and fValA
  fValU[0] = 0.0;
  fValA[0] = 0.0;
  for (Int_t i = 1; i < fDimA; i++)
    {
      Double_t x = GetFfgX(i);
      fValU[i] = x * x;
      fValA[i] = fValA[i-1] + fPolyA0[i-1] * (fValU[i]-fValU[i-1]);
    }
}

//______________________________________________________________________________
void TVpMaterial::PrintATable()
{
  // Print A table.  Use for testing purposes.  The default size of A table is
  // 1024 rows.
  //
  // Example:
  // matWater->PrintATable();
  // 0.000000e+00    0.000000e+00
  // 1.000000e+12    3.446864e+11
  // 1.014971e+12    3.498467e+11
  // ...
  // 3.882869e+18    6.045622e+15
  // 3.940999e+18    nan            // FIX the NaN

  for (Int_t i = 0; i < fDimA; i++)
    printf("%e\t%e\n", fValU[i], fValA[i]);
}

//______________________________________________________________________________
void TVpMaterial::PrintStatusGridData(std::ostream& out)
{
  // Print status.  The listing is very long; complete Csg and Iub arrays are
  // printed.
  //
  // Example:
  // ofstream out("tables.dat");
  // matWater->PrintStatusGridData(out);
  // out.close();

  out.setf(std::ios_base::scientific);
  out << "<TVpMaterialGridData::PrintStatus>\n"
      << "$Id: TVpMaterial.cc 62 2009-06-27 10:54:08Z malusek $\n";
  out << std::setprecision(4);
  
  // Iub arrays 
  out << "<Iub arrays>\n"
      << "fDimIub: " << fDimIub << '\n'
      << "Energy min: " << fEnergyMinIUB << '\n'
      << "Energy max: " << fEnergyMaxIUB << '\n'
      << "    i  fEnergyIUB     fInIub0     fInIub1     fPhIub0     fPhIub1     fCoIub0     fCoIub1\n";
  for (Int_t i = 0; i < fDimIub; i++)
    out  << std::setw(5) << i << ' '
	 << std::setw(11) << fEnergyIUB[i]<< ' '
	 << std::setw(11) << fInIub0[i] << ' '
	 << std::setw(11) << fInIub1[i] << ' '
	 << std::setw(11) << fPhIub0[i] << ' '
	 << std::setw(11) << fPhIub1[i] << ' '
	 << std::setw(11) << fCoIub0[i] << ' '
	 << std::setw(11) << fCoIub1[i] << '\n';
  out << "</Iub arrays>\n";
  
  // Csg arrays
  out << "<Csg arrays>\n"
      << "fDimIub: " << fDimIub << '\n'
      << "Energy min: " << fEnergyMinIUB << '\n'
      << "Energy max: " << fEnergyMaxIUB << '\n'
      << "    i  fEnergyIUB     fInCsg0     fInCsg1     fPhCsg0     fPhCsg1     fCoCsg0     fCoCsg1     fToCsg0     fToCsg1\n";
  for (Int_t i = 0; i < fDimIub; i++)
    out  << std::setw(5) << i << ' '
	 << std::setw(11) << fEnergyIUB[i]<< ' '
	 << std::setw(11) << fInCsg0[i] << ' '
	 << std::setw(11) << fInCsg1[i] << ' '
	 << std::setw(11) << fPhCsg0[i] << ' '
	 << std::setw(11) << fPhCsg1[i] << ' '
	 << std::setw(11) << fCoCsg0[i] << ' '
	 << std::setw(11) << fCoCsg1[i] << ' '
	 << std::setw(11) << fToCsg0[i] << ' '
	 << std::setw(11) << fToCsg1[i] << '\n';
  out << "</Csg arrays>\n";

  out << std::setprecision(6);  // Back to defaults

  // Iub values
  out << "<Iub values>\n"
      << "fDimIub: " << fDimIub << '\n'
      << "Energy min: " << fEnergyMinIUB << '\n'
      << "Energy max: " << fEnergyMaxIUB << '\n'
      << "    i    fEnergyIUB         InIub         PhIub         CoIub\n";
  for (Int_t i = 0; i < fDimIub; i++)
    out  << std::setw(5) << i << ' '
	 << std::setw(13) << fEnergyIUB[i]<< ' '
	 << std::setw(13) << GetInIub(fEnergyIUB[i]) << ' '
	 << std::setw(13) << GetPhIub(fEnergyIUB[i]) << ' '
	 << std::setw(13) << GetCoIub(fEnergyIUB[i]) << '\n';
  out << "</Iub values>\n";

  // Csg values
  out << "<Csg values>\n"
      << "fDimIub: " << fDimIub << '\n'
      << "Energy min: " << fEnergyMinIUB << '\n'
      << "Energy max: " << fEnergyMaxIUB << '\n'
      << "    i    fEnergyIUB         InCsg         CoCsg         PhCsg         ToCsg\n"
      << "              [keV]      [cm^2/g]      [cm^2/g]      [cm^2/g]        [1/cm]\n";
  for (Int_t i = 0; i < fDimIub; i++)
    out  << std::setw(5) << i << ' '
	 << std::setw(13) << fEnergyIUB[i]<< ' '
	 << std::setw(13) << GetInCsg(fEnergyIUB[i]) << ' '
	 << std::setw(13) << GetPhCsg(fEnergyIUB[i]) << ' '
	 << std::setw(13) << GetCoCsg(fEnergyIUB[i]) << ' '
	 << std::setw(13) << GetToCsg(fEnergyIUB[i]) << '\n';
  out << "</Csg values>\n";

  // Sfg arrays
  out << "<Sfg arrays>\n"
      << "fUseSf: " << fUseSf << '\n'
      << "fDimLnSfg: " << fDimLnSfg << '\n'
      << "fLnSfgXMin: " << fLnSfgXMin << '\n'
      << "fLnSfgXMax: " << fLnSfgXMax << '\n';
  if (fUseSf != 0)
    {
      for (Int_t i = 0; i < fDimLnSfg; i++)
	out  << std::setw(5) << i << ' '
	     << std::setw(13) << fLnSfg0[i] << ' '
	     << std::setw(13) << fLnSfg1[i] << '\n';
    }
  out << "</Sfg arrays>\n";
  
  // Sfg values
  out << "<Sfg values>\n"
      << "fUseSf: " << fUseSf << '\n'
      << "fDimLnSfg: " << fDimLnSfg << '\n'
      << "InXGMin: " << exp(fLnSfgXMin) << '\n'
      << "InXGMax: " << exp(fLnSfgXMax) << '\n';
  
  if (fUseSf != 0)
    {
      if (fDimLnSfg > 0)
	out  << std::setw(5) << 0 << ' '
	     << std::setw(13) << 0 << ' '
	 << std::setw(13) << GetSfg(0) << '\n';
      for (Int_t i = 1; i <= fDimLnSfg; i++)
	{
	  Double_t x = exp(fLnSfgXMin + (i-1)*fDLnSfgX);
	  out  << std::setw(5) << i << ' '
	       << std::setw(13) << x << ' '
	       << std::setw(13) << GetSfg(x) << '\n';
	}
    }
  out << "</Sfg values>\n";
  
  // Ffg arrays
  out << "<Ffg arrays>\n"
      << "fUseFf: " << fUseFf << '\n'
      << "fDimLnFfg: " << fDimLnFfg << '\n'
      << "fLnFfgXMin: " << fLnFfgXMin << '\n'
      << "fLnFfgXMax: " << fLnFfgXMax << '\n';
  if (fUseFf != 0)
    {
      for (Int_t i = 0; i < fDimLnFfg; i++)
	out  << std::setw(5) << i << ' '
	     << std::setw(13) << fLnFfg0[i] << ' '
	     << std::setw(13) << fLnFfg1[i] << '\n';
    }
  out << "</Ffg arrays>\n";
  
  // Ffg values
  out << "<Ffg values>\n"
      << "fUseFf: " << fUseFf << '\n'
      << "fDimLnFfg: " << fDimLnFfg << '\n'
      << "InXGMin: " << exp(fLnFfgXMin) << '\n'
      << "InXGMax: " << exp(fLnFfgXMax) << '\n';
  if (fUseFf != 0)
    {
      if (fDimLnFfg > 0)
	out  << std::setw(5) << 0 << ' '
	     << std::setw(13) << 0 << ' '
	     << std::setw(13) << GetFfg(0) << '\n';
      for (Int_t i = 1; i <= fDimLnFfg; i++)
	{
	  Double_t x = exp(fLnFfgXMin + (i-1)*fDLnFfgX);
	  out  << std::setw(5) << i << ' '
	       << std::setw(13) << x << ' '
	       << std::setw(13) << GetFfg(x) << '\n';
	}
    }
  out << "</Ffg values>\n";
  
  out << "</TVpMaterialGridData::PrintStatus>" << std::endl;
}

//______________________________________________________________________________
Double_t TVpMaterial::EvaluateIntegrand(Int_t selector, Double_t par, Double_t x)
{
  // Provide value to the TVpIntegral.
  //
  // selector = 0 ... InPdfTheta
  // selector = 1 ... CoPdfTheta
  // selector = 2 ... InPdfOmega
  // selector = 3 ... CoPdfOmega

  switch (selector)
    {
    case 0:
      return GetInPdfTheta(par, cos(x));
    case 1:
      return GetCoPdfTheta(par, cos(x));
    case 2:
      return GetInPdfOmega(par, cos(x)) * 2 * M_PI * sin(x);
    case 3:
      return GetCoPdfOmega(par, cos(x)) * 2 * M_PI * sin(x);

    default:
      std::cerr << "Error: TVpMaterial::EvaluateIntegrand: "
		<< "Invalid selector\n";
    }
  return 0.0;
}


//______________________________________________________________________________
TH2F *TVpMaterial::GetHistCoCSG(Int_t numEnergy, Int_t numXi)
{
  // Make a new histogram

  Double_t dXi = 2.0 / numXi;
  Double_t denergy = (fEnergyMaxCsf - fEnergyMinCsf) / numEnergy;

  TH2F *h2 = new TH2F("CoCSG", "CoCSG", numEnergy, fEnergyMinCsf, fEnergyMaxCsf, numXi, -1, 1);

  for (Int_t i = 0; i <  numEnergy; i++)  
    for (Int_t j = 0; j < numXi; j++)
      {
	Double_t energy = fEnergyMinCsf + i * denergy;
	Double_t xi = -1.0 + j * dXi;
	h2->SetCellContent(i+1, j+1, GetCoCsTheta(energy, xi));
      }
  h2->SetXTitle("Energy [keV]");
  h2->SetYTitle("#xi [1]");
  h2->SetZTitle("#sigma [cm^{2} g^{-1}]");
  h2->GetXaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetTitleOffset(1.5);
  return h2;
}

//______________________________________________________________________________
TH2F *TVpMaterial::GetHistInCSG(Int_t numEnergy, Int_t numTheta)
{
  // Make a new histogram

  Double_t dtheta = M_PI / numTheta;
  Double_t denergy = (fEnergyMaxCsf - fEnergyMinCsf) / numEnergy;

  TH2F *h2 = new TH2F("InCSG", "InCSG", numEnergy, fEnergyMinCsf, fEnergyMaxCsf, numTheta, 0, M_PI);

  for (Int_t i = 0; i <  numEnergy; i++)  
    for (Int_t j = 0; j < numTheta; j++)
      {
	Double_t energy = fEnergyMinCsf + (i + 0.5) * denergy;
	Double_t theta = (j + 0.5) * dtheta;
	h2->SetCellContent(i+1, j+1, GetInCsTheta(energy, cos(theta)));
      }
  h2->SetXTitle("Energy [keV]");
  h2->SetYTitle("Theta [rad]");
  h2->SetZTitle("#sigma [cm^{2} g^{-1}]");
  return h2;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistKnCsI(Int_t numEnergy = 128)
{
  // Klein-Nishina
  
  Double_t denergy = (fEnergyMaxCsf - fEnergyMinCsf) / numEnergy;
  TH1F *h1 = new TH1F("KNI", "KNI", numEnergy, fEnergyMinCsf, fEnergyMaxCsf);
  for (Int_t i = 0; i <  numEnergy; i++)
    {
      Double_t energy = fEnergyMinCsf + denergy * i;
      h1->SetBinContent(i+1, GetKnCsI(energy));
    }
  h1->SetXTitle("Energy [keV]");
  h1->SetYTitle("#sigma");
  return h1;
}

//______________________________________________________________________________
TH2F *TVpMaterial::GetHistKnCsTheta2D(Int_t numEnergy, Int_t numTheta)
{
  // Make a new histogram

  Double_t dtheta = M_PI / numTheta;
  Double_t denergy = (fEnergyMaxCsf - fEnergyMinCsf) / numEnergy;

  TH2F *h2 = new TH2F("KND", "KND", numEnergy, fEnergyMinCsf, fEnergyMaxCsf, numTheta, 0, M_PI);

  for (Int_t i = 0; i <  numEnergy; i++)  
    for (Int_t j = 0; j < numTheta; j++)
      {
	Double_t energy = fEnergyMinCsf + i * denergy;
	Double_t theta = j * dtheta;
	h2->SetCellContent(i+1, j+1, GetKnCsTheta(energy, theta));
      }
  h2->SetXTitle("Energy [keV]");
  h2->SetYTitle("Theta [rad]");
  h2->SetZTitle("#sigma [cm^2]");
  return h2;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistA()
{
  // Make a new histogram

  TH1F *h1 = new TH1F("A", "A", fDimA-1, fValU);
  for (Int_t i = 0; i <  fDimA; i++)
    h1->SetBinContent(i+1, fValA[i]);
  h1->SetXTitle("u = x^2");
  h1->SetYTitle("A");
  return h1;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistCoCSG(Double_t energy, Int_t numXi)
{
  // Make a new histogram

  Double_t dXi = 2.0 / numXi;

  TH1F *h = new TH1F("CoCsG", "CoCsG", numXi, -1, 1);

  for (Int_t i = 0; i < numXi; i++)
    {
      Double_t xi = -1.0 + i * dXi;
      h->SetBinContent(i+1, GetCoCsTheta(energy, xi));
    }
  h->SetXTitle("#xi [1]");
  h->SetYTitle("#sigma [cm^{2} g^{-1}]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistCoPdfXi(Double_t energy, Int_t numXi)
{
  // Make a new histogram

  Double_t dXi = 2.0 / numXi;

  TH1F *h = new TH1F("CoCsXi", "CoCsXi", numXi, -1, 1);

  for (Int_t i = 0; i < numXi; i++)
    {
      Double_t xi = -1.0 + i * dXi;
      h->SetBinContent(i+1, GetCoPdfTheta(energy, xi));
    }
  h->SetXTitle("#xi [1]");
  h->SetYTitle("pdf");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistCoPdfTheta(Double_t energy, Int_t numTheta)
{
  // Return PDF of the "per-theta" angular distribution of coherently
  // scattered photons with a given energy as a function of the scattering
  // angle theta.  The PDF is based on tabulated data and evaluated in the
  // center of each bin.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("CoPdfTheta", "CoPdfTheta", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetCoPdfTheta(energy, cos(theta)));
    }
  h->SetXTitle("#theta [1]");
  h->SetYTitle("(1/#sigma) d #sigma / d #theta  [1/rad]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistInPdfTheta(Double_t energy, Int_t numTheta)
{
  // Return PDF of the "per-theta" angular distribution of incoherently
  // scattered photons with a given energy as a function of the scattering
  // angle theta.  The PDF is based on tabulated data and evaluated in the
  // center of each bin.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("InPdfTheta", "InPdfTheta", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetInPdfTheta(energy, cos(theta)));
    }
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("(1/#sigma) d #sigma / d #theta  [1/rad]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistCoPdfOmega(Double_t energy, Int_t numTheta)
{
  // Return PDF of the "per-omega" angular distribution of coherently
  // scattered photons with a given energy as a function of the scattering
  // angle theta.  The PDF is based on tabulated data and evaluated in the
  // center of each bin.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("CoPdfOmega", "CoPdfOmega", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetCoPdfOmega(energy, cos(theta)));
    }
  h->SetXTitle("#theta [1]");
  h->SetYTitle("(1/#sigma) d #sigma / d #Omega  [1/sr]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistInPdfOmega(Double_t energy, Int_t numTheta)
{
  // Return PDF of the "per-omega" angular distribution of incoherently
  // scattered photons with a given energy as a function of the scattering
  // angle theta.  The PDF is based on tabulated data and evaluated in the
  // center of each bin.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("InPdfOmega", "InPdfOmega", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetInPdfOmega(energy, cos(theta)));
    }
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("(1/#sigma) d #sigma / d #Omega  [1/sr]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistKnCsTheta(Double_t energy, Int_t numTheta)
{
  // Return "per-theta" Klein-Nishina differential cross section for a
  // scattering of a photon on an electron in barn/radian as a function of the
  // scattering angle theta.  The cross section is based on an analytical
  // formula and evaluated in the center of each bin.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("KnCsTheta", "KnCsTheta", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetKnCsTheta(energy, cos(theta)));
    }
  h->Scale(1.0 / TVpConstant::fBarn);  // convert from cm^2 to barn
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("d #sigma / d #theta  [barn/rad]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistInCsTheta(Double_t energy, Int_t numTheta)
{
  // Return the histogram corresonding to the function GetInCsTheta(), i.e.
  // either the "per-theta" differential cross section for incoherent
  // scattering of a photon on an atom in barn/radian if the incoherent
  // scattering function is used or the per-theta Klein-Nishina differential
  // cross section in barn/radian otherwise.  NOTE: If relative values of the
  // scattering function are used then absolute values of the differential
  // cross section are not correct.  Values are not renormalized to match the
  // cross section.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("InCsTheta", "InCsTheta", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetInCsTheta(energy, cos(theta)));
    }
  h->Scale(1.0 / TVpConstant::fBarn);  // convert from cm^2 to barn
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("d #sigma / d #theta  [barn/rad]");
  std::cout << "Warning: ABSOLUTE values of the differential cross section are correct "
	    << "only if absolute values of the scattering function are correct.  For "
	    << "compounds and mixtures only RELATIVE values are often available."
	    << std::endl;
  return h;
}

//______________________________________________________________________________
TH1F *TVpMaterial::GetHistCoCsTheta(Double_t energy, Int_t numTheta)
{
  // Return the histogram corresonding to the function GetCoCsTheta(), i.e.
  // either the "per-theta" differential cross section for incoherent
  // scattering of a photon on an atom or molecule in barn/radian if the form
  // factor function is used or the per-theta Thompson differential cross
  // section in barn/radian otherwise.  NOTE: If relative values of the form
  // factor are used then absolute values of the differential cross section
  // are not correct.  Values are not renormalized to match the cross section.
  // 
  // Parameters:
  // energy ..... energy of the photon in keV
  // numTheta ... number of bins of the TH1F histogram 

  Double_t dTheta = M_PI / numTheta;

  TH1F *h = new TH1F("CoCsTheta", "CoCsTheta", numTheta, 0, M_PI);

  for (Int_t i = 0; i < numTheta; i++)
    {
      Double_t theta = (i + 0.5) * dTheta;
      h->SetBinContent(i+1, GetCoCsTheta(energy, cos(theta)));
    }
  h->Scale(1.0 / TVpConstant::fBarn);  // convert from cm^2 to barn
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("d #sigma / d #theta  [barn/rad]");
  std::cout << "Warning: ABSOLUTE values of the differential cross section are correct "
	    << "only if absolute values of the form factor are correct.  For compounds "
	    << "and mixtures only RELATIVE values are often available."  << std::endl;
  return h;
}
