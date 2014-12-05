//______________________________________________________________________________
//
// TVpMaterialFileData stores material data read from data files.  These data
// are used by the class TVpMaterial to create grid data tables which are used
// in MC simulations.
//
// Naming convention:
// Csf ... cross sections based on file data
// Sff ... scattering function (incoherent scattering), based on file data
// Fff ... form factor (coherent scattering), based on file data

#include <istream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "misc.h"
#include "TVpMath.h"
#include "TVpMaterialFileData.h"

ClassImp(TVpMaterialFileData)

//______________________________________________________________________________
TVpMaterialFileData::TVpMaterialFileData()
{
  // Default constructor. Initialize all data members to 0

  fEnergyCsf = fInCsf = fPhCsf = fCoCsf = fCoXf = fCoFff = fInXf = fInSff = 0;
  fName = fFileNameMAT = fFileNameCFF = fFileNameISF;
  fUseFf = fUseSf = 0;
  fDensity = 0.0;
  fDimCsf = fDimFff = fDimSff = 0;
}

//______________________________________________________________________________
TVpMaterialFileData::~TVpMaterialFileData()
{
  // Destructor. Delete all allocated arrays.

  delete [] fEnergyCsf;
  delete [] fInCsf;
  delete [] fPhCsf;
  delete [] fCoCsf;
  delete [] fCoXf;
  delete [] fCoFff;
  delete [] fInXf;
  delete [] fInSff;
  delete [] fName;
  delete [] fFileNameMAT;
  delete [] fFileNameCFF;
  delete [] fFileNameISF;
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::GetIndexCsf(Double_t energy) const
{
  // Return index. If energy is out of range, return -1 (under minimum) or -2
  // (over maximum) and write a warning to cerr.
  //
  // -1        0         1                 fDimCsf-2           -2
  // -----|---------|-------............|----------------|-------------
  // fEnergyCsf[0] fEnergyCsf[1] fEnergyCsf[fDimCsf-2] fEnergyCsf[fDimCsf-1]

  Int_t index = TVpMath::FindIndexByBinarySearch(fDimCsf, fEnergyCsf, energy);
  if (index < 0)
    std::cerr << "Warning: TVpMaterialFileData::GetIndexCsf: energy = "
	      << energy << " out of range "
	      << "(" << fEnergyCsf[0] << ", " << fEnergyCsf[fDimCsf-1] << ")\n";
  return index;
}

//______________________________________________________________________________
Double_t TVpMaterialFileData::GetInCsf(Double_t energy) const
{
  // Return incoherent scattering cross section
  // Source: log-log interpolated material table

  Int_t ind;

  if ((ind = GetIndexCsf(energy)) == -1)
    return 0.0;       // out of range

  return TVpMath::GetLogLogInterpolation(fEnergyCsf[ind], fInCsf[ind],
					 fEnergyCsf[ind+1], fInCsf[ind+1], energy);
}

//______________________________________________________________________________
Double_t TVpMaterialFileData::GetPhCsf(Double_t energy) const
{
  // Return photoeffect cross section
  // Source: log-log interpolated material table

  Int_t ind;

  if ((ind = GetIndexCsf(energy)) == -1)
    return 0.0;       // out of range

  return TVpMath::GetLogLogInterpolation(fEnergyCsf[ind], fPhCsf[ind],
					 fEnergyCsf[ind+1], fPhCsf[ind+1], energy);
}

//______________________________________________________________________________
Double_t TVpMaterialFileData::GetCoCsf(Double_t energy) const
{
  // Return the total coherent CS by log-log interpolating MAT file values.

  Int_t ind;

  if ((ind = GetIndexCsf(energy)) == -1)
    return 0.0;       // out of range

  return TVpMath::GetLogLogInterpolation(fEnergyCsf[ind], fCoCsf[ind],
					 fEnergyCsf[ind+1], fCoCsf[ind+1], energy);
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::ReadMatFile(const Char_t *fileName)
{
  // Read cross section data. The format is apparent from the example.
  // Energy interval must be 1 keV in current implementation.
  // Material ... Descriptive name up to 80 characters, the rest is ignored
  // Density .... Density of the material
  // Energy ..... Photon energy. Currently, this field is not read.
  // Co_CS ...... Coherent scattering cross section
  // In_CS ...... Incoherent scattering cross section
  // Ph_CS ...... Photoelectric scattering cross section
  //
  // === Beginning of the example =================================
  // # Material: Brain (whole), adult
  // # Density: 1.04 g/cm^3
  // # Dimension: 213
  // # Energy      Co_CS     In_CS     Ph_CS
  // # (MeV)       (cm2/g)   (cm2/g)   (cm2/g)
  //  1.00000e-03  1.33e+00  1.37e-02  3.70e+03
  //  1.07200e-03  1.31e+00  1.55e-02  3.07e+03
  // ...
  //  2.00000e-01  1.35e-03  1.35e-01  3.17e-04
  // == End of the example ========================================
  //
  // RETURN VALUE:
  // 0 ... success
  // 1 ... openning of the file failed

  // Debug info
  fprintf(stderr, "<Info TVpMaterialFileData::ReadMatFile(%s)>\n", fileName);

  FILE *fp;
  Char_t matName[80];
  Int_t strLen;
  const Double_t  cMeV2keV = 1.0e+3;             // conversion coefficient
  
  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpMaterial::ReadMATFile: Cannot open the file %s: ",
	      fileName);
      perror("");
      return 1;
    }
   
  // Read the file
  fscanf(fp, "# Material: %[^\n]\n", matName);       // Read line 1
  fscanf(fp, "# Density: %lf%*[^\n]\n", &fDensity);  // Read line 2
  fscanf(fp, "# Dimension: %d\n", &fDimCsf);          // Read line 3
  fscanf(fp, "%*[^\n]\n");                           // Skip line 4
  fscanf(fp, "%*[^\n]\n");                           // Skip line 5

  // Copy name
  strLen = strlen(matName) + 1;
  fName = new Char_t[strLen];
  strncpy(fName, matName, strLen);

  // Allocate CS data arrays
  fEnergyCsf = new Double_t[fDimCsf];
  fInCsf = new Double_t[fDimCsf];
  fPhCsf = new Double_t[fDimCsf];
  fCoCsf = new Double_t[fDimCsf];

  // Read cross section data
  for (Int_t i = 0; i < fDimCsf; i++)
    {
      fscanf(fp, "%lf %lf %lf %lf", &fEnergyCsf[i], &fCoCsf[i], &fInCsf[i], &fPhCsf[i]);
      fEnergyCsf[i] *= cMeV2keV;
    }

  // Initialize some variables
  fEnergyMinCsf = fEnergyCsf[0];
  fEnergyMaxCsf = fEnergyCsf[fDimCsf-1];

  Int_t fnameLength = strlen(fileName) + 1;
  fFileNameMAT = new Char_t[fnameLength];
  strncpy(fFileNameMAT, fileName, fnameLength);

  // Close the file
  fclose(fp);
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::ReadCffFile(const Char_t *fileName)
{
  // Read coherent scattering form factor data.
  //
  // Format example:
  // # Form factor, format 2.0
  // # Material: Oxygen
  // # Dimension: 123
  // # x          Fun
  // # [1/cm]     [1]
  // 0.000000e+00 8.000000e+00
  // ...
  // 1.562500e+01 1.078010e-03

  // Debug info
  fprintf(stderr, "<Info TVpMaterialFileData::ReadCffFile(%s)>\n", fileName);

  FILE *fp;
  Char_t matName[128];
  Int_t strLen;
  
  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpMaterial::ReadCffFile: Cannot open the file %s: ",
	      fileName);
      perror("");
      return 1;
    }
   
  // Read the file
  fscanf(fp, "%*[^\n]\n");                           // Skip line 1
  fscanf(fp, "# Material: %[^\n]\n", matName);       // Read line 2
  fscanf(fp, "# Dimension: %d\n", &fDimFff);         // Read line 3
  fscanf(fp, "%*[^\n]\n");                           // Skip line 4
  fscanf(fp, "%*[^\n]\n");                           // Skip line 5

  // Copy name
  strLen = strlen(matName) + 1;
  fNameCFF = new Char_t[strLen];
  strncpy(fNameCFF, matName, strLen);

  // Allocate data arrays
  fCoXf = new Double_t[fDimFff];
  fCoFff = new Double_t[fDimFff];
  
  // Read data
  for (Int_t i = 0; i < fDimFff; i++)
    fscanf(fp, "%lf %lf", &fCoXf[i], &fCoFff[i]);
  
  Int_t fnameLength = strlen(fileName) + 1;
  fFileNameCFF = new Char_t[fnameLength];
  strncpy(fFileNameCFF, fileName, fnameLength);

  // Close the file
  fclose(fp);

  fUseFf = 1;
  return 0;
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::ReadIsfFile(const Char_t *fileName)
{
  // Read incoherent scattering scatering function data
  //
  // Format example:
  // # Incoherent scattering function, format 2.0
  // # Material:  Water, ICRU-46
  // # Density: 1.000000e+00
  // # Dimension: 369
  // # x          Sun
  // # [1/cm]     [1]
  // 0.000000e+00 0.000000e+00
  // ...
  // 1.665000e+09 3.331789e+00

  // Debug info
  fprintf(stderr, "<Info TVpMaterialFileData::ReadIsfFile(%s)>\n", fileName);

  FILE *fp;
  Char_t matName[128];
  Int_t strLen;
  
  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpMaterial::ReadIsfFile: Cannot open the file %s: ",
	      fileName);
      perror("");
      return 1;
    }
   
  // Read the file
  fscanf(fp, "%*[^\n]\n");                           // Skip line 1
  fscanf(fp, "# Material: %[^\n]\n", matName);       // Read line 2
  fscanf(fp, "# Dimension: %d\n", &fDimSff);         // Read line 3
  fscanf(fp, "%*[^\n]\n");                           // Skip line 4
  fscanf(fp, "%*[^\n]\n");                           // Skip line 5

  // Copy name
  strLen = strlen(matName) + 1;
  fNameISF = new Char_t[strLen];
  strncpy(fNameCFF, matName, strLen);

  // Allocate data arrays
  fInXf = new Double_t[fDimSff];
  fInSff = new Double_t[fDimSff];
  
  // Read data
  for (Int_t i = 0; i < fDimSff; i++)
    fscanf(fp, "%lf %lf", &fInXf[i], &fInSff[i]);
  
  Int_t fnameLength = strlen(fileName) + 1;
  fFileNameISF = new Char_t[fnameLength];
  strncpy(fFileNameISF, fileName, fnameLength);

  // Close the file
  fclose(fp);
  
  fUseSf = 1;
  return 0;
}

//______________________________________________________________________________
Double_t TVpMaterialFileData::GetSff(Double_t x) const
{
  // Return S(x) from InSF

  Int_t index;
  if ((index = GetIndexSff(x)) < 0)
    return 0.0;  // out of range
  
  Double_t S;  
  if (index == 0)  // lin-lin interpolation
    S = TVpMath::GetLinLinInterpolation(fInXf[index], fInSff[index],
					fInXf[index+1], fInSff[index+1], x);
  else  // log-log interpolation
    S = TVpMath::GetLogLogInterpolation(fInXf[index], fInSff[index],
					fInXf[index+1], fInSff[index+1], x);
  return S;
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::GetIndexFff(Double_t x) const
{
  // Return index. If x is out of range, return -1 (under minimum) or -2 (over
  // maximum) and write a warning to cerr.
  //
  // -1        0         1             fDimFff-2      -2
  // -----|---------|-------.........|------------|-------------
  //  fCoXf[0]   fCoXf[1]   fCoXf[fDimFff-2]  fCoXf[fDimFff-1]

  Int_t index = TVpMath::FindIndexByBinarySearch(fDimFff, fCoXf, x);
  if (index < 0)
    std::cerr << "Warning: TVpMaterialFileData::GetIndexFff: x = "
	      << x << " out of range "
	      << "(" << fCoXf[0] << ", " << fCoXf[fDimFff-1] << ")\n";
  return index;
}

//______________________________________________________________________________
Int_t TVpMaterialFileData::GetIndexSff(Double_t x) const
{
  // Return index. If x is out of range, return -1 (under minimum) or -2 (over
  // maximum) and write a warning to cerr.
  //
  // -1        0         1             fDimSff-2      -2
  // -----|---------|-------.........|------------|-------------
  //  fInXf[0]   fInXf[1]   fInXf[fDimSff-2]  fInXf[fDimSff-1]

  Int_t index = TVpMath::FindIndexByBinarySearch(fDimSff, fInXf, x);
  if (index < 0)
    std::cerr << "Warning: TVpMaterialFileData::GetIndexSff: x = "
	      << x << " out of range "
	      << "(" << fInXf[0] << ", " << fInXf[fDimSff-1] << ")\n";
  return index;
}

//______________________________________________________________________________
Double_t TVpMaterialFileData::GetFff(Double_t x) const
{
  // Return CS

  Int_t ind;

  if ((ind = GetIndexFff(x)) < 0)
    return 0.0;  // out of range
  
  Double_t value = TVpMath::GetLinLinInterpolation(fCoXf[ind], fCoFff[ind],
						   fCoXf[ind+1], fCoFff[ind+1], x);
  return value;
}

//______________________________________________________________________________
void TVpMaterialFileData::PrintStatus()
{
  // Print status
  std::cout.setf(std::ios_base::scientific);
  std::cout << "<TVpMaterialFileData::PrintStatus>\n"
	    << "$Id: TVpMaterialFileData.cc 62 2009-06-27 10:54:08Z malusek $\n"
	    << "File name MAT: " << fFileNameMAT << '\n'
	    << "File name CFF: " << fFileNameCFF << '\n'
	    << "File name ISF: " << fFileNameISF << '\n'
	    << "Material name MAT: " << fName << '\n'
	    << "Material name CFF: " << fNameCFF << '\n'
	    << "Material name ISF: " << fNameISF << '\n'
	    << "Density: " << fDensity << '\n'
	    << "\n*** Cross sections\n"
	    << "Dimension: " << fDimCsf << '\n'
	    << "E [kev]      CoCS [cm2/g] InCS [cm2/g] PhCS [cm2/g]\n"
	    << fEnergyCsf[0] <<' '<< fCoCsf[0] <<' '<< fInCsf[0] <<' '<< fPhCsf[0] <<'\n'
	    << fEnergyCsf[1] <<' '<< fCoCsf[1] <<' '<< fInCsf[1] <<' '<< fPhCsf[1] <<'\n'
	    << "...\n"
	    << fEnergyCsf[fDimCsf-1] << ' '
	    << fCoCsf[fDimCsf-1] << ' '
	    << fInCsf[fDimCsf-1] << ' '
	    << fPhCsf[fDimCsf-1] << '\n'
	    << "\n*** Coherent scattering form factor\n"
	    << "Dimension: " << fDimFff << '\n'
	    << "CoX [1/cm]   CoFF [1]\n"
	    << fCoXf[0] << ' ' << fCoFff[0] << '\n'
	    << fCoXf[1] << ' ' << fCoFff[1] << '\n'
	    << "...\n"
	    << fCoXf[fDimFff-1] << ' ' << fCoFff[fDimFff-1] << '\n'
   	    << "\n*** Incoherent scattering scattering function\n"
	    << "Dimension: " << fDimSff << '\n'
	    << "InX [1/cm]   InSF [1]\n"
	    << fInXf[0] << ' ' << fInSff[0] << '\n'
	    << fInXf[1] << ' ' << fInSff[1] << '\n'
	    << "...\n"
	    << fInXf[fDimSff-1] << ' ' << fInSff[fDimSff-1] << '\n';
  
  
  std::cout << "<TVpMaterialFileData::PrintStatus>" << std::endl;
}

//______________________________________________________________________________
TGraph *TVpMaterialFileData::GetGraphFileData(EFileData fileData) const
{
  // Return a graph.  Skip the point (x[0] = 0, y[0]) which causes
  // problems in log-log plots for fCoFff and fInSff.

  switch (fileData)
    {
    case kCoCsf:
      return new TGraph(fDimCsf, fEnergyCsf, fCoCsf);
    case kInCsf:
      return new TGraph(fDimCsf, fEnergyCsf, fInCsf);
    case kPhCsf:
      return new TGraph(fDimCsf, fEnergyCsf, fPhCsf);
    case kFff:
      return new TGraph(fDimFff-1, fCoXf+1, fCoFff+1);
    case kSff:
      return new TGraph(fDimSff-1, fInXf+1, fInSff+1);
    default:
      std::cerr << "Error: TVpMaterialFileData::GetGraphFileData: "
		<< "incorrect file data\n";
    }
  return 0;
}

//______________________________________________________________________________
void TVpMaterialFileData::DrawGraphFileData(EFileData fileData, const Char_t *opt) const
{
  // Draw graph

  const Char_t *titleX, *titleY, *title;

  switch (fileData)
    {
    case kCoCsf:
      title = "CoCsf";
      titleX = "Energy [keV]";
      titleY = "#Sigma [cm^{2} g^{-1}]";
      break;
    case kInCsf:
      title = "InCsf";
      titleX = "Energy [keV]";
      titleY = "#Sigma [cm^{2} g^{-1}]";
      break;
    case kPhCsf:
      title = "PhCsf";
      titleX = "Energy [keV]";
      titleY = "#Sigma [cm^{2} g^{-1}]";
      break;
    case kFff:
      title = "Fff";
      titleX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      titleY = "F [1]";
      break;
    case kSff:
      title = "Sff";
      titleX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      titleY = "S [1]";
      break;
        default:
      std::cerr << "Error: TVpMaterialFileData::GetGraphFileData: "
		<< "incorrect file data\n";
      break;
    }
  TGraph *graph = GetGraphFileData(fileData);
  // graph->Draw("ALP");
  graph->Draw(opt);
  graph->GetXaxis()->SetTitle(titleX);
  graph->GetYaxis()->SetTitle(titleY);
  graph->SetTitle(title);
  graph->Draw(opt);
}
