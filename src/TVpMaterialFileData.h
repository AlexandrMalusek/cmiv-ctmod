#ifndef TVpMaterialFileData_h
#define TVpMaterialFileData_h

#include "TObject.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

class TVpMaterialFileData
{
 public:
  enum EFileData {kCoCsf, kInCsf, kPhCsf, kFff, kSff};
  
  Char_t     *fFileNameMAT;      //! Cross section data file name
  Char_t     *fFileNameCFF;      //! Coherent scattering Form Factor data file name
  Char_t     *fFileNameISF;      //! Scatering function data file name

  Char_t     *fName;             //! material name taken from the MAT file
  Char_t     *fNameCFF;          //! material name taken from the CFF file
  Char_t     *fNameISF;          //! material name taken from the ISF file

  Double_t    fDensity;          //  material density

  // Cross section data
  Int_t       fDimCsf;           //  Dimension of Csf arrays
  Double_t   *fEnergyCsf;        //! Energy grid of cross section data [keV]
  Double_t   *fInCsf;            //! Incoherent scattering cross section [g/cm2]
  Double_t   *fPhCsf;            //! Photoefect cross section [g/cm2]
  Double_t   *fCoCsf;            //! Coherent scattering cross section [g/cm2]
  Double_t    fEnergyMinCsf;     //  Min Csf energy [keV]
  Double_t    fEnergyMaxCsf;     //  Max Csf energy [keV]

  // Coherent scattering form factor data
  Int_t       fUseFf;            //  Use-CFF flag
  Int_t       fDimFff;           //  Dimension of the fCoXf and fCoFff arrays
  Double_t   *fCoXf;             //! x grid of CoFF data [1/cm]
  Double_t   *fCoFff;            //! Coherent scattering form factor data

  // Incoherent scattering scattering function data
  Int_t       fUseSf;            //  Use-ISF flag
  Int_t       fDimSff;           //  Dimension of the fInXf and fInSff arrays
  Double_t   *fInXf;             //! x grid of InSF data [1/cm]
  Double_t   *fInSff;            //! Incoherent scattering scattering function

  TVpMaterialFileData();
  virtual ~TVpMaterialFileData();
  inline Double_t    GetDensity() const;
  Int_t              ReadMatFile(const Char_t *filename);
  Int_t              ReadCffFile(const Char_t *filename);
  Int_t              ReadIsfFile(const Char_t *filename);
  void               PrintStatus();

  // Cross section
  Int_t              GetIndexCsf(Double_t energy) const;
  Double_t           GetInCsf(Double_t energy) const;
  Double_t           GetPhCsf(Double_t energy) const;
  Double_t           GetCoCsf(Double_t energy) const;

  // Form factor
  Int_t              GetIndexFff(Double_t x) const;
  inline Int_t       GetUseFf() const;
  Double_t           GetFff(Double_t x) const;
  inline void        SetUseFf(Int_t useCFF = 1);

  // Scatering function
  Int_t              GetIndexSff(Double_t x) const;
  inline Int_t       GetUseSf() const;
  Double_t           GetSff(Double_t x) const;
  inline void        SetUseSf(Int_t useISF = 1);

  TGraph *GetGraphFileData(EFileData fileData) const;
  void DrawGraphFileData(EFileData fileData, const Char_t *opt = "AP") const;

  ClassDef(TVpMaterialFileData,1) // Material data read from files, implementation class
};

//______________________________________________________________________________
inline Double_t TVpMaterialFileData::GetDensity() const
{
  // Return mass density in g/cm^3.

  return fDensity;
}

//______________________________________________________________________________
inline Int_t TVpMaterialFileData::GetUseFf() const
{
  // If form factors are switched off, return 0.  Otherwise return a non zero
  // value.

  return fUseFf;
}

//______________________________________________________________________________
inline Int_t TVpMaterialFileData::GetUseSf() const
{
  // If icoherent scattering functions are switched off, return 0.  Otherwise
  // return a non zero value.

  return fUseSf;
}

//______________________________________________________________________________
inline void TVpMaterialFileData::SetUseFf(Int_t useCFF)
{
  // Set the form factors switch.  useCFF=1 is the default.
  //
  // Input parameters:
  // - useCFF - value of the switch (0=CFF turned off, non-zero=CFF turned on)

  fUseFf = useCFF;
}

//______________________________________________________________________________
inline void TVpMaterialFileData::SetUseSf(Int_t useISF)
{
  // Set the icoherent scattering functions switch.  useISF=1 is the default.
  //
  // Input parameters:
  // - useISF - value of the switch (0=ISF turned off, non-zero=ISF turned on)

  fUseSf = useISF;
}

#endif // TVpMaterialFileData
