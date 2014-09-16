#ifndef TVpMaterial_h
#define TVpMaterial_h
#include "TObject.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include <iosfwd>
#include <cmath>
#include "TVpConstant.h"
#include "TVpMaterialGridData.h"
#include "TVpMaterialFileData.h"
#include "TVpMaterialDefaults.h"
#include "TVpIntegral.h"



class TVpMaterial : public TVpMaterialGridData, public TVpMaterialFileData,
		    public TVpIntegral
{
 public:
  enum ECrossSection {kInCS = 1, kPhCS = 2, kCoCS = 4};

  Int_t       fDimA;             //  Dimension of the A function (see the theory)
  Double_t   *fValU;             //! Values of the u variable
  Double_t   *fValA;             //! Values of the A function
  Double_t   *fPolyA0;           //! Polynomial approx. of A

 private:
  Int_t   InitializeCsgIub();
  Int_t   InitializeSfg();
  Int_t   InitializeFfg();
  Int_t   InitializeNog();
  void    InitializeATable();

 public:
  TVpMaterial();
  TVpMaterial(const Char_t *fileNameMat, const Char_t *fileNameCff = 0,
	      const Char_t *fileNameIsf = 0);
  virtual ~TVpMaterial();
  inline Char_t     *GetName();
  inline Char_t     *GetFileName();
  Int_t              GetBin(Double_t, Int_t edgeDim, Double_t *lowEdge);
  Double_t           GetA(Double_t u);
  Double_t           GetInvA(Double_t au);
  Double_t           GetCoXSample(Double_t x);

  Double_t           GetCoPdfTheta(Double_t energy, Double_t cosTheta);
  Double_t           GetInPdfTheta(Double_t energy, Double_t cosTheta);
  Double_t           GetCoPdfOmega(Double_t energy, Double_t cosTheta);
  Double_t           GetInPdfOmega(Double_t energy, Double_t cosTheta);
  Double_t           GetCoPdfIntegral(Double_t energy);
  Double_t           GetInPdfIntegral(Double_t energy);
  Double_t           GetCoPdfOmegaIntegral(Double_t energy);
  Double_t           GetInPdfOmegaIntegral(Double_t energy);

  Double_t           GetCoCsTheta(Double_t energy, Double_t cosTheta);
  Double_t           GetInCsTheta(Double_t energy, Double_t cosTheta);
  Double_t           GetKnCsTheta(Double_t energy, Double_t cosTheta);
  Double_t           GetKnCsOmega(Double_t energy, Double_t cosTheta);
  Double_t           GetKnCsI(Double_t energy);
  Double_t           GetCSSum(Double_t energy, Int_t csSum);
  inline Double_t    GetMCSSum(Double_t energy, Int_t csSum);

  Int_t              Initialize();
  Int_t              Initialize(TVpMaterialDefaults& materialDefaults);

  void               SetDensity(Double_t density);
  void               PrintATable();
  void               PrintStatusGridData(std::ostream& out = std::cout);
  inline Double_t    GetMomentumTransferX(Double_t energy, Double_t cosTheta);
  Double_t           EvaluateIntegrand(Int_t selector, Double_t par, Double_t x);

  TH2F *GetHistCoCSG(Int_t numEnery, Int_t numXi);
  TH2F *GetHistInCSG(Int_t numEnery, Int_t numTheta);
  TH1F *GetHistCoCSG(Double_t energy, Int_t numXi);

  TH1F *GetHistCoPdfXi(Double_t energy, Int_t numXi);

  TH1F *GetHistCoPdfTheta(Double_t energy, Int_t numTheta);
  TH1F *GetHistInPdfTheta(Double_t energy, Int_t numTheta);
  TH1F *GetHistCoPdfOmega(Double_t energy, Int_t numTheta);
  TH1F *GetHistInPdfOmega(Double_t energy, Int_t numTheta);

  TH1F *GetHistKnCsTheta(Double_t energy, Int_t numTheta);
  TH1F *GetHistInCsTheta(Double_t energy, Int_t numTheta);
  TH1F *GetHistCoCsTheta(Double_t energy, Int_t numTheta);

  TH1F *GetHistKnCsI(Int_t numEnergy);
  TH2F *GetHistKnCsTheta2D(Int_t numEnergy, Int_t numTheta);

  TH1F *GetHistA();

  ClassDef(TVpMaterial,1) // Material properties, user class
};

typedef TVpMaterial* TVpMaterialPtr;

//______________________________________________________________________________
inline Char_t *TVpMaterial::GetName()
{
  return fName;
}

//______________________________________________________________________________
inline Char_t *TVpMaterial::GetFileName()
{
  return fFileNameMAT;
}

//______________________________________________________________________________
inline Double_t TVpMaterial::GetMCSSum(Double_t energy, Int_t csSum)
{
  return GetCSSum(energy, csSum) * fDensity;
}

//______________________________________________________________________________
inline Double_t TVpMaterial::GetMomentumTransferX(Double_t energy, Double_t cosTheta)
{
  // Return the momentum transfer.  Here, x is in [1/cm]

  return energy * sqrt(0.5 * (1.0 - cosTheta)) / TVpConstant::hc;
}

#endif  // TVpMaterial_h
