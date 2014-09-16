#ifndef TVpBeamHardeningCorrection_h
#define TVpBeamHardeningCorrection_h

#include "TObject.h"
#include "TGraph.h"
#include "TVpMaterial.h"
#include "TVpSpectrum.h"
#include "TVpDetectorResponse.h"

class TVpBeamHardeningCorrection
{
 public:
  TVpMaterial        *fMaterialPtr;      //! Material, water by default
  TVpSpectrum         *fSpectrumPtr;     //! Spectrum
  TVpDetectorResponse *fDetResponsePtr;  //! Detector response
  Int_t                fDim;             //  Grid dimension
  Double_t            *fX;               //! Material thickness in cm
  Double_t            *fP;               //! Radiological path
  Double_t             fS;               //  Effective source intensity
  Double_t             fDx;              //  bin width
  Double_t             fEffLac;          //  Effective Lac in 1/cm

  TVpBeamHardeningCorrection(TVpMaterial *materialPtr,
			     TVpSpectrum *spectrumPtr,
			     TVpDetectorResponse *detResponsePtr,
			     Int_t dim = 1000,
			     Double_t maxThickness = 100);
  virtual ~TVpBeamHardeningCorrection();

  void            SetEffectiveLac(Double_t effLac);
  Double_t        GetGridP(Double_t x) const;
  Double_t        GetGridX(Double_t P) const;
  Double_t        CalculateP(Double_t x) const;
  Double_t        CalculateS() const;
  Double_t        Correct(Double_t I1overI0) const;
  Double_t        Correct(Double_t I2overI1, Double_t I1overI0) const;
  Double_t        Correct(Double_t I2, Double_t I1, Double_t S) const;
  inline Double_t GetS() const;
  Double_t        GetI0(Double_t distance) const;

  TGraph *GetGraphGridP() const;
  TGraph *GetGraphMonoP() const;

  ClassDef(TVpBeamHardeningCorrection,1) // Beam hardening correction
};
//______________________________________________________________________________
inline Double_t TVpBeamHardeningCorrection::GetS() const
{
  return fS;
}

#endif  // TVpBeamHardeningCorrection_h
