#ifndef TVpDimensionsTomograph_h
#define TVpDimensionsTomograph_h

#include "TObject.h"
#include <iostream>

class TVpDimensionsTomograph
{
 public:
  Int_t    fNdl;  // Number of detector elements in the perpendicular direction
  Int_t    fNdw;  // Number of detector elements in the axial direction
  Double_t fSad;  // Source - axis distance in cm
  Double_t fSdd;  // Source - detector distance in cm
  Double_t fDal;  // Detector array length in the  perpendicular direction in cm
  Double_t fDaw;  // Detector array width (in the axial direction) in cm
  
  TVpDimensionsTomograph();
  virtual ~TVpDimensionsTomograph() {};

  void Read(std::istream &in = std::cin);
  void Write(std::ostream &out = std::cout) const;

  Int_t    GetNdl() const {return fNdl;}
  Int_t    GetNdw() const {return fNdw;}
  Double_t GetSad() const {return fSad;}
  Double_t GetSdd() const {return fSdd;}
  Double_t GetDal() const {return fDal;}
  Double_t GetDaw() const {return fDaw;}

  void SetNdl(Int_t ndl);
  void SetNdw(Int_t ndw);
  void SetSad(Double_t sad);
  void SetSdd(Double_t sdd);
  void SetDal(Double_t dal);
  void SetDaw(Double_t daw);

  ClassDef(TVpDimensionsTomograph,1) // Dimensions of a CT scanner
};

#endif  // TVpDimensionsTomograph_h
