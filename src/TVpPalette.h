#ifndef TVpPalette_h
#define TVpPalette_h

#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>

class TVpPalette
{
 public:
  TVpPalette() {};
  virtual ~TVpPalette() {};

  void SetGrey(const Int_t colNum = 128, Int_t startIndex = 300);
  void SetW2B(const Int_t colNum = 128, Int_t startIndex = 300);
  void SetB2W(const Int_t colNum = 128, Int_t startIndex = 300);

  ClassDef(TVpPalette,1) // A color palette, obsolete
};

#endif
