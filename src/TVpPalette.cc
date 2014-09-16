//______________________________________________________________________________
//
// TVpPalette sets the color or black and white palette in ROOT. The
// implementation of TVpPalette is not quite good and may be changed
// in the future.
//
// Example:
//
// TVpPalette *palette = new TVpPalette();
// palette->SetW2B(128, 256);

#include "TVpPalette.h"

ClassImp(TVpPalette)

//______________________________________________________________________________
void TVpPalette::SetGrey(const Int_t colNum, Int_t startIndex)
{
  // Set the white->black palette.
  // colNum ....... Number of colors in the palette
  // startIndex ... starting index of allocated colors
  // This function is deprecated, use SetW2B() instead. 

  SetW2B(colNum, startIndex);
}

//______________________________________________________________________________
void TVpPalette::SetW2B(const Int_t colNum, Int_t startIndex)
{
  // Set the white->black palette.
  // colNum ....... Number of colors in the palette
  // startIndex ... starting index of allocated colors
  
  Int_t palette[colNum];

  for (Int_t i=0; i < colNum; i++)
    {
      Float_t val = 1 - i/(Float_t)colNum;
      if (! gROOT->GetColor(startIndex+i) )
	TColor *color = new TColor(startIndex+i, val, val, val, "");
      else
	{
	  TColor *color = gROOT->GetColor(startIndex+i);
	  color->SetRGB(val, val, val);
	}
      palette[i] = startIndex + i;
    }
  gStyle->SetPalette(colNum, palette);
}

//______________________________________________________________________________
void TVpPalette::SetB2W(const Int_t colNum, Int_t startIndex)
{
  // Set the black->white palette.
  // colNum ....... Number of colors in the palette
  // startIndex ... starting index of allocated colors
  
  Int_t palette[colNum];

  for (Int_t i=0; i < colNum; i++)
    {
      Float_t val = i/(Float_t)colNum;
      if (! gROOT->GetColor(startIndex+i) )
	TColor *color = new TColor(startIndex+i, val, val, val, "");
      else
	{
	  TColor *color = gROOT->GetColor(startIndex+i);
	  color->SetRGB(val, val, val);
	}
      palette[i] = startIndex + i;
    }
  gStyle->SetPalette(colNum, palette);
}
