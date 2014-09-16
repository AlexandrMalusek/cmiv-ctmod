#ifndef TVpTrack_h
#define TVpTrack_h

#include "TObject.h"

class TVpTrack
{
 public:
  Int_t        fNumEnergyBoundary;    //  Number of energy boundaries
  Int_t        *fColor;               //! Array of colors
  Float_t      *fEnergyBoundary;      //! Array of energy boundaries
  Int_t        fNumPolyLines;         //  Number of polylines
  TPolyline3D  *fPolyLine;            //! Array of polylines
  
  TVpTrack();
  ~TVpTrack();
  void SetColors(Int_t *color);
  void SetEnergyBoundaries(Float_t *energyBoundary);
  AddEvent(Event *event);
  Draw();
  
  ClassDef(TVpTrack,1) // Visualization of particle tracks
};

#endif  // TVpTrack_h
