#ifndef TVpRectangleCut_h
#define TVpRectangleCut_h

#include "TObject.h"
#include "TObject.h"
#include "TH2.h"
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpGeometry.h"
#include "TVpSource.h"

class TVpRectangleCut
{
 public:
  TVpVector3  fCenter;  // Center of the rectangle
  TVpVector3  fOrigin;  // Origin of the mesh of central points.
  TVpVector3  fBaseX;   // Step between cental points in the x direction
  TVpVector3  fBaseY;   // Step between cental points in the y direction
  Double_t    fSizeX;   // Size of the rectangle in the x direction
  Double_t    fSizeY;   // Size of the rectangle in the y direction
  Int_t       fNx;      // Number of tiles in the x direction
  Int_t       fNy;      // Number of tiles in the y direction

  TVpRectangleCut(TVpVector3 center, TVpVector3 baseX, TVpVector3 baseY,
		  Double_t sizeX, Double_t sizeY, Int_t dimX, Int_t dimY);
  virtual ~TVpRectangleCut();

  TH2S *GetIndex(TVpGeometry *geometry);
  TH2F *GetLac(TVpGeometry *geometry, Double_t energy, Int_t ndivision = 1);
  TH2F *GetMCSSum(TVpGeometry *geometry, Double_t energy,
		  Int_t csSum, Int_t ndivision = 1);
  TH2F *RegisterSourceParticle(TVpSource *source, Long_t nparticles);
  Int_t WriteLac(Char_t *filename, TVpGeometry *geometry, Double_t energy,
		 Int_t ndivision = 1);
  Int_t WriteMCSSum(Char_t *filename, TVpGeometry *geometry, Double_t energy,
		    Int_t csSum, Int_t ndivision = 1);
  void  Draw();
  void  DrawSourceParticle(TVpSource *source, Long_t nparticles);

  ClassDef(TVpRectangleCut,1) // Visualisation of a rectangular cut of a geometry 
};

#endif // TVpRectangleCut_h
