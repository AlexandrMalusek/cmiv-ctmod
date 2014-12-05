#ifndef TVpPointDetectorArray_h
#define TVpPointDetectorArray_h

////////////////////////////////////////////////////////////////////////////////
//
// TVpPointDetectorArray
//
// 2D array of point detectors.
// 
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include "TVpVector3.h"
#include "TVpPointDetector.h"
#include "TVpGeometry.h"
#include "TVpDetectorResponse.h"

#include "TObject.h"
#include "TTree.h"
#include "fomEvent.h"
#include "TH1.h"
#include "TH2.h"

class TVpPointDetectorArray : public TVpObjectLocation
{
 public:
  Int_t             fNx;            // number of detectors in the X-direction
  Int_t             fNy;            // number of detectors in the Y-direction
  
  Int_t             fNumOfDets;     // Number of detectors in the array
  Int_t             fNumOfChannels; // Number of channels of point detectors
  Double_t          fMinEnergy;     // Lower energy cutoff
  Double_t          fMaxEnergy;     // Upper energy cutoff

  TVpPointDetector *fPointDetector; //! Array of point detectors

  TVpPointDetectorArray();
  TVpPointDetectorArray
    (Int_t nx, Int_t ny, Int_t numOfChannels, Double_t minEnergy, Double_t maxEnergy,
     TVpDetectorResponse *detectorResponse = 0);
  TVpPointDetectorArray(const Char_t *fileName);
  virtual ~TVpPointDetectorArray();
  
  inline Int_t      GetIndex(Int_t iy, Int_t ix) const;
  inline Int_t      GetNx() const;
  inline Int_t      GetNy() const;
  inline Double_t   GetMeanEstimate(Int_t iy, Int_t ix) const;
  inline Double_t   GetVarianceEstimate(Int_t iy, Int_t ix) const;
  inline Double_t   GetMeanEstimate(Int_t index) const;
  inline Double_t   GetVarianceEstimate(Int_t index) const;
  inline void       GetPosition(Int_t iy, Int_t ix, TVpVector3 *positionPtr) const;
  inline void       GetPosition(Int_t index, TVpVector3 *positionPtr) const;
  inline void       RegisterParticle(Int_t index, TVpParticle *particle);
  inline void       RegisterVirtualPhoton(TVpGeometry *geometry, TVpPhoton *photonPtr);
  void              MultiplyChannelContentBy(Double_t factor);
  inline Int_t      GetNumOfDetectors() const;
  virtual void      PrintStatus(std::ostream &out = std::cout) const;
  virtual Int_t     GetType() const {return 0;}
  virtual Int_t     WritePdaFile(const Char_t *fileName, Int_t withErrors = 0);
  virtual Int_t     WritePdaFile(FILE *fp, Int_t withErrors = 0);
  virtual Int_t     ReadPdaFile(const Char_t *fileName);
  virtual Int_t     ReadPdaFile(FILE *fp);
  Int_t             WriteValuesToTabFile(const Char_t *fileName);
  Int_t             WriteErrorsToTabFile(const Char_t *fileName);
  void              SetScoredQuantity(TVpPointDetector::EScoredQuantity scoredQuantity);
  virtual TVpPointDetectorArray* Clone(Int_t numDivX = 1, Int_t numDivY = 1) const { return 0;};
  void              Zero();
  void              UpdateCountersAtEndOfHistory(ULong_t numOfHistories);
  void              PrintYProfile(Int_t sliceX, std::ostream &out = std::cout) const;
  virtual void      SetDefaultPosition() {};
  void              SetActiveTranslation(TVpVector3 *traVec);
  void              SetActiveRotation(TVpMatrix3x3 *rotMat);
  void     SetActiveTransformation(TVpMatrix3x3 *rotMat, TVpVector3 *traVec);

  enum EDrawOption {kDrawElementPoints = 1, kDrawElementAxes = 2, kDrawPdaAxes = 4};
  void            FillFomEvent(UInt_t time, TTree *treePtr, FomEvent *fomEvent) const;
  virtual void    Draw(Int_t option = kDrawElementPoints) const;
  virtual TH2F   *GetImage(Int_t isRotated = 0) const;
  virtual TH1F   *GetYProfile(Int_t slice, Int_t isCm = 0, const Char_t *name = 0) const;
  virtual TGraph *GetYProfileGraph(Int_t slice, Int_t isCm = 0, const Char_t *name = 0) const;
  inline void     DedEventInitialize(Int_t ix, Int_t iy, const Char_t *fileName);
  inline void     DedEventClose(Int_t ix, Int_t iy);
  TH1F           *GetHistEnergySpectrum(Int_t ix, Int_t iy) const;

  ClassDef(TVpPointDetectorArray,1) // Point detector array, 2D
};

//______________________________________________________________________________
inline Int_t TVpPointDetectorArray::GetIndex(Int_t iy, Int_t ix) const
{
  // Return the detector index

  if (iy < 0 || iy >= fNy || ix < 0 || ix >= fNx)
    {
      fprintf(stderr,
	      "TVpPointDetectorArray::GetIndex: indexes out of range: iy = %d ix = %d\n",
	      iy, ix);
      exit(1);
    }
  return iy * fNx + ix;
}

//______________________________________________________________________________
inline Int_t TVpPointDetectorArray::GetNx() const
{
  return fNx;
}

//______________________________________________________________________________
inline Int_t TVpPointDetectorArray::GetNy() const
{
  return fNy;
}

//______________________________________________________________________________
inline Double_t TVpPointDetectorArray::GetMeanEstimate(Int_t iy, Int_t ix) const
{
  // Return channel content sum for a given point detector

  return fPointDetector[GetIndex(iy, ix)].GetMeanEstimate();
}

//______________________________________________________________________________
inline Double_t TVpPointDetectorArray::GetVarianceEstimate(Int_t iy, Int_t ix) const
{
  // Return channel content2 sum for a given point detector

  return fPointDetector[GetIndex(iy, ix)].GetVarianceEstimate();
}

//______________________________________________________________________________
inline Double_t TVpPointDetectorArray::GetMeanEstimate(Int_t index) const
{
  // Return channel content sum for a given point detector

  return fPointDetector[index].GetMeanEstimate();
}

//______________________________________________________________________________
inline Double_t TVpPointDetectorArray::GetVarianceEstimate(Int_t index) const
{
  // Return channel content2 sum for a given point detector

  return fPointDetector[index].GetVarianceEstimate();
}

//______________________________________________________________________________
inline void TVpPointDetectorArray::GetPosition(Int_t iy, Int_t ix,
					       TVpVector3 *positionPtr) const
{
  // Get position of the detector element

  fPointDetector[GetIndex(iy, ix)].GetPosition(positionPtr);
}

//______________________________________________________________________________
inline void TVpPointDetectorArray::GetPosition(Int_t index, TVpVector3 *positionPtr) const
{
  // Get position of the detector element

  fPointDetector[index].GetPosition(positionPtr);
}

//______________________________________________________________________________
inline Int_t TVpPointDetectorArray::GetNumOfDetectors() const
{
  return fNumOfDets;
}

//______________________________________________________________________________
inline void TVpPointDetectorArray::RegisterParticle(Int_t index, TVpParticle *particle)
{
  // Register contribution from the particle to a point detector.  This
  // routine is used by the analytic projection.

  fPointDetector[index].RegisterParticle(particle);
}

//______________________________________________________________________________
inline void TVpPointDetectorArray::RegisterVirtualPhoton(TVpGeometry *geometry,
							 TVpPhoton *photonPtr)
{
  // Register possible contribution from a particle to all detector elements
  
  for (Int_t index = 0; index < fNumOfDets; index++)
    fPointDetector[index].RegisterVirtualPhoton(geometry, photonPtr);
}

//______________________________________________________________________________
void TVpPointDetectorArray::DedEventInitialize(Int_t ix, Int_t iy, const Char_t *fileName)
{
  // Initialize collection of data for the point detector (ix, iy)

  fPointDetector[GetIndex(iy, ix)].DedEventInitialize(fileName);
}

//______________________________________________________________________________
void TVpPointDetectorArray::DedEventClose(Int_t ix, Int_t iy)
{
  // Close the event file associated with the point detector (ix, iy)
  
  fPointDetector[GetIndex(iy, ix)].DedEventClose();
}

#endif // TVpPointDetectorArray
