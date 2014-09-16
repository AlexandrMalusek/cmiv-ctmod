#ifndef TVpVoxelArray_h
#define TVpVoxelArray_h

#include "TObject.h"
#include "TH2.h"
#include "TArrayC.h"
#include "TVectorF.h"
#include "TTree.h"
#include "TVpVoxelArray.h"
#include "TVpSpectrum.h"
#include "TVpDetectorResponse.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iosfwd>
#include "TVpMaterial.h"
#include "TVpParticle.h"
#include "TVpBox.h"
#include "TVpMaterialManager.h"

class TVpVoxelArray : public TVpBox
{
 private:
  static const int fMaxNumOfTissues = 256;// Max number of tissues
  Double_t         fDeltaEdge;           // Delta vicinity of edges
  Double_t         fDeltaArray;          // Delta vicinity of the voxel array
  Double_t         fInfinity;            // A number bigger than the size of the voxel array

  void            Initialize1();
  void            Initialize2(Int_t nx, Int_t ny, Int_t nz, Double_t arraySizeX,
			     Double_t arraySizeY, Double_t arraySizeZ);
  void            SetMinMaxTissue();
  Double_t        GetEnergyBehindVoxel(TVpSpectrum *spectrum, Int_t iz, Int_t iy, Int_t ix,
				       const Char_t *direction);

 public:
  Double_t        fEnergyMinIUB;    // Global material setting
  Double_t        fEnergyMaxIUB;    // Global material setting
  Double_t        fEnergyStepIUB;   // Global material setting

  std::string     fFileNameVAM; // VAM data file name
  std::string     fFileNameVAD; // VAD data file name
  std::string     fFileNameVAN; // VAN data file name
  std::string     fFileNameMMG; // MMG data file name
  Double_t        fDpX;         // Distance between planes perpendicular to the X-direction
  Double_t        fDpY;         // Distance between planes perpendicular to the Y-direction
  Double_t        fDpZ;         // Distance between planes perpendicular to the Z-direction

  Int_t           fNx;          // Number of voxels in the X-direction
  Int_t           fNy;          // Number of voxels in the Y-direction
  Int_t           fNz;          // Number of voxels in the Z-direction
  
  UChar_t        *fTissue;      //! Array of tissue index numbers (dim = fNx*fNy*fNz)
  TVpMaterialPtr fMaterial[fMaxNumOfTissues];    //! Array of material pointers
  std::string    fTissueName[fMaxNumOfTissues];  //! Array of tissue names
  TVpMaterialManager *fMaterialManager;   //! Material manager

  Float_t        *fVoxDensity;  //! Array of mass densities  (dim = fNx*fNy*fNz)
  Float_t        *fVoxEimp;     //! Sum of imparted energies per voxel (dim = fNx*fNy*fNz)
  Float_t        *fVoxEimp2;    //! Sum of squares of imp. en. per voxel (dim = fNx*fNy*fNz)
  Float_t        *fTisEimp;     //! Sum of imp. en. per tissue (dim = fMaxNumOfTissues)
  Float_t        *fTisEimp2;    //! Sum of sq. of imp. en. per tissue (dim = fMaxNumOfTissues)

  TVpMaterial    *fWater;       //! Pointer to water material
 private:
  UChar_t         fMinTissue;   //  minimum tissue number
  UChar_t         fMaxTissue;   //  maximum tissue number
  Int_t           fUnit;        //  TH2 units; 0 = voxels, 1 = cm
  Double_t        fWatchTime;   //  Time interval between reporting (sec)
  
 public:
  TVpVoxelArray();
  TVpVoxelArray(Char_t *name, Int_t index, Int_t nx, Int_t ny, Int_t nz,
		Double_t sizeX, Double_t sizeY, Double_t sizeZ);
  TVpVoxelArray(Char_t *fileNameVAM, Char_t *name = "voxel array", Int_t index = 0);
  virtual           ~TVpVoxelArray();
  inline Int_t      GetNX() const;
  inline Int_t      GetNY() const;
  inline Int_t      GetNZ() const;
  inline Int_t      GetIndex(Int_t iz, Int_t iy, Int_t ix) const;
  inline void       GetVoxelIndices(Int_t& ai, Int_t& aj, Int_t& ak, Double_t rLoc[]) const;
  inline Double_t   GetVoxelVolume() const;
  inline Double_t   GetVoxelDensity(Int_t iz, Int_t iy, Int_t ix) const;
  Double_t          GetTissueMass(Int_t index) const;
  inline Double_t   GetVoxelLac(Int_t iz, Int_t iy, Int_t ix, Double_t energy) const;
  Double_t          GetLac(Int_t iz, Int_t iy, Int_t ix, TVpSpectrum *spectrum) const;
  TVpMaterial*      GetMaterial(Double_t rLoc[]) const;
  Int_t             GetTissueIndex(Double_t rLoc[]) const;
  virtual Double_t  GetOpticalPathInside(Double_t r[], Double_t u[], Double_t energy, 
					 Double_t distance) const;
  virtual Double_t  GetPathLength(Double_t r[], Double_t u[], Double_t energy,
				  Double_t opticalPath) const;
  Int_t             GetSubIndex(Double_t r[]) const;
  Int_t             ReadVADFile(const Char_t *fileName);
  Int_t             ReadVAMFile(const Char_t *fileName);
  Int_t             ReadVANFile(const Char_t *fileName);
  Int_t             WriteVADFile(const Char_t *fileName) const;
  Int_t             WriteVADFile(std::ostream &out = std::cout) const;
  Int_t             WriteVMDFile(const Char_t *fileName) const;
  Int_t             WriteVMDFile(std::ostream &out = std::cout) const;
  void              PrintStatus(std::ostream &out = std::cout) const;
  void              PrintTissueMassTable(std::ostream &out = std::cout) const;
  void              PrintMatTable() const;
  void              SetTissue(UChar_t tissueIndex, Int_t minX, Int_t maxX,
			      Int_t minY, Int_t maxY,
			      Int_t minZ, Int_t maxZ);
  Int_t             SetVoxelDensity(Int_t iz, Int_t iy, Int_t ix, Double_t density);
  void              SetScoreEnergyImpartedOn(Int_t option);
  void              SetScoreEnergyImpartedOff(Int_t option);

  Int_t         ReadVPFile(Char_t *filename, Int_t verbose = 0);
  Int_t         ReadPGMFile(Char_t *filePrefix, Int_t verbose = 0);
  void          WritePGMZ(Char_t *filePrefix);
	 
  TH2S         *GetMatIndSliceX(Int_t slice) const;
  TH2S         *GetMatIndSliceY(Int_t slice) const;
  TH2S         *GetMatIndSliceZ(Int_t slice) const;

  TH2F         *GetLacSliceX(Int_t slice, TVpSpectrum *spectrum) const;
  TH2F         *GetLacSliceY(Int_t slice, TVpSpectrum *spectrum) const;
  TH2F         *GetLacSliceZ(Int_t slice, TVpSpectrum *spectrum) const;
	 
  TH2F         *GetCtnSliceX(Int_t slice, TVpSpectrum *spectrum) const;
  TH2F         *GetCtnSliceY(Int_t slice, TVpSpectrum *spectrum) const;
  TH2F         *GetCtnSliceZ(Int_t slice, TVpSpectrum *spectrum) const;

  TH2F         *GetRadiographX(TVpSpectrum *spectrum, Int_t verbose = kTRUE,
			       Int_t quantity = 0,
			       TVpDetectorResponse *detectorResponse = 0) const;
  TH2F         *GetRadiographY(TVpSpectrum *spectrum, Int_t verbose = kTRUE,
			       Int_t quantity = 0,
			       TVpDetectorResponse *detectorResponse = 0) const;
  TH2F         *GetRadiographZ(TVpSpectrum *spectrum, Int_t verbose = kTRUE,
			       Int_t quantity = 0,
			       TVpDetectorResponse *detectorResponse = 0) const;

  TH2F         *GetAbsorbedDoseSliceX(Int_t slice) const;
  TH2F         *GetAbsorbedDoseSliceY(Int_t slice) const;
  TH2F         *GetAbsorbedDoseSliceZ(Int_t slice) const;

  TH2F         *GetEnergyImpartedSliceX(Int_t slice) const;
  TH2F         *GetEnergyImpartedSliceY(Int_t slice) const;
  TH2F         *GetEnergyImpartedSliceZ(Int_t slice) const;

  TH2F         *GetRandRadiographY(Int_t numPhotons) const;
  TH2F         *GetRayRadiographZ(Double_t energy) const;

  void          CalculateEnergyImpartedProjectionX(TVpSpectrum *spectrum,
						   Int_t positiveDir = 1);
  void          CalculateEnergyImpartedProjectionY(TVpSpectrum *spectrum,
						   Int_t positiveDir = 1);
  void          CalculateEnergyImpartedProjectionZ(TVpSpectrum *spectrum,
						   Int_t positiveDir = 1);

  inline void   SetWater(TVpMaterial *water);
  inline void   SetUnit(Int_t unit);
  void          PrintVoi(Char_t tissueNumber) const;
  void          PrintVoiTable() const;
  void          PrintTStatTable() const;
  void          SetSingleTissue(Int_t t);
  void          ReduceTissues(Int_t tn, TArrayC& zero, TArrayC& keep);
  void          ChangeTissue(TArrayC& oldTN, TArrayC& newTN);
  void          Segment(TVectorF &ctNumber, TVectorF &threshold);
  void          Segment(TTree *treeWithCtNumbers, TVectorF &threshold,
			Int_t verbose = 0);
  void          RebinAndSegment(TTree *treeWithCtNumbers, TVectorF &threshold,
				Int_t mx, Int_t my, Int_t mz, Int_t method = 0,
				Int_t verbose = 0);
  void          SetWatchTime(Double_t watchTime) { fWatchTime = watchTime; }
  
  void          RotateX180();
  void          RotateY180();
  void          RotateZ180();
  void          RotateX90();
  void          RotateY90();
  void          RotateZ90();

  Int_t         WriteVAMFile(Char_t *fileName);
  Int_t         WriteVANFile(Char_t *fileName);

  ClassDef(TVpVoxelArray,1) // Geometry solid: a voxel array
};

//______________________________________________________________________________
inline Int_t TVpVoxelArray::GetNX() const
{
  // Return voxel array dimension in x-direction

  return fNx;
}

//______________________________________________________________________________
inline Int_t TVpVoxelArray::GetNY() const
{
  // Return voxel array dimension in y-direction

  return fNy;
}

//______________________________________________________________________________
inline Int_t TVpVoxelArray::GetNZ() const
{
  // Return voxel array dimension in z-direction

  return fNz;
}

//______________________________________________________________________________
inline Int_t TVpVoxelArray::GetIndex(Int_t iz, Int_t iy, Int_t ix) const
{
  // Return the index of the fTissue array for a voxel (ix, iy, iz) or -1 if
  // indices are out of range.

  if (iz < 0 || iz >= fNz || iy < 0 || iy >= fNy || ix < 0 || ix >= fNx)
    {
      fprintf(stderr,
	      "Error: TVpVoxelArray::GetIndexVoxel: Index out of range: iz = %d iy = %d ix = %d\n",
	      iz, iy, ix);
      return -1;
    }
  return iz*fNy*fNx + iy*fNx + ix;
}

//______________________________________________________________________________
inline Double_t TVpVoxelArray::GetVoxelVolume() const
{
  // Return the volume of a single voxel in cm^3.  All voxel have the same
  // volume.

  return fDpX*fDpY*fDpZ;
}

//______________________________________________________________________________
inline Double_t TVpVoxelArray::GetVoxelDensity(Int_t iz, Int_t iy, Int_t ix) const
{
  // Return the mass density (in g/cm^3) of a voxel with indices (ix, iy, iz).

  return fMaterial[fTissue[GetIndex(iz,iy,ix)]]->fDensity;
}

//______________________________________________________________________________
inline void TVpVoxelArray::GetVoxelIndices(Int_t& ai, Int_t& aj, Int_t& ak,
					   Double_t r[]) const
{
  // Set voxel indices corresponding to a point r[3].  The points is specified
  // in local coordinate system (LCS) of the solid; it must be inside the
  // voxel array.
  //
  // Method:
  // Since the ray tracing routine may position r outside the voxel array
  // owing to rounding errors, a repositioning that guarantees the point r is
  // inside the voxel array is performed.  The shift given by the epsVicinity
  // parameter (= 100 nm) should be quite OK for 1 mm voxels which are common
  // in CT.
  //
  // Input:
  // - r[3] - position of the point in (LCS) in cm
  //
  // Output:
  // - ai - voxel index in x-direction
  // - aj - voxel index in y-direction
  // - ak - voxel index in z-direction

  const Double_t epsVicinity = 1.0e-5;  // cm = 100 nm
  Double_t d;

  // Move the particle behind the epsVicinity if needed
  if (r[0] < epsVicinity)
    r[0] = epsVicinity;
  else if (r[0] > (d = fSizeX - epsVicinity))
    r[0] = d;

  if (r[1] < epsVicinity)
    r[1] = epsVicinity;
  else if (r[1] > (d = fSizeY - epsVicinity))
    r[1] = d;

  if (r[2] < epsVicinity)
    r[2] = epsVicinity;
  else if (r[2] > (d = fSizeZ - epsVicinity))
    r[2] = d;

  // Calculate indices.  They shouldn't be out of range.
  ai = (Int_t) (r[0] / fDpX);
  aj = (Int_t) (r[1] / fDpY);
  ak = (Int_t) (r[2] / fDpZ);
}

inline Double_t TVpVoxelArray::GetVoxelLac(Int_t iz, Int_t iy, Int_t ix,
					   Double_t energy) const
{
  // Return the linear attenuation coefficient of a voxel.
  //
  // Input:
  // - iz, iy, ix - voxel indices
  // - energy - photon energy in keV

  // This algotithm will be redesigned to avoid repetitive calls to GetToCsg().
  Int_t d_index = GetIndex(iz, iy, ix);  // get 1d-array index
  if (d_index == -1)
    return 0;
  Int_t d_tn = fTissue[d_index];
  Double_t lac = fMaterial[d_tn]->GetToCsg(energy);
  return lac;
}

//______________________________________________________________________________
inline void TVpVoxelArray::SetUnit(Int_t unit)
{
  // Set the unit used in graphs.
  //
  // Input:
  // - unit - (0 = voxel, 1 = cm);

  fUnit = unit;
}

//______________________________________________________________________________
inline void TVpVoxelArray::SetWater(TVpMaterial *water)
{
  // Set the pointer to the material defining water.  Water is needed to
  // calculate CT numbers; it is not needed for MC simulations.
  //
  // Input:
  // - water - pointer to a material defining water

  fWater = water;
}

#endif  // TVpVoxelArray_h
