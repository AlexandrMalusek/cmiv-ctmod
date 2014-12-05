//______________________________________________________________________________
//
// TVpVoxelArray defines a voxel array with equally sized voxels. Ray tracing
// routines treat the voxel array as a box with an internal structure.
//
// Local coordinate system: The origin is in the corner, edges coincide
// with the coordinate system axes; see the figure:
//Begin_Html
/*
<img src="png/voxel_array_1.png">
*/
//End_Html
// 
// See TVpSolid for the description of the translation and rotation matrices.
//
// TVpVoxelArray contains sizes and dimensions of the voxel array and for each
// voxel in the voxel array it contains a material index number. It also
// contains an array of pointers to TVpMaterial objects which describe
// material properties of each material contained in the voxel array.
//
// Member functions are useful for ray tracing and voxel array
// transformations.
//
// Limitations: 
// - A material manager which initializes voxel array materials is
// automatically created when a VAM file is read.  In the future, a global
// material manager may be used instead.
// - Copy and assignement constructors haven't been implemented yet.
//
// ************************************************************************
// Rotations:
// CTmod uses the following patient coordinate system:
//  - x-axis: toe-to-head direction 
//  - y-axis: left-hand-to-right-hand 
//  - z-axis: back-to-front
// It corresponds to the head-first-into-gantry patient orientation.  Voxel
// arrays in this system are labeled with the suffix _so. Alternatively, CTmod
// can use the Voxman patient coordinate system (suffix _vo):
//  - x-axis: head-to-toe direction 
//  - y-axis: right-hand-to-left-hand 
//  - z-axis: back-to-front
// It corresponds to the feet-first-to-gantry patient orientation.  Rotation
// routines:
// axis  90 deg       180 deg
// ----  ----------   ------------
//    x  RotateX90()  RotateX180()
//    y  RotateY90()  RotateY180()
//    z  RotateZ90()  RotateZ180()
//
// Example:
// gSystem->Load("libRCTmod.so");
// TVpVoxelArray *va = new TVpVoxelArray();
// va->ReadVADFile("vct_cd01.vad");
// va->RotateY90();
// va->RotateX180();
// va->RotateZ180();
// va->WriteVADFile("phantom_so.vad");
//
// ************************************************************************
// Visualization:
// axis  material index     Lac             CT number       Radiograph
// ----  -----------------  --------------  --------------  ----------------
//    x  GetMatIndSliceX()  GetLacSliceX()  GetCtnSliceX()  GetRadiographX()
//    y  GetMatIndSliceY()  GetLacSliceY()  GetCtnSliceY()  GetRadiographY()
//    z  GetMatIndSliceZ()  GetLacSliceZ()  GetCtnSliceZ()  GetRadiographZ()
//
// axis  energy imparted            absorbed dose
// ----  -------------------------  -----------------------
//    x  GetEnergyImpartedSliceX()  GetAbsorbedDoseSliceX()
//    x  GetEnergyImpartedSliceY()  GetAbsorbedDoseSliceY()
//    x  GetEnergyImpartedSliceZ()  GetAbsorbedDoseSliceZ()
//
// Each function has an example attached.  The canvas size was X:400x250,
// Y:600x250, and Z:500x250.
// 
// ************************************************************************
// Segmentation:
//
// ************************************************************************
// Mass density:
// (1) mass density of each voxel is given by the corresponding material
// density.  This concept is used in threshold segmentation with step density
// functions.  This is the default setting.
// (2) each voxel has its own mass density.  This concept is used in threshold
// segmentation with linear density functions.  To use it, read the mass
// densities via ReadVMDFile().  The feature can be turned on or off via
// SetUseVoxelDensityOn() and SetUseVoxelDensityOff(), respectively.
//
// ************************************************************************
// Energy imparted:
// The following options are provided:
// 0 - score energy imparted to the box and corresponding variance 
// 1 - score energy imparted to each tissue and corresponding variance
// 2 - score energy imparted to each voxel and corresponding variance
// Scoring can be turned on and off using SetScoreEnergyImpartedOn() and
// SetScoreEnergyImpartedOff().  Option 2 may require large amount of memory:
// a 1 GiB of RAM will be allocated to store energy imparted per voxel for a
// 512x512x512 voxel array.
//
// Example:
// gSystem->Load("libRCTmod.so");
// TVpVoxelArray *va = new TVpVoxelArray("phantom/va_bg/va_bg.vam", "va_1", 1);

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdio>
#include "TVpVoxelArray.h"

#include "TBranch.h"
#include "TVpMath.h"
#include "ctNumber.h"

ClassImp(TVpVoxelArray)

using namespace std;

//______________________________________________________________________________
TVpVoxelArray::TVpVoxelArray()
  : TVpBox()
{
  // Default constructor.
  
  Initialize1();
}

//______________________________________________________________________________
TVpVoxelArray::TVpVoxelArray(const Char_t *name, Int_t index,Int_t nx, Int_t ny, Int_t nz,
			     Double_t sizeX, Double_t sizeY, Double_t sizeZ)
  : TVpBox(name, index, sizeX, sizeY, sizeZ)
{
  // Constructor of an empty voxel array.  It is used e.g. to voxelize a
  // geometry.  Otherwise an initialization of a voxel array from a VAM file
  // is used.
  //
  // Input:
  // - name - user-defined descriptive name of the solid
  // - index - unique solid index
  // - nx - dimension of the voxel array in x-direction
  // - ny - dimension of the voxel array in y-direction
  // - nz - dimension of the voxel array in z-direction
  // - sizeX - size of voxel array in x-direction in cm
  // - sizeY - size of voxel array in y-direction in cm
  // - sizeZ - size of voxel array in z-direction in cm
  //
  // Example:
  // root [] TVpVoxelArray *vaPtr = new TVpVoxelArray("va_1", 1, 10, 15, 20,
  //         10.0, 15.0, 20.0);

  Initialize1();
  Initialize2(nx, ny, nz, sizeX, sizeY, sizeZ);
}
  
//______________________________________________________________________________
TVpVoxelArray::TVpVoxelArray(const Char_t *fileNameVAM,
			     const Char_t *name, Int_t index)
  : TVpBox(name, index, 0.0, 0.0, 0.0)
{
  // Constructor with full initialization from the VAM file.
  
  Initialize1();
  ReadVAMFile(fileNameVAM);
  SetMinMaxTissue();
  fUnit = 0;            // voxels
  fWatchTime = 5;       // seconds
}

//______________________________________________________________________________
TVpVoxelArray::~TVpVoxelArray()
{
  // Default destructor.  fWater is not owned by the TVpVoxelArray.

  delete[] fTissue;
  delete[] fVoxDensity;
  delete[] fVoxEimp;
  delete[] fVoxEimp2;
  delete[] fTisEimp;
  delete[] fTisEimp2;
  delete fMaterialManager;
}

//______________________________________________________________________________
void TVpVoxelArray::Initialize1()
{
  // Initialization of constants.

  fDeltaEdge = 1e-15;
  fDeltaArray = 0.0;    // was 1e-10;
  fInfinity = 1e30;

  // Pointers
  fTissue = 0;
  fMaterialManager = 0;
  fWater = 0;
  fVoxDensity = fVoxEimp = fVoxEimp2 = fTisEimp = fTisEimp2 = 0;

  fMinTissue = 0;       // Initial value, real value not yet known
  fMaxTissue = 255;     // Initial value, real value not yet known
  fUnit = 0;            // voxels
  fWatchTime = 5;       // seconds
}

//______________________________________________________________________________
void TVpVoxelArray::Initialize2(Int_t nx, Int_t ny, Int_t nz, Double_t arraySizeX,
			  Double_t arraySizeY, Double_t arraySizeZ)
{
  // Initialization of dimensions and arrays.
  //
  // Input:
  // - nx - dimension of the voxel array in x-direction
  // - ny - dimension of the voxel array in y-direction
  // - nz - dimension of the voxel array in z-direction
  // - arraySizeX - size of voxel array in x-direction in cm
  // - arraySizeY - size of voxel array in y-direction in cm
  // - arraySizeZ - size of voxel array in z-direction in cm

  if (nx < 0 || ny < 0 || nz < 0 ||
      arraySizeX < 0.0 || arraySizeY < 0.0 || arraySizeZ < 0.0)
    {
      cerr << "Error: TVpVoxelArray::Initialize2: A dimension or size is lower "
		<< "than zero.  Initialization is aborted.\n";
      return;
    }
  fNx = nx; fNy = ny; fNz = nz;
  fSizeX = arraySizeX;
  fSizeY = arraySizeY;
  fSizeZ = arraySizeZ;
  fDpX = arraySizeX / nx;
  fDpY = arraySizeY / ny;
  fDpZ = arraySizeZ / nz;
  
  // Delete and allocate arrays
  delete[] fTissue;
  fTissue = new UChar_t[fNx * fNy * fNz];
  if (fTissue == 0)
    cerr << "Error: TVpVoxelArray::Initialize2: Allocation of tissue indices ("
	 << fNx * fNy * fNz << " bytes) failed.\n";
}

//______________________________________________________________________________
Int_t TVpVoxelArray::ReadVAMFile(const Char_t *fileName)
{
  // Read voxel array data from a VAM (Voxel Array Master) file.  New material
  // manager is created and material data are initialized.
  //
  // Input:
  // - fileName - name of the VAM file
  //
  // Format of the VAM file:
  // # Format: VAM 2.0
  // # MMG file: <name_of_the_MMG_file>
  // # VAN file: <name_of_the_VAN_file>
  // # VAD file: <name_of_the_VAD_file>
  //
  // Format example:
  // # Format: VAM 2.0
  // # MMG file: phantom/phantom2/phantom2.mmg
  // # VAN file: phantom/phantom2/phantom2.van
  // # VAD file: phantom/phantom2/phantom2.vad
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray();
  // root [] va->ReadVAMFile("phantom/phantom2/phantom2.vam");
  // <Info TVpMaterialFileData::ReadMatFile(material/air.mat)>
  // <Info TVpMaterialFileData::ReadCffFile(material/air.cff)>
  // <Info TVpMaterialFileData::ReadIsfFile(material/air.isf)>
  // ...
  // <Info TVpMaterialFileData::ReadIsfFile(material/breast.isf)>

  std::string recName;                 // The record name e.g. "# Size"
  std::string recRest;                 // The rest of the record
  std::string name;

  std::ifstream vam(fileName);
  if (!vam)
    {
      std::cerr << "Error: VoxelArray::ReadVAMFile: Cannot open file: " 
		<< fileName << '\n';
      return 1;
    }

  // Read and process the file header
  // # Format:
  getline(vam, recName, ':');
  if (recName != std::string("# Format"))
    {
      cerr << "Error: VoxelArray::ReadVAMFile: Record \"# Format:\" is missing.\n";
      return 2;
    }
  getline(vam, recRest);
  if (recRest != std::string(" VAM 2.0"))
    {
      cerr << "Error: VoxelArray::ReadVAMFile: Format version should be:\n"
	   << "\" VAM 2.0\"" << " but is:\n"
	   << '"' << recRest << "\"\n";
      return 3;
    }
  
  // # MMG file:
  getline(vam, recName, ':');
  if (recName != std::string("# MMG file"))
    {
      cerr << "Error: VoxelArray::ReadVAMFile: Record \"# MMG File:\" is missing.\n";
      return 4;
    }
  vam >> name;
  getline(vam, recRest);
  fMaterialManager = new TVpMaterialManager(name.c_str());
  fFileNameMMG = name;

  // # VAN file:
  getline(vam, recName, ':');
  if (recName != std::string("# VAN file"))
    {
      cerr << "Error: VoxelArray::ReadVAMFile: Record \"# VAN File:\" is missing.\n";
      return 5;
    }
  vam >> name;
  getline(vam, recRest);
  ReadVANFile(name.c_str());

  // # VAD file:
  getline(vam, recName, ':');
  if (recName != std::string("# VAD file"))
    {
      cerr << "Error: VoxelArray::ReadVAMFile: Record \"# VAD File:\" is missing.\n";
      return 6;
    }
  vam >> name;
  getline(vam, recRest);
  ReadVADFile(name.c_str());

  return 0;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::ReadVANFile(const Char_t *fileName)
{
  // Read VAN file.  It contains a list of material-data file-names.  If no
  // error is detected then return 0.  See the source code for other values.
  // This function is called by ReadVAMFile().
  //
  // Input:
  // - fileName - name of the VAN file
  //
  // Format of the VAN file:
  // # Format: VAN 2.0
  // 0  <file_name_1.mat> <description_1>
  // 1  <file_name_2.mat> <description_2>
  // ...
  //
  // Format example:
  // # Format: VAN 2.0
  // 0       material/skeleton_spongiosa.mat         bone
  // 1       material/water.mat                      water

  std::string recName;                 // The record name e.g. "# Format"
  std::string recRest;                 // The rest of the record
  int tissue;                          // Tissue index
  std::string materialFileName;        // Material file name
  std::string tissueName;              // Descriptive tissue name

  ifstream van(fileName);
  if (!van)
    {
      cerr << "Error: VoxelArray::ReadVANFile: Cannot open file: " << fileName << '\n';
      return 1;
    }

  // Read and process the file header
  // # Format:
  getline(van, recName, ':');
  if (recName != std::string("# Format"))
    {
      cerr << "Error: VoxelArray::ReadVANFile: Record \"# Format:\" is missing.\n";
      return 2;
    }
  getline(van, recRest);
  if (recRest != std::string(" VAN 2.0"))
    {
      cerr << "Error: VoxelArray::ReadVANFile: Format version should be:\n"
	   << "\" VAN 2.0\"" << " but is:\n"
	   << '"' << recRest << "\"\n";
      return 3;
    }

  // Read and process the data
  while (van >> tissue >> materialFileName && getline(van, recRest))
    {
      if (tissue < 0 || tissue >= fMaxNumOfTissues)
	{
	  cerr << "Error: VoxelArray::ReadVANFile: Tissue index "
	       << tissue << " is out of range [0,255].\n";
	  return 4;
	}
      if (recRest.length() == 0)
	{
	  cerr << "Error: VoxelArray::ReadVANFile: Tissue name is missing for tissue "
	       << tissue << ".\n";
	  return 5;
	}

      fTissueName[tissue] = recRest.substr(recRest.find_first_not_of(" \t"),
					   recRest.find_last_not_of(" \t"));
      fMaterial[tissue] = fMaterialManager->RegisterMaterial(materialFileName.c_str());
    }

  fFileNameVAN = std::string(fileName);
  return 0;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::GetSubIndex(Double_t r[]) const
{
  // Return tissue index for position r which is inside the voxel array.  If
  // the position r is outside the voxel array then return -1.  See
  // GetVoxelIndices for more information.
  //
  // Input:
  // - r[3] - position in local coordinate system of the solid in cm
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/va_bg/va_bg.vam", "va_1", 1);
  // root [] Double_t r[3] = {1.5, 2.5, 3.5};
  // root [] cout << va->GetSubIndex(r) << endl;
  // 0

  Int_t ai, aj, ak, index;
  GetVoxelIndices(ai, aj, ak, r);
  index = GetIndex(ak,aj,ai);
  if (index == -1)
    return -1;
  return fTissue[index]; 
}

//______________________________________________________________________________
Double_t TVpVoxelArray::GetOpticalPathInside(Double_t r[], Double_t u[],
					     Double_t energy, Double_t distance) const
{
  // Return optical path of a line segment specified by position r[3],
  // direction u[3], and length.  The position and direction are in local
  // coordinate system of the solid (LCS). The point r must be inside the
  // voxel array.
  //
  // Input:
  // - r[3] - ray's position in LCS in cm
  // - u[3] - ray's direction in LCS.  Must be a unit vector.
  // - energy - photon energy in keV
  // - distance - path length in cm
  //
  // Note: Optical path (sometimes also called radiological path) is a line
  // integral of the linear attenuation coefficient.  It is dimensionless.
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/va_bg/va_bg.vam", "va_1", 1);
  // root [] Double_t r[3] = {0.5, 0.5, 0.5};
  // root [] Double_t u[3] = {1, 0, 0};
  // root [] cout << va->GetOpticalPathInside(r, u, 60, 9) << endl;
  // 2.52793

  Int_t plane;                      // voxel entrance side index
  Int_t aiCur, ajCur, akCur;        // current voxel indexes
  Int_t aiNew, ajNew, akNew;        // next voxel indexes
  Int_t di, dj, dk;                 // voxel index step (-1,0,1)
  Double_t dx, dy, dz;
  Double_t tOld, tOldX, tOldY, tOldZ;
  Double_t tNew, tNewX, tNewY, tNewZ;
  Double_t opticalPath = 0.0;

  // Calculate voxel indices aiNew, ajNew, and akNew (ranging from 0)
  // corresponding to the point r.
  GetVoxelIndices(aiNew, ajNew, akNew, r);

  // Find:
  // 1. voxel index increments di, dj, and dk (possible values = -1, 0, 1)
  // 2. lengths dx, dy, and dz of line segments delimited by adjacent
  // voxel-defining planes.
  // 3. distances tOldX, tOldY, and tOldZ to nearest intersections (in the
  // opposite direction) with voxel-defining planes.  These are negative
  // numbers or 0.

  // X component
  if (u[0] == 0.0)
    {
      di = 0;
      dx = fInfinity;
      tOldX = -fInfinity/2;
    }
  else if (u[0] > 0.0)
    {
      di = 1;
      dx = fDpX / u[0];
      tOldX = (fDpX * aiNew - r[0]) / u[0];
    }
  else  // (u[0] < 0.0)
    {
      di = -1;
      dx = -fDpX / u[0];
      tOldX = (fDpX * (aiNew + 1) - r[0]) / u[0];
    }

  // Y component
  if (u[1] == 0.0)
    {
      dj = 0;
      dy = fInfinity;
      tOldY = -fInfinity/2;
    }
  else if (u[1] > 0.0)
    {
      dj = 1;
      dy = fDpY / u[1];
      tOldY = (fDpY * ajNew - r[1]) / u[1];
    }
  else  // (u[1] < 0.0)
    {
      dj = -1;
      dy = -fDpY / u[1];
      tOldY = (fDpY * (ajNew + 1) - r[1]) / u[1];
    }

  // Z component
  if (u[2] == 0.0)
    {
      dk = 0;
      dz = fInfinity;
      tOldZ = -fInfinity/2;
    }
  else if (u[2] > 0.0)
    {
      dk = 1;
      dz = fDpZ / u[2];
      tOldZ = (fDpZ * akNew - r[2]) / u[2];
    }
  else  // (u[2] < 0.0)
    {
      dk = -1;
      dz = -fDpZ / u[2];
      tOldZ = (fDpZ * (akNew + 1) - r[2]) / u[2];
    }  

  tOld = 0.0;  // the path starts in the point r, thus tOld must be 0.0

  do
    {
      // Here is the main loop. First, calculate tNewX, tNewY, and tNewZ (new
      // distances between voxel-defining planes and the point r).  Select the
      // shortest one (tNew).  Test for intersections near edges and update
      // the current voxel indexes (aiOld, ajOld, akOld) accordingly.
      // a?Cur ... current voxel
      // a?New ... next voxel

      // Switch to the new voxel
      aiCur = aiNew;
      ajCur = ajNew;
      akCur = akNew;
      
      // Calculate new distances to intersections of the ray with
      // voxel-defining planes.
      tNewX = tOldX + dx;
      tNewY = tOldY + dy;
      tNewZ = tOldZ + dz;
      
      // Select shortest distance tNew = min(tNewX, tNewY, tNewZ).  The
      // variable "plane" defines which plane is crossed. (X = 1, Y = 2, Z =
      // 3).
      tNew = tNewX;
      plane = 0;
      if (tNewY < tNew)
	{
	  tNew = tNewY;
	  plane = 1;
	}
      if (tNewZ < tNew)
	{
	  tNew = tNewZ;
	  plane = 2;
	}
      
      // Calculate new voxel indices and check if the corresponding tNew? is
      // in the delta-vicinity of tNew.  If yes then change the appropriate
      // voxel index.
      //
      // As a result, the algorithm "jumps" over edges and does not calculate
      // the contribution to optical path which originates from the
      // delta-vicinity.  It is assumed that the resulting error is small.
      // The benefit is that the algorithm will not get into troubles due to
      // rounding errors close to edges and vertices.
      switch (plane)
	{
	case 0:                            // Crossing X plane
	  tOldX = tNewX;
	  aiNew += di;
	  if (tNewY - tNewX < fDeltaEdge)  // We already know (tnewX <= tNewY)
	    {                              // therefore, no need for fabs()
	      tOldY = tNewY;
	      ajNew += dj;
	    }
	  if (tNewZ - tNewX < fDeltaEdge)
	    {
	      tOldZ = tNewZ;
	      akNew += dk;
	    }
	  break;

	case 1:                            // Crossing Y plane
	  tOldY = tNewY;
	  ajNew += dj;
	  if (tNewX - tNewY < fDeltaEdge)
	    {
	      tOldX = tNewX;
	      aiNew += di;
	    }
	  if (tNewZ - tNewY < fDeltaEdge)
	    {
	      tOldZ = tNewZ;
	      akNew += dk;
	    }
	  break;
	  
	case 2:                       // Crossing Z plane
	  tOldZ = tNewZ;
	  akNew += dk;
	  if (tNewX - tNewZ < fDeltaEdge)
	    {
	      tOldX = tNewX;
	      aiNew += di;
	    }
	  if (tNewY - tNewZ < fDeltaEdge)
	    {
	      tOldY = tNewY;
	      ajNew += dj;
	    }
	  break;
	}

#if 0
      Int_t d_index = GetIndex(akCur, ajCur, aiCur);  // get 1d-array index

      // Something is terribly wrong if the current voxel is outside the voxel
      // array.  Print an error message and return.
      if (d_index == -1)
	return opticalPath;

      // Find the linear attenuation coefficient
      Int_t d_tn = fTissue[d_index];
      Double_t lac = fMaterial[d_tn]->GetToCsg(energy);
#endif
      Double_t lac = GetVoxelLac(akCur, ajCur, aiCur, energy);

      if (tNew >= distance)
	{
	  // the particle is at the end
	  opticalPath += (distance - tOld) * lac;
	  return opticalPath;
	}
      else
	opticalPath += (tNew - tOld) * lac;

      tOld = tNew;
    }
  // Repeat only when the new voxel is inside the voxel array.
  while (aiNew >= 0 && aiNew < fNx &&
	 ajNew >= 0 && ajNew < fNy &&
	 akNew >= 0 && akNew < fNz);
    
  return opticalPath;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::ReadVADFile(const Char_t *fileName)
{
  // Read VAD file.  It contains voxel array size, dimensions, and material
  // index of each voxel.  If no error is detected then return 0.  See the
  // source code for other values.  This function is called by ReadVAMFile().
  //
  // Input:
  // - fileName - name of the VAD file
  // Format of the VAN file:
  // # Format: VAD 2.0
  // # Size: <sizeX> <sizeY> <sizeZ> cm
  // # Dimension: <nx> <ny> <nz>
  // # Index size: 1 byte
  // <h(1,1,1)>...<h(1,1,nz)>
  // ...
  // <h(1,ny,1)>...<h(1,ny,nz)>
  //
  // <h(2,1,1)>...<h(2,1,nz)>
  // ...
  // <h(2,ny,1)>...<h(2,ny,nz)>
  // ...
  // <h(nx,ny,1)>...<h(nx,ny,nz)>
  //
  // Notes: (1) hexidecimal numbers are used to represent bytes, (2) ny x nz
  // matrices are separated by a newline.
  //
  // Format example:
  // # Format: VAD 2.0
  // # Size: 10 10 10 cm
  // # Dimension: 10 10 10
  // # Index size: 1 byte
  // 001F0000000000000000
  // 00000000000000000000
  // 00000000000000000000
  // 
  // 00000000000000000000
  // 00000000000000000000

  char c;
  std::string recName;                 // The record name e.g. "# Size"
  std::string recRest;                 // The rest of the record

  ifstream vad(fileName);
  if (!vad)
    {
      cerr << "Error: VoxelArray::ReadVADFile: Cannot open file: " << fileName << '\n';
      return 1;
    }

  // Read and process the file header
  // # Format:
  getline(vad, recName, ':');
  if (recName != std::string("# Format"))
    {
      cerr << "Error: VoxelArray::ReadVADFile: Record \"# Format:\" is missing.\n";
      return 2;
    }
  getline(vad, recRest);
  if (recRest != std::string(" VAD 2.0"))
    {
      cerr << "Error: VoxelArray::ReadVADFile: Format version should be:\n"
	   << "\" VAD 2.0\"" << " but is:\n"
	   << '"' << recRest << "\"\n";
      return 3;
    }

  // # Size:
  getline(vad, recName, ':');
  if (recName != std::string("# Size"))
    {
      cerr << "Error: VoxelArray::ReadVADFile: Record \"# Size:\" is missing.\n";
      return 4;
    }
  vad >> fSizeX >> fSizeY >> fSizeZ;
  getline(vad, recRest);
  
  // # Dimension:
  getline(vad, recName, ':');
  if (recName != std::string("# Dimension"))
    {
      cerr << "Error: VoxelArray::ReadVADFile: Record \"# Dimension:\" is missing.\n";
      return 5;
    }
  vad >> fNx >> fNy >> fNz;
  getline(vad, recRest);
  
  // # Index size:
  getline(vad, recName, ':');
  if (recName != std::string("# Index size"))
    {
      cerr << "Error: VoxelArray::ReadVADFile: Record \"# Index size:\" is missing.\n";
      return 6;
    }
  getline(vad, recRest);

  // Alocate the arrays, the re-asignement of values is harmless
  Initialize2(fNx, fNy, fNz, fSizeX, fSizeY, fSizeZ);

  // Read the data
  Int_t tissue;
  Char_t num[3];
  Char_t *endptr;
  num[2] = '\0';
  for (Int_t i = 0; i < fNx; i++)
    {
      for (Int_t j = 0; j < fNy; j++)
	{
	  for (Int_t k = 0; k < fNz; k++)
	    {
	      vad.read(num, 2);
	      tissue = strtol(num, &endptr, 16);
	      if (*endptr != '\0')
		{
		  cerr << "Error: VoxelArray::ReadVADFile: \"" << num
		       << "\" is not a hexadecimal number.\n";
		  return 8;
		}
	      fTissue[GetIndex(k,j,i)] = tissue;
	    }
	  vad.get(c);                                   // Read '\n'
	  if (c != '\n')
	    {
	      cerr << "Error: VoxelArray::ReadVADFile: \"" << c
		   << "\" read instead of '\\n', file line = "
		   << i*(fNy+1) + j + 5 << "\n.";
	      return 9;
	    }
	}
      vad.get(c);                                       // Read '\n'
      if (c != '\n')
	{
	  cerr << "Error: VoxelArray::ReadVADFile: \"" << c
	       << "\" read instead of '\\n', file line = "
	       << (i+1)*(fNy+1) + 4 << "\n.";
	  return 10;
	} 
    }
  
  if (!vad.good())
    {
      cerr << "Error: VoxelArray::ReadVADFile: An error occured while reading data.\n";
      return 9;
    }
  
  vad.close();

  fFileNameVAD = std::string(fileName);

  return 0;
} 

//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVADFile(const Char_t *fileName) const
{
  // Write a VAD file.  If no error is detected then return 0.  See the source
  // code for other values and ReadVADFile for the VAD format description.
  // This function may be used for saving a voxelized geometry to a file.
  //
  // Input:
  // - fileName - name of the VAD file
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/va_bg/va_bg.vam", "va_1", 1);
  // root [] va->WriteVADFile("va.vad");

  std::ofstream out(fileName);
  if (!out)
    {
      std::cerr << "Error: TVpVoxelArray::WriteVADFile: Cannot open file: "
                << fileName << '\n';
      return 1;
    }
  return WriteVADFile(out);
}

//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVADFile(std::ostream &out) const
{
  // Write VAD file.  If no error is detected then return 0.  See
  // WriteVADFile(const Char_t *fileName) for more information.
  //
  // Input:
  // - out - output stream (default = cout)

  Int_t tissue;

  out << "# Format: VAD 2.0\n"
      << "# Size: " << fSizeX << ' ' << fSizeY << ' ' << fSizeZ << " cm\n"
      << "# Dimension: " << fNx << ' ' << fNy << ' ' << fNz << '\n'
      << "# Index size: 1 byte\n";
  
  out.setf(ios_base::hex, ios_base::basefield);
  out.setf(ios_base::uppercase);
  for (Int_t ix = 0; ix < fNx; ix++)
    {
      for (Int_t iy = 0; iy < fNy; iy++)
	{
	  for (Int_t iz = 0; iz < fNz; iz++)
	    {
	      tissue = fTissue[GetIndex(iz, iy, ix)];
	      out << std::setw(2) << std::setfill('0') << tissue;
	    }
	  out << '\n';
	}
      out << '\n';
    }
  
  return 0;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVMDFile(const Char_t *fileName) const
{
  // Write a VMD file.  If no error is detected then return 0.  See the source
  // code for other values and ReadVMDFile for the VMD format description.
  // This function may be used for saving a voxelized geometry to a file.
  //
  // Input:
  // - fileName - name of the VMD file

  std::ofstream out(fileName);
  if (!out)
    {
      std::cerr << "Error: TVpVoxelArray::WriteVMDFile: Cannot open file: "
                << fileName << '\n';
      return 1;
    }
  return WriteVMDFile(out);
}

//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVMDFile(std::ostream &out) const
{
  // Write VMD file.  If no error is detected then return 0.  See
  // WriteVMDFile(const Char_t *fileName) for more information.
  //
  // Input:
  // - out - output stream (default = cout)

  out << scientific;
  for (Int_t ix = 0; ix < fNx; ix++)
    for (Int_t iy = 0; iy < fNy; iy++)
      for (Int_t iz = 0; iz < fNz; iz++)
	out << fVoxDensity[GetIndex(iz, iy, ix)] << '\n';
  return 0;
}

//______________________________________________________________________________
Double_t TVpVoxelArray::GetPathLength(Double_t r[], Double_t u[], Double_t energy,
				      Double_t opticalPath) const
{
  // Return path length of a line segment specified by position r[3],
  // direction u[3], and optical path.  The position and direction are in local
  // coordinate system of the solid (LCS). The point r must be inside the
  // voxel array.
  //
  // Input:
  // - r[3] - ray's position in LCS in cm
  // - u[3] - ray's direction in LCS.  Must be a unit vector.
  // - energy - photon energy in keV
  // - opticalPath - optical path (dimensionless)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/va_bg/va_bg.vam", "va_1", 1);
  // root [] Double_t r[3] = {0.5, 0.5, 0.5};
  // root [] Double_t u[3] = {1, 0, 0};
  // root [] cout << va->GetPathLength(r, u, 60, 2.52793) << endl;
  // 9.00001

  Int_t plane;                      // voxel entrance side index
  Int_t aiCur, ajCur, akCur;        // current voxel indexes
  Int_t aiNew, ajNew, akNew;        // next voxel indexes
  Int_t di, dj, dk;                 // voxel index step (-1,0,1)
  Double_t dx, dy, dz;
  Double_t tOld, tOldX, tOldY, tOldZ;
  Double_t tNew, tNewX, tNewY, tNewZ;
  Double_t oldOpticalPath, newOpticalPath;
  Double_t realPath;
  Double_t tmcs;                    // Total macroscopic cross section

  oldOpticalPath = newOpticalPath = 0.0;

  GetVoxelIndices(aiNew, ajNew, akNew, r);

  // Find:
  // 1. voxel index increments di, dj, and dk (possible values = -1, 0, 1)
  // 2. lengths dx, dy, and dz of line segments delimited by adjacent
  // voxel-defining planes.
  // 3. distances tOldX, tOldY, and tOldZ to nearest intersections (in the
  // opposite direction) with voxel-defining planes.  These are negative
  // numbers or 0.

  // X component
  if (u[0] == 0.0)
    {
      di = 0;
      dx = fInfinity;
      tOldX = -fInfinity/2;
    }
  else if (u[0] > 0.0)
    {
      di = 1;
      dx = fDpX / u[0];
      tOldX = (fDpX * aiNew - r[0]) / u[0];
    }
  else  // (u[0] < 0.0)
    {
      di = -1;
      dx = -fDpX / u[0];
      tOldX = (fDpX * (aiNew + 1) - r[0]) / u[0];
    }

  // Y component
  if (u[1] == 0.0)
    {
      dj = 0;
      dy = fInfinity;
      tOldY = -fInfinity/2;
    }
  else if (u[1] > 0.0)
    {
      dj = 1;
      dy = fDpY / u[1];
      tOldY = (fDpY * ajNew - r[1]) / u[1];
    }
  else  // (u[1] < 0.0)
    {
      dj = -1;
      dy = -fDpY / u[1];
      tOldY = (fDpY * (ajNew + 1) - r[1]) / u[1];
    }

  // Z component
  if (u[2] == 0.0)
    {
      dk = 0;
      dz = fInfinity;
      tOldZ = -fInfinity/2;
    }
  else if (u[2] > 0.0)
    {
      dk = 1;
      dz = fDpZ / u[2];
      tOldZ = (fDpZ * akNew - r[2]) / u[2];
    }
  else  // (u[2] < 0.0)
    {
      dk = -1;
      dz = -fDpZ / u[2];
      tOldZ = (fDpZ * (akNew + 1) - r[2]) / u[2];
    }  

  tOld = 0.0;  // the path starts in the point r, thus tOld must be 0.0
  
  do
    {
      // Here is the main loop. First, we calculate tNewX, tNewY, and
      // tNewZ (new distances between X, Y, and Z planes and the r
      // point), and then we select the shortest one. We also test for
      // intersections near edges and vertexes and update the current
      // voxel indexes (aiOld, ajOld, akOld) properly.
      // a?Cur ... current voxel
      // a?New ... next voxel

      // Switch to new voxel
      aiCur = aiNew;
      ajCur = ajNew;
      akCur = akNew;
      
      // Calculate new distances
      tNewX = tOldX + dx;
      tNewY = tOldY + dy;
      tNewZ = tOldZ + dz;
      
      // Select minimal distance tNew = min(tNewX, tNewY, tNewZ)
      // The variable "plane" tells which plane we are going to cross. 
      tNew = tNewX;
      plane = 0;
      if (tNewY < tNew)
	{
	  tNew = tNewY;
	  plane = 1;
	}
      if (tNewZ < tNew)
	{
	  tNew = tNewZ;
	  plane = 2;
	}
      
      // Check if tNew? is in a delta vicinity of tNew and calculate new indexes
      switch (plane)
	{
	case 0:                            // Crossing X plane
	  tOldX = tNewX;
	  aiNew += di;
	  if (tNewY - tNewX < fDeltaEdge)  // We already know (tnewX <= tNewY)
	    {                              // therefore, no need for fabs()
	      tOldY = tNewY;
	      ajNew += dj;
	    }
	  if (tNewZ - tNewX < fDeltaEdge)
	    {
	      tOldZ = tNewZ;
	      akNew += dk;
	    }
	  break;

	case 1:                            // Crossing Y plane
	  tOldY = tNewY;
	  ajNew += dj;
	  if (tNewX - tNewY < fDeltaEdge)
	    {
	      tOldX = tNewX;
	      aiNew += di;
	    }
	  if (tNewZ - tNewY < fDeltaEdge)
	    {
	      tOldZ = tNewZ;
	      akNew += dk;
	    }
	  break;
	  
	case 2:                       // Crossing Z plane
	  tOldZ = tNewZ;
	  akNew += dk;
	  if (tNewX - tNewZ < fDeltaEdge)
	    {
	      tOldX = tNewX;
	      aiNew += di;
	    }
	  if (tNewY - tNewZ < fDeltaEdge)
	    {
	      tOldY = tNewY;
	      ajNew += dj;
	    }
	  break;
	}

      Int_t d_index = GetIndex(akCur, ajCur, aiCur);
      if (d_index == -1)
	return 1.0e30;   // This shouldn't happen.
      Int_t d_tn = fTissue[d_index];

      // How much will the optical path increase if we cross the voxel?
      tmcs = fMaterial[d_tn]->GetToCsg(energy);
      oldOpticalPath = newOpticalPath;
      newOpticalPath += (tNew - tOld) * tmcs;

      // If the increment is too big then the end voxel has been reached.
      // Adjust the real path in this case.
      if (newOpticalPath > opticalPath)
	{
	  realPath = tOld + (opticalPath - oldOpticalPath) / tmcs;
	  return realPath;
	}

      tOld = tNew;
    }
  while (aiNew >= 0 && aiNew < fNx &&
	 ajNew >= 0 && ajNew < fNy &&
	 akNew >= 0 && akNew < fNz);

  // The particle leaves the voxel array.  Return a big number.
  return 1.0e30;
}

//______________________________________________________________________________
TVpMaterial* TVpVoxelArray::GetMaterial(Double_t rLoc[]) const
{
  // Return the material.

  Int_t tissueIndex = GetTissueIndex(rLoc);
  return fMaterial[tissueIndex];
}

//______________________________________________________________________________
Int_t TVpVoxelArray::GetTissueIndex(Double_t rLoc[]) const
{
  // Return the tissue index 

  Int_t ai, aj, ak;
  GetVoxelIndices(ai, aj, ak, rLoc);
  Int_t globalIndex = GetIndex(ak, aj, ai);
  Int_t tissueIndex = fTissue[globalIndex];
  return tissueIndex;
}

//______________________________________________________________________________
void TVpVoxelArray::PrintStatus(std::ostream &out) const
{
  // Print the object status
  //
  // Input:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpVoxelArray *vaPtr = new TVpVoxelArray("va_1", 1, 10, 15, 20,
  //         10.0, 15.0, 20.0);
  // root [] vaPtr->PrintStatus()
  // <TVpVoxelArray>
  // Dimension: 10   15      20
  // <TVpBox>
  // Size [cm]:      10      15      20
  // <TVpSolid>
  // Name:   va_1
  // Index:  1
  // Material:       Undefined or voxel array
  // Position of origin in Universe coordinates [cm]
  //    0.00000000e+00   0.00000000e+00   0.00000000e+00
  // Rotation matrix:
  //    1.00000000e+00   0.00000000e+00   0.00000000e+00
  //    0.00000000e+00   1.00000000e+00   0.00000000e+00
  //    0.00000000e+00   0.00000000e+00   1.00000000e+00
  // Overlap solid indices:  None
  // Base solid indices:     None
  // </TVpSolid>
  // </TVpBox>
  // </TVpVoxelArray>

  out << "<TVpVoxelArray>\n"
      << "$Id: TVpVoxelArray.cc 63 2009-06-27 19:20:01Z malusek $\n"
      << "Dimension: " << fNx << '\t' << fNy << '\t' << fNz << '\n';
  TVpBox::PrintStatus(out);
  out << "</TVpVoxelArray>" << endl;
}

//______________________________________________________________________________
void TVpVoxelArray::SetMinMaxTissue()
{
  // Set the fMinTissue and fMaxTissue data memebers according to
  // voxel phantom tissue data

  Int_t nxyz = fNx * fNy * fNz;

  fMinTissue = 255;    // I suppose UChar_t is in the range 0,...,255
  fMaxTissue = 0;

  for (Int_t i = 0; i < nxyz; i++)
    {
      if (fTissue[i] < fMinTissue)
	fMinTissue = fTissue[i];
      if (fTissue[i] > fMaxTissue)
	fMaxTissue = fTissue[i];
    }
}  


//______________________________________________________________________________
Double_t TVpVoxelArray::GetTissueMass(Int_t index) const
{
  // Return the total tissue mass.  Used for organ dose calculation.

  Int_t dim = fNx*fNy*fNz;
  Int_t nVoxels = 0;
  // Calculate the number of voxels containing this tissue
  for (Int_t i = 0; i < dim; i++)
    if (fTissue[i] == index)
      nVoxels++;

  Double_t density = (fMaterial[index] == 0) ? -1.0 : fMaterial[index]->GetDensity();
  Double_t mass = nVoxels * GetVoxelVolume() * density;
  return mass;
}

//______________________________________________________________________________
Double_t TVpVoxelArray::GetLac(Int_t iz, Int_t iy, Int_t ix, TVpSpectrum *spectrum) const
{
  // Return the linear attenuation coefficient (in 1/cm) of a voxel for a
  // given spectrum or 0.0 if indices are out of range.
  //
  // Input parameters:
  // - iz - voxel index, z-axis (range=0,..,GetNZ()-1)
  // - iy - voxel index, y-axis (range=0,..,GetNY()-1)
  // - ix - voxel index, x-axis (range=0,..,GetNX()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");  
  // root [] Double_t lac = va->GetLac(30, 20, 10, new TVpSpectrum(70.0));
  // root [] cout << lac << endl;
  // 2.067869e-04
  
  Int_t index = GetIndex(iz, iy, ix);
  if (index < 0)
    return 0.0;
  Int_t matIndex = fTissue[index];
  Double_t lac = fMaterial[matIndex]->GetToCsg(spectrum);
  return lac;
}

//______________________________________________________________________________
void TVpVoxelArray::PrintTissueMassTable(std::ostream &out) const
{
  // Print masses (in g) of all tissues.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] va->PrintTissueMassTable()
  // 0       4.742865e+01    outside phantom
  // 1       1.904698e+03    skin
  // ...
  // 57      9.019977e+00    teeth

  out << std::scientific;
  for (Int_t i = fMinTissue; i <= fMaxTissue; i++)
    if (fMaterial[i] != 0)
      out << i << '\t' << GetTissueMass(i) << '\t' << fTissueName[i].c_str() << endl;
}

//______________________________________________________________________________
void TVpVoxelArray::PrintMatTable() const
{
  // List the material name for each material index.
  //
  // Columns:
  // MI = material index
  //  C = coherent scattering form factor (0=not used, 1=used)
  //  I = incoherent scattering function (0=not used, 1=used)
  // Tissue name comes from the .van file, material name comes from .mat files
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");   
  // root [] va->PrintMatTable()
  //  MI C I Tissue name          Material name
  //   0 1 1 outside phantom      Air, VOXMAN composition
  //   1 1 1 skin                 Skin, adult, ICRU-46
  // ...
  //  57 1 1 teeth                Skeleton-spongiosa, adult, ICRU-46

  printf(" MI C I Tissue name          Material name\n");
  for (Int_t i = fMinTissue; i <= fMaxTissue; i++)
    if (fMaterial[i] != 0)
      printf("%3d %d %d %-20s %s\n", i,
	     fMaterial[i]->GetUseFf(), fMaterial[i]->GetUseSf(),
	     fTissueName[i].c_str(), fMaterial[i]->GetName());
}

//______________________________________________________________________________
void TVpVoxelArray::SetTissue(UChar_t tissueIndex, Int_t minX, Int_t maxX,
			      Int_t minY, Int_t maxY,
			      Int_t minZ, Int_t maxZ)
{
  // Set the tissue index in a box inside the voxel array.  This
  // function is useful for debugging purposes only.

  if (minX < 0 || minX > maxX || maxX > fNx-1
      || minY < 0 || minY > maxY || maxY > fNy-1
      || minZ < 0 || minZ > maxZ || maxZ > fNz-1)
    {
      std::cerr << "Error: TVpVoxelArray::SetTissue: incorrect min and max values\n";
      return;
    }
  for (Int_t ix = minX; ix <= maxX; ix++)
    for (Int_t iy = minY; iy <= maxY; iy++)
      for (Int_t iz = minZ; iz <= maxZ; iz++)
	{
	  Int_t index = GetIndex(iz, iy, ix);
	  fTissue[index] = tissueIndex;
	}
}

//______________________________________________________________________________
void TVpVoxelArray::SetScoreEnergyImpartedOn(Int_t option)
{
  // Turn on the scoring of energy imparted. The option is the sum of
  // 0 - score energy imparted to the box and corresponding variance 
  // 1 - score energy imparted to each tissue and corresponding variance
  // 2 - score energy imparted to each voxel and corresponding variance
  //
  // Input parameters:
  // - option - scoring option, see above (range=0,1,2)
  //
  // Method:
  // Reallocate arrays.  If an array is allocated then it is automatically used
  // by corresponding routines.  As a consequence, this routine cannot be
  // called before the dimension of the voxel array is known for option==2.
  // Allocated arrays are zeroed.

  if (option == 1 || option == 3)  // per tissue
    {
      delete[] fTisEimp;
      delete[] fTisEimp2;
      fTisEimp = new Float_t[fMaxNumOfTissues];
      fTisEimp2 = new Float_t[fMaxNumOfTissues];
      for (Int_t i = 0; i < fMaxNumOfTissues; i++)
	fTisEimp[i] = fTisEimp2[i] = 0.0;
    }
  if (option == 2)  // per voxel
    {
      delete[] fVoxEimp;
      delete[] fVoxEimp2;
      Int_t dim = fNx*fNy*fNz;
      fVoxEimp = new Float_t[dim];
      fVoxEimp2 = new Float_t[dim];
      for (Int_t i = 0; i < dim; i++)
	fVoxEimp[i] = fVoxEimp2[i] = 0.0;
    }
}

//______________________________________________________________________________
void TVpVoxelArray::SetScoreEnergyImpartedOff(Int_t option)
{
  // Turn on the scoring of energy imparted. The option is the sum of
  // 1 - score energy imparted to each tissue and corresponding variance
  // 2 - score energy imparted to each voxel and corresponding variance
  //
  // Input parameters:
  // - option - scoring option, see above (range=1,2)
  // 
  // Method:
  // Deallocate arrays.

  if (option == 1 || option == 3)  // per tissue
    {
      delete[] fTisEimp;
      delete[] fTisEimp2;
      fTisEimp = fTisEimp2 = 0;
    }
  if (option == 2)  // per voxel
    {
      delete[] fVoxEimp;
      delete[] fVoxEimp2;
      fVoxEimp = fVoxEimp2 = 0;
    }
}


#include "TVpParticle.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TStopwatch.h"

void swap(Double_t& d1, Double_t& d2)
{
  Double_t td = d1;
  d1 = d2;
  d2 = td;
}

void swap(Int_t& i1, Int_t& i2)
{
  Int_t ti = i1;
  i1 = i2;
  i2 = ti;
}


//______________________________________________________________________________
void TVpVoxelArray::PrintVoi(Char_t tissueNumber) const
{
  // Print the Volume Of Interest for the organ labeled by tissueNumber.  If
  // the tissueNumber is out of range then the output has the min and max
  // voxel indicies swapped, see the example.
  //
  // Input parameters:
  // - tissueNumber - tissue index (range=0,..)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam"); 
  // root [] va->PrintVoi(2)  // brain
  // 3 31 49 83 8 53          // minx maxx miny maxy minz maxz
  // root [] va->PrintVoi(0)  // out of range
  // 0 237 0 127 2 78

  Int_t minX, maxX, minY, maxY, minZ, maxZ;

  minX = fNx-1; maxX = 0;
  minY = fNy-1; maxY = 0;
  minZ = fNz-1; maxZ = 0;

  for (Int_t i = 0; i < fNz; i++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t k = 0; k < fNx; k++)
	if (fTissue[GetIndex(i, j, k)] == tissueNumber)
	  {
	    if (k < minX)
	      minX = k;
	    if (k > maxX)
	      maxX = k;

	    if (j < minY)
	      minY = j;
	    if (j > maxY)
	      maxY = j;

	    if (i < minZ)
	      minZ = i;
	    if (i > maxZ)
	      maxZ = i;
	  }
  printf("%d %d %d %d %d %d\n", minX, maxX, minY, maxY, minZ, maxZ);
}

//______________________________________________________________________________
void TVpVoxelArray::PrintVoiTable() const
{
  // Print the table of Volumes Of Interest.
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam"); 
  // root [] va->PrintVoiTable()
  // #MatInd minX    maxX    minY    maxY    minZ    maxZ    Tissue name
  // 0       0       237     0       127     2       78      outside phantom
  // 1       1       236     0       127     2       78      skin
  // ...
  // 57      31      37      58      73      45      58      teeth

  Int_t tn;
  Int_t dim = fMaxTissue - fMinTissue + 1;  // Tissue Range

  // Allocate "min" and "max" arrays
  Int_t *minX = new Int_t[dim];
  Int_t *maxX = new Int_t[dim];
  Int_t *minY = new Int_t[dim];
  Int_t *maxY = new Int_t[dim];
  Int_t *minZ = new Int_t[dim];
  Int_t *maxZ = new Int_t[dim];

  // Initialize "min" and "max" arrays to border values
  for (Int_t i = 0; i < dim; i++)
    {
      minX[i] = fNx - 1;
      minY[i] = fNy - 1;
      minZ[i] = fNz - 1;
      maxX[i] = maxY[i] = maxZ[i] = 0;
    }

  for (Int_t i = 0; i < fNz; i++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t k = 0; k < fNx; k++)
	{
	  tn = fTissue[GetIndex(i, j, k)] - fMinTissue; 
	  
	  if (k < minX[tn])
	    minX[tn] = k;
	  if (k > maxX[tn])
	    maxX[tn] = k;
	  
	  if (j < minY[tn])
	    minY[tn] = j;
	  if (j > maxY[tn])
	    maxY[tn] = j;
	  
	  if (i < minZ[tn])
	    minZ[tn] = i;
	  if (i > maxZ[tn])
	    maxZ[tn] = i;
	}

  // Print the table
  // cout is not used because ROOT cannot redirect cout output to a file.
  printf("#MatInd\tminX\tmaxX\tminY\tmaxY\tminZ\tmaxZ\tTissue name\n");
  for (Int_t i = 0; i < dim; i++)
    if (minX[i] <= maxX[i])        // the same should be true for Y and Z
      printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", i+fMinTissue, minX[i], maxX[i],
	     minY[i], maxY[i], minZ[i], maxZ[i], fTissueName[i].c_str());
  
  // Delete "min" and "max" arrays
  delete [] minX;
  delete [] minY;
  delete [] minZ;
  delete [] maxX;
  delete [] maxY;
  delete [] maxZ;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::ReadVPFile(const Char_t *filename, Int_t verbose)
{
  // Read Voxel Phantom File. Read tissue value for each voxel from
  // an ASCII data file filename.

  FILE *fp;
  const Int_t stringWidth = 4; // VP format uses fixed length 3-digits fields
  Char_t sn[stringWidth];
  Int_t index = 0;

  // Open the file
  if ((fp = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "TVpVoxelArray::ReadVPFile: Cannot open the file %s: ",
	      filename);
      perror("");
      return 1;
    }
  
  // Read the file
  for (Int_t i = 0; i < fNz; i++)
    {
      for (Int_t j = 0; j < fNy; j++)
	{
	  for (Int_t k = 0; k < fNx; k++)
	    {
	      if (fgets(sn, stringWidth, fp) == NULL)
		{
		  fprintf(stderr, "TVpVoxelArray::ReadVPFile: Cannot read from file %s: ix = %d, iy = %d, iz = %d\n", filename, k, j, i);
		  return 2;
		}
	      fTissue[index++] = atoi(sn);
	    }
	  fgets(sn, 2, fp);     // skip the \n character
	}
      if (verbose)
	printf("iz = %d\n", i);
    }

  // Close the file
  fclose(fp);

  // Set the fMinTissue and fMaxTissue data members
  SetMinMaxTissue();

  return 0;
}

//______________________________________________________________________________
TH2S *TVpVoxelArray::GetMatIndSliceX(Int_t slice) const
{
  // Return a slice perpendicular to the x-axis containing material indices.
  // If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNX()-1)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // TH2S *h = va->GetMatIndSliceX(100);
  // h->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetMatIndSliceX.png">
  */
  //End_Html

  Char_t name[128];                   // histogram's name

  if ( slice < 0 || fNx <= slice )   // check the range
    {
      std::cerr << "Error: TVpVoxelArray::GetMatIndSliceX: slice out of range.\n";
      return 0;
    }
  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "Material indices, slice: X%04d", slice);
  TH2S *h = new TH2S(name, name, fNy, 0, rangeY, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      h->SetCellContent(j+1, k+1, fTissue[GetIndex(k, j, slice)]);

  h->SetMinimum(fMinTissue);
  h->SetMaximum(fMaxTissue);

  if (fUnit == 0)
    {
      h->SetXTitle("Y [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("Y / cm");
      h->SetYTitle("Z / cm");
    }

  return h;
}

//______________________________________________________________________________
TH2S *TVpVoxelArray::GetMatIndSliceY(Int_t slice) const
{
  // Return a slice perpendicular to the y-axis containing material indices.
  // If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNY()-1)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // TH2S *h = va->GetMatIndSliceY(64);
  // h->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetMatIndSliceY.png">
  */
  //End_Html

  Char_t name[128];                   // histogram's name
 
  if ( slice < 0 || fNy <= slice )   // check the range
    {
      std::cerr << "Error: TVpVoxelArray::GetMatIndSliceY: slice out of range.\n";
      return 0;
    }

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "Material indices, slice: Y%04d", slice);
  TH2S *h = new TH2S(name, name, fNx, 0, rangeX, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t i = 0; i < fNx; i++)
      h->SetCellContent(i+1 , k+1, fTissue[GetIndex(k, slice, i)]);

  h->SetMinimum(fMinTissue);
  h->SetMaximum(fMaxTissue);

  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("X / cm");
      h->SetYTitle("Z / cm");
    }

  return h;
}

//______________________________________________________________________________
TH2S *TVpVoxelArray::GetMatIndSliceZ(Int_t slice) const
{
  // Return a slice perpendicular to the z-axis containing material indices.
  // If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNZ()-1)
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // TH2S *h = va->GetMatIndSliceZ(30);
  // h->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetMatIndSliceZ.png">
  */
  //End_Html

  Char_t name[128];                   // histogram's name

  if ( slice < 0 || fNz <= slice )   // check the range
    {
      std::cerr << "Error: TVpVoxelArray::GetMatIndSliceZ: slice out of range.\n";
      return 0;
    }

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;

  sprintf(name, "Material indices, slice: Z%04d", slice);
  TH2S *h = new TH2S(name, name, fNx, 0, rangeX, fNy, 0, rangeY);
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      h->SetCellContent(i+1, j+1, fTissue[GetIndex(slice, j, i)]);

  h->SetMinimum(fMinTissue);
  h->SetMaximum(fMaxTissue);

  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Y [voxel]");
    }
  else
    {
      h->SetXTitle("X / cm");
      h->SetYTitle("Y / cm");
    }

  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetLacSliceX(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the x-axis containing values of the
  // linear attenuation coefficient (in 1/cm) corresponding to a given
  // spectrum.  If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNX()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *lac = va->GetLacSliceX(100, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] lac->SetContour(256);
  // root [] lac->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetLacSliceX.png">
  */
  //End_Html

  Char_t name[128];                 // histogram's name
  Int_t matIndex;                   // material index
  Double_t totalMCS;

  if ( slice < 0 || fNx <= slice )   // check the range
    return 0;

  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "Linear attenuation cefficient [1/cm], slice: X%04d", slice);
  TH2F *h = new TH2F(name, name, fNy, 0, rangeY, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      {
	matIndex = fTissue[GetIndex(k, j, slice)];
	totalMCS = fMaterial[matIndex]->GetToCsg(spectrum);
	h->SetCellContent(j+1, k+1, totalMCS);
      }
  
  if (fUnit == 0)
    {
      h->SetXTitle("Y [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("Y [cm]");
      h->SetYTitle("Z [cm]");
    } 
  
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetLacSliceY(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the y-axis containing values of the
  // linear attenuation coefficient (in 1/cm) corresponding to a given
  // spectrum.  If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNY()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *lac = va->GetLacSliceY(64, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] lac->SetContour(256);
  // root [] lac->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetLacSliceY.png">
  */
  //End_Html

  Char_t name[128];                 // histogram's name
  Int_t matIndex;                   // material index
  Double_t totalMCS;
 
  if ( slice < 0 || fNy <= slice )   // check the range
    return 0;

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "Linear attenuation cefficient [1/cm], slice: Y%04d", slice);
  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t i = 0; i < fNx; i++)
      {
	matIndex = fTissue[GetIndex(k, slice, i)];
	totalMCS = fMaterial[matIndex]->GetToCsg(spectrum);
	h->SetCellContent(i+1, k+1, totalMCS);
      }
  
  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("X [cm]");
      h->SetYTitle("Z [cm]");
    }

  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetLacSliceZ(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the Z-axis containing values of the
  // linear attenuation coefficient (in 1/cm) corresponding to a given
  // spectrum.  If the slice is out of range, return 0.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNZ()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *lac = va->GetLacSliceZ(30, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] lac->SetContour(256);
  // root [] lac->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetLacSliceZ.png">
  */
  //End_Html

  Char_t name[128];                 // histogram's name
  Int_t matIndex;                   // material index
  Double_t totalMCS;

  if ( slice < 0 || fNz <= slice )   // check the range
    return 0;

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;

  sprintf(name, "Linear attenuation cefficient [1/cm], slice: Z%04d", slice);
  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNy, 0, rangeY);
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      {
	matIndex = fTissue[GetIndex(slice, j, i)];
	totalMCS = fMaterial[matIndex]->GetToCsg(spectrum);
	h->SetCellContent(i+1, j+1, totalMCS);
      }

  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Y [voxel]");
    }
  else
    {
      h->SetXTitle("X [cm]");
      h->SetYTitle("Y [cm]");
    }

  return h;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetCtnSliceX(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the x-axis containing CT numbers
  // corresponding to a given spectrum.  If the slice is out of range, return
  // 0.  The material water must be set via SetWater() before this function is
  // called.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNX()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TVpMaterial *water = new TVpMaterial("material/water.mat");
  // root [] water->Initialize();
  // root [] va->SetWater(water);
  // root [] TH2F *ctn = va->GetCtnSliceX(100, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] ctn->SetContour(256);
  // root [] ctn->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetCtnSliceX.png">
  */
  //End_Html

  Char_t name[128];                         // histogram's name
  Int_t matIndex;                           // material index
  Double_t totalMCSTissue, totalMCSWater;

  if ( slice < 0 || fNx <= slice )   // check the range
    return 0;

  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "CT number, slice: X%04d", slice);
  TH2F *h = new TH2F(name, name, fNy, 0, rangeY, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      {
	matIndex = fTissue[GetIndex(k, j, slice)];
	totalMCSTissue = fMaterial[matIndex]->GetToCsg(spectrum);
	totalMCSWater = fWater->GetToCsg(spectrum);
	h->SetCellContent(j+1, k+1, 1000*(totalMCSTissue-totalMCSWater)/totalMCSWater);
      }
  
  if (fUnit == 0)
    {
      h->SetXTitle("Y [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("Y [cm]");
      h->SetYTitle("Z [cm]");
    } 
  
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetCtnSliceY(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the y-axis containing CT numbers
  // corresponding to a given spectrum.  If the slice is out of range, return
  // 0.  The material water must be set via SetWater() before this function is
  // called.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNY()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TVpMaterial *water = new TVpMaterial("material/water.mat");
  // root [] water->Initialize();
  // root [] va->SetWater(water);
  // root [] TH2F *ctn = va->GetCtnSliceY(64, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] ctn->SetContour(256);
  // root [] ctn->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetCtnSliceY.png">
  */
  //End_Html

  Char_t name[128];                         // histogram's name
  Int_t matIndex;                           // material index
  Double_t totalMCSTissue, totalMCSWater;

  if ( slice < 0 || fNy <= slice )   // check the range
    return 0;

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "CT number, slice: Y%04d", slice);
  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t i = 0; i < fNx; i++)
      {
	matIndex = fTissue[GetIndex(k, slice, i)];
	totalMCSTissue = fMaterial[matIndex]->GetToCsg(spectrum);
	totalMCSWater = fWater->GetToCsg(spectrum);
	h->SetCellContent(i+1, k+1, 1000*(totalMCSTissue-totalMCSWater)/totalMCSWater);
      }
  
  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("X [cm]");
      h->SetYTitle("Z [cm]");
    } 
  
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetCtnSliceZ(Int_t slice, TVpSpectrum *spectrum) const
{
  // Return a slice perpendicular to the z-axis containing CT numbers
  // corresponding to a given spectrum.  If the slice is out of range, return
  // 0.  The material water must be set via SetWater() before this function is
  // called.
  //
  // Input parameters:
  // - slice - voxel index of the slice (range=0,..,GetNZ()-1)
  // - spectrum - X-ray spectrum
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TVpMaterial *water = new TVpMaterial("material/water.mat");
  // root [] water->Initialize();
  // root [] va->SetWater(water);
  // root [] TH2F *ctn = va->GetCtnSliceZ(30, new TVpSpectrum(70.0));
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] ctn->SetContour(256);
  // root [] ctn->Draw("COLZ");
  //Begin_Html
  /*
    <img src="png/fig_GetCtnSliceZ.png">
  */
  //End_Html

  Char_t name[128];                         // histogram's name
  Int_t matIndex;                           // material index
  Double_t totalMCSTissue, totalMCSWater;

  if ( slice < 0 || fNz <= slice )   // check the range
    return 0;

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;

  sprintf(name, "CT number, slice: Z%04d", slice);
  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNy, 0, rangeY);
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      {
	matIndex = fTissue[GetIndex(slice, j, i)];
	totalMCSTissue = fMaterial[matIndex]->GetToCsg(spectrum);
	totalMCSWater = fWater->GetToCsg(spectrum);
	h->SetCellContent(i+1, j+1, 1000*(totalMCSTissue-totalMCSWater)/totalMCSWater);
      }
  
  if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Y [voxel]");
    }
  else
    {
      h->SetXTitle("X [cm]");
      h->SetYTitle("Y [cm]");
    } 
  
  return h;
}

//______________________________________________________________________________
void TVpVoxelArray::WritePGMZ(const Char_t *filePrefix)
{
  // Write the TVpVoxelArray in the PGM Data Format

  Char_t fn[1025];     // I don't check for buffer overflow !
  FILE *fp;
  Int_t index = 0;
  
  for (Int_t k = 0; k < fNz; k++)
    {
      // Open the file
      sprintf(fn, "%s.%d", filePrefix, k+1);
      if ((fp = fopen(fn, "w")) == NULL)
	{
	  exit(1);
	}

      // Write the file's header
      fprintf(fp, "P5\n");
      fprintf(fp, "%d %d\n", fNx, fNy);
      fprintf(fp, "255\n");

      // Write numbers in binary format
      for (Int_t i = 0; i < fNx; i++)
	for (Int_t j = 0; j < fNy; j++)
	  fputc(fTissue[index++], fp);

      // fputc(fTissue[GetIndex(k, j, i)], fp);
      
      // Close the file
      fclose(fp);
    }
}

//______________________________________________________________________________
Int_t TVpVoxelArray::ReadPGMFile(const Char_t *filePrefix, Int_t verbose)
{
  FILE *fp;
  Int_t c;
  Int_t index = 0;
  Char_t line[256];     // one line
  Char_t fn[1025];      // file name


  for (Int_t k = 0; k < fNz; k++)
    {
      sprintf(fn, "%s.%d", filePrefix, k+1);
    
      // Open the file
      if ((fp = fopen(fn, "r")) == NULL)
	{
	  fprintf(stderr, "TVpVoxelArray::ReadPGMFile: Cannot open the file %s: ",
		  fn);
	  perror("");
	  return 1;
	}
      
      /* skip 3 lines */
      fgets(line, 256, fp);
      fgets(line, 256, fp);
      fgets(line, 256, fp);
      
      // Read the data
      for (Int_t j = 0; j < fNy; j++)
	for (Int_t i = 0; i < fNx; i++)
	  if (( c = getc(fp)) == EOF)
	    {
	      fTissue[index++] = c;
	      fprintf(stderr, "TVpVoxelArray::ReadPGMFile: Cannot read from file %s: ix = %d, iy = %d\n", fn, i, j);
	      return 2;
	    } 

      if (verbose)
	printf("file = %s\n", fn);
      
      // Close the file
      fclose(fp);
    }

  // Set the fMinTissue and fMaxTissue data members
  SetMinMaxTissue();

  return 0;
}

//______________________________________________________________________________
void TVpVoxelArray::PrintTStatTable() const
{
  // Print tissue statistics.
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");  
  // root [] va->PrintTStatTable();
  // #MatInd  NumVoxels Mass       Volume     Density    Tissue name
  // #                  [g]        [cm^3]     [g/cm^3]
  // #------ ---------- ---------- ---------- ---------- --------------------
  // 0          1344341 4.743e+01  3.952e+04  1.200e-03  outside phantom
  // 1            59436 1.905e+03  1.747e+03  1.090e+00  skin
  // ...
  // 57             260 9.020e+00  7.644e+00  1.180e+00  teeth

  Double_t density, volume, mass;
  const Int_t dim = 256;             // Number of UChar_t values
  Int_t *numVox = new Int_t[dim];    // Number of voxels for each tissue
  Int_t nxyz = fNx * fNy *fNz;
  Double_t voxelVolume = GetVoxelVolume();
  
  // Set initial values of the nv array to 0
  for (Int_t i = 0; i < dim; i++)
    numVox[i] = 0;
  
  for (Int_t i = 0; i < nxyz; i++)
    numVox[fTissue[i]]++;
  
  printf("%s\t%-10s %-10s %-10s %-10s %-s\n", "#MatInd", " NumVoxels", "Mass",
	 "Volume", "Density", "Tissue name");
  printf("%s\t%-10s %-10s %-10s %-10s %-s\n", "#", "", "[g]",
	 "[cm^3]", "[g/cm^3]", "");
  printf("%s\t%-10s %-10s %-10s %-10s %-s\n", "#------", "----------",
	 "----------", "----------", "----------", "--------------------");
  
  for (Int_t i = fMinTissue; i <= fMaxTissue; i++)
    if (numVox[i] != 0)
      {
	density = fMaterial[i]->GetDensity();
	volume = numVox[i] * voxelVolume;
	mass = volume * density;
	printf("%d\t%10d %-10.3e %-10.3e %-10.3e %s\n", i, numVox[i], mass,
	       volume, density, fTissueName[i].c_str());
      }

  delete [] numVox;
}


//______________________________________________________________________________
void TVpVoxelArray::SetSingleTissue(Int_t t)
{
  for (Int_t i = 0; i < fNx * fNy * fNz; i++)
    if (fTissue[i] != 0)
      fTissue[i] = t;
}

//______________________________________________________________________________
void TVpVoxelArray::ReduceTissues(Int_t tn, TArrayC& zero, TArrayC& keep)
{
  // All tissue numbers listed in the "zero" array are set to 0.
  // All tissue numbers listed in the "keep" array are kept unchanged.
  // Remaining tissues are set to "tn".
  
  Int_t j, njz, njk, nxyz;

  njz = zero.GetSize();
  njk = keep.GetSize();
  printf("njz = %d, njk = %d\n", njz, njk);
  nxyz = fNx * fNy * fNz;

  for (Int_t i = 0; i < nxyz; i++)
    {
      // Process "zero" array
      for (j = 0; j < njz; j++)
	if (fTissue[i] == zero[j])
	  {
	    fTissue[i] = 0;
	    goto next_voxel;
	  }
      
      // Process "keep" array
      for (j = 0; j < njk; j++)
      	if (fTissue[i] == keep[j])
      	  goto next_voxel;
      
      // Set the tissue number
      fTissue[i] = tn;
      
    next_voxel:
      continue;
    }
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetRadiographX(TVpSpectrum *spectrum, Int_t verbose,
				    Int_t quantity,
				    TVpDetectorResponse *detectorResponse) const
{
  // Return a parallel projection radiograph corresponding to a given spectrum
  // and scored quantity.  Photons impinge in the x-axis direction.
  //
  // Input parameters:
  // - spectrum - X-ray spectrum
  // - verbose - verbosity level (0=silent, 1=verbose)
  // - quantity - scored quantity (0=fluence, 1=energy fluence, 2=energy imparted)
  // - detectorResponse - energy absorption efficiency function
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *phi = va->GetRadiographX(new TVpSpectrum(70.0));  // fluence
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetW2B();
  // root [] phi->SetContour(256);
  // root [] phi->Draw("colz");
  // // Detector response
  // root [] TVpDetectorResponse *dr = new TVpDetectorResponse();
  // root [] dr->ReadAerFile("detector/BaFBrI.aer");
  // root [] TH2F *eps = va->GetRadiographX(new TVpSpectrum(70.0), 1, 2, dr);
  // root [] eps->SetContour(256);
  // root [] eps->Draw("colz");  // shown below
  //Begin_Html
  /*
    <img src="png/fig_GetRadiographX.png">
  */
  //End_Html

  const Char_t *name = "Parallel projection radiograph X-axis";  // histogram's name
  Double_t sum, eSum;
  Double_t energy, weight;
  TStopwatch stopwatch;
  TStopwatch stopwatchTotal;
  
  Int_t watchTotal = fNy * fNz;
  Int_t watchCurrent = 0;

  if (quantity < 0 || quantity > 2)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographX: quantity out of range.\n";
      return 0;
    }
  if (quantity == 2 && detectorResponse == 0)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographX: detectorResponse missing.\n";
      return 0;
    } 

  if (verbose == kTRUE)
    {
      stopwatch.Start();
      stopwatchTotal.Start();
    }

  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  TH2F *h = new TH2F(name, name, fNy, 0, rangeY, fNz, 0, rangeZ);
  for (Int_t i = 0; i < fNy; i++)
    {
      for (Int_t j = 0; j < fNz; j++)
	{
	  // Integrate over the spectrum (both cont. and disc.)
	  Int_t maxChannel = spectrum->GetNumOfChannels();
	  eSum = 0.0;
	  for (Int_t channel = 0; channel < maxChannel; channel++)
	    {
	      spectrum->GetEnergyAndWeight(channel, energy, weight);
	      sum = 0.0;
	      for (Int_t k = 0; k < fNx; k++)
		sum += fMaterial[fTissue[GetIndex(j, i, k)]]->GetToCsg(energy);

	      switch (quantity)
		{
		  case 0:  // fluence
		    eSum += weight * exp(-sum * fDpX);
		    break;
		  case 1:  // energy fluence
		    eSum += energy * weight * exp(-sum * fDpX);
		    break;
		  case 2:  // detector response
		    eSum += energy * weight * exp(-sum * fDpX) *
		      detectorResponse->GetResponse(1.0, energy);
		    break;
		}
	    }
	  h->SetCellContent(i+1, j+1, eSum);
	  
	  // Check the time
	  watchCurrent++;
	  if (stopwatch.RealTime() > fWatchTime)
	    {
	      printf(" %3.0f%%", (100.0 * watchCurrent)/watchTotal);
	      stopwatch.Start(1);
	    }
	  else
	    stopwatch.Continue();
	}
    }

   if (fUnit == 0)
    {
      h->SetXTitle("Y [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("Y / cm");
      h->SetYTitle("Z / cm");
    }
  stopwatchTotal.Print();
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetRadiographY(TVpSpectrum *spectrum, Int_t verbose,
				    Int_t quantity,
				    TVpDetectorResponse *detectorResponse) const
{
  // Return a parallel projection radiograph corresponding to a given spectrum
  // and scored quantity.  Photons impinge in the y-axis direction.
  //
  // Input parameters:
  // - spectrum - X-ray spectrum
  // - verbose - verbosity level (0=silent, 1=verbose)
  // - quantity - scored quantity (0=fluence, 1=energy fluence, 2=energy imparted)
  // - detectorResponse - energy absorption efficiency function
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *phi = va->GetRadiographY(new TVpSpectrum(70.0));  // fluence
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetW2B();
  // root [] phi->SetContour(256);
  // root [] phi->Draw("colz");
  // // Detector response
  // root [] TVpDetectorResponse *dr = new TVpDetectorResponse();
  // root [] dr->ReadAerFile("detector/BaFBrI.aer");
  // root [] TH2F *eps = va->GetRadiographY(new TVpSpectrum(70.0), 1, 2, dr);
  // root [] eps->SetContour(256);
  // root [] eps->Draw("colz");  // shown below
  //Begin_Html
  /*
    <img src="png/fig_GetRadiographY.png">
  */
  //End_Html

  const Char_t *name = "Parallel projection radiograph Y-axis";  // histogram's name
  Double_t sum, eSum;
  Double_t energy, weight;
  TStopwatch stopwatch;
  TStopwatch stopwatchTotal;
  
  Int_t watchTotal = fNx * fNz;
  Int_t watchCurrent = 0;

  if (quantity < 0 || quantity > 2)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographY: quantity out of range.\n";
      return 0;
    }
  if (quantity == 2 && detectorResponse == 0)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographY: detectorResponse missing.\n";
      return 0;
    } 

  if (verbose == kTRUE)
    {
      stopwatch.Start();
      stopwatchTotal.Start();
    }

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNz, 0, rangeZ);
  for (Int_t i = 0; i < fNx; i++)
    {
      for (Int_t j = 0; j < fNz; j++)
	{
	  // Integrate over the spectrum (both cont. and disc.)
	  Int_t maxChannel = spectrum->GetNumOfChannels();
	  eSum = 0.0;
	  for (Int_t channel = 0; channel < maxChannel; channel++)
	    {
	      spectrum->GetEnergyAndWeight(channel, energy, weight);
	      sum = 0.0;
	      for (Int_t k = 0; k < fNy; k++)
		sum += fMaterial[fTissue[GetIndex(j, k, i)]]->GetToCsg(energy);

	      switch (quantity)
		{
		  case 0:  // fluence
		    eSum += weight * exp(-sum * fDpY);
		    break;
		  case 1:  // energy fluence
		    eSum += energy * weight * exp(-sum * fDpY);
		    break;
		  case 2:  // detector response
		    eSum += energy * weight * exp(-sum * fDpY) *
		      detectorResponse->GetResponse(1.0, energy);
		}
	    }
	  h->SetCellContent(i+1, j+1, eSum);
	  
	  // Check the time
	  watchCurrent++;
	  if (stopwatch.RealTime() > fWatchTime)
	    {
	      printf(" %3.0f%%", (100.0 * watchCurrent)/watchTotal);
	      stopwatch.Start(1);
	    }
	  else
	    stopwatch.Continue();
	}
    }

   if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("X / cm");
      h->SetYTitle("Z / cm");
    }
  stopwatchTotal.Print();
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetRadiographZ(TVpSpectrum *spectrum, Int_t verbose,
				    Int_t quantity,
				    TVpDetectorResponse *detectorResponse) const
{
  // Return a parallel projection radiograph corresponding to a given spectrum
  // and scored quantity.  Photons impinge in the y-axis direction.
  //
  // Input parameters:
  // - spectrum - X-ray spectrum
  // - verbose - verbosity level (0=silent, 1=verbose)
  // - quantity - scored quantity (0=fluence, 1=energy fluence, 2=energy imparted)
  // - detectorResponse - energy absorption efficiency function
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] TH2F *phi = va->GetRadiographZ(new TVpSpectrum(70.0));  // fluence
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetW2B();
  // root [] phi->SetContour(256);
  // root [] phi->Draw("colz");
  // // Detector response
  // root [] TVpDetectorResponse *dr = new TVpDetectorResponse();
  // root [] dr->ReadAerFile("detector/BaFBrI.aer");
  // root [] TH2F *eps = va->GetRadiographZ(new TVpSpectrum(70.0), 1, 2, dr);
  // root [] eps->SetContour(256);
  // root [] eps->Draw("colz");  // shown below
  //Begin_Html
  /*
    <img src="png/fig_GetRadiographZ.png">
  */
  //End_Html

  const Char_t *name = "Parallel projection radiograph Z-axis";  // histogram's name
  Double_t sum, eSum;
  Double_t energy, weight;
  TStopwatch stopwatch;
  TStopwatch stopwatchTotal;
  
  Int_t watchTotal = fNx * fNy;
  Int_t watchCurrent = 0;

  if (quantity < 0 || quantity > 2)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographZ: quantity out of range.\n";
      return 0;
    }
  if (quantity == 2 && detectorResponse == 0)
    {
      std::cerr << "Error: TVpVoxelArray::GetRadiographZ: detectorResponse missing.\n";
      return 0;
    } 

  if (verbose == kTRUE)
    {
      stopwatch.Start();
      stopwatchTotal.Start();
    }

  Float_t rangeX = (fUnit == 0)? fNx : fSizeX;
  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;

  TH2F *h = new TH2F(name, name, fNx, 0, rangeX, fNy, 0, rangeY);
  for (Int_t i = 0; i < fNx; i++)
    {
      for (Int_t j = 0; j < fNy; j++)
	{
	  // Integrate over the spectrum (both cont. and disc.)
	  Int_t maxChannel = spectrum->GetNumOfChannels();
	  eSum = 0.0;
	  for (Int_t channel = 0; channel < maxChannel; channel++)
	    {
	      spectrum->GetEnergyAndWeight(channel, energy, weight);
	      sum = 0.0;
	      for (Int_t k = 0; k < fNz; k++)
		sum += fMaterial[fTissue[GetIndex(k, j, i)]]->GetToCsg(energy);
	      
	      switch (quantity)
		{
		  case 0:  // fluence
		    eSum += weight * exp(-sum * fDpZ);
		    break;
		  case 1:  // energy fluence
		    eSum += energy * weight * exp(-sum * fDpZ);
		    break;
		  case 2:  // detector response
		    eSum += energy * weight * exp(-sum * fDpZ) *
		      detectorResponse->GetResponse(1.0, energy);
		}
	    }
	  h->SetCellContent(i+1, j+1, eSum);
	  
	  // Check the time
	  watchCurrent++;
	  if (stopwatch.RealTime() > fWatchTime)
	    {
	      printf(" %3.0f%%", (100.0 * watchCurrent)/watchTotal);
	      stopwatch.Start(1);
	    }
	  else
	    stopwatch.Continue();
	}
    }

   if (fUnit == 0)
    {
      h->SetXTitle("X [voxel]");
      h->SetYTitle("Y [voxel]");
    }
  else
    {
      h->SetXTitle("X / cm");
      h->SetYTitle("Y / cm");
    }
  stopwatchTotal.Print();
  return h;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetAbsorbedDoseSliceX(Int_t slice) const
{
  // Return
  
  return 0;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetAbsorbedDoseSliceY(Int_t slice) const
{
  // Return

  return 0;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetAbsorbedDoseSliceZ(Int_t slice) const
{
  // Return
  
  return 0;
} 


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetEnergyImpartedSliceX(Int_t slice) const
{
  // Return energy imparted.
  //
  // Example:
  // root [] TVpVoxelArray *va = new TVpVoxelArray("phantom/phantom2/phantom2.vam");
  // root [] va->SetScoreEnergyImpartedOn(2);
  // root [] va->CalculateEnergyImpartedProjectionX();
  // root [] TH2F *h = va->GetEnergyImpartedSliceX(100);
  // root [] TVpPalette *pal = new TVpPalette();
  // root [] pal->SetB2W();
  // root [] h->SetContour(256);
  // root [] h->Draw("COLZ");

  Char_t name[128];                   // histogram's name

  if ( slice < 0 || fNx <= slice )   // check the range
    {
      std::cerr << "Error: TVpVoxelArray::GetEnergyImpartedSliceX: "
		<< "slice out of range.\n";
      return 0;
    }
  if (fVoxEimp == 0)
    {
      std::cerr << "Error: TVpVoxelArray::GetEnergyImpartedSliceX: "
		<< "scoring was not turned on.\n";
      return 0;
    }

  Float_t rangeY = (fUnit == 0)? fNy : fSizeY;
  Float_t rangeZ = (fUnit == 0)? fNz : fSizeZ;

  sprintf(name, "Energy imparted, slice: X%04d", slice);
  TH2F *h = new TH2F(name, name, fNy, 0, rangeY, fNz, 0, rangeZ);
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      h->SetCellContent(j+1, k+1, fVoxEimp[GetIndex(k, j, slice)]);
  
  if (fUnit == 0)
    {
      h->SetXTitle("Y [voxel]");
      h->SetYTitle("Z [voxel]");
    }
  else
    {
      h->SetXTitle("Y / cm");
      h->SetYTitle("Z / cm");
    }

  return h;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetEnergyImpartedSliceY(Int_t slice) const
{
  // Return 
  
  return 0;
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetEnergyImpartedSliceZ(Int_t slice) const
{
  // Return

  return 0;
}


//______________________________________________________________________________
Double_t TVpVoxelArray::GetEnergyBehindVoxel(TVpSpectrum *spectrum,
					     Int_t iz, Int_t iy, Int_t ix,
					     const Char_t *direction)
{
  // Return energy behind a given voxel

  Double_t psi = 0.0;
  Double_t energy, weight, radPath;
  Int_t nChannels = spectrum->GetNumOfChannels();


  if (!strcmp(direction, "+z"))
    {
      for (Int_t channel = 0; channel < nChannels; channel++)
	{
	  spectrum->GetEnergyAndWeight(channel, energy, weight);
	  // Calculate radiological path for the given energy
	  radPath = 0.0;
	  for (Int_t jz = 0; jz <= iz; jz++)
	    radPath += GetVoxelLac(jz, iy, ix, energy);
	  radPath *= fDpZ;
	  psi += energy * weight * exp(-radPath);
	}
    }

  return radPath;
}

//______________________________________________________________________________
void TVpVoxelArray::CalculateEnergyImpartedProjectionX(TVpSpectrum *spectrum,
						       Int_t positiveDir)
{
  // Calculate
}


//______________________________________________________________________________
void TVpVoxelArray::CalculateEnergyImpartedProjectionY(TVpSpectrum *spectrum,
						       Int_t positiveDir)
{
  // Calculate
}


//______________________________________________________________________________
void TVpVoxelArray::CalculateEnergyImpartedProjectionZ(TVpSpectrum *spectrum,
						       Int_t positiveDir)
{
  // Calculate energy imparted

  Double_t ebv[fNz+1];  // energy behind a voxel
  Int_t index;
  for (Int_t ix = 0; ix < fNx; ix++)
    for (Int_t iy = 0; iy < fNy; iy++)
      {
	for (Int_t iz = 0; iz < fNz; iz++)
	  ebv[iz] = GetEnergyBehindVoxel(spectrum, iz, iy, ix,
					 (positiveDir == 1) ? "+z" : "-z");
	for (Int_t iz = 0; iz < fNz; iz++)
	  {
	    if (positiveDir == 1)
	      {
		index = GetIndex(iz, iy, ix);
		fVoxEimp[index] = (iz == 0) ? ebv[0] : ebv[iz] - ebv[iz-1];
	      }
	  }
      }
  
}


//______________________________________________________________________________
void TVpVoxelArray::ChangeTissue(TArrayC& oldTN, TArrayC& newTN)
{
  // Change tissue numbers listed in the "oldTN" array to tissue numbers
  // listed in the "newTN" array. 
  
  Int_t dimOld, dimNew, dim, nxyz;

  // Get the lower array dimension
  dimOld = oldTN.GetSize();
  dimNew = newTN.GetSize();
  dim = (dimOld < dimNew) ? dimOld : dimNew;

  nxyz = fNx * fNy * fNz;

  for (Int_t i = 0; i < nxyz; i++)
    for (Int_t j = 0; j < dim; j++)
      if (fTissue[i] == oldTN[j])
	fTissue[i] = newTN[j];

  // Keep the fMinTissue and fMaxTissue data members up to date
  SetMinMaxTissue();
}


//______________________________________________________________________________
TH2F *TVpVoxelArray::GetRandRadiographY(Int_t numPhotons) const
{
  // Create a TH2F histogram containing a noisy Roentgen scan
  // TH2F *h = vp.GetRandRadiographZ(10); h->Draw("COLZ");

  const Char_t *name = "RandRadiograph Z-axis";  // histogram's name
  const Int_t maxRand = 32767;
  Int_t i, k;
  Float_t sum, rndX, rndZ;

  TH2F *h = new TH2F(name, name, fNx, 0, fNx, fNz, 0, fNz);
  for (Int_t m = 0; m < numPhotons; m++)
    {
      rndX = rand()/(float)maxRand * fNx;
      rndZ = rand()/(double)maxRand * fNz;
      i = (Int_t) rndX;
      k = (Int_t) rndZ;
      sum = 0.0;
      for (Int_t j = 0; j < fNy; j++)
	sum += fMaterial[fTissue[GetIndex(k, j, i)]]->GetToCsg(100);
      h->Fill(rndX, rndZ, sum);
    }
      
  return h;
}

//______________________________________________________________________________
TH2F *TVpVoxelArray::GetRayRadiographZ(Double_t energy) const
{
  // Create a TH2F histogram containing a Roentgen scan
  // TH2F *h = vp.GetRadiographZ(10); h->Draw("COLZ");

  const Char_t *name = "Radiograph Z-axis";  // histogram's name
  Float_t eps = 0.01;
  Double_t r[3], u[3];

  u[0] = 0.0; u[1] = 0.0; u[2] = 1.0; 
  r[2] = -1.0;
  TH2F *h = new TH2F(name, name, fNx, 0, fNx, fNy, 0, fNy);
  for (Int_t ix = 0; ix < fNx-1; ix++)
    for (Int_t iy = 0; iy < fNy-1; iy++)
      {
	r[0] = ix * fDpX + eps;
	r[1] = iy * fDpY + eps;
	Double_t opticalPath = GetOpticalPathInside(r, u, energy, 1e30);
	h->SetBinContent(ix+1, iy+1, exp(-opticalPath));
      }

  return h;
}

//______________________________________________________________________________
void TVpVoxelArray::RotateX180()
{
  // Rotate the voxel array 180 deg around the X axis.  No temporary array is
  // generated.  See the class description for an example.
  
  UChar_t tis;
  Int_t ind1, ind2, ny2, nz2;

  nz2 = fNz / 2;
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t k = 0; k < nz2; k++)
	{
	  ind1 = GetIndex(k, j, i);
	  ind2 = GetIndex(fNz-1-k, fNy-1-j, i);
	  tis = fTissue[ind1];
	  fTissue[ind1] = fTissue[ind2];
	  fTissue[ind2] = tis;
	}

  if ( fNz % 2 != 0)           // fNz is odd
    {
      ny2 = fNy/2;
      for (Int_t i = 0; i < fNx; i++)
	for (Int_t j = 0; j < ny2; j++)
	  {
	    ind1 = GetIndex(nz2, j, i);
	    ind2 = GetIndex(nz2, fNy-1-j, i);
	    tis = fTissue[ind1];
	    fTissue[ind1] = fTissue[ind2];
	    fTissue[ind2] = tis;
	  }
    }
}

//______________________________________________________________________________
void TVpVoxelArray::RotateY180()
{
  // Rotate the voxel array 180 deg around the Y axis.  No temporary array is
  // generated.  See the class description for an example.
  
  UChar_t tis;
  Int_t ind1, ind2, nx2, nz2;

  nz2 = fNz / 2;
  for (Int_t j = 0; j < fNy; j++)
    for (Int_t i = 0; i < fNx; i++)
      for (Int_t k = 0; k < nz2; k++)
	{
	  ind1 = GetIndex(k, j, i);
	  ind2 = GetIndex(fNz-1-k, j, fNx-1-i);
	  tis = fTissue[ind1];
	  fTissue[ind1] = fTissue[ind2];
	  fTissue[ind2] = tis;
	}

  if ( fNz % 2 != 0)           // fNz is odd
    {
      nx2 = fNx/2;
      for (Int_t j = 0; j < fNy; j++)
	for (Int_t i = 0; i < nx2; i++)
	  {
	    ind1 = GetIndex(nz2, j, i);
	    ind2 = GetIndex(nz2, j, fNx-1-i);
	    tis = fTissue[ind1];
	    fTissue[ind1] = fTissue[ind2];
	    fTissue[ind2] = tis;
	  }
    }
}

//______________________________________________________________________________
void TVpVoxelArray::RotateZ180()
{
  // Rotate the voxel array 180 deg around the Z axis.  No temporary array is
  // generated.  See the class description for an example.
  
  UChar_t tis;
  Int_t ind1, ind2, nx2, ny2;

  ny2 = fNy / 2;
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t i = 0; i < fNx; i++)
      for (Int_t j = 0; j < ny2; j++)
	{
	  ind1 = GetIndex(k, j, i);
	  ind2 = GetIndex(k, fNy-1-j, fNx-1-i);
	  tis = fTissue[ind1];
	  fTissue[ind1] = fTissue[ind2];
	  fTissue[ind2] = tis;
	}

  if ( fNy % 2 != 0)           // fNy is odd
    {
      nx2 = fNx/2;
      for (Int_t k = 0; k < fNz; k++)
	for (Int_t i = 0; i < nx2; i++)
	  {
	    ind1 = GetIndex(k, ny2, i);
	    ind2 = GetIndex(k, ny2, fNx-1-i);
	    tis = fTissue[ind1];
	    fTissue[ind1] = fTissue[ind2];
	    fTissue[ind2] = tis;
	  }
    }
}

//______________________________________________________________________________
void TVpVoxelArray::RotateX90()
{
  // Rotate the voxel array 90 deg around the X axis. (Swap Y and Z
  // coordinates.)  A temporay array is generated.  See the class description
  // for an example.
  
  Int_t oldInd, newInd;
  Int_t newNx = fNx;
  Int_t newNy = fNz;
  Int_t nxyz = fNz * fNy * fNx;

  // Allocate the new array
  UChar_t * tmpTissue = new UChar_t[nxyz];
  if (tmpTissue == 0)
    {
      fprintf(stderr, "Error: TVpVoxelArray::RotateY90: failed to allocate array.\n");
      return;
    }
  
  // Copy the old array to the new array
  for (Int_t i = 0; i < nxyz; i++)
    tmpTissue[i] = fTissue[i];

  // Rotate the array
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t i = 0; i < fNx; i++)
	{
	  oldInd = GetIndex(k, j, i);
	  newInd = j*newNy*newNx + k*newNx + i;
	  fTissue[newInd] = tmpTissue[oldInd];
	}

  // Delete the temporary array
  delete tmpTissue;

  // Swap coordinate - related variables
  swap(fNy, fNz);
  swap(fSizeY, fSizeZ);
  swap(fDpY, fDpZ);
}

//______________________________________________________________________________
void TVpVoxelArray::RotateY90()
{
  // Rotate the voxel array 90 deg around the Y axis. (Swap X and Z
  // coordinates.)  A temporay array is generated.  See the class description
  // for an example.
  
  Int_t oldInd, newInd;
  Int_t newNx = fNz;
  Int_t newNy = fNy;
  Int_t nxyz = fNz * fNy * fNx;

  // Allocate the new array
  UChar_t * tmpTissue = new UChar_t[nxyz];
  if (tmpTissue == 0)
    {
      fprintf(stderr, "Error: TVpVoxelArray::RotateY90: failed to allocate array.\n");
      return;
    }
  
  // Copy the old array to the new array
  for (Int_t i = 0; i < nxyz; i++)
    tmpTissue[i] = fTissue[i];

  // Rotate the array
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t i = 0; i < fNx; i++)
	{
	  oldInd = GetIndex(k, j, i);
	  newInd = i*newNy*newNx + j*newNx + k;
	  fTissue[newInd] = tmpTissue[oldInd];
	}

  // Delete the temporary array
  delete tmpTissue;

  // Swap coordinate - related variables
  swap(fNx, fNz);
  swap(fSizeX, fSizeZ);
  swap(fDpX, fDpZ);
}

//______________________________________________________________________________
void TVpVoxelArray::RotateZ90()
{
  // Rotate the voxel array 90 deg around the Z axis. (Swap X and Y
  // coordinates.)  A temporay array is generated.  See the class description
  // for an example.
  
  Int_t oldInd, newInd;
  Int_t newNx = fNy;
  Int_t newNy = fNx;

  Int_t nxyz = fNz * fNy * fNx;

  // Allocate the new array
  UChar_t * tmpTissue = new UChar_t[nxyz];
  if (tmpTissue == 0)
    {
      fprintf(stderr, "Error: TVpVoxelArray::RotateZ90: failed to allocate array.\n");
      return;
    }
  
  // Copy the old array to the new array
  for (Int_t i = 0; i < nxyz; i++)
    tmpTissue[i] = fTissue[i];

  // Rotate the array
  for (Int_t k = 0; k < fNz; k++)
    for (Int_t j = 0; j < fNy; j++)
      for (Int_t i = 0; i < fNx; i++)
	{
	  oldInd = GetIndex(k, j, i);
	  newInd = k*newNy*newNx + i*newNx + j;
	  fTissue[newInd] = tmpTissue[oldInd];
	}

  // Delete the temporary array
  delete tmpTissue;

  // Swap coordinate - related variables
  swap(fNx, fNy);
  swap(fSizeX, fSizeY);
  swap(fDpX, fDpY);
}


//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVAMFile(const Char_t *fileName)
{
  // Write the VAM file.  See ReadVAMFile() for the format description.
  //
  // Input parameters:
  // - fileName - name of the VAM file
  
  ofstream vam(fileName);
  if (!vam)
    {
      cerr << "TVpVoxelArray::WriteVAMFile: Cannot open file: " << fileName << '\n';
      return 1;
    }
  
  vam << "# Format: VAM 2.0\n";
  vam << "# MMG file: " << fFileNameMMG << '\n';
  vam << "# VAN file: " << fFileNameVAN << '\n';
  vam << "# VAD file: " << fFileNameVAD << '\n';

  if (!vam.good())
    {
      cerr << "Error: TVpVoxelArray::WriteVAMFile: An output error occured.\n";
      return 2;
    }

  return 0;
}

//______________________________________________________________________________
Int_t TVpVoxelArray::WriteVANFile(const Char_t *fileName)
{
  // Write the VAN file.  See ReadVANFile for the format description.
  //
  // Input parameters:
  // - fileName - name of the VAN file
  
  ofstream van(fileName);
  if (!van)
    {
      cerr << "TVpVoxelArray::WriteVANFile: Cannot open file: " << fileName << '\n';
      return 1;
    }
  
  van << "# Format: VAN 2.0\n";
  for (Int_t i = 0; i < fMaxNumOfTissues; i++)
    {
      if (fMaterial[i] != 0)
	van << i << '\t' << fMaterial[i]->GetFileName() << '\t'
	    << fMaterial[i]->GetName() << '\n';
    }

  if (!van.good())
    {
      cerr << "Error: TVpVoxelArray::WriteVANFile: An output error occured.\n";
      return 2;
    }

  return 0;
}


//______________________________________________________________________________
void TVpVoxelArray::Segment(TVectorF &ctNumber, TVectorF &threshold)
{
  // Segment the CT dataset given by the array of CT numbers according to
  // thresholds.

  // Check that dimensions match
  Int_t nT = threshold.GetNoElements();
  if (fNx*fNy*fNz != ctNumber.GetNoElements())
    {
      std::cerr << "Error: TVpVoxelArray::Segment: Dimensions do not match"
		<< ", fNx*fNy*fNz = " << fNx*fNy*fNz
		<< ", ctNumber.GetNoElements() = " << ctNumber.GetNoElements() << std::endl;
      return;
    }

  Int_t index;
  Double_t ctn;
  for (Int_t iz = 0; iz < fNz; iz++)
    for (Int_t iy = 0; iy < fNy; iy++)
      for (Int_t ix = 0; ix < fNx; ix++)
	{
	  index = GetIndex(iz, iy, ix);
	  ctn = ctNumber[index];

	  // Find the material index
	  Int_t im;
	  bool found = 0;
	  for (im = 0; im < nT; im++)
	    if (ctn < threshold[im])
	      {
		found = 1;
		fTissue[index] = (UChar_t) im;
		break;
	      }
	  if (! found)
	    fTissue[index] = (UChar_t) nT;
	}
}

//______________________________________________________________________________
void TVpVoxelArray::Segment(TTree *treeWithCtNumbers, TVectorF &threshold,
			    Int_t verbose)
{
  // Segment the CT dataset in a ROOT's Tree to the array of CT numbers
  // according to thresholds.
  //
  // Input parameters:
  // - treeWithCtNumbers - tree with a ctn branch
  // - threshold - array of CT thresholds

  Int_t nT = threshold.GetNoElements();  // Number of thresholds

  // Get the branch
  CtNumber ctn;
  TBranch *branch = treeWithCtNumbers->GetBranch("ctn");
  branch->SetAddress(&ctn);
  Int_t nEntries = treeWithCtNumbers->GetEntries();

  // Check dimensions match
  if (fNx*fNy*fNz != nEntries)
    {
      std::cerr << "Error: TVpVoxelArray::Segment: Dimensions do not match"
		<< ", fNx*fNy*fNz = " << fNx*fNy*fNz
		<< ", nEntries = " << nEntries << std::endl;
      return;
    }

  Int_t index;       // to index the data member fTissue
  Int_t iEntry = 0;  // to index the Tree

  for (Int_t iz = 0; iz < fNz; iz++)
    for (Int_t iy = 0; iy < fNy; iy++)
      for (Int_t ix = 0; ix < fNx; ix++)
	{
	  index = GetIndex(iz, iy, ix);
	  branch->GetEntry(iEntry++);
	  if (verbose > 0 && iEntry <= 10)
	    std::cerr << ctn.v << ' ';
	  
	  // Find corresponding material index
	  Int_t im;  // material index
	  bool found = 0;
	  for (im = 0; im < nT; im++)
	    if (ctn.v < threshold[im])
	      {
		found = 1;
		fTissue[index] = (UChar_t) im;
		break;
	      }
	  if (! found)
	    fTissue[index] = (UChar_t) nT;
	}
}

//______________________________________________________________________________
void TVpVoxelArray::RebinAndSegment(TTree *treeWithCtNumbers, TVectorF &threshold,
				    Int_t mx, Int_t my, Int_t mz, Int_t method,
				    Int_t verbose)
{
  // Rebin and segment the CT dataset in a ROOT's Tree to the array of CT
  // numbers according to thresholds.
  //
  // Input parameters:
  // - treeWithCtNumbers - tree with a ctn branch
  // - threshold - array of CT thresholds
  // - mx, my, mz - the number of merged bins in x, y, and z direction
  // - method - 0=faster, with temporary array, 1=slower, no temporary array
  // - verbose - verbosity level (default = 0)
  //
  // Method:
  // This routine combines TVpMath::RebinCtArray() and
  // TVpVoxelArray::Segment().  Data in
  // the Tree are read via random access, segmented and stored in the voxel
  // array.
  //
  // Example:
  // const Int_t nT = 16;
  // Float_t T[nT] = {-975, -925, -885, -750, -600, -560, -400, -200, -130,
  //                  -25, 50, 130, 300, 500, 700, 1000};
  // TVectorF threshold(nT, T);
  // TFile file("vct_cd01.root");
  // TTree *tree = (TTree *) file.Get("Tctn");
  // TVpVoxelArray *va = new TVpVoxelArray("va", 0, nxn, nyn, nzn, sizeX, sizeY, sizeZ);
  // va->RebinAndSegment(tree, threshold, mx, my, mz);

  Int_t nT = threshold.GetNoElements();  // Number of thresholds
  Int_t nxo = fNx * mx;  // number of bins in the original array in x direction
  Int_t nyo = fNy * my;  // number of bins in the original array in y direction
  Int_t nzo = fNz * mz;  // number of bins in the original array in z direction
  Int_t nm = mx * my * mz;

  // Get the branch
  CtNumber ctn;
  TBranch *branch = treeWithCtNumbers->GetBranch("ctn");
  branch->SetAddress(&ctn);
  Int_t nEntries = treeWithCtNumbers->GetEntries();

  // Check that dimensions match
  if (nxo*nyo*nzo != nEntries)
    {
      std::cerr << "Error: TVpVoxelArray::RebinAndSegment: Dimensions do not match"
		<< ", nxo*nyo*nzo = " << nxo*nyo*nzo
		<< ", nEntries = " << nEntries
		<< ".  Segmentation aborted.\n";
      return;
    }

  switch (method)
    {
      case 0:  // faster, with temporary arrays
	// Read CT numbers into a vector
	{
	  TVectorF vo(nxo * nyo * nzo);  // original vector
	  std::cerr << "Reading... ";
	  for (Int_t i = 0; i < nEntries; i++)
	    {
	      branch->GetEntry(i);
	      vo[i] = ctn.v;
	    }
	  std::cerr << "Rebinning... ";
	  TVectorF vn;  // rebinned vector
	  Int_t nxn, nyn, nzn;
	  TVpMath::RebinCtArray(vo, nxo, nyo, nzo, mx, my, mz,
				vn, nxn, nyn, nzn);
	  std::cerr << "Segmenting... ";
	  Segment(vn, threshold);
	  std::cerr << "Done.\n";
	}
	break;
	
      case 1:  // slower, no temporary array
	Int_t index;   // to index the data member fTissue
	Int_t iEntry;  // to index the Tree
	Int_t cta;     // average CT number
	
	for (Int_t iz = 0; iz < fNz; iz++)
	  {
	    for (Int_t iy = 0; iy < fNy; iy++)
	      for (Int_t ix = 0; ix < fNx; ix++)
		{
		  index = GetIndex(iz, iy, ix);
		  Double_t sum = 0.0;
		  for (Int_t izb = 0; izb < mz; izb++)  // small box indices
		    for (Int_t iyb = 0; iyb < my; iyb++)
		      for (Int_t ixb = 0; ixb < mx; ixb++)
			{
			  Int_t izo = iz * mz + izb;
			  Int_t iyo = iy * my + iyb;
			  Int_t ixo = ix * mx + ixb;
			  iEntry = TVpMath::Get3dArrayIndex
			    (ixo, iyo, izo, nxo, nyo, nzo);
			  branch->GetEntry(iEntry);
			  sum += ctn.v;
			}
		  cta = Int_t (sum / nm);  // The average CT number
		  
		  // Find corresponding material index
		  Int_t im;  // material index
		  bool found = 0;
		  for (im = 0; im < nT; im++)
		    if (cta < threshold[im])
		      {
			found = 1;
			fTissue[index] = (UChar_t) im;
			break;
		      }
		  if (! found)
		    fTissue[index] = (UChar_t) nT;
		}
	    std::cerr << 'y'; 
	  }
	std::cerr << 'z';
	break;
      default:
	break;
    }
}

//______________________________________________________________________________
Int_t TVpVoxelArray::SetVoxelDensity(Int_t iz, Int_t iy, Int_t ix, Double_t density)
{
  // Set voxel density.  If the fVoxDensity array has not been allocated then
  // allocate it.

  if (fVoxDensity == 0)
    { // Allocate the array
      fVoxDensity = new Float_t[fNx * fNy * fNz];   
      if (fVoxDensity == 0)
	cerr << "Error: TVpVoxelArray::SetVoxelDensity: Allocation of density array ("
	     << fNx * fNy * fNz << " bytes) failed.\n";
    }
  Int_t index = GetIndex(iz, iy, ix);
  if (index < 0)
    return -1;
  fVoxDensity[index] = density;
  return 0;
}
