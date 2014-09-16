#include "AVpAbout.h"

ClassImp(AVpAbout)

//______________________________________________________________________________
//
// CTmod is a C++ toolkit for the calculation of projections in computed
// tomography (CT).  It can also be used for the calculation of projections in
// planar radiography and for threshold segmentation of CT phantoms.
//
// Design concepts:
//
// - Containers (e.g. TVpSetupTomograph) contain pointers to objects
// (e.g. TVpSource).


//______________________________________________________________________________
AVpAbout::AVpAbout()
{
  // Default constructor.
}

//______________________________________________________________________________
AVpAbout::~AVpAbout()
{
  // Default destructor.
}

//______________________________________________________________________________
void AVpAbout::Classes()
{
  // ******************************************************************************
  // Run configuration
  // ******************************************************************************
  // TVpSetupDetector
  // TVpSetupTomograph
  // TVpPitch
  // TVpPitchHelical
  // TVpPitchSequential
  //
  // TVpSetupTomograph is a container which references individual CT components
  // and can be used for calculation of analytical projections.
  // TVpSetupDetector has not been fully implemented yet.  Both classes derive
  // from the abstract class AvSetup.
  //
  // ******************************************************************************
  // Run management
  // ******************************************************************************
  // TVpRunManager
  // TVpRunManagerDetector
  // TVpRunManagerTomograph
  //
  // TVpRunManager controls the Monte Carlo run.  TVpRunManagerTomograph adjusts
  // some routines for the CT configuration.  TVpRunManagerDetector has not been
  // fully implemented yet.
  //
  // ******************************************************************************
  // Detector
  // ******************************************************************************
  // TVpPointDetector
  // TVpPointDetectorArray
  // TVpPointDetectorArrayCylindrical
  // TVpPointDetectorArrayPlanar
  // TVpDetectorResponse
  // TVpAntiscatterGrid
  // 
  // TVpPointDetector scores selected quantities at a given point in space.
  // TVpPointDetectorArray arranges point detectors into two dimensional mesh.
  // TVpPointDetectorArrayCylindrical and TVpPointDetectorArrayPlanar inherit
  // from TVpPointDetectorArray and assume that the mesh is cylindrical or
  // planar, respectively.
  //
  // ******************************************************************************
  // Particle source
  // ******************************************************************************
  // TVpSource
  // TVpSourceIso
  // TVpSourceParallelX
  // TVpSourceCyliderFanBeam
  // TVpSourceCylinderFanBeam2
  // TVpSourcePlaneFanBeam
  // TVpSourcePlaneFanBeam2
  // TVpSpectrum
  //
  // ******************************************************************************
  // Geometry
  // ******************************************************************************
  // TVpGeometry
  // TVpSolid
  // TVpSolidNeighbor
  // TVpSolidStatistics
  //
  // ******************************************************************************
  // Solids
  // ******************************************************************************
  // TVpBox
  // TVpCylinder
  // TVpEllipsoid
  // TVpSphere
  // TVpVoxelArray
  //
  // ******************************************************************************
  // Materials
  // ******************************************************************************
  // TVpMaterial
  // TVpMaterialDefaults
  // TVpMaterialFileData
  // TVpMaterialGridData
  // TVpMaterialManager
  //
  // ******************************************************************************
  // Particles
  // ******************************************************************************
  // TVpParticle
  // TVpPhoton
  //
  // ******************************************************************************
  // Mathematics
  // ******************************************************************************
  // TVpMath
  // TVpRandom
  // TVpRandom2
  // TVpRandom3
  // TVpVector3
  // TVpMatrix3x3
  // TVpIntegral
  // TVpVolumeIntegrator
  //
  // ******************************************************************************
  // Image
  // ******************************************************************************
  // TVpSinogram
  // 
  // ******************************************************************************
  // Visualization
  // ******************************************************************************
  // TVpRectangleCut
  //
  // ******************************************************************************
  // Miscelaneous
  // ******************************************************************************
  // TVpConstant
  // TVpPalette
  //
  // ******************************************************************************
  // Testing
  // ******************************************************************************
  // TVpSolidTester
}
