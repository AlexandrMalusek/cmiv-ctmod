//______________________________________________________________________________
//
// The TVpSetupTomograph class corresponds to a real CT scanner with a phantom
// inside.  It serves as:
//
// 1. a container referencing individual CT components
// 2. a tool for calculating analytical projections
//
// The source of photons (X-ray tube) is referenced by the fSourcePtr, the
// phantom by the fGeometryPtr, and the detector array by the fDetectorPtr
// pointers.
//
// The position of the source and detector array is described via a single
// parameter called angle here.  The member function SetPosition(angle) uses
// the TVpPitch class to calculate the corresponding rotation matrix and
// translation vector and then changes the source and detector coordinates.
// Angle == 0 correspods to the default position of the source and detector
// array.
//
// Note: SetPosition() modifies the source and detector array coordinates.
// Keep it in mind when these components are shared by several
// TVpSetupTomograph classes.
//
// The geometry describes the phantom only.  In principle, it may also include
// mechanical components of the CT scanner but, in this case, the validity of
// the collision density estimator variance reduction technique must be
// verified since contributions from interactions in the close vicinity of a
// point detector deteriorate convergency speed of the variance reduction
// method.  If only a phantom in vacuum is used then no interactions can occur
// in the vicinity of the point detector array.
//
// Analytical projections use photon spectrum provided by the TVpSpectrum
// class. The resulting energy planar fluence (optionally weighted by the
// detector energy absorption efficiency function) is normalized per one
// source particle.
//
// An example:
//
// // Define individual components
// TVpSource *sourcePtr = ...
// TVpGeometry *geometryPtr = ...
// TVpPointDetectorArray *pdaPtr = ...
// TVpPitch *pitchPtr = ...
//
// // Define the tomograph
// TVpSetupTomograph *tomographPtr = new TVpSetupTomograph(
//    sourcePtr,    // The X-ray source 
//    geometryPtr,  // The phantom geometry
//    pdaPtr,       // The point detector array
//    pitchPtr);    // The movement trajectory
//
// // Set default position of the source and the PDA
// tomographPtr->SetSrcMatVec(matI, TVpVector3(0, 0, srcDist));
// tomographPtr->SetPdaMatVec(matI, TVpVector3(0, 0, srcDist-pdaSrcDist));
//
// // Move the source and the PDA
// tomographPtr->SetPosition(0);
//______________________________________________________________________________

#include <math.h>
#include "TVpSetupTomograph.h"
#include "TVpRunManager.h"

ClassImp(TVpSetupTomograph)

//______________________________________________________________________________
TVpSetupTomograph::TVpSetupTomograph()
{
  // Create a tomograph and initialize all data memebers to zero.

  fSourcePtr = fSource0Ptr = fSource1Ptr = 0;
  fGeometryPtr = 0;
  fDetectorPtr = 0;
  fPitchPtr = 0;
  fDimensionsPtr = 0;
  fPdaType = kPdaUnknown;
  fSourceType = kSourceUnknown;
}

//______________________________________________________________________________
TVpSetupTomograph::TVpSetupTomograph
(TVpSource *sourcePtr, TVpGeometry *geometryPtr, TVpPointDetectorArray *detectorPtr,
 TVpPitch *pitchPtr, TVpDimensionsTomograph *dimensionsPtr)
{
  // Create a tomograph and setup pointers to its components.  Components
  // retain their original positions.
  //
  // Input parameters:
  // - sourcePtr - Pointer to an x-ray source
  // - geometryPtr - Pointer to a geometry defining the phantom
  // - detectorPtr - Pointer to a point detector array
  // - pitchPtr - Pointer to a source and detector trajectory
  // - dimensionsPtr - Pointer to a list of dimensions

  fSourcePtr = fSource0Ptr = sourcePtr;
  fSource1Ptr = 0;
  fGeometryPtr = geometryPtr;
  fDetectorPtr = detectorPtr;
  fPitchPtr = pitchPtr;
  fDimensionsPtr = dimensionsPtr;
  fPdaType = (EPdaType) fDetectorPtr->GetType();
  fSourceType = (ESourceType) fSourcePtr->GetType();
}

//______________________________________________________________________________
TVpSetupTomograph::~TVpSetupTomograph()
{
  // Delete the tomograph.  Do not call destructors of individual
  // components.  The concept is that there may be e.g. two tomographs
  // sharing the same geometry.
}

//______________________________________________________________________________
void TVpSetupTomograph::SetPosition(Double_t rotationAngle, Int_t sliceIndex)
{
  // Set the source and detector position.
  //
  // Input parameters:
  // - angle - rotation angle in rad

  // Transformation of the source
  // Default position
  TVpVector3 srcTraVecDef(0, 0, fDimensionsPtr->GetSad());
  TVpMatrix3x3 srcRotMatDef(1,0,0, 0,1,0, 0,0,1);
  // New position
  TVpVector3 srcTraVecNew;
  TVpMatrix3x3 srcRotMatNew;
  fPitchPtr->GetActiveTransformation
    (rotationAngle, sliceIndex, srcRotMatDef, srcTraVecDef, srcRotMatNew, srcTraVecNew);
  fSource0Ptr->SetActiveTransformation(&srcRotMatNew, &srcTraVecNew);
  if (fSource1Ptr != 0)
    fSource1Ptr->SetActiveTransformation(&srcRotMatNew, &srcTraVecNew);
  
  // Transformation of the PDA, FIX for PDAplanar
  // Default position
  TVpVector3 pdaTraVecDef;
  switch (fPdaType)
    {
      case kPdaPlanar:
	pdaTraVecDef = TVpVector3(0, 0, fDimensionsPtr->GetSad() - fDimensionsPtr->GetSdd());
	break;
      case kPdaCylindrical:
	pdaTraVecDef = TVpVector3(0, 0, fDimensionsPtr->GetSad());
	break;
      case kPdaUnknown:
	std::cerr << "Error: TVpSetupTomograph::SetPosition: Unimplemented for kPdaUnknown.\n";
	break;
    }
  TVpMatrix3x3 pdaRotMatDef(1,0,0, 0,1,0, 0,0,1);
  // New position
  TVpVector3 pdaTraVecNew;
  TVpMatrix3x3 pdaRotMatNew;
  fPitchPtr->GetActiveTransformation
    (rotationAngle, sliceIndex, pdaRotMatDef, pdaTraVecDef, pdaRotMatNew, pdaTraVecNew);
  fDetectorPtr->SetActiveTransformation(&pdaRotMatNew, &pdaTraVecNew);
}

//______________________________________________________________________________
void TVpSetupTomograph::AnalyticProjection(TVpPointDetectorArray *pda, 
					   Int_t numOfPoints,
					   Int_t withPhantom)
{
  // Calculate an analaytic projection.  The scored quantity is defined via
  // TVpPointDetectorArray::fScoredQuantity.  Use SetPosition() to rotate the
  // scanner.
  //
  // Limitation: Particles are emitted from the center of the source only,
  // volume distribution hasn't been implemented yet.
  //
  // Both the X-ray source and the point detector array (PDA) use the Universe
  // Coordinate System (UCS).  It is supposed here that the particle is
  // initially located in the "universe" solid.
  //
  // Input parameters:
  // - pda - Pointer to an alternative PDA which may use higher resolution for
  //   the primary projection.  If not specified then the default PDA referenced
  //   by the TVpSetupTomograph is used.
  // - numOfPoints - A contribution to each point detector is averaged over
  //   numOfPoints points in the source.  This feature is useful only for
  //   non-point sources.
  // - withPhantom - 0 = no phantom, 1 = with phantom

  TVpVector3 endPoint;      // location of the point detector, UCS
  TVpVector3 sourcePoint;   // location of the X-ray source, UCS
  TVpParticle particle;     // default history number is 0
  Double_t opticalPath;     // = \int \mu dx
  Double_t spacePath;       // = \int dx
  
  if (pda == 0)
    pda = fDetectorPtr;

  pda->Zero();
  // Set some useful variables
  TVpSpectrum *spectrum = fSourcePtr->fSpectrumPtr;
  Int_t numOfChannels = spectrum->GetNumOfChannels();
  
  // Calculate contribution to each point detector
  for (Int_t i = 0; i < pda->GetNumOfDetectors(); ++i)
    {
      // Sum over source points
      for (Int_t spIdx = 0; spIdx < numOfPoints; ++spIdx)
	{
	  pda->GetPosition(i, &endPoint);
	  sourcePoint = fSourcePtr->fTraVecL2u;
	  spacePath = norm(endPoint - sourcePoint);
	  
	  // Sum over all energy channels in the spectrum
	  for (Int_t channel = 0; channel < numOfChannels; ++channel)
	    {
	      // Get a particle from the source. Particle's weight
	      // corresponds to the channel intensity. Each spectrum is
	      // normalized so the sum of weights over all channels is 1.
	      fSourcePtr->GetParticleHeadedToPoint(&particle, &endPoint, channel);
	      
	      if (particle.fWeight > 0.0)  // Ignore zero-weight particles
		{
		  // Optical path depends on particle's energy
		  opticalPath = withPhantom ? fGeometryPtr->GetOpticalPath(particle, spacePath)
		    : 0.0;
		  
		  // Note that the source emits 1 photon into 4Pi.
		  particle.fWeight *= exp(-opticalPath) 
		    / (4 * M_PI * spacePath * spacePath * numOfPoints);

		  // Register the virtual particle by the point detector
		  pda->RegisterParticle(i, &particle);
		}
	    }
	}
    }
  pda->UpdateCountersAtEndOfHistory(1);
}

//______________________________________________________________________________
void TVpSetupTomograph::AnalyticSingleScatterProjection(TVpVolumeIntegrator *viPtr,
							TVpPointDetectorArray *pdaPtr,
							Int_t numOfPoints)
{
  // Calculate single scatter projection.
  //
  // Input parameters:
  // - viPtr - volume integrator pointer
  // - pdaPtr - point detector array pointer
  
  TVpVector3 endPointUcs;      // location of the point detector, UCS
  TVpVector3 sourcePointUcs;   // location of the X-ray source, UCS
  TVpParticle particle;        // default history number is 0
  Double_t opticalPath;        // = \int \mu dx
  Double_t spacePath;          // = \int dx
  
  if (pdaPtr == 0)
    pdaPtr = fDetectorPtr;

  pdaPtr->Zero();
  // Set some useful variables
  TVpSpectrum *spectrum = fSourcePtr->fSpectrumPtr;
  Int_t numOfChannels = spectrum->GetNumOfChannels();
  
  // Calculate contribution to each point detector
  for (Int_t i = 0; i < pdaPtr->GetNumOfDetectors(); ++i)
    {
      // Sum over source points
      for (Int_t spIdx = 0; spIdx < numOfPoints; ++spIdx)
	{
	  pdaPtr->GetPosition(i, &endPointUcs);
	  sourcePointUcs = fSourcePtr->fTraVecL2u;
	  spacePath = norm(endPointUcs - sourcePointUcs);
	  
	  // Sum over all energy channels in the spectrum
	  for (Int_t channel = 0; channel < numOfChannels; ++channel)
	    {
	      // Get a particle from the source. Particle's weight
	      // corresponds to the channel intensity. Each spectrum is
	      // normalized so the sum of weights over all channels is 1.
	      fSourcePtr->GetParticleHeadedToPoint(&particle, &endPointUcs, channel);
	      
	      if (particle.fWeight > 0.0)  // Ignore zero-weight particles
		{
		  // Optical path depends on particle's energy
		  opticalPath = fGeometryPtr->GetOpticalPath(particle, spacePath);
		  
		  // Note that the source emits 1 photon into 4Pi.
		  particle.fWeight *= exp(-opticalPath) 
		    / (4 * M_PI * spacePath * spacePath * numOfPoints);

		  // Register the virtual particle by the point detector
		  pdaPtr->RegisterParticle(i, &particle);
		}
	    }
	}
    }
  pdaPtr->UpdateCountersAtEndOfHistory(1);
}

//______________________________________________________________________________
void TVpSetupTomograph::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream (default = cout)

  out << "<TVpSetupTomograph>\n"
      << "$Id: TVpSetupTomograph.cc 62 2009-06-27 10:54:08Z malusek $\n";
  fSourcePtr->PrintStatus(out);
  if (fSource1Ptr != 0)
    {
      out << "Alternative source:\n";
      fSource1Ptr->PrintStatus(out);
    }
  fGeometryPtr->PrintStatus(out);
  fDetectorPtr->PrintStatus(out);
  fPitchPtr->PrintStatus(out);
  out << "</TVpSetupTomograph>\n";
}

//______________________________________________________________________________
void TVpSetupTomograph::ActivateSource(Int_t index)
{
  // Activate source with index "index".  The second source can implement a cone
  // beam whose size is reduced so that it covers the phantom only.  Primary
  // projection must be calculated with the full cone beam but the simulation of
  // photons which do not hit the phantom is a wasting of time.  The secondary
  // source may thus be used for the scatter projection.
  //
  // Input parameters:
  // - index - source index (range = 0, 1)

  switch (index)
    {
    case 0 :
      fSourcePtr = fSource0Ptr;
      break;
    case 1 :
      fSourcePtr = fSource1Ptr;
      break;
    default:
      std::cerr << "Error: TVpSetupTomograph::ActivateSource: Index out of range" << std::endl;
    }
}

//______________________________________________________________________________
void TVpSetupTomograph::SetSource1(TVpSource *source1)
{
  // Set alternative source.  Set the pointer to an alternative X-ray source.
  //
  // Input parameters:
  // - source1 - Pointer to an alternative source

  fSource1Ptr = source1;
}
