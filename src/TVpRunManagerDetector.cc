#include "TVpRunManagerDetector.h"
#include "misc.h"

ClassImp(TVpRunManagerDetector)

//______________________________________________________________________________
TVpRunManagerDetector::TVpRunManagerDetector()
{
  fGeometryPtr = 0;
  fPhotonWeight = 0.0;
}

//______________________________________________________________________________
TVpRunManagerDetector::TVpRunManagerDetector(TVpSetupDetector *setupDetectorPtr)
  : TVpRunManager(setupDetectorPtr)
{
  fGeometryPtr = setupDetectorPtr->fGeometry;
  fSetupDetectorPtr = setupDetectorPtr;
  fPhotonWeight = 0.0;
}

//______________________________________________________________________________
TVpRunManagerDetector::~TVpRunManagerDetector()
{
  // A destructor
}

//______________________________________________________________________________
Int_t TVpRunManagerDetector::PhotonLife(TVpPhoton *photonPtr)
{
  TVpGeometry::EStepResult stepResult;
  Int_t noTrackEnd;
  Double_t t;

  do
    {
      // Sample optical free path
      Double_t opticalPath = -log(getRand());
      
      // Calculate corresponding free path
      Double_t path = photonPtr->fSolid->GetPathLength(photonPtr->GetPos(),
						       photonPtr->GetDir(),
						       photonPtr->GetEnergy(),
						       opticalPath);

      // Move the particle
      switch (stepResult= fGeometryPtr->Step(photonPtr, path, t))
	{
	case TVpGeometry::kTerminated:
	  photonPtr->fLastInteraction = TVpParticle::kKilled;
	  FillPedEvent1(photonPtr);
	  FillPedEvent2(photonPtr);
	  return 0;           

	case TVpGeometry::kNoChange:
	  // An interaction will occur
	  switch (photonPtr->Interaction(fSurvivalBiasing))
	    {
	    case TVpParticle::kPhotoefect:
	      FillPedEvent1(photonPtr);
	      photonPtr->Photoeffect();
	      FillPedEvent2(photonPtr); 
	      return 0;

	    case TVpParticle::kComptonScattering:
	      FillPedEvent1(photonPtr);
	      noTrackEnd = photonPtr->ComptonScattering(fEnergyCutoff);
	      FillPedEvent2(photonPtr);
	      if (noTrackEnd == 0)
		{ // Track end, the life is finished.  Produce a virtual interaction.
		  photonPtr->fLastInteraction = TVpParticle::kTrackEnd;
		  FillPedEvent1(photonPtr);
		  photonPtr->fEnergy = 0.0;  // Deposit the energy locally.
		  FillPedEvent2(photonPtr);
		  return 0;
		}	      
	      break;
	      
	    case TVpParticle::kCoherentScattering:
	      FillPedEvent1(photonPtr);
	      photonPtr->CoherentScattering();
	      FillPedEvent2(photonPtr);
	      break;
	      std::cerr << "Error: PhotonLife: Unknown interaction type sampled.\n";	  
	    default:
	      
	      return -1;
	    }
	default:
	  // We ignore other cases
	  break;
	}

#ifdef USE_HISTORY     
      // Store the particle in its history
      if (fDebugOptions & kKeepHistory)
	photonPtr->AddToHistory();
#endif // USE_HISTORY

    }
  while (photonPtr->GetWeight() > fPhotonWeight);

  return 0;
}
