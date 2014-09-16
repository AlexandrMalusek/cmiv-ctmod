//______________________________________________________________________________
//
// TVpRunManagerTomograph defines the tomograph-setup-specific memeber
// function PhotonLife() and performs normalization of scored quantities via
// Run().

#include "TVpRunManagerTomograph.h"
#include "misc.h"
#include "fomEvent.h"

ClassImp(TVpRunManagerTomograph)

//______________________________________________________________________________
TVpRunManagerTomograph::TVpRunManagerTomograph()
  : TVpRunManager()
{
  // Default constructor.  Deprecated function, do not use it.

  fGeometryPtr = 0;
}

//______________________________________________________________________________
TVpRunManagerTomograph::TVpRunManagerTomograph(TVpSetupTomograph *setupTomographPtr)
  : TVpRunManager(setupTomographPtr)
{
  // Constructor.

  fGeometryPtr = setupTomographPtr->fGeometryPtr;
  fSetupTomographPtr = setupTomographPtr;
}

//______________________________________________________________________________
TVpRunManagerTomograph::~TVpRunManagerTomograph()
{
  // A destructor
}

//______________________________________________________________________________
Int_t TVpRunManagerTomograph::PhotonLife(TVpPhoton *photonPtr)
{
  TVpGeometry::EStepResult stepResult;
  Int_t noTrackEnd;
  Double_t t;

  do
    {
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
		  // Register the particle's contribution to the PDA
		  (fSetupTomographPtr->GetPointDetectorArrayPtr())->
		    RegisterVirtualPhoton(fGeometryPtr, photonPtr);
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
		  // Register the particle's contribution to the PDA
		  (fSetupTomographPtr->GetPointDetectorArrayPtr())->
		    RegisterVirtualPhoton(fGeometryPtr, photonPtr);
		  FillPedEvent1(photonPtr);
		  photonPtr->CoherentScattering();
		  FillPedEvent2(photonPtr);
		  break;
		  
		default:
		  std::cerr << "Error: PhotonLife: Unknown interaction type sampled.\n";
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
      while (photonPtr->GetWeight() > fWeightCutoff);
    }
  while (RussianRoulette(photonPtr));
  
  return 0;
}

//______________________________________________________________________________
void TVpRunManagerTomograph::Run(Long_t numOfHistories)
{
  // TVpTomograph specific Run() function

  // Run the generic function
  TVpRunManager::Run(numOfHistories);

  // Update PDA counters
  (fSetupTomographPtr->GetPointDetectorArrayPtr())->
    UpdateCountersAtEndOfHistory(numOfHistories);
}

//______________________________________________________________________________
void TVpRunManagerTomograph::FillFomEvent(UInt_t time) const
{
  // Pre-fill the FOM event.  To avoid many #ifded ... #endif in the code,
  // this routine is called in every case but its body may be empty.

  if (fFomEventFilePtr == 0)
    return;
  
  // Get current time!!!
  fSetupTomographPtr->fDetectorPtr->FillFomEvent(time, fFomEventTreePtr, &fomEvent);
}
