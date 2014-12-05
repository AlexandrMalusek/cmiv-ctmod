//______________________________________________________________________________
//
// TVpRunManager is a base class which controls the Monte Carlo run.  It
// provides the generic member function Run() and the virtual function
// PhotonLife().  Run() pushes photons on the stack and starts their transport
// via PhotonLife().  PhotonLife() depends on a particular setup and therefore
// must be defined in a derived class, e.g. TVpRunManagerTomograph or
// TVpRunManagerDetector.
//
// Scoring is performed using the function PhotonLife() defined in a derived
// class.  As a consequence, the detector response must be normalized in the
// Run() function of the derived class.
//
// *******************************
// PED (Particle Event Data) files
// *******************************
//
// The PED file structure is defined in pedEvent.h as follows:
//
// typedef struct {
//  Int_t   c;  //  1, Particle classification index, e.g. 0 photon, 1 electron
//  Int_t   i;  //  2, Interaction type
//  Float_t x0; //  3, x[0]
//  Float_t x1; //  4, x[1]
//  Float_t x2; //  5, x[2]
//  Float_t u0; //  6, u[0]
//  Float_t u1; //  7, u[1]
//  Float_t u2; //  8, u[2]
//  Float_t e;  //  9, Energy
//  Float_t w;  // 10, Weight
//  Float_t d;  // 11, Energy imparted (= E2 - E1)
//  Int_t b;    // 12, Body index
//  Int_t v;    // 13, Voxel index (usually 0 for non voxel geometry)
//  Int_t t;    // 14, Track number
//  Int_t p;    // 15, Parent track number
//  Int_t s;    // 16, Shower number, (0, 1, ...)
// } PedEvent;
//
// Interaction types are defined in TVpParticle.h:
// 0 - kKilled                  4 - kCoherentScattering
// 1 - kTakenFromStack          5 - kPairProduction
// 2 - kPhotoefect              6 - kTrackEnd
// 3 - kComptonScattering
//
// To produce the PED file ctmod_ped.root:
// TVpRunManagerTomograph *trmPtr = new TVpRunManagerTomograph(setupTomographPtr);
// trmPtr->PedEventInitialize("ctmod_ped.root");  // open PED file
// trmPtr->Run(10000);                            // run simulation
// trmPtr->PedEventClose();                       // close PED file    
//
// To analyze the file ctmod_ped.root in a ROOT script:
// #include "pedEvent.h"
// PedEvent event;                                 // a structure describing the event
// TFile *file = new TFile("ctmod_ped.root");
// TTree *tree = (TTree *)file->Get("Tped");
// tree->SetBranchAddress("PedEvent", &event);
// for (Int_t i = 0; i < tree->GetEntries(); i++)  // loop over all events
//   {
//     tree->GetEntry(i);                          // set the structure event
//     cout << event.s << endl;                    // access a member
//   }

#include <cmath>
#include <iostream>
#include <csignal>
#include <iomanip>
#include <unistd.h>
#include "TVpRunManager.h"

ClassImp(TVpRunManager)

void alarm_signal_handler (int signum);

PedEvent pedEvent;
FomEvent fomEvent;
struct sigaction new_alarm_signal_action;    // used by the alarm handler
TVpRunManager *gTVpRunManagerPtr;            // used by the alarm handler

//______________________________________________________________________________
TVpRunManager::TVpRunManager()
{
  // Default constructor

  fSetupPtr = 0;
  constructorInit();
}

//______________________________________________________________________________
TVpRunManager::TVpRunManager(AVpSetup *setupPtr)
{
  // Constructor with initialization.
  //
  // Input parameters:
  // - setupPtr - pointer to the setup

  fSetupPtr = setupPtr;
  constructorInit();
}

//______________________________________________________________________________
void TVpRunManager::constructorInit()
{
  // Init function used by constructors.
  
  fDebugOptions = 0;
  fCurentHistory = fNumOfHistories = 0L;
  fEnergyCutoff = 10.0;  // keV
  fWeightCutoff = 0.0;
  fRussianRouletteWeight = 1.0;
  fSurvivalBiasing = 0;
  fProgRepTimeInterval = 60;

  // The alarm signal
  new_alarm_signal_action.sa_handler = alarm_signal_handler;
  sigaction(SIGALRM, &new_alarm_signal_action, NULL);
  alarm(5);  // first report after 5 s then 60 s
  gTVpRunManagerPtr = this;

  // Open the CTmod run file stream
  fCtmodRunFilePtr = new std::ofstream("ctmod_run.xml");
  if (! (*fCtmodRunFilePtr))
    std::cerr << "Error: Cannot open ctmod_run.xml." << std::endl;

  fPedEventTreePtr = 0;
  fPedEventFilePtr = 0;
  fFomEventTreePtr = 0;
  fFomEventFilePtr = 0;
}

//______________________________________________________________________________
TVpRunManager::~TVpRunManager()
{
  // Default destructor


}

//______________________________________________________________________________
void TVpRunManager::Run(Long_t numOfHistories)
{
  // Simulate numOfHistories photon histories.  A photon is generated in the
  // source, pushed on the stack, and its life is simulated.  All secondary
  // photons are processed before a new history starts.  Scored quantities are
  // NOT normalized here.
  //
  // Input parameters:
  // - numOfHistories - number of simulated histories

  TVpPhoton photon;
  fNumOfHistories = numOfHistories;
  time(&fStartTime);
  PrintProgressReport(1);

  for (fCurentHistory = 0; fCurentHistory < numOfHistories; fCurentHistory++)
    {
      // Setup photon's parameters and push it on the stack
      fSetupPtr->SetSourceParticle(photon);
      photon.fLastInteraction = TVpParticle::kTakenFromStack;
      photon.SetHistoryNumber(fCurentHistory);
      fPhotonStack.push(photon);
      
#ifdef USE_HISTORY
      if (fDebugOptions & kKeepHistory)
	photon.AddToHistory();
#endif // USE_HISTORY

      do
	{
	  // Get photon's parameters from the stack
	  photon = fPhotonStack.top();
	  fPhotonStack.pop();
	  FillPedEvent1(&photon);
	  FillPedEvent2(&photon);
	  // Simulate photon's life
	  PhotonLife(&photon);

#ifdef USE_HISTORY	  
	  if (fDebugOptions & kDrawTrajectories)
	    photon.DrawTrajectory();
	  if (fDebugOptions & kKeepHistory)
	    photon.ClearHistory();
#endif // USE_HISTORY
	}
      while (!fPhotonStack.empty());
    }
  PrintProgressReport(2);
  PrintCtmodRunFile();
}

//______________________________________________________________________________
void TVpRunManager::SetEnergyCutoff(Double_t energyCutoff)
{
  // Set the energy cutoff.  If a photon energy is lower than this value then
  // the photon is absorbed locally as a track end.  Typical value is 10 keV.
  //
  // Input parameters:
  // - energyCutoff - energy cutoff in keV

  if (energyCutoff < 0.0)
    {
      std::cerr << "Error: TVpRunManager::SetEnergyCutoff: energyCutoff = "
		<< energyCutoff << " must be >= 0."
		<< "Setting energyCutoff = 10 keV instead.\n";
      energyCutoff = 10.0;
    }
  fEnergyCutoff = energyCutoff;
}

//______________________________________________________________________________
void TVpRunManager::SetWeightCutoff(Double_t weightCutoff)
{
  // Set the weight cutoff.  If photon's weight is lower than this value then
  // Russian roulette is played with the photon.  The value must be in [0,1].
  // If 0 then Russian Roulette is not played (a photon cannot achieve zero
  // weight in a finite number of steps in a MC simulation).  If 1 then
  // Russian roulette is always played.  See also SetRussianRouletteWeight().
  //
  // Input parameters:
  // - weightCutoff - weight cutoff (range = [0,1])

  if (weightCutoff < 0.0 || weightCutoff > 1.0)
    {
      std::cerr << "Error: TVpRunManager::SetWeightCutoff: weightCutoff = "
		<< weightCutoff << " is out of range [0,1]."
		<< "Setting weightCutoff = 0.1 instead.\n";
      weightCutoff = 0.1;
    }
  fWeightCutoff = weightCutoff;
}

//______________________________________________________________________________
void TVpRunManager::SetRussianRouletteWeight(Double_t russianRouletteWeight)
{
  // Set the weight of the photon which survives Russian roulette.  The number
  // must be greater than 0.  The higher the number is the more likely it is
  // the photon will not survive the Russian roulette.
  //
  // Input parameters:
  // - russianRouletteWeight - russian roulette weight

  if (fRussianRouletteWeight < 0.0)
    {
      std::cerr << "Error: TVpRunManager::SetRussianRouletteWeight: russianRouletteWeight = "
		<<  russianRouletteWeight << " must be greater than 0.0."
		<< "Setting russianRouletteWeight = 0.5 instead.\n";
      russianRouletteWeight = 0.5;
    }
  fRussianRouletteWeight = russianRouletteWeight;
}

//______________________________________________________________________________
void TVpRunManager::SetSurvivalBiasing(Int_t survivalBiasing)
{
  // Turn on the survival biasing variance reduction technique.  If 0 then
  // analog simulation is used.  If 1 then survival biasing is used.  Other
  // value results in an error message and no change.
  //
  // Input parameters:
  // - survivalBiasing - 0 = analog simulation, 1 = survival biasing.

  if (survivalBiasing != 0 && survivalBiasing != 1)
    {
      std::cerr << "Error: TVpRunManager::SetSurvivalBiasing: survivalBiasing = "
		<< survivalBiasing << " must be 0 or 1."
		<< "The new value is ignored.\n";
      return;
    }
  fSurvivalBiasing = survivalBiasing;
}

//______________________________________________________________________________
void TVpRunManager::SetFomHistoryStep(Long_t fomHistoryStep)
{
  // Set the figure of merit ... 

  fFomHistoryStep = fomHistoryStep;
}

//______________________________________________________________________________
void TVpRunManager::SetProgressReportTimeInterval(time_t interval)
{
  // Set the progress report time interval in seconds.
  //
  // Input parameters:
  // - interval - time interval in s

  fProgRepTimeInterval = interval;
}

//______________________________________________________________________________
time_t TVpRunManager::GetProgressReportTimeInterval() const
{
  // Get the progress report time interval in seconds.

  return fProgRepTimeInterval;
}

//______________________________________________________________________________
time_t TVpRunManager::GetElapsedTime() const
{
  // Get elapsed time in second.  The time counter is initialized to zero when
  // the constructor is called.
  
  return time(0) - fStartTime;
}

//______________________________________________________________________________
void TVpRunManager::PrintProgressReport(Int_t option) const
{
  // Print progress report.
  //
  // Input parameters:
  // - option - 0 = normal report, 1 = run start, 2 = run end
  //
  // Processing history x of x, x% done. Estimated time to finish: x.

  Int_t pct;
  const Int_t maxTs = 20;
  Char_t stmString[maxTs];  // start time string
  Char_t ctmString[maxTs];  // current time string
  Char_t etmString[maxTs];  // end time string

  // Start time
  struct tm *stm = localtime(&fStartTime);
  strftime(stmString, maxTs, "%Y-%m-%d %H:%M:%S", stm);

  // Current time
  time_t current_time = time(0);
  struct tm *ctm = localtime(&current_time);
  strftime(ctmString, maxTs, "%Y-%m-%d %H:%M:%S", ctm);

  Double_t elapsed_time = difftime(current_time, fStartTime);
  Double_t total_time = (fCurentHistory == 0) ? 0.0 : 
    (fNumOfHistories * elapsed_time) / fCurentHistory;

  // End time
  time_t end_time = (time_t) (fStartTime + total_time);
  struct tm *etm = localtime(&end_time);
  strftime(etmString, maxTs, "%Y-%m-%d %H:%M:%S", etm);

  switch (option)
    {
    case 0:  // run progress
      pct = (Int_t) ((fCurentHistory * 100.0) / fNumOfHistories);
      std::cerr << "<Progress: "
		<< ctmString
		<< " history " << fCurentHistory << " of "
		<< fNumOfHistories << ", " << pct << " % done.\n"
		<< "  Elapsed wall time = "
		<< elapsed_time / 3600.0
		<< " h, Estimated wall time of run = "
		<< total_time / 3600 
		<< " h.\n"
		<< "  Estimated end: "
		<< etmString
		<< ">\n";
      break;

    case 1:  // run start
      std::cerr << "<Progress: "
		<< ctmString
		<< " the run started.>\n";
      break;

    case 2:  // run end
      std::cerr << "<Progress: "
		<< ctmString
		<< " the run finished.>\n"
		<< "<Performance: Histories per wall time = "
		<< fNumOfHistories / total_time
		<< " 1/s.>\n";
      break;
    }
}

//______________________________________________________________________________
void TVpRunManager::PrintCtmodRunFile()
{
  // Print the CTmod run file

  const Int_t maxTs = 20;
  Char_t stmString[maxTs];  // start time string
  Char_t ctmString[maxTs];  // current time string
  Char_t etmString[maxTs];  // end time string

  // Start time
  struct tm *stm = localtime(&fStartTime);
  strftime(stmString, maxTs, "%Y-%m-%d %H:%M:%S", stm);

  // Current time
  time_t current_time = time(0);
  struct tm *ctm = localtime(&current_time);
  strftime(ctmString, maxTs, "%Y-%m-%d %H:%M:%S", ctm);

  Double_t elapsed_time = difftime(current_time, fStartTime);
  Double_t total_time = (fNumOfHistories * elapsed_time) / fCurentHistory;

  // End time
  time_t end_time = (time_t) (fStartTime + total_time);
  struct tm *etm = localtime(&end_time);
  strftime(etmString, maxTs, "%Y-%m-%d %H:%M:%S", etm);

  (*fCtmodRunFilePtr)
    << "<statistics_simulation_time>\n"
    << "  Start: " << stmString << '\n'
    << "  End:   " << etmString << '\n'
    << "  Total elapsed wall time: " << elapsed_time << " s\n"
    << "  Total elapsed wall time: " << elapsed_time / 60 << " m\n"
    << "  Total elapsed wall time: " << elapsed_time / 3600 << " h\n"
    << "  Elapsed wall time per history: "
    << elapsed_time / fNumOfHistories << " s\n"
    << "  Number of histories per wall second: "
    << (Int_t) (fNumOfHistories / elapsed_time) << " 1/s\n"
    << "</statistics_simulation_time>"  << std::endl;
  
  // Print geometry-related statistics
  fSetupPtr->GetGeometryPtr()->PrintStatistics(*fCtmodRunFilePtr, fNumOfHistories);

  fCtmodRunFilePtr->close();
}


//______________________________________________________________________________
void TVpRunManager::PedEventInitialize(const Char_t *fileName)
{
  // Create a ROOT tree.  The branch address is initialized to an external
  // object.
  //
  // Input aparameters:
  // - fileName - name of the PED file in ROOT format
  
  fPedEventFilePtr = new TFile(fileName, "RECREATE");
  fPedEventTreePtr = new TTree("Tped", "PedEvent");
  fPedEventTreePtr->Branch("PedEvent", &pedEvent,
			   "c/I:i:x0/F:x1:x2:u0:u1:u2:e:w:d:b/I:v:t:p:s");
  fPedEventTreePtr->Print();
}

//______________________________________________________________________________
void TVpRunManager::PedEventClose()
{
  // Write and close the file
  fPedEventFilePtr->Write();
  fPedEventTreePtr->Print();
  fPedEventFilePtr->Close();
}

//______________________________________________________________________________
void TVpRunManager::FomEventInitialize(const Char_t *fileName)
{
  // Create a ROOT tree.  The branch address is initialized to an external
  // object.
  
  fFomEventFilePtr = new TFile(fileName, "RECREATE");
  fFomEventTreePtr = new TTree("Tfom", "FOM event");
  fFomEventTreePtr->Branch("FOM event", &fomEvent, "i/I:c:h/i:t:m/F:v");
  fFomEventTreePtr->Print();
}

//______________________________________________________________________________
void TVpRunManager::FomEventClose()
{
  // Write and close the file
  fFomEventFilePtr->Write();
  fFomEventTreePtr->Print();
  fFomEventFilePtr->Close();
}

void alarm_signal_handler (int signum)
{
    sigaction (SIGALRM, &new_alarm_signal_action, NULL);
    alarm(gTVpRunManagerPtr->GetProgressReportTimeInterval());
    gTVpRunManagerPtr->PrintProgressReport();    
}
