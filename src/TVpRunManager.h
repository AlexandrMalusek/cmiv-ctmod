#ifndef TVpRunManager_h
#define TVpRunManager_h

#include <stack>
#include <ctime>
#include <fstream>
#include "TVpPhoton.h"
#include "TVpGeometry.h"
#include "AVpSetup.h"
#include "TVpPhoton.h"
#include "pedEvent.h"
#include "fomEvent.h"
#include "misc.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

extern PedEvent pedEvent;

class TVpRunManager
{
 private:
  void            constructorInit();
  
 protected:
  std::ofstream   *fCtmodRunFilePtr;            //! CTmod run file

  inline void  FillPedEvent1(TVpPhoton *photonPtr) const;
  inline void  FillPedEvent2(TVpPhoton *photonPtr) const;
  virtual void FillFomEvent(UInt_t time) const = 0;

 public:
  enum EDebugOption {kWritePed = 1, kKeepHistory = 2, kDrawTrajectories = 4};

  Int_t                 fDebugOptions;          //  Debug options
  std::stack<TVpPhoton> fPhotonStack;           //! Photon stack      
  Long_t                fCurentHistory;         //! Current history
  Long_t                fNumOfHistories;        //! Total number of histories
  AVpSetup             *fSetupPtr;              //! The setup
  Double_t              fEnergyCutoff;          //! Global photon energy cutoff
  Double_t              fWeightCutoff;          //! Global photon weight cutoff
  Double_t              fRussianRouletteWeight; //! Photon weight after Russian roulette
  Int_t                 fSurvivalBiasing;       //! Use survival biasing
  time_t                fStartTime;             //! Start time in Unix seconds
  Int_t                 fProgRepTimeInterval;   //! Progress report time interval
  Long_t                fFomHistoryStep;        //! Number of histories between FOM events

  TTree                *fPedEventTreePtr;       //! Pointer to the PED tree
  TFile                *fPedEventFilePtr;       //! Pointer to the PED file
  TTree                *fFomEventTreePtr;       //! Pointer to the FOM tree
  TFile                *fFomEventFilePtr;       //! Pointer to the FOM file

  TVpRunManager();
  TVpRunManager(AVpSetup *setupPtr);
  virtual ~TVpRunManager();

  inline void       SetDebugOptions(Int_t debug);
  inline Int_t      RussianRoulette(TVpPhoton *photonPtr);
  virtual Int_t     PhotonLife(TVpPhoton *photonPtr) = 0;
  virtual void      Run(Long_t numOfHistories);
  void              SetEnergyCutoff(Double_t weightCutoff);
  void              SetWeightCutoff(Double_t weightCutoff);
  void              SetRussianRouletteWeight(Double_t russianRouletteWeight);
  void              SetSurvivalBiasing(Int_t survivalBiasing);
  void              SetProgressReportTimeInterval(time_t interval);
  time_t            GetProgressReportTimeInterval() const;
  time_t            GetElapsedTime() const;
  void              PrintProgressReport(Int_t option = 0) const;
  void              PrintCtmodRunFile();
  void              SetFomHistoryStep(Long_t fomHistoryStep);

  void              PedEventInitialize(const Char_t *fileName);
  void              PedEventClose();
  void              FomEventInitialize(const Char_t *fileName);
  void              FomEventClose();

  ClassDef(TVpRunManager,1) // Run manager
};

//______________________________________________________________________________
inline void TVpRunManager::SetDebugOptions(Int_t debugOptions)
{
  fDebugOptions = debugOptions;

  // The history must be kept in order to draw trajectories
  if (fDebugOptions & kDrawTrajectories)
    fDebugOptions |= kKeepHistory;
}

//______________________________________________________________________________
inline Int_t TVpRunManager::RussianRoulette(TVpPhoton *photonPtr)
{
  // Play Russian roulette with the photon.  Return 0 if killed, 0 if survived.

  Double_t ws = photonPtr->GetWeight() / fRussianRouletteWeight;
  return (getRand() < ws) ? 1 : 0;
}

//______________________________________________________________________________
void TVpRunManager::FillPedEvent1(TVpPhoton *photonPtr) const
{
  // Pre-fill the PED event.  To avoid many #ifded ... #endif in the code,
  // this routine is called in every case but its body may be empty.

  if (fPedEventFilePtr == 0)
    return;

  TVpVector3 pos = photonPtr->fSolid->PosLocToUni(photonPtr->fPosition);
  TVpVector3 dir = photonPtr->fSolid->DirLocToUni(photonPtr->fDirection);
  pedEvent.c = 0;  // Photon
  pedEvent.i = photonPtr->fLastInteraction;
  pedEvent.x0 = pos.GetX();
  pedEvent.x1 = pos.GetY();
  pedEvent.x2 = pos.GetZ();
  pedEvent.u0 = dir.GetX();
  pedEvent.u1 = dir.GetY();
  pedEvent.u2 = dir.GetZ();
  pedEvent.e = photonPtr->fEnergy;
  pedEvent.w = photonPtr->fWeight;
  pedEvent.b = photonPtr->fSolid->GetIndex();
  pedEvent.v = photonPtr->fSolid->GetSubIndex(photonPtr->GetPos());
  pedEvent.t = fCurentHistory;
  pedEvent.p = 0;
  pedEvent.s = fCurentHistory;
}

//______________________________________________________________________________
void TVpRunManager::FillPedEvent2(TVpPhoton *photonPtr) const
{
  // Fill the PED event.  To avoid many #ifded ... #endif in the code, this
  // routine is called in every case but its body may be empty.

  if (fPedEventFilePtr == 0)
    return;
  pedEvent.d = pedEvent.e - photonPtr->fEnergy;
  fPedEventTreePtr->Fill();
}

#endif


