#ifndef TVpPhoton_h
#define TVpPhoton_h

#include "TObject.h"
#include "TH1.h"
#include "TVpParticle.h"

class TVpPhoton : public TVpParticle
{
 public:
  inline TVpPhoton();
  inline ~TVpPhoton();

  Int_t                 ComptonScattering(Double_t energyCutoff);
  void                  CoherentScattering();
  void                  Photoeffect();
  EParticleInteraction  Interaction(Int_t survivalBiasing);
  void                  VirtualScatterToDirection(const TVpVector3& dirLocal);

  TH1F *GetHistCoScXi(Int_t nevents, Int_t nbins = 128);
  TH1F *GetHistCoScTheta(Int_t nevents, Int_t nbins = 128);
  TH1F *GetHistCoScPhi(Int_t nevents, Int_t nbins = 128);
  TH1F *GetHistInScTheta(Int_t nevents, Int_t nbins = 128);
  TH1F *GetHistInScPhi(Int_t nevents, Int_t nbins = 128);

  ClassDef(TVpPhoton,1) // Photon interaction sampling routines
};

inline TVpPhoton::TVpPhoton()
{
  // Default constructor
}

inline TVpPhoton::~TVpPhoton()
{
  // Destructor
}

#endif
