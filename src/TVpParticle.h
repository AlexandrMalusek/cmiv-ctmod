#ifndef TVpParticle_h
#define TVpParticle_h

#include "TObject.h"
#include <iosfwd>
//#include <iomanip>
#include <vector>
#include "TVpVector3.h"
#include "TVpSolid.h"

class TVpGeometry;

class TVpParticle
{
 protected:
  Double_t     fSinTheta;    // Sinus of azimuthal scattering angle 
  Double_t     fCosTheta;    // Cosinus of azimuthal scattering angle 
  Double_t     fSinPhi;      // Sinus of polar scattering angle 
  Double_t     fCosPhi;      // Cosinus of polar scattering angle 

 public:
  enum EParticleInteraction {kKilled = 0, kTakenFromStack = 1,  kPhotoefect = 2,
			     kComptonScattering = 3, kCoherentScattering = 4,
			     kPairProduction = 5, kTrackEnd = 6};

  Int_t                fType;            //  Particle type
  Double_t             fEnergy;          //  Particle energy in keV
  Double_t             fWeight;          //  Statistical weight, 1 is the default
  TVpVector3           fPosition;        //  Position in local coordinates, cm
  TVpVector3           fDirection;       //  Direction in local coordinates
  TVpSolid            *fSolid;           //! Pointer to the current solid
  EParticleInteraction fLastInteraction; //  Last interaction
  ULong_t              fHistoryNumber;   //  History number
#ifdef USE_HISTORY
  std::vector<TVpParticle> fHistory;// History of particle (after-interaction vals)
#endif  

  inline TVpParticle();
  inline virtual     ~TVpParticle();
  inline void        Set(Int_t type, Double_t energy, Double_t weight, Double_t x,
		       Double_t y, Double_t z, Double_t u, Double_t v, Double_t w,
		       TVpSolid *solid, ULong_t historyNumber = 0);
  inline void        SetPosition(Double_t x, Double_t y, Double_t z);
  inline void        SetPosition(const TVpVector3& pos);
  inline void        SetDirection(Double_t u, Double_t v, Double_t w);
  inline void        SetDirection(const TVpVector3& dir);
  inline void        SetEnergy(Double_t energy);
  TVpVector3         GetUniversalPosition() const;
  TVpVector3         GetUniversalDirection() const;
  inline void        SetSolid(TVpSolid *solid);
  inline TVpSolid   *GetSolidPtr();
  inline void        SetHistoryNumber(ULong_t historyNumber);
  inline Double_t    GetWeight() const;
  inline void        SetWeight(Double_t weight);
  inline Double_t    GetEnergy() const;
  void               NewDirection();
  void               SetIsotropicSinPhiCosPhi();
  inline Double_t   *GetPos();
  inline Double_t   *GetDir();
  inline TVpVector3  GetLocPos() const;
  inline TVpVector3  GetLocDir() const;

#ifdef USE_HISTORY
  void               AddToHistory();
  void               ClearHistory();
#endif // USE_HISTORY

  friend std::ostream&  operator<<(std::ostream& s, TVpParticle& p);

#ifdef USE_HISTORY
  void               DrawTrajectory();
#endif // USE_HISTORY
  ClassDef(TVpParticle,1) // Kinematic parameters of a particle
};


//______________________________________________________________________________
inline TVpParticle::TVpParticle()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpParticle::~TVpParticle()
{
  // Destructor
}

//______________________________________________________________________________
inline void TVpParticle::Set(Int_t type, Double_t energy, Double_t weight, Double_t x,
	       Double_t y, Double_t z, Double_t u, Double_t v, Double_t w, TVpSolid *solid,
	       ULong_t historyNumber)
{
  fType = type; fEnergy = energy; fWeight = weight;
  fPosition.Set(x, y, z);
  fDirection.Set(u, v, w);
  fSolid = solid;
  fHistoryNumber = historyNumber;
}

//______________________________________________________________________________
inline void TVpParticle::SetPosition(Double_t x, Double_t y, Double_t z)
{
  fPosition.Set(x, y, z);
}

//______________________________________________________________________________
inline void TVpParticle::SetPosition(const TVpVector3& pos)
{
  fPosition = pos;
}

//______________________________________________________________________________
inline void TVpParticle::SetDirection(Double_t u, Double_t v, Double_t w)
{
  fDirection.Set(u, v, w);
}

//______________________________________________________________________________
inline void TVpParticle::SetDirection(const TVpVector3& dir)
{
  fDirection = dir;
}

//______________________________________________________________________________
inline void TVpParticle::SetEnergy(Double_t energy)
{
  fEnergy = energy;
}

//______________________________________________________________________________
inline void TVpParticle::SetSolid(TVpSolid *solid)
{
  fSolid = solid;
}

//______________________________________________________________________________
inline void TVpParticle::SetHistoryNumber(ULong_t historyNumber)
{
  fHistoryNumber = historyNumber;
}

//______________________________________________________________________________
inline Double_t TVpParticle::GetWeight() const
{
  return fWeight;
}

//______________________________________________________________________________
inline void TVpParticle::SetWeight(Double_t weight)
{
  fWeight = weight;
}

//______________________________________________________________________________
inline Double_t TVpParticle::GetEnergy() const
{
  return fEnergy;
}

//______________________________________________________________________________
inline Double_t *TVpParticle::GetPos()
{
  return fPosition.fR;
}

//______________________________________________________________________________
inline Double_t *TVpParticle::GetDir()
{
  return fDirection.fR;
}

//______________________________________________________________________________
inline TVpVector3 TVpParticle::GetLocPos() const
{
  // Return position in the current solid coordinate system.

  return fPosition;
}

//______________________________________________________________________________
inline TVpVector3 TVpParticle::GetLocDir() const
{
  // Return direction in the current solid coordinate system.

  return fDirection;
}

//______________________________________________________________________________
inline TVpSolid *TVpParticle::GetSolidPtr()
{
  // Return the pointer to the current solid.

  return fSolid;
}

#endif
