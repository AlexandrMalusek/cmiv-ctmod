//______________________________________________________________________________
//
// TVpPhoton implements functions which simulate photon interactions:
// ComptonScattering(), CoherentScattering(), and Photoeffect().  Interaction
// type is sampled by the function Interaction().

#include <math.h>
#include "TVpPhoton.h"
#include "misc.h"
#include "TVpMaterial.h"
#include "TVpConstant.h"

ClassImp(TVpPhoton)

//______________________________________________________________________________
Int_t TVpPhoton::ComptonScattering(Double_t energyCutoff)
{
  // Simulate incoherent (Compton) scattering.  Based on EGS4 and Persliden.
  // Incoherent scattering function is used if material->GetUseSf() returns a
  // non-zero value.
  //
  // The energy of the scttered photon is sampled and the corresponding
  // scattering (polar) angle is calculated.  Azimuthal angle is sampled
  // isotropicaly.  The new photon direction is then calculated.  Compton
  // profiles are not used - the electron is at rest.
  //
  // If the photon's energy falls below the energyCutoff then fill in the
  // track end statistics and return 0.  Otherwise return 1.

  Double_t alpha1, alpha2, sumalpha;
  Double_t e0, e, k0, k0rme, t;
  Double_t ere = 511.0;           // Do it some other way
  TVpMaterial *material = fSolid->GetMaterial(GetPos());

  k0 = fEnergy;
  k0rme = fEnergy / ere;
  e0 = 1.0 / (1.0 + 2.0 * k0rme);
  alpha1 = -log(e0);
  alpha2 = (1.0 - e0*e0) / 2.0;
  sumalpha = alpha1 + alpha2;

  if (!material->GetUseSf())
    {
      do
	{ 
	  if (alpha1 >= sumalpha * getRand())
	    e = exp(-alpha1 * getRand());
	  else
	    if (k0rme >= (k0rme + 1.0) * getRand())
	      e = e0 +(1.0 - e0) * getRandTriang();
	    else
	      e = e0 + (1.0 - e0) * getRand();
	  
	  t = (1.0 - e) / (e * k0rme);
	  fSinTheta = t * (2.0 - t);	// r-value is equal to fSinTheta**2
	}
      while (e * fSinTheta / (1.0 + e*e) > getRand());
      fSinTheta = sqrt(fSinTheta);
      fCosTheta = 1.0 - t;
    }
  else
    {
      Double_t xMax = fEnergy / TVpConstant::hc; // [1/cm]
      Double_t Smax = material->GetSfg(xMax);  // [1]
      Double_t S;
      do
	{
	  do
	    {
	      if (alpha1 >= sumalpha * getRand())
		e = exp(-alpha1 * getRand());
	      else
		if (k0rme >= (k0rme + 1.0) * getRand())
		  e = e0 +(1.0 - e0) * getRandTriang();
		else
		  e = e0 + (1.0 - e0) * getRand();
	      
	      t = (1.0 - e) / (e * k0rme);
	      fSinTheta = t * (2.0 - t);	// r-value is equal to fSinTheta**2
	    }
	  while (e * fSinTheta / (1.0 + e*e) > getRand());
	  fSinTheta = sqrt(fSinTheta);
	  fCosTheta = 1.0 - t;
	  Double_t sinThetaHalf = sqrt(0.5*(1.0 - fCosTheta));
	  Double_t x = fEnergy * sinThetaHalf / TVpConstant::hc;
	  S = material->GetSfg(x);
	}
      while (getRand() > S / Smax);
    }

  SetIsotropicSinPhiCosPhi();
  fEnergy *= e;
  NewDirection();

  // Score statistics
  fSolid->fEnergyImpartedIn += fWeight * (k0 - fEnergy);  // Score energy imparted
  fSolid->fNumOfInterIn++;

  // Handle track ends
  if (fEnergy <= energyCutoff)
    {
      fSolid->fEnergyImpartedTr += fWeight * fEnergy;  // all energy imparted locally
      fSolid->fNumOfInterTr++;
      return 0;
    }
  return 1;
}

//______________________________________________________________________________
void TVpPhoton::CoherentScattering()
{
  // Simulate coherent (Rayleigh) scattering.  Based on Persliden.  Form
  // factors are used if material->GetUseFf() returns a non-zero value.
  //
  // Scattering (polar) angle is sampled directly if the form factor is not
  // used (from the Thompson cross section) otherwise a rejection technique is
  // used.  Azimuthal angle is sampled isotropically.  New photon direction is
  // calculated and no energy is imparted.

  Double_t q, ash;
  Double_t u, umax, Aumax;
  Double_t fxi;
  TVpMaterial *material = fSolid->GetMaterial(GetPos());

  // Get theta
  if (material->GetUseFf())
    {
      // Use coherent scattering form factors
        umax = fEnergy / TVpConstant::hc;
	umax = umax * umax;
	Aumax = material->GetA(umax);
	do
	  {
	    u = material->GetInvA(getRand() * Aumax);
	    fCosTheta = 1.0 - 2.0 * u / umax;
	    fxi = 0.5 * (1.0 + fCosTheta * fCosTheta);
	  }
	while (getRand() > fxi);
	fSinTheta = sqrt(1.0 - fCosTheta*fCosTheta);
    }
  else
    {
      // Sample according to differential Thompson cross section
      q = 2.0 - 4.0 * getRand();
      ash = log(q + sqrt(1.0 + q*q)) / 3.0;
      fCosTheta = exp(-ash) - exp(ash);
      fSinTheta = sqrt(1.0 - fCosTheta*fCosTheta);
    }

  SetIsotropicSinPhiCosPhi();
  NewDirection();
  // Score statistics
  fSolid->fNumOfInterCo++;
}

//______________________________________________________________________________
void TVpPhoton::Photoeffect()
{
  // Simulate photoefect.  The photon is absorbed and no characteristic
  // radiation is generated.  
  
  // Score statistics
  fSolid->fEnergyImpartedPh += fWeight * fEnergy;  // all energy imparted locally
  fSolid->fNumOfInterPh++;
}

//______________________________________________________________________________
TVpPhoton::EParticleInteraction TVpPhoton::Interaction(Int_t survivalBiasing)
{
  // Sample the interaction type.  Return the sampled value and set the
  // fLastInteraction.
  
  Double_t rnd = getRand();
  TVpMaterial *material = fSolid->GetMaterial(GetPos());
  
  if (survivalBiasing == 0)
    { // analog simulation
      if (rnd < material->GetInIub(fEnergy))
	fLastInteraction = kComptonScattering;
      else if (rnd < material->GetPhIub(fEnergy))
	fLastInteraction = kPhotoefect;
      else
	fLastInteraction = kCoherentScattering;
    }
  else
    { // survival biasing
      Double_t inIub = material->GetInIub(fEnergy);
      Double_t phIub = material->GetPhIub(fEnergy);
      Double_t ps = inIub + 1.0 - phIub;  // probability of scatter
      // Select scattering interaction
      if (rnd*ps < inIub)
	fLastInteraction = kComptonScattering;
      else
	fLastInteraction = kCoherentScattering;
      // Deposit energy corresponding to the absorbed weight
      fSolid->fEnergyImpartedPh += fWeight * fEnergy * (1.0 - ps);  
      // Compensate the weight
      fWeight *= ps;
    }
  return fLastInteraction;
}

//______________________________________________________________________________
void TVpPhoton::VirtualScatterToDirection(const TVpVector3& dirLocal)
{
  // Scatter the particle into a specified angle. Modify its (1)
  // energy, (2) direction, (3) and weight, and return distance to the
  // point. fLastInteraction must cointain the interaction type.

  // Calculate distance and define new direction
  Double_t cosTheta = fDirection | dirLocal;
  fDirection = dirLocal;
  TVpMaterial *materialPtr = fSolid->GetMaterial(GetPos());

  // Modify the particle's weight according to the angular distribution p.d.f.
  Double_t k, rk, pdf;
  switch (fLastInteraction)
    {
    case kComptonScattering:
      k = fEnergy / TVpConstant::fElectronRestEnergy;
      rk = 1.0 / (1.0 + k * (1.0 - cosTheta));
      pdf = materialPtr->GetInPdfOmega(fEnergy, cosTheta);
      fWeight *= pdf;
      fEnergy *= rk;
      break;

    case kCoherentScattering:
      pdf = materialPtr->GetCoPdfOmega(fEnergy, cosTheta);
      fWeight *= pdf;
      break;
      
#if 0
      // will be deleted
      Double_t k, k2, k3, u, u2, pi, sigTot, rk, f1;
    case kComptonScattering:
      k = fEnergy / 511.0;
      k2 = k * k;
      k3 = k2 * k;
      u = (1.0 + 2.0 * k);
      u2 = u * u;
      pi = 3.141592654;
      sigTot = 2*pi * (2*k*(2+8*k+9*k2+k3) + u2*(k2-2-2*k)*log(u)) / (k3*u2);
      rk = 1.0 / (1.0 + k * (1.0 - cosTheta));
      f1 = ((rk * rk) * (rk + 1.0 / rk - 1.0 + cosTheta*cosTheta) ) / sigTot;
      fWeight *= f1;
      fEnergy *= rk;
      break;
      
    case kCoherentScattering:
      // the integral ower Omega gives 2Pi
      fWeight *= 3.0/8.0 * (1.0 + cosTheta*cosTheta);
      break;
#endif

    default:
      std::cerr << "Error: TVpParticle::VirtualScatterToPoint: incorrect interaction type: "
		<< fLastInteraction << std::endl;
    }
}


//______________________________________________________________________________
TH1F *TVpPhoton::GetHistCoScXi(Int_t nevents, Int_t nbins)
{
  // Sample values of xi=cos(Theta), where Theta is the scattering angle, for
  // coherent scattering.  Fill a TH1F histogram with these values and return
  // a pointer to the histogram.  The original photon's parameters are not
  // changed.

  TVpPhoton orgPhoton = *this;  // keep original values

  TH1F *h = new TH1F("CoScXi", "CoScXi", nbins, -1.0, 1.0);  
  for (Int_t i = 0; i < nevents; i++)
    {
      CoherentScattering();
      h->Fill(fCosTheta);
      *this = orgPhoton;        // restore original values
    }
  h->SetXTitle("#xi = cos(#theta) [1]");
  h->SetYTitle("N [1]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpPhoton::GetHistCoScTheta(Int_t nevents, Int_t nbins)
{
  // Sample values of the scattering angle Theta for coherent scattering.
  // Fill a TH1F histogram with these values and return a pointer to the
  // histogram.  The original photon's parameters are not changed.

  Double_t theta;
  TVpPhoton orgPhoton = *this;  // keep original values

  TH1F *h = new TH1F("CoScTheta", "CoScTheta", nbins, 0, M_PI);  
  for (Int_t i = 0; i < nevents; i++)
    {
      CoherentScattering();
      theta = acos(fCosTheta);
      h->Fill(theta);
      *this = orgPhoton;        // restore original values
    }
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("N [1]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpPhoton::GetHistCoScPhi(Int_t nevents, Int_t nbins)
{
  // Sample values of the azimuthal angle Phi for coherent scattering.  Fill a
  // TH1F histogram with these values and return a pointer to the histogram.
  // In theory, the sampled Phi has uniform distribution U(0, 2Pi).  The
  // original photon's parameters are not changed.

  Double_t phi;
  TVpPhoton orgPhoton = *this;  // keep original values

  TH1F *h = new TH1F("CoScPhi", "CoScPhi", nbins, 0, 2*M_PI);  
  for (Int_t i = 0; i < nevents; i++)
    {
      CoherentScattering();
      phi = (fSinPhi >= 0.0)? acos(fCosPhi) : acos(fCosPhi) + M_PI ;
      h->Fill(phi);
      *this = orgPhoton;        // restore original values
    }
  h->SetXTitle("#phi [rad]");
  h->SetYTitle("N [1]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpPhoton::GetHistInScTheta(Int_t nevents, Int_t nbins)
{
  // Sample values of the scattering angle Theta for incoherent scattering.
  // Fill a TH1F histogram with these values and return a pointer to the
  // histogram.  The original photon's parameters are not changed.

  Double_t theta;
  TVpPhoton orgPhoton = *this;  // keep original values

  TH1F *h = new TH1F("InScTheta", "InScTheta", nbins, 0, M_PI);  
  for (Int_t i = 0; i < nevents; i++)
    {
      ComptonScattering(0.0);
      theta = acos(fCosTheta);
      h->Fill(theta);
      *this = orgPhoton;        // restore original values
    }
  h->SetXTitle("#theta [rad]");
  h->SetYTitle("N [1]");
  return h;
}

//______________________________________________________________________________
TH1F *TVpPhoton::GetHistInScPhi(Int_t nevents, Int_t nbins)
{ 
  // Sample values of the azimuthal angle Phi for incoherent scattering.  Fill
  // a TH1F histogram with these values and return a pointer to the histogram.
  // In theory, the sampled Phi has uniform distribution U(0, 2Pi).  The
  // original photon's parameters are not changed.

  Double_t phi;
  TVpPhoton orgPhoton = *this;  // keep original values

  TH1F *h = new TH1F("InScPhi", "InScPhi", nbins, 0, 2*M_PI);  
  for (Int_t i = 0; i < nevents; i++)
    {
      ComptonScattering(0.0);
      phi = (fSinPhi >= 0.0)? acos(fCosPhi) : acos(fCosPhi) + M_PI ;
      h->Fill(phi);
      *this = orgPhoton;         // restore original values
    }
  h->SetXTitle("#phi [rad]");
  h->SetYTitle("N [1]");
  return h;
}
