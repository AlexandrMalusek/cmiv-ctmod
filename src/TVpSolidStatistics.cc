//______________________________________________________________________________
//
// TVpSolidStatistics contains data members which are used for dosimetry
// related scoring.  Each solid scores the energy imparted, the number of
// interactions and the number of particles enetering and leaving the solid.
// Coherent scattering, incoherent scattering, photoeffect and track ends are
// scored separately.
//
// Scoring is performed by TVpPhoton::ComptonScattering(),
// TVpPhoton::CoherentScattering(), and TVpPhoton::Photoeffect().

#include "TVpSolidStatistics.h"

ClassImp(TVpSolidStatistics)

//______________________________________________________________________________
TVpSolidStatistics::TVpSolidStatistics()
{
  // Default constructor. Zero all data members.

  fEnergyImpartedPh = fEnergyImpartedIn = fEnergyImpartedTr = 0.0;
  fNumOfParticlesIn = fNumOfParticlesOut = 0.0;
  fNumOfInterPh = fNumOfInterIn = fNumOfInterCo = fNumOfInterTr = 0L;
}

//______________________________________________________________________________
TVpSolidStatistics::~TVpSolidStatistics()
{
  // Destructor
}

//______________________________________________________________________________
void TVpSolidStatistics::PrintStatistics(std::ostream &out, Long_t numOfHist) const
{
  // Print statistics

  Double_t energyImpartedTot = fEnergyImpartedPh + fEnergyImpartedIn + fEnergyImpartedTr;
  Double_t numOfInterTot = fNumOfInterPh + fNumOfInterIn + fNumOfInterCo;
  Double_t numOfAbsorbedParticles = fNumOfParticlesIn - fNumOfParticlesOut;

  out << "  <photoeffect>\n"
      << "    Number of interactions = " << fNumOfInterPh << '\n'
      << "    Energy imparted per history = "
      << fEnergyImpartedPh / numOfHist  << " keV\n"
      << "    Energy imparted per interaction = ";
  if (fNumOfInterPh == 0)
    out << "NA";
  else
    out << fEnergyImpartedPh / fNumOfInterPh;
  out << " keV\n"
      << "  </photoeffect>\n";

  out << "  <incoherent_scattering>\n"
      << "    Number of interactions = " << fNumOfInterIn << '\n'
      << "    Energy imparted per history = "
      << fEnergyImpartedIn / numOfHist << " keV\n"
      << "    Energy imparted per interaction = ";
  if (fNumOfInterIn == 0)
    out << "NA";
  else
    out << fEnergyImpartedIn / fNumOfInterIn;
  out << " keV\n"
      << "  </incoherent_scattering>\n";
    
  out << "  <coherent_scattering>\n"
      << "    Number of interactions = " << fNumOfInterCo << '\n'
      << "  </coherent_scattering>\n";
    
  out << "  <track_ends>\n"
      << "    Number = " << fNumOfInterTr << '\n'
      << "    Energy imparted per history = "
      << fEnergyImpartedTr / numOfHist << " keV\n"
      << "    Energy imparted per track end = ";
  if (fNumOfInterTr == 0)
    out << "NA";
  else
    out << fEnergyImpartedTr / fNumOfInterTr;
  out << " keV\n"
      << "  </track_end>\n";
  
  out << "  <all_interactions>\n"
      << "    Number of interactions = " << numOfInterTot << '\n'
      << "    Energy imparted per history = "
      << energyImpartedTot / numOfHist << " keV\n"
      << "    Energy imparted per interaction = ";
  if (numOfInterTot == 0)
    out << "NA";
  else
    out << energyImpartedTot / numOfInterTot;
  out << " keV\n"
      << "  </all_interactions>\n";
  
  out << "  Number of particles in = " << fNumOfParticlesIn << '\n'
      << "  Number of particles out = " << fNumOfParticlesOut << '\n'
      << "  Number of absorbed particles = " << numOfAbsorbedParticles << '\n';
}

//______________________________________________________________________________
Double_t TVpSolidStatistics::GetEnergyImparted(Long_t numOfHist) const
{
  // Return energy imparted per one history.  Assume numOfHist histories were
  // simulated.

  return (fEnergyImpartedPh + fEnergyImpartedIn + fEnergyImpartedTr) / numOfHist;
}
