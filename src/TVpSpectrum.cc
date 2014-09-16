//______________________________________________________________________________
//
// TVpSpectrum implements an X-ray spectrum which consists of continuous and
// discrete parts. Sum of weights of all channels (continuous and discrete) is
// set to 1.0 when the spectrum is read.
//
// Continuous part consists of GetNumOfContChannels(), channel width is 1.0
// keV, and channel centers are 1.0 keV, 2.0 keV, ... This implementation may
// change in the future.
//
// Discrete part consists of GetNumOfDiscChannels().  Each channel is has an
// energy and weight.
//
// Currently, the rejection sampling and the cumulative distribution function
// are used for the continuous and discrete parts, respectively.  These
// methods are not effective but the time spent in the GetRandEnergy() is
// relatively small:
// TStopwatch swatch; 
// swatch.Start(); TH1F *h = spectrum->GetRandHist(10000000); swatch.Print();
// Real time 0:00:09, CP time 9.710
// i.e. about 10 s for 1e7 samples on a 2.4 GHz P4.
//
// Examples:
//
// // Draw the spectrum
// gSystem->Load("libRCTmod.so");
// TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
// TH1F *h = spectrum->GetHist();
// h->Draw("");
//Begin_Html
/*
<img src="png/spectrum.png">
*/
//End_Html
//
// // Sample energy in MC simulations
// for (Int_t i = 0; i < numOfEnergies; ++i)
//   h->Fill(spectrum->GetRandEnergy())
//
// // Loop in analytical calculations
// for (Int_t channel = 0; channel < spectrum->GetNumOfChannels(); ++channel)
//   {
//     spectrum->GetEnergyAndWeight(channel, energy, weight)
//     // use the energy and weight
//   }

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iomanip>
#include "TVpSpectrum.h"
#include "misc.h"

ClassImp(TVpSpectrum)

//______________________________________________________________________________
TVpSpectrum::TVpSpectrum()
{
  // Default constructor.  Create a spectrum and initialize data members to
  // zero.

  fName = 0;
  fNumOfChannels = fNumOfContChannels = fNumOfDiscChannels = 0;
  fContChannelWeight = fDiscChannelWeight = fDiscChannelEnergy = 0;
}

//______________________________________________________________________________
TVpSpectrum::TVpSpectrum(const Char_t *fileName)
{
  // Create a spectrum from data in the SPE file. Problems are reported to
  // stderr.
  //
  // Input parameters:
  // - fileName - name of the SPE file

  TVpSpectrum::ReadSpeFile(fileName);
}

//______________________________________________________________________________
TVpSpectrum::TVpSpectrum(Double_t energy)
{
  // Create a spectrum with a single energy line.
  //
  // Input parameters:
  // - energy - photon energy in keV

  // Use the energy level as the name of the spectrum
  fName = new Char_t[40];
  sprintf(fName, "Monoenergetic %10.3e keV", energy);
  
  fNumOfContChannels = 0;
  fNumOfChannels = fNumOfDiscChannels = 1;
  

  // Allocate 1 discrete channel and no continuous channels
  fDiscChannelEnergy = new Double_t[1];
  fDiscChannelWeight = new Double_t[1];
  fDiscDist = new Double_t[1];
  fContChannelWeight = 0;
  
  // Set the energy and weight. 
  fDiscChannelEnergy[0] = energy;
  fDiscChannelWeight[0] = 1;

  // Make sure additional data members are properly initialized.
  Normalize();
}

//______________________________________________________________________________
TVpSpectrum::~TVpSpectrum()
{
  // Destructor.

  delete [] fName;
  delete [] fContChannelWeight;
  delete [] fDiscChannelEnergy;
  delete [] fDiscChannelWeight;
}

//______________________________________________________________________________
Int_t TVpSpectrum::ReadSpeFile(const Char_t *fileName)
{ 
  // Read spectrum in SPE format.
  //
  // Input parameters:
  // - fileName - name of the SPE file
  //
  // Limitation:
  // Only 1 keV interval is allowed now.  This may change in the future.
  //
  // Format:
  // # Name: W, 70.0 kV, 15 deg, 2.5 mm Al
  // # Continuos spectrum
  // # Dimension: 70
  // 1.0000e+00 0
  // 2.0000e+00 0
  // ...  
  // 6.9000e+01 90
  // 7.0000e+01 7
  // # Discrete spectrum
  // # Dimension: 4
  // 5.7980e+01 5
  // 5.9320e+01 8
  // 6.7200e+01 2
  // 6.9100e+01 0
  //
  // Energy in the continuous part is given for the bin center and starts at
  // 1.0 keV.  Bin width is 1 keV, energy values are ignored.

  FILE *fp;
  const Int_t stringWidth = 80;         // max name length
  Char_t line[stringWidth+1];           // line
  Int_t strLen;

  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr,
	      "TVpSpectrum::ReadSpeFile: Cannot open the file %s: ",
	      fileName);
      perror("");
      return 1;
    }

  // Read the first line
  fscanf(fp, "# Name: %[^\n]\n", line);  // unsafe
  strLen = strlen(line);
  fName = new Char_t[strLen+1];
  strcpy(fName, line);

  // Skip the second line
  fscanf(fp, "%*[^\n]\n");

  // Read dimension of the continuous part and allocate the array
  fscanf(fp, "# Dimension:%d", &fNumOfContChannels);
  fContChannelWeight = new Double_t[fNumOfContChannels];
  
  // Read the spectrum (energy field is not stored)
  for (Int_t i = 0; i < fNumOfContChannels; i++)
    fscanf(fp, "%*e%le", &fContChannelWeight[i]);
  
  // Skip one line (# Discrete part)
  fscanf(fp, "\n%*[^\n]\n");
  
  // Read dimension of the discrete part and allocate the arrays
  fscanf(fp, "# Dimension:%d", &fNumOfDiscChannels);
  fDiscChannelEnergy = new Double_t[fNumOfDiscChannels];
  fDiscChannelWeight = new Double_t[fNumOfDiscChannels];
  fDiscDist = new Double_t[fNumOfDiscChannels];

  // Read the discrete part
  for (Int_t i = 0; i < fNumOfDiscChannels; i++)
    fscanf(fp, "%le%le", &fDiscChannelEnergy[i], &fDiscChannelWeight[i]);
  
  // Close the file
  fclose(fp);

  // Set the total number of channels
  fNumOfChannels = fNumOfContChannels + fNumOfDiscChannels;

  // Normalize the spectrum
  Normalize();

  return 0;
}

//______________________________________________________________________________
void TVpSpectrum::Normalize()
{
  // Normalize the spectrum so that the sum of weights of all channels
  // (continuous and discrete) equals 1.0.
  
  Double_t sum = 0.0;
  Double_t sumCont = 0.0;
  Double_t sumDisc = 0.0;

  // Continuous part
  for (Int_t i = 0; i < fNumOfContChannels; i++)
    sumCont += fContChannelWeight[i];

  // Discrete part
  for (Int_t i = 0; i < fNumOfDiscChannels; i++)
    sumDisc += fDiscChannelWeight[i];

  sum = sumCont + sumDisc;
  fContFraction = sumCont / sum;

  // Normalize continuous part, also find the maximum weight
  fContMaximum = 0.0;                  // We suppose weight is >= 0 !! 
  for (Int_t i = 0; i < fNumOfContChannels; i++)
    {
      fContChannelWeight[i] /= sum;
      if (fContChannelWeight[i] > fContMaximum)
	fContMaximum = fContChannelWeight[i];
    }

  // Normalize discrete part
  for (Int_t i = 0; i < fNumOfDiscChannels; i++)
    fDiscChannelWeight[i] /= sum;

  // Calculate the distribution function of the discrete distribution
  Double_t discFraction = 1.0 - fContFraction;
  fDiscDist[0] = fDiscChannelWeight[0] / discFraction;
  for (Int_t i = 1; i < fNumOfDiscChannels; i++)
    fDiscDist[i] = fDiscDist[i-1] + fDiscChannelWeight[i] / discFraction;
}

//______________________________________________________________________________
void TVpSpectrum::GetEnergyAndWeight(Int_t channel, Double_t& energy,
				     Double_t& weight) const
{
  // Get energy and weight of a generalized channel.  Generalized channels
  // cover both the continuous and discrete parts.
  //
  // Input parameters:
  // - channel - generalized channel number (range=0,..,GetNumOfChannels()-1)
  //
  // Output parameters:
  // - energy - channel energy in keV
  // - weight - channel weight
  //
  // Example:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe"); 
  // root [] Double_t energy, weight; 
  // root [] spectrum->GetEnergyAndWeight(80, energy, weight);
  // root [] cout << energy << ' ' << weight << endl;
  // 8.10000000e+01 1.16863106e-02

  if (channel < fNumOfContChannels)
    GetContEnergyAndWeight(channel, energy, weight);
  else
    GetDiscEnergyAndWeight(channel-fNumOfContChannels, energy, weight);
}

//______________________________________________________________________________
void TVpSpectrum::GetContEnergyAndWeight(Int_t channel, Double_t& energy,
					 Double_t& weight) const
{
  // Get energy and weight of a continuous channel.  Continuous channels cover
  // the continuous part of the spectrum.
  //
  // Input parameters:
  // - channel - continuous channel number (range=0,..,GetNumOfContChannels()-1)
  //
  // Output parameters:
  // - energy - channel energy in keV
  // - weight - channel weight

#ifndef TVp_NO_RANGE_CHECKING
  if (channel < 0 || channel >= fNumOfContChannels)
    {
      fprintf(stderr, "TVpSpectrum::GetContEnergyAndWeight: channel %d out of range\n",
	      channel);
      return;
    }
#endif

  // i-th channel corresponds to energy (i+0.5, i+1.5).
  energy = channel + 1.0;
  weight = fContChannelWeight[channel];
}


//______________________________________________________________________________
void TVpSpectrum::GetDiscEnergyAndWeight(Int_t channel, Double_t& energy,
					 Double_t& weight) const
{
  // Get energy and weight of a discrete channel.  Discrete channels cover the
  // discrete part of the spectrum.
  //
  // Input parameters:
  // - channel - discrete channel number (range=0,..,GetNumOfDiscChannels()-1)
  //
  // Output parameters:
  // - energy - channel energy in keV
  // - weight - channel weight

#ifndef TVp_NO_RANGE_CHECKING
  if (channel < 0 || channel >= fNumOfDiscChannels)
    {
      fprintf(stderr, "TVpSpectrum::GetDiscEnergyAndWeight: channel %d out of range\n",
	      channel);
      return;
    }
#endif

  energy = fDiscChannelEnergy[channel];
  weight = fDiscChannelWeight[channel];
}

//______________________________________________________________________________
Double_t TVpSpectrum::GetRandEnergy() const
{
  // Return sampled photon's energy from the spectrum.
  //
  // Method:
  // Sample energy from the continuous or discrete part according to the
  // fContFraction.  See GetContRandEnergy() and GetDiscRandEnergy().
  //
  // Example:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");  
  // root [] cout << spectrum->GetRandEnergy() << endl;
  // 5.79869953e+01

  if (getRand() < fContFraction)
    return GetContRandEnergy();
  return GetDiscRandEnergy();
}


//______________________________________________________________________________
Double_t TVpSpectrum::GetContRandEnergy() const
{
  // Return sampled photon's energy from the continuous part of the spectrum.
  //
  // Method:
  // The rejection sampling is used.  The continuous part is approximated by a
  // step function with bin width of 1 keV. This method is not efficient but
  // see the note in class description.

  Double_t energy;

  do
    energy = getRand() * fNumOfContChannels;
  while (GetContWeight(energy) < getRand() * fContMaximum);

  return energy + 0.5;
}

//______________________________________________________________________________
Double_t TVpSpectrum::GetDiscRandEnergy() const
{
  // Return sampled photon's energy from the discrete part of the spectrum.
  //
  // Method:
  // Sampling based on cumulative distribution function with linear search is
  // used.  This method is not efficient but see the note in class
  // description.

  Double_t rnd = getRand();

  for (Int_t i = 0; i < fNumOfDiscChannels; i++)
    if (rnd < fDiscDist[i])
      return fDiscChannelEnergy[i];

  return fDiscChannelEnergy[fNumOfDiscChannels-1]; // shouldn't happen
}

//______________________________________________________________________________
Double_t TVpSpectrum::GetPhi() const
{
  // Return the sum of all weights.  It may be interpreted as the photon
  // fluence or the number of photons.
  //
  // Example:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  // root [] cout << spectrum->GetPhi() << endl;
  // 1

  Double_t energy, weight;
  Double_t sum = 0.0;
  for (Int_t i = 0; i < fNumOfChannels; i++)
    {
      GetEnergyAndWeight(i, energy, weight);
      sum += weight;
    }
  return sum;
}

//______________________________________________________________________________
Double_t TVpSpectrum::GetPsi() const
{
  // Return the sum of weights multiplied by corresponding energies.  It may
  // be interpreted as the energy fluence.  The ratio GetPsi()/GetPhi() may be
  // interpreted as the mean energy of the spectrum.
  //
  // Examples:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  // root [] cout << "Energy fluence = " << spectrum->GetPsi() << " keV/m^2" << endl;
  // Energy fluence = 67.8092 keV/m^2
  // root [] cout << "Mean energy = " << spectrum->GetPsi()/spectrum->GetPhi() << " keV" << endl;
  // Mean energy = 67.8092 keV

  Double_t energy, weight;
  Double_t sum = 0.0;
  for (Int_t i = 0; i < fNumOfChannels; i++)
    {
      GetEnergyAndWeight(i, energy, weight);
      sum += energy * weight;
    }
  return sum;
}

//______________________________________________________________________________
void TVpSpectrum::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream (default = cout)
  //
  // Example:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  // root [2] spectrum->PrintStatus();
  // <TVpSpectrum>
  // Name: 120.0 kV
  // Number of continuous channels: 120
  // 1.00000000e+00 0.00000000e+00
  // 2.00000000e+00 0.00000000e+00
  // ...
  // 1.20000000e+02 2.00108058e-05
  // Number of discrete channels: 4
  // 5.79800000e+01 3.26876513e-02
  // 5.93200000e+01 5.94020771e-02
  // 6.72000000e+01 2.46933344e-02
  // 6.91000000e+01 6.78366318e-03
  // Weight sum = 1.00000000e+00
  //</TVpSpectrum>

  Double_t energy, weight;
  Double_t weightSum = 0;
  
  out << "<TVpSpectrum>\n"
      << "$Id: TVpSpectrum.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Name: " << fName << '\n'
      << "Number of continuous channels: " << fNumOfContChannels << '\n'
      << std::scientific << std::setprecision(8);
  // Continuous channels  
  for (Int_t channel = 0; channel < fNumOfContChannels; ++channel)
    {
      GetContEnergyAndWeight(channel, energy, weight);
      out << energy << ' ' << weight << '\n';
      weightSum += weight;
    }

  // Discrete Channels
  out << "Number of discrete channels: " << fNumOfDiscChannels << '\n';
  for (Int_t channel = 0; channel < fNumOfDiscChannels; ++channel)
    {
      GetDiscEnergyAndWeight(channel, energy, weight);
      out << energy << ' ' << weight << '\n';
      weightSum += weight;
    }
  
  out << "Weight sum = " << weightSum << '\n'
      << "</TVpSpectrum>\n";
}

//______________________________________________________________________________
TH1F *TVpSpectrum::GetHist()
{
  // Return a histogram containing the spectrum.  Discrete part is added to
  // the continuous part.  The number of 1 keV channels is set so that all
  // discrete channels fit in.  It is assumed that discrete channels are
  // ordered.
  //
  // Example:
  // root[] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  // root[] TH1F *h = spectrum->GetHist();
  // root[] h->Draw("");

  Double_t energy, weight;

  // Set the number of histogram bins
  Int_t numOfBins = fNumOfContChannels;
  if (fNumOfDiscChannels > 0)
    {
      Double_t maxDiscEnergy = fDiscChannelEnergy[fNumOfDiscChannels-1];    
      if (maxDiscEnergy > fNumOfContChannels + 0.5)
	numOfBins = Int_t (maxDiscEnergy + 0.5);
    }

  TH1F *h = new TH1F("allParts", fName, numOfBins, 0.5, numOfBins + 0.5);
  // Fill in the spectrum
  for (Int_t i = 0; i < GetNumOfChannels(); i++)
    {
      GetEnergyAndWeight(i, energy, weight);
      h->Fill(energy, weight);
    }

  h->SetXTitle("Energy / keV");
  h->SetYTitle("Weight per channel");

  return h;
}


//______________________________________________________________________________
TH1F *TVpSpectrum::GetRandHist(Int_t nEvent)
{
  // Return a histogram filled with random samples of photon energies from
  // both the discrete and continuous parts of the spectrum.  Use to check
  // that GetRandEnergy() works properly.
  //
  // Input parameters:
  // - nEvent - number of samples
  //
  // Example:
  // root [] TVpSpectrum *spectrum = new TVpSpectrum("spectra/w120_16_Cu05.spe");
  // root [] TH1F *h = spectrum->GetRandHist(10000);
  // root [] h->Draw();

  // Set the number of histogram bins (the same code as in GetHist() is used)
  Int_t numOfBins = fNumOfContChannels;
  if (fNumOfDiscChannels > 0)
    {
      Double_t maxDiscEnergy = fDiscChannelEnergy[fNumOfDiscChannels-1];    
      if (maxDiscEnergy > fNumOfContChannels + 0.5)
	numOfBins = Int_t (maxDiscEnergy + 0.5);
    }

  TH1F *h = new TH1F("randAll", fName, numOfBins, 0.5, numOfBins + 0.5);
  for (Int_t i = 0; i < nEvent; i++)
    h->Fill(GetRandEnergy());
  
  h->SetXTitle("Energy / keV");
  h->SetYTitle("Number of events per channel");

  return h;
}
