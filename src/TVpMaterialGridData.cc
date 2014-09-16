//______________________________________________________________________________
//
// The TVpMaterialGridData class provides cross section tables used
// during the MC simulation run. The design prefers fast access to
// data.
//
// Naming convention:
// Csg ... cross sections based on grid data
// Sfg ... scattering function, incoherent scattering, based on grid data
// Ffg ... form factor, coherent scattering, based on grid data
//
// Macroscopic cross section data in an equidistant energy grid:
// InCsg ... incoherent scatttering, [cm2/g]
// PhCsg ... photoeffect, [cm2/g]
// CoCsg ... coherent scattering, [cm2/g]
// 
// Data for interaction type sampling in an equidistant energy grid.
// fInIub ... fInCsg/fToCsg, [1]
// fPhIub ... (fInCsg+fPhCsg)/fToCsg, [1]
// fCoCsg0 ... (fInCsg+fPhCsg+fCoCsg)/fToCsg = 1, [1]
//
// Linear interpolation is based on:
// value = coef0[channel] + coef1[channel] * energy
// e.g.
// inCsg = fInCsg0[channel] + fInCsg1[channel] * energy
//
// Iub, Csg, Nog:
// ************************************************************************
// Lin-lin interpolation is performed, for instance for InIub
//
// inIub = fInIub0[i] +  fInIub1[i] * E,
//
// where i = (E - fEnergyMinIUB) / energyStepIUB.
//
//                         energyStepIUB
// fEnergyMinIUB           <-------->       fEnergyMaxIUB
//       |--------|--------|--------|...--------|
//      E[0]     E[1]     E[2]     E[3]     E[fDimIub]
// index:     0        1         2    fDimIub-1
//
// E[i] = fEnergyMinIUB + i * energyStepIUB.
//
// Sfg:
// ************************************************************************
// The function is nice looking in log-log scale but x[0] is 0.  Thus a mixed
// approach is used: first bin (x[0], x[1]) uses lin-lin interpolation
//
// S = fLnSfg0[0] + fLnSfg1[0] * x,
//
// and the rest (x[1], x[fDimLnSfg]) uses log-log interpolation
//
// S = exp(fLnSfg0[i] + fLnSfg1[i] * log(x)),
//
// where i = 1 + (int) (log(x) - fLnSfgXMin) / fDLnSfgX.  The bin width
// fDLnSfgX = (fLnSfgXMax - fLnSfgXMin)/(fDimLnSfg - 1) is fixed in the log
// scale.
//      
//                     fDLnSfgX
//       fLnSfgXMin   <-------->       fLnSfgXMax
//  |--------|--------|--------|...--------|
// x[0]     x[1]     x[2]     x[3]     x[fDimLnSfg]
// index: 0      1         2    fDimLnSfg-1
//
// x[0] = 0
// x[i] = exp(fLnSfgXMin + (i-1)*fDLnSfgX), 1 <= i <= fDimLnSfg
// 
// Ffg:
// ************************************************************************
// The same concept as for Sfg is applied.  The function is nice looking in
// log-log scale but x[0] is the atomic number.  Thus a mixed approach is
// used: first bin (x[0], x[1]) uses lin-lin interpolation
//
// F = fLnFfg0[0] + fLnFfg1[0] * x,
//
// and the rest (x[1], x[fDimLnFfg]) uses log-log interpolation
//
// F = exp(fLnFfg0[i] + fLnFfg1[i] * log(x)),
//
// where i = 1 + (int) (log(x) - fLnFfgXMin) / fDLnFfgX.  The bin width
// fDLnFfgX = (fLnFfgXMax - fLnFfgXMin)/(fDimLnFfg - 1) is fixed in the log
// scale.
//      
//                     fDLnFfgX
//       fLnFfgXMin   <-------->       fLnFfgXMax
//  |--------|--------|--------|...--------|
// x[0]     x[1]     x[2]     x[3]     x[fDimLnFfg]
// index: 0      1         2    fDimLnFfg-1
//
// x[0] = 0
// x[i] = exp(fLnFfgXMin + (i-1)*fDLnFfgX), 1 <= i <= fDimLnFfg

#include <iostream>
#include <iomanip>
#include "TVpMaterialGridData.h"

ClassImp(TVpMaterialGridData)

//______________________________________________________________________________
TVpMaterialGridData::TVpMaterialGridData()
{
  // Default constructor. Initialize all data members to 0
  
  fDimIub = fDimLnSfg = 0;
  fInIub0 = fInIub1 = fPhIub0 = fPhIub1 = fCoIub0 = fCoIub1 = 0;
  fInCsg0 = fInCsg1 = fPhCsg0 = fPhCsg1 = fCoCsg0 = fCoCsg1 = 0;
  fToCsg0 = fToCsg1 = fEnergyIUB = 0;
  fInNtg0 = fInNtg1 = fCoNtg0 = fCoNtg1 = 0;
  fInNog0 = fInNog1 = fCoNog0 = fCoNog1 = 0;
  fLnSfg0 = fLnSfg1 = 0;
  fLnFfg0 = fLnFfg1 = 0;
  fLnSfgXMin = fLnSfgXMax = fDLnSfgX = 0.0;
  fLnFfgXMin = fLnFfgXMax = fDLnFfgX = 0.0;
}

//______________________________________________________________________________
TVpMaterialGridData::~TVpMaterialGridData()
{
  // Destructor.

  delete [] fEnergyIUB;
  delete [] fToCsg0; delete [] fToCsg1;
  delete [] fPhIub0; delete [] fPhIub1; 
  delete [] fCoIub0; delete [] fCoIub1;
  delete [] fInIub0; delete [] fInIub1;
  delete [] fInCsg0; delete [] fInCsg1;
  delete [] fPhCsg0; delete [] fPhCsg1;
  delete [] fCoCsg0; delete [] fCoCsg1;
  delete [] fInNtg0; delete [] fInNtg1;
  delete [] fCoNtg0; delete [] fCoNtg1;
  delete [] fLnSfg0; delete [] fLnSfg1;
  delete [] fLnFfg0; delete [] fLnFfg1;
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetToCsg(TVpSpectrum *spectrum)
{
  // Return mean MCS.

  return GetLac(spectrum);
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetLac(TVpSpectrum *spectrum, Int_t weighting)
{
  // Return mean linear attenuation coefficient.
  //
  // Input:
  // - spectrum - the spectrum
  // - weighting - 0 = fluence, 1= energy fluence wighting

  Double_t energy, weight;
  Double_t sum = 0;

  // Integrate over the spectrum (both cont. and disc.)
  Int_t maxChannel = spectrum->GetNumOfChannels();
  if (weighting == 0)  // fluence
    {
      for (Int_t channel = 0; channel < maxChannel; channel++)
	{
	  spectrum->GetEnergyAndWeight(channel, energy, weight);
	  sum += weight * GetToCsg(energy);
	}
      sum /= spectrum->GetPhi();  // normalize to unit fluence
    }
  else  // energy fluence
      {
      for (Int_t channel = 0; channel < maxChannel; channel++)
	{
	  spectrum->GetEnergyAndWeight(channel, energy, weight);
	  sum += weight * energy * GetToCsg(energy);
	}
      sum /= spectrum->GetPsi();  // normalize to unit energy fluence
      }  
  
  return sum;
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetSfg(Double_t x) const
{
  // Return the value of the scattering function (incoherent scattering).  Use
  // lin-lin interpolation in the first bin (i = 0) and log-log interpolation
  // in other bins.  Bins i = 1, 2, ... have the same width in log scale.

  // bin 0 uses lin-lin interpolation
  // fLnSfg0[0] is sfg0[0], and fLnSfg1[0] is sfg1[0]; the logarithm is not
  // used
  if (x < 0) // this should never happen
    {
      std::cerr << "Error: TVpMaterialGridData::GetSfg: "
		<< "x < 0 cm^{-1}, "
		<< "x = " << x << " cm^{-1}\n";
      return 0.0;
    }
  if (x == 0.0)   // The special case x = 0
    return fLnSfg0[0];

  Double_t lnX = log(x);
  if (lnX <= fLnSfgXMin)
    {
      Double_t S =  fLnSfg0[0] + fLnSfg1[0] * x;
      return S;
    }

  if (lnX > fLnSfgXMax) // this should never happen
    {
      std::cerr << "Error: TVpMaterialGridData::GetSfg: "
		<< "x out of range, "
		<< "x = " << x << " cm^{-1}, Xmax = " << exp(fLnSfgXMax) << " cm^{-1}\n";
      return 0.0;
    }

  // bins 1, 2, ... use log-log interpolation
  Int_t index = 1 + (int) ((lnX - fLnSfgXMin) / fDLnSfgX);
  Double_t lnS = fLnSfg0[index] + fLnSfg1[index] * lnX;
  Double_t S = exp(lnS);
  return S;
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetFfg(Double_t x) const
{
  // Return the value of the form factor (coherent scattering).  Use lin-lin
  // interpolation in the first bin (i = 0) and log-log interpolation in other
  // bins.  Bins i = 1, 2, ... have the same width in log scale.

  Double_t lnX = log(x);

  // bin 0 uses lin-lin interpolation
  // fLnFfg0[0] is sfg0[0], and fLnFfg1[0] is sfg1[0]; the logarithm is not
  // used
  if (lnX <= fLnFfgXMin)
    {
      Double_t F =  fLnFfg0[0] + fLnFfg1[0] * x;
      return F;
    }

  if (lnX > fLnFfgXMax) // this should never happen
    {
      std::cerr << "Error: TVpMaterialGridData::GetFfg: "
		<< "x out of range, "
		<< "x = " << x << " cm^{-1}, Xmax = " << exp(fLnFfgXMax) << " cm^{-1}\n";
      return 0.0;
    }

  // bins 1, 2, ... use log-log interpolation
  Int_t index = 1 + (int) ((lnX - fLnFfgXMin) / fDLnFfgX);
  Double_t lnF = fLnFfg0[index] + fLnFfg1[index] * lnX;
  Double_t F = exp(lnF);
  return F;
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetFfgValue(Int_t index) const
{
  // Return

  if (index == 0)
    return fLnFfg0[0];
  Double_t logxi = fLnFfgXMin + (index-1)*fDLnFfgX;
  Double_t F = exp(fLnFfg0[index] + fLnFfg1[index] * logxi);
  return F;
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetFfgX(Int_t index) const
{
  // Return
  
  if (index == 0)
    return 0.0;
  Double_t logxi = fLnFfgXMin + (index-1)*fDLnFfgX;
  return exp(logxi);
}

//______________________________________________________________________________
Double_t TVpMaterialGridData::GetGridDataValue(EGridGraphType graphType, Double_t energy)
{
  switch (graphType)
    {
    case kCoCsg:
      return GetCoCsg(energy);
    case kInCsg:
      return GetInCsg(energy);
    case kPhCsg:
      return GetPhCsg(energy);
    case kToCsg:
      return GetToCsg(energy);
    case kCoIub:
      return GetCoIub(energy);
    case kInIub:
      return GetInIub(energy);
    case kPhIub:
      return GetPhIub(energy);
    default:
      std::cerr << "Error:  TVpMaterialGridData::GetGridDataValue\n";
    }
  return 0.0;
}


//______________________________________________________________________________
TH1F *TVpMaterialGridData::GetHistGridData(EGridGraphType graphType)
{
  // Return a histogram of grid data values.
  // 
  // Csg and Iub histograms: The number of histogram bins is equal to the
  // fDimIub, fix bin size is used.  The value is calculated in the center of
  // each bin.
  //
  // Sfg histogram: The number of histogram bins is equal to fDimLnSfg.  The
  // grid is equidistant in log scale, therefore variable bin width is used.
  // The value is evaluated in the centre of the log-scale bin.

  const Char_t *labelX, *labelY, *title;
  TH1F *h;

  switch (graphType)
    {
    case kCoIub:
      title = "CoIub";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{In}+#Sigma_{Ph}+#Sigma_{Co}/#Sigma_{t} [1]";
      break;
    case kInIub:
      title = "InIub";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{In}/#Sigma_{t} [1]";      
      break;
    case kPhIub:
      title = "PhIub";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{In}+#Sigma_{Ph}/#Sigma_{t} [1]";      
      break;
    case kCoCsg:
      title = "CoCsg";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{Co} [cm^{2}/g]";      
      break;
    case kInCsg:
      title = "InCsg";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{In} [cm^{2}/g]";      
      break;
    case kPhCsg:
      title = "PhCsg";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{Ph} [cm^{2}/g]";      
      break;
    case kToCsg:
      title = "ToCsg";
      labelX = "Energy [keV]";
      labelY = "#Sigma_{t} [cm^{-1}]";      
      break;
    case kSfg:
      title = "Sfg";
      labelX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      labelY = "S [1]";
      break;
    case kFfg:
      title = "Ffg";
      labelX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      labelY = "F [1]";
      break;
    default:
      std::cerr << "Error: TVpMaterialGridData::GetHistGridData"
		<< "Invalid graph type\n"; 
      return 0;
    }

  if (graphType == kSfg)
    { // *** Sfg ***
      if (fDimLnSfg == 0)
	{
	  std::cerr << "Error: TVpMaterialGridData::GetHistGridData: "
		    << "fDimLnSfg == 0\n";
	  return 0;
	}
      // Bin edges.  The number of points is the number of intervals plus 1.
      Double_t *x = new Double_t[fDimLnSfg+1];
      x[0] = 0;
      for (Int_t i = 1; i <= fDimLnSfg; i++)
	x[i] = exp(fLnSfgXMin + (i-1) * fDLnSfgX );
      h = new TH1F(title, title, fDimLnSfg, x);
      Double_t xm;  // x in the middle
      Double_t S;   // Sfg value
      // Bin 0
      xm = exp(0.5 * fLnSfgXMin);
      S = GetSfg(xm);
      h->SetBinContent(1, S);
      // Other bins
      for (Int_t i = 1; i < fDimLnSfg; i++)
	{
	  xm = exp(fLnSfgXMin + (i-0.5) * fDLnSfgX); // ((i-1) + (i))/2 = i-0.5
	  S = GetSfg(xm);
	  h->SetBinContent(i+1, S);
	}
      delete [] x;
    }
  else if (graphType == kFfg)
    { // *** Ffg ***, the algorithm is the same as for Sfg only variable names
      // are changed
      if (fDimLnFfg == 0)
	{
	  std::cerr << "Error: TVpMaterialGridData::GetHistGridData: "
		    << "fDimLnFfg == 0\n";
	  return 0;
	}
      // Bin edges.  The number of points is the number of intervals plus 1.
      Double_t *x = new Double_t[fDimLnFfg+1];
      x[0] = 0;
      for (Int_t i = 1; i <= fDimLnFfg; i++)
	x[i] = exp(fLnFfgXMin + (i-1) * fDLnFfgX );
      h = new TH1F(title, title, fDimLnFfg, x);
      Double_t xm;  // x in the middle
      Double_t S;   // Ffg value
      // Bin 0
      xm = exp(0.5 * fLnFfgXMin);
      S = GetFfg(xm);
      h->SetBinContent(1, S);
      // Other bins
      for (Int_t i = 1; i < fDimLnFfg; i++)
	{
	  xm = exp(fLnFfgXMin + (i-0.5) * fDLnFfgX); // ((i-1) + (i))/2 = i-0.5
	  S = GetFfg(xm);
	  h->SetBinContent(i+1, S);
	}
      delete [] x;
    }
  else 
    { // *** Csg and Iub ***
      if (fDimIub == 0)
	{
	  std::cerr << "Error: TVpMaterialGridData::GetHistGridData: "
		    << "fDimIub == 0\n";
	  return 0;
	}
      h = new TH1F(title, title, fDimIub, fEnergyMinIUB, fEnergyMaxIUB);
      Double_t de = (fEnergyMaxIUB - fEnergyMinIUB) / fDimIub;
      Double_t energy;
      for (Int_t i = 0; i < fDimIub; i++)
	{
	  energy = fEnergyMinIUB + (i + 0.5) * de;
	  h->SetBinContent(i+1, GetGridDataValue(graphType, energy));
	}
    }
  h->SetXTitle(labelX);
  h->SetYTitle(labelY);
  h->SetTitle(title);
  return h;
}

//______________________________________________________________________________
TGraph *TVpMaterialGridData::GetGraphGridData(EGridGraphType graphType) const
{
  // Return a graph.  Skip the point (x[0] = 0, y[0]) which causes
  // problems in log-log plots for fCoFff and fInSff.
  //
  // Limitation:  implemented only for kSfg and kFfg
  
  Double_t *x, *F, *S;

  switch (graphType)
    {
    case kFfg:
      // x[0] == 0 is skipped and thus i -> i+1
      x = new Double_t[fDimLnFfg];
      F = new Double_t[fDimLnFfg];
      for (Int_t i = 0; i < fDimLnFfg; i++)
	{
	  x[i] = GetFfgX(i+1);
	  F[i] = GetFfgValue(i+1);
	}
      return new TGraph(fDimLnFfg, x, F);
    case kSfg:
      // x[0] == 0 is skipped and thus i-1 -> i
      x = new Double_t[fDimLnSfg];
      S = new Double_t[fDimLnSfg];
      for (Int_t i = 0; i < fDimLnSfg; i++)
	{
	  Double_t lnx = fLnSfgXMin + i*fDLnSfgX;
	  x[i] = exp(lnx);
	  S[i] = exp(fLnSfg0[i] + fLnSfg1[i] * lnx);
	}
      return new TGraph(fDimLnSfg, x, S);
    default:
      std::cerr << "Error: TVpMaterialGridData::GetGraphGridData: "
		<< "incorrect graphType\n";
    }
  return 0;
}

//______________________________________________________________________________
void TVpMaterialGridData::DrawGraphGridData(EGridGraphType graphType, Char_t *opt) const
{
  // Draw graph

  const Char_t *titleX, *titleY, *title;

  switch (graphType)
    {
    case kFfg:
      title = "Ffg";
      titleX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      titleY = "F [1]";
      break;
    case kSfg:
      title = "Sfg";
      titleX = "x = sin(#theta/2)/#lambda [cm^{-1}]";
      titleY = "S [1]";
      break;
        default:
      std::cerr << "Error: TVpMaterialGridData::GetGraphGridData: "
		<< "incorrect graphType\n";
      break;
    }
  TGraph *graph = GetGraphGridData(graphType);
  graph->Draw(opt);
  graph->GetXaxis()->SetTitle(titleX);
  graph->GetYaxis()->SetTitle(titleY);
  graph->SetTitle(title);
  graph->Draw(opt);
}
