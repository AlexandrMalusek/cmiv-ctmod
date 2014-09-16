//______________________________________________________________________________
//
// TVpMath provides mathematical functions which can be used by other classes.
//
// Some interpolation methods do not perform range checking.  It may be
// changed in the future.

#include <cmath>
#include <iostream>
#include "TVpMath.h"

ClassImp(TVpMath)

//______________________________________________________________________________
TVpMath::TVpMath()
{
  // Default constructor.
}

//______________________________________________________________________________
Double_t TVpMath::GetLinLinInterpolation(Double_t x1, Double_t y1,
					 Double_t x2, Double_t y2, Double_t x)
{
  // Lin-Lin interpolation: a linear interpolation in lin-lin scale, x1 < x2.
  //
  // Input parameters:
  // - (x1, y1), (x2, y2) - interpolation points
  // - x - independent variable
  //
  // Output:
  // Returns interpolated value.
  //
  // References:
  // - Derived by A.M., searching for a book reference ...

  Double_t y = y1 + (y2 - y1)/(x2 - x1)*(x - x1);
  return y;
}

//______________________________________________________________________________
Double_t TVpMath::GetLogLogInterpolation(Double_t x1, Double_t y1,
					 Double_t x2, Double_t y2, Double_t x)
{
  // Log-log interpolation: a linear inrepolation in log-log scale, x1 < x2.
  //
  // Input parameters:
  // - (x1, y1), (x2, y2) - interpolation points
  // - x - independent variable
  //
  // Output:
  // Returns interpolated value.
  //
  // References:
  // - Derived by A.M., searching for a book reference ...
  // - ENDF/B-VI recommends log-log interpolation for some cross section data

  // Log scale
  Double_t lx1 = log(x1);  Double_t ly1 = log(y1);
  Double_t lx2 = log(x2);  Double_t ly2 = log(y2);
  Double_t lx = log(x);
  // Linear interpolation
  Double_t ly = (ly2 - ly1)/(lx2 - lx1)*(lx - lx1) + ly1;
  // Back to normal scale
  return exp(ly);
}

//______________________________________________________________________________
Double_t TVpMath::GetLinLinInterpolation(Double_t x, Int_t n, Double_t x0, Double_t dx,
					 Double_t *y)
{
  // Lin-Lin interpolation: a linear interpolation in lin-lin scale,
  // equidistant grid.  If x is out of range then extrapolate.
  //
  // Input parameters:
  // - x - independent variable
  // - n - number of grid points
  // - x0 - starting point of the grid
  // - dx - step between points, dx = x_{i+1} - x_i
  // - y - array of function values
  //
  // Output:
  // Returns an interpolated or extrapolated value.
  //
  // References:
  // - Derived by A.M., searching for a book reference ...

  // Checks
  if (n < 2)
    {
      std::cerr << "Error: TVpMath::GetLinLinInterpolation: n must be at least 2, n = "
		<< n << ". Returning 0.0\n";
      return 0.0;
    }

  // Find i so that x_i <= x < x_{i+1}
  Int_t i = (int) ((x - x0) / dx);
  if (i < 0)
    i = 0;   // Extrapolate, use points (x_0, y[0]) and (x_1, y[1])
  else if (i >= n - 1)
    i = n - 2;  // Extrapolate, use points (x_{n-2}, y[n-2]) and (x_{n-1}, y[n-1])

  Double_t xi = x0 + i * dx;  // corresponds to y[i]
  Double_t val = y[i] + (y[i+1] - y[i]) / dx * (x - xi);
  return val;
}

//______________________________________________________________________________
void TVpMath::SetPoly1(Double_t& p0, Double_t& p1, Double_t x1, Double_t y1,
		       Double_t x2, Double_t y2)
{
  // Set coefficients p0, p1 of the 1-st order polynomial which crosses points
  // (x1,y1) and (x2,y2).
  //
  // Input parameters:
  // - (x1, y1), (x2, y2) - interpolation points
  //
  // Output:
  // - p0, p1 - polynomial coefficients
  //
  // Method:
  // Find the solution of
  // y1 = p0 + p1 * x1
  // y2 = p0 + p1 * x2
  //
  // Reference:
  // - Derived by A.M.

  p1 = (y2 - y1) / (x2 - x1);
  p0 = y1 - p1*x1;
}

//______________________________________________________________________________
Int_t TVpMath:: FindIndexByBinarySearch(Int_t n, Double_t *g, Double_t x)
{
  // Return index of the interval where x lies: if index > 0 is the returned
  // value then x lies in [g[index], g[index+1]).  Intervals are given by the
  // grid points g[i], i=0,...,n-1.  If x is out of range, return -1 (under
  // minimum) or -2 (over maximum).
  //
  // If x == g[i] then return i but if x == g[n-1] then return -2.  This way,
  // it is possible to use both g[i] and g[i+1] in an interpolation routine.
  //
  // -1        0         1                  n         -2
  // -----|---------|-------.........|------------|-------------
  //     g[0]      g[1]            g[n-2]       g[n-1]

  // test range
  if (x < g[0])
    return -1;
  if (x >= g[n-1])
    return -2;

  Int_t mid;  // middle value
  Int_t min = 0;
  Int_t max = n-1;

  do
    {
      mid = (min + max) / 2;
      if (x < g[mid])
	max = mid;
      else
	min = mid;
    }
  while (max - min > 1);
  return min;
}

//______________________________________________________________________________
Int_t TVpMath::FindIndexByLinearSearch(Int_t n, Double_t *g, Double_t x)
{
  // Return index of the interval where x lies: if index is the returned value
  // then x lies in [g[index], g[index+1]). Intervals are given by the grid
  // points g[i], i=0,...,n-1.  If x is out of range, return -1 (under
  // minimum) or -2 (over maximum).  If x == g[i] then return i.
  //
  // -1        0         1                  n         -2
  // -----|---------|-------.........|------------|-------------
  //     g[0]      g[1]            g[n-2]       g[n-1]
  
  if (x < g[0])
    return -1;
  
  // Slow but easier to program now
  for (Int_t i = 1; i < n; i++)
    if (x < g[i])
      return i-1;
  
  return -2;
}


//______________________________________________________________________________
void TVpMath::RebinCtArray(TVectorF &cto, Int_t nxo, Int_t nyo, Int_t nzo,
			   Int_t mx, Int_t my, Int_t mz,
			   TVectorF &ctn, Int_t& nxn, Int_t& nyn, Int_t& nzn)
{
  // Returned a rebinned array.  The vector "cto" is interpreted as a nxo x
  // nyo x nzo array.  A new nxn x nyn x nzn array is created and each element
  // in this array is an average of mx x my x mz numbers of the original
  // array.  The function is used to rebin voxel arrays containing CT numbers
  // because CT numbers (and linear attenuation coefficients) can be averadged
  // this way.
  //
  // Input parameters:
  // - cto - original vector
  // - nxo, nyo, nzo - dimension of the original array
  // - mx, my, mz - number of elements that will be merged in each direction
  //
  // Output parameters:
  // - nxn, nyn, nzn - dimension of the new array
  //
  // Example:

  nxn = nxo / mx;
  nyn = nyo / my;
  nzn = nzo / mz;
  Int_t nm = mx * my * mz;
  ctn.ResizeTo(nxn * nyn * nzn);
  for (Int_t izn = 0; izn < nzn; izn++)
    for (Int_t iyn = 0; iyn < nyn; iyn++)
      for (Int_t ixn = 0; ixn < nxn; ixn++)
	{
	  Int_t indexn = Get3dArrayIndex(ixn, iyn, izn, nxn, nyn, nzn);
	  Double_t sum = 0.0;
	  for (Int_t iz = 0; iz < mz; iz++)
	    for (Int_t iy = 0; iy < my; iy++)
	      for (Int_t ix = 0; ix < mx; ix++)
		{
		  Int_t izo = izn * mz + iz;
		  Int_t iyo = iyn * my + iy;
		  Int_t ixo = ixn * mx + ix;
		  Int_t indexo = Get3dArrayIndex(ixo, iyo, izo, nxo, nyo, nzo);
		  sum += cto[indexo];
		}
	  ctn[indexn] = sum / nm;
	}
}
