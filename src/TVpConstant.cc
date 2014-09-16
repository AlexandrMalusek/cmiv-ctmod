//______________________________________________________________________________
//
// TVpConstant defines values of physical constants and conversion factors
// used in CTmod.  All CTmod quantities use the following internal system of
// units:
//
// Basic units:
// length [cm]
// energy [keV]
// mass [g]
//
// According to the "Guide for the Use of the International System of Unites",
// (B.N. Taylor, NIST Special publication 811, 1995 edition), a value of a
// quantity A can be written as a product of a number and a unit, symbolicaly
// A = {A}[A], where {A} is the numerical value and [A] is the unit.
//
// To convert e.g. a cross section from barns to cm^2, multiply it by the
// corresponding conversion factor, symbolicaly
//
// cross_section / cm^2 = (cross_section / barn) * TVpConstant::fBarn
//
// or
//
// {cross_section}_{cm^2} = {cross_section}_{barn} * TVpConstant::fBarn
//
// In words: the cross section in cm^2 can be obtained from the cross section
// in barns by multiplying it with the factor TVpConstant::fBarn.

#include "TVpConstant.h"

ClassImp(TVpConstant)

// Plank's constant in [keV cm]
Double_t TVpConstant::hc = 1.2398e-7;

// Electron rest energy in [keV]
Double_t TVpConstant::fElectronRestEnergy =  511.0034;

// Classical electron radius in [cm]
Double_t TVpConstant::fClassicalElectronRadius =  2.8179380e-13;

// Classical electron radius devided by 2.0 in [cm]
Double_t TVpConstant::fClassicalElectronRadius22 =
TVpConstant::fClassicalElectronRadius * TVpConstant::fClassicalElectronRadius / 2.0;

// Barn to cm^2 
Double_t TVpConstant::fBarn = 1e-24;

// Micrometer to cm
Double_t TVpConstant::fMicroMeter = 1e-4;
