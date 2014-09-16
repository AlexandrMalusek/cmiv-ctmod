// Figure Of Merit data provide information about convergence speed and
// stability. FOM is defined as
//
// FOM = v/(t m^2)

#ifndef fomEvent_h
#define fomEvent_h

typedef struct {
  Int_t   i;  // point detector index
  Int_t   c;  // point detector channel index (-1 for a single channel)
  UInt_t  h;  // history number
  UInt_t  t;  // wall time in seconds
  Float_t m;  // mean estimate
  Float_t v;  // variance estimate
} FomEvent;

extern FomEvent fomEvent;

#endif // fomEvent_h
