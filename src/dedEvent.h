// Detector Event Data - parameters of interactions (coherent and incoherent
// scatter) which contribute to a specified point detector.  DED events can be
// stored in a ROOT Tree and used to estimate PDFs of scored quantities.

#ifndef dedEvent_h
#define dedEvent_h

typedef struct {
  Int_t   i;    //  Interaction type
  Float_t x0;   //  x[0]
  Float_t x1;   //  x[1]
  Float_t x2;   //  x[2]
  Float_t u0;   //  u[0]
  Float_t u1;   //  u[1]
  Float_t u2;   //  u[2]
  Float_t e;    //  Energy
  Float_t w;    //  Weight
  Float_t c;    //  Contribution
  Float_t a;    //  xi = cos(theta), where theta is the incidence angle
  Int_t s;      //  Solid index
  Int_t t;      //  Solid SubIndex (e.g. tissue number)
  Int_t h;      //  History number
} DedEvent;

#endif // dedEvent_h
