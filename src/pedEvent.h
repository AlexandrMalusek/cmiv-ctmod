// Particle Event Data - parameters of particle interactions (coherent and
// incoherent scatter, photoeffect).  PED events can be stored in a ROOT Tree
// and used to visualize particle tracks, to calculate spatial distribution of
// absorbed dose, ...

#ifndef pedEvent_h
#define pedEvent_h

typedef struct {
  Int_t   c;  //  1, Particle classification index, e.g. 0 photon, 1 electron
  Int_t   i;  //  2, Interaction type, e.g. 0 Compton scattering, 1 photoeffect
  Float_t x0; //  3, x[0]
  Float_t x1; //  4, x[1]
  Float_t x2; //  5, x[2]
  Float_t u0; //  6, u[0]
  Float_t u1; //  7, u[1]
  Float_t u2; //  8, u[2]
  Float_t e;  //  9, Energy
  Float_t w;  // 10, Weight
  Float_t d;  // 11, Energy imparted (= E2 - E1)
  Int_t b;    // 12, Body index
  Int_t v;    // 13, Voxel index (usually 0 for non voxel geometry)
  Int_t t;    // 14, Track number
  Int_t p;    // 15, Parent track number
  Int_t s;    // 16, Shower number
} PedEvent;

#endif // pedEvent_h
