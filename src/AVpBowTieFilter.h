#ifndef AVpBowTieFilter_h
#define AVpBowTieFilter_h

#include "TObject.h"
#include <iosfwd>
#include "TVpVector3.h"

class AVpBowTieFilter
{
 public:
  virtual ~AVpBowTieFilter() {};

  virtual Double_t GetTransmission(Double_t energy, const TVpVector3 &directionL)
    const = 0;
  virtual void     PrintStatus(std::ostream &out = std::cout) const = 0;
};

#endif  // AVpBowTieFilter_h
