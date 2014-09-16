#ifndef AVpSetup_h
#define AVpSetup_h

#include "TObject.h"
#include <iosfwd>
#include "TVpGeometry.h"

class AVpSetup
{
 public:
  virtual ~AVpSetup() {};

  virtual void          SetSourceParticle(TVpParticle& particle) = 0;
  virtual void          PrintStatus(std::ostream &out = std::cout) const = 0;
  virtual TVpGeometry  *GetGeometryPtr() const = 0;
};

#endif  // AVpSetup_h
