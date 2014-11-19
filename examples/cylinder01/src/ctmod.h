#ifndef ctmod01_h
#define ctmod01_h

#include "TVpGeometry.h"
#include "TVpSetupTomograph.h"
#include "TVpMaterialManager.h"

TVpMaterialManager *getMaterialManager();
TVpGeometry *getGeometry(TVpMaterialManager *mmPtr);
TVpSetupTomograph *getTomograph(TVpMaterialManager *mmPtr);

#endif 
