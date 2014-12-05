//______________________________________________________________________________
//
// TVpRectangleCut defines a fNx x fNy pavement of a sizeX x sizeY
// rectangle in a geometry. The cut can be used to get maps (2D
// histograms) of solid indicies or linear attenuation coefficients.
// Values are determined in the center of each tile or if ndivision>1
// the average value of ndivision x ndivision points is returned.
//
// The pavement:
//   ----------------  sizeY/2
//   |  |  |  |  |  |
//   ----------------
//   |  |  |  |  |  |
//   ----------------
//   |  |  |  |  |  |
//   ---------------- -sizeY/2
//  -sizeX/2  sizeX/2
//
// The example above defines the 5 x 3 pavement of the rectangle.
//
// One tile:
// -----------
// | x  x  x |
// |         |
// | x  x  x |
// |         |
// | x  x  x |
// -----------
//
// Points inside one tile are located so that they create a regular
// mesh when more tiles are joined together. In the example above,
// ndivisions = 3 and the average value of the linear attenuation
// coefficient is calculated from 9 points.
//
// Notes:
//
// 1. ndivisions = 1 for solid indicies. Larger value does not make
// sence.
//
// 2. Linear attenuation coefficient can be calculated for a single
// energy or a spectrum. No beam attenuation inside the geometry is
// considered so values obtained for a spectrum are not realistic.

#include "TVpRectangleCut.h"
#include "TPolyLine3D.h"

ClassImp(TVpRectangleCut)

//______________________________________________________________________________
TVpRectangleCut::TVpRectangleCut
(TVpVector3 center, TVpVector3 baseX, TVpVector3 baseY,
 Double_t sizeX, Double_t sizeY, Int_t nx, Int_t ny)
{
  // RectangleCut constructor. The pavement is specified by its
  // center, two orthonormal base vectors, the size of the rectangle
  // and the number of tiles in each direction.

  fCenter = center;
  fBaseX = (sizeX/nx) * normalize(baseX);
  fBaseY = (sizeY/ny) * normalize(baseY);
  fSizeX = sizeX;
  fSizeY = sizeY;
  fNx = nx;
  fNy = ny;
  fOrigin = fCenter - ((nx-1)/2.0) * fBaseX - ((ny-1)/2.0) * fBaseY;
}

//______________________________________________________________________________
TVpRectangleCut::~TVpRectangleCut()
{
  // Destructor
}

//______________________________________________________________________________
TH2S *TVpRectangleCut::GetIndex(TVpGeometry *geometry)
{
  // Return a 2D histogram containing the solid index in the central
  // point of each tile. Index for points outside the universe is set
  // to -1.

  Int_t index;
  TVpSolid *solid;
  TVpVector3 position;
  Char_t title[100];
  TVpVector3 u = normalize(fBaseX * fBaseY);
  
  sprintf(title, "Indices, r=(%.3e,%.3e,%.3e) u=(%.3e,%.3e,%.3e)",
	  fCenter.GetX(), fCenter.GetY(), fCenter.GetZ(),
	  u.GetX(), u.GetY(), u.GetZ());
  TH2S *h = new TH2S("index", title, fNx, -fSizeX/2.0, fSizeX/2.0,
		     fNy, -fSizeY/2.0, fSizeY/2.0);
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      {
	position = fOrigin + i * fBaseX + j * fBaseY;
	solid = geometry->GetSolid(position);
	index = (solid == 0) ? -1 : solid->GetIndex();
	h->SetCellContent(i+1, j+1, index);
      }

  return h;
}

//______________________________________________________________________________
TH2F *TVpRectangleCut::GetLac(TVpGeometry *geometry, Double_t energy, Int_t ndivision)
{
  // Display Lac
  
  return GetMCSSum(geometry, energy, 0, ndivision);
}

//______________________________________________________________________________
TH2F *TVpRectangleCut::GetMCSSum(TVpGeometry *geometry, Double_t energy,
				 Int_t csSum, Int_t ndivision)
{
  // Return 2D histogram of the Macroscopic Cross Section.
  // geometry   ... pointer to the geometry
  // energy     ... energy [keV]
  // csSum      ... include only specified cross sections, see later
  // ndivisions ... number of averaging points in x and y directions
  //
  // Cross section types:
  // kAllCS = 0 ... all cross sections
  // kInCS = 1 ... incoherent scattering
  // kPhCS = 2 ... photoeffect
  // kCoCS = 4 ... coherent scattering
  // Example: kInCS + kCoCS = 5, incoherent and coherent scattering is included

  TVpVector3 position;
  Char_t title[100];
  TVpVector3 u = normalize(fBaseX * fBaseY);
  Int_t npoints = ndivision * ndivision;
  Double_t a1 = 1.0 / (2.0*ndivision);
  Double_t a2 = 1.0 / (Double_t) ndivision;
  
  sprintf(title, "LAC, r=(%.3e,%.3e,%.3e) u=(%.3e,%.3e,%.3e)",
	  fCenter.GetX(), fCenter.GetY(), fCenter.GetZ(),
	  u.GetX(), u.GetY(), u.GetZ());
  TH2F *h = new TH2F("lac", title, fNx, -fSizeX/2.0, fSizeX/2.0,
		     fNy, -fSizeY/2.0, fSizeY/2.0);
  for (Int_t i = 0; i < fNx; i++)
    for (Int_t j = 0; j < fNy; j++)
      {
	Double_t sum = 0.0;
	for (Int_t id = 0; id < ndivision; id++)
	  for (Int_t jd = 0; jd < ndivision; jd++)
	    {
	      position = fOrigin + (a1+i+a2*id) * fBaseX + (a1+j+a2*jd) * fBaseY;
	      TVpSolid *solid = geometry->GetSolid(position);
	      sum += (csSum == 0) ? solid->GetMaterial(position.Get())->GetLac(energy):
		solid->GetMaterial(position.Get())->GetMCSSum(energy, csSum);
	    }
	h->SetCellContent(i+1, j+1, sum / npoints);
      }

  return h;
}

//______________________________________________________________________________
TH2F *TVpRectangleCut::RegisterSourceParticle(TVpSource *source, Long_t nparticles)
{
  // Register

  Char_t title[100];
  TVpVector3 u = normalize(fBaseX * fBaseY);
  TVpVector3 b1 = normalize(fBaseX);
  TVpVector3 b2 = normalize(fBaseY);

  sprintf(title, "cross, r=(%.3e,%.3e,%.3e) u=(%.3e,%.3e,%.3e)",
	  fCenter.GetX(), fCenter.GetY(), fCenter.GetZ(),
	  u.GetX(), u.GetY(), u.GetZ());
  TH2F *h = new TH2F("cross", title, fNx, -fSizeX/2.0, fSizeX/2.0,
		     fNy, -fSizeY/2.0, fSizeY/2.0);
  
  TVpParticle particle;
  TVpVector3 cross;  // crossing point of the particle ray and the cut plane
  for (Double_t l = 0; l < nparticles; l++)
    {
      source->GetParticle(&particle);
      Double_t t = (1 / (u | particle.fDirection)) *
	((fCenter - particle.fPosition) | particle.fDirection);
      if (t < 0)  // Don't consider particles travelling in the opposite direction
	continue;
      cross = particle.fPosition + t * particle.fDirection;
      Double_t t1 = (cross - fCenter)|b1;
      Double_t t2 = (cross - fCenter)|b2;
      h->Fill(t1, t2);
    }
  return h;
}

//______________________________________________________________________________
Int_t TVpRectangleCut::WriteLac(const Char_t *filename, TVpGeometry *geometry,
				Double_t energy, Int_t ndivision)
{
  // Write LAC

  return WriteMCSSum(filename, geometry, energy, 0, ndivision);
}

//______________________________________________________________________________
Int_t TVpRectangleCut::WriteMCSSum(const Char_t *filename, TVpGeometry *geometry,
				   Double_t energy, Int_t csSum, Int_t ndivision)
{
  // Write Macroscopic Cross Section (MCS) Sum

  Int_t npoints = ndivision * ndivision;
  Double_t a1 = 1.0 / (2.0*ndivision);
  Double_t a2 = 1.0 / (Double_t) ndivision;
  FILE *fp;

  if ((fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "Error: TVpRectangleCut::WriteLac: Cannot open file %s:", filename);
      perror("");
      return 1;
    }
  TVpVector3 position;

  for (Int_t j = 0; j < fNy; j++)  
    for (Int_t i = 0; i < fNx; i++)
      {
	Double_t sum = 0.0;
	for (Int_t id = 0; id < ndivision; id++)
	  for (Int_t jd = 0; jd < ndivision; jd++)
	    {
	      position = fOrigin + (a1+i+a2*id) * fBaseX + (a1+j+a2*jd) * fBaseY;
	      TVpSolid *solid = geometry->GetSolid(position);
	      sum += (csSum == 0) ? solid->GetMaterial(position.Get())->GetLac(energy):
		solid->GetMaterial(position.Get())->GetMCSSum(energy, csSum);

	    }
	fprintf(fp, "%e\n", sum / npoints);
      }
  fclose(fp);
  return 0;
}



//______________________________________________________________________________
void TVpRectangleCut::Draw()
{
  // Draw the rectangle (cut plane) in 3D

  TVpVector3 p0 = fOrigin;
  TVpVector3 p1 = fOrigin + fNx * fBaseX;
  TVpVector3 p2 = fOrigin + fNx * fBaseX + fNy * fBaseY;
  TVpVector3 p3 = fOrigin + fNy * fBaseY;
  TPolyLine3D *border = new TPolyLine3D(5);
  border->SetPoint(0, p0.GetX(), p0.GetY(), p0.GetZ());
  border->SetPoint(1, p1.GetX(), p1.GetY(), p1.GetZ());
  border->SetPoint(2, p2.GetX(), p2.GetY(), p2.GetZ());
  border->SetPoint(3, p3.GetX(), p3.GetY(), p3.GetZ());
  border->SetPoint(4, p0.GetX(), p0.GetY(), p0.GetZ());
  border->SetLineColor(1);
  border->Draw(); 
}

//______________________________________________________________________________
void TVpRectangleCut::DrawSourceParticle(TVpSource *source, Long_t nparticles)
{
  // Draw the rectangle (cut plane) in 3D

  TVpVector3 n = normalize(fBaseX * fBaseY);  // normal to the surface
  TVpParticle particle;
  TVpVector3 cross;  // crossing point of the particle ray and the cut plane
  for (Double_t l = 0; l < nparticles; l++)
    {
      source->GetParticle(&particle);
      // std::cout << particle << std::endl;
      
      Double_t t = (1 / (n | particle.fDirection)) * ((fCenter - particle.fPosition) | n);
      if (t < 0)  // Don't consider particles travelling in the opposite direction
	continue;
      cross = particle.fPosition + t * particle.fDirection;

      TPolyLine3D *line = new TPolyLine3D(2);
      line->SetPoint(0, particle.fPosition.GetX(), particle.fPosition.GetY(),
		     particle.fPosition.GetZ());
      line->SetPoint(1, cross.GetX(), cross.GetY(), cross.GetZ());
      line->SetLineColor(1);
      line->Draw();
    }
}


