#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

int main()
{
  polyscope::init();
  
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
  auto h=0.5; //gridstep
  params( "polynomial", "goursat" )( "gridstep", h );
  auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  auto K               = SH3::getKSpace( params );
  auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
  auto surface         = SH3::makeDigitalSurface( binary_image, K, params );
  auto embedder        = SH3::getCellEmbedder( K );
  SH3::Cell2Index c2i;
  auto surfels         = SH3::getSurfelRange( surface, params );
  auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
  //Need to convert the faces
  std::vector<std::vector<std::size_t>> faces;
  for(auto &face: primalSurface->allFaces())
    faces.push_back(primalSurface->verticesAroundFace( face ));
  
  polyscope::registerPointCloud("Pointels", primalSurface->positions());
  
  auto digsurf = polyscope::registerSurfaceMesh("Primal surface", primalSurface->positions(), faces);
  digsurf->setEdgeWidth(1.0);
  digsurf->setEdgeColor({1.,1.,1.});
  
  //Computing some differential quantities
  params("r-radius", 5.0);
  auto normals   = SHG3::getCTrivialNormalVectors(surface, surfels, params);
  auto normalsTrivial   = SHG3::getTrivialNormalVectors(K,surfels);
  auto normalsII = SHG3::getIINormalVectors(binary_image, surfels, params);
  auto Mcurv     = SHG3::getIIMeanCurvatures(binary_image, surfels, params);
  auto Gcurv     = SHG3::getIIGaussianCurvatures(binary_image, surfels, params);

  //Surfel area measure
  std::vector<double> areaMeasure(surfels.size());
  for(auto i=0; i < areaMeasure.size(); ++i)
    areaMeasure[i] = normalsTrivial[i].dot(normalsII[i]);
   
  
  //Attaching quantities
  digsurf->addFaceVectorQuantity("Trivial normal vectors", normalsTrivial);
  digsurf->addFaceVectorQuantity("CTrivial normal vectors", normals);
  digsurf->addFaceVectorQuantity("II normal vectors", normalsII);
  digsurf->addFaceScalarQuantity("II mean curvature", Mcurv);
  digsurf->addFaceScalarQuantity("II Gaussian curvature", Gcurv);
  digsurf->addFaceScalarQuantity("Surfel area measure", areaMeasure);
  
  
  polyscope::show();
  return 0;
}
