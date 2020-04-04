#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"


#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "args/args.hxx"

#include "CorrectedNormalCurrentFormula.h"
#include "NormalCycleFormula.h"
#include "RusinkiewiczCurvatureFormula.h"

using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;


using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

float clampM=10.0;

// Computing principal curvatures k1 and k2
static
std::pair<RealVector, RealVector>
curvDirFromTensor(const CorrectedNormalCurrentFormula<RealVector,RealVector>::RealTensor &tensor,
                  const double area,
                  const RealVector &N)
{
  CorrectedNormalCurrentFormula<RealVector,RealVector>::RealTensor M = tensor;
  M += M.transpose();
  M *= 0.5;
  const double   coef_N = 1000.0 * area;
  // Adding 1000 area n x n to anisotropic measure
  for ( int j = 0; j < 3; j++ )
    for ( int k = 0; k < 3; k++ )
      M( j, k ) += coef_N * N[ j ] * N[ k ];
  auto V = M;
  RealVector L;
  std::cout<<N.norm()<<std::endl;
  EigenDecomposition< 3, double>::getEigenDecomposition( M, V, L );
  return std::pair<RealVector,RealVector>(V.column( 1 ),  V.column( 0 ));
}



// Example computation function -- this one computes and registers a scalar
// quantity
void doWork()
{
  geometry->requireVertexGaussianCurvatures();
  geometry->requireFaceNormals();
  geometry->requireVertexNormals();
  geometry->requireVertexPositions();
  psMesh->addVertexScalarQuantity("curvature polyscope",
                                  geometry->vertexGaussianCurvatures,
                                  polyscope::DataType::SYMMETRIC);
  // Set vertex tangent spaces
  geometry->requireFaceTangentBasis();
  FaceData<Vector3> vBasisX(*mesh);
  for(Face f : mesh->faces()) {
    vBasisX[f] = geometry->faceTangentBasis[f][0];
  }
  psMesh->setFaceTangentBasisX(vBasisX);

  
  FaceData<double> cncMean(*mesh);
  FaceData<double> cncGauss(*mesh);
  FaceData<double> rusMean(*mesh);
  FaceData<double> rusGauss(*mesh);
  FaceData<double> m0(*mesh);
  FaceData<double> m1(*mesh);
  FaceData<double> m2(*mesh);
  VertexData<RealPoint> normV(*mesh);

  FaceData<RealVector> d1(*mesh);
  FaceData<RealVector> d2(*mesh);

  FaceData<std::complex<double>> intd1(*mesh);
  FaceData<std::complex<double>> intd2(*mesh);

  VertexData<double> ncGauss(*mesh);
  EdgeData<double> ncMean(*mesh);
  
  auto clamp= [](double v){ return (v< -clampM)? -clampM: (v>clampM)? clampM: v; };
  
  
  //Default polyscope GC does NC
  for(auto vert: mesh->vertices())
    ncGauss[vert] = clamp(geometry->vertexGaussianCurvatures[vert]);
  
  for(auto face: mesh->faces())
  {
    auto vb = face.adjacentVertices().begin();
    auto A =  geometry->vertexPositions[*vb];
    auto n =  geometry->vertexNormals[*vb];
    RealPoint pA(A.x,A.y, A.z);
    RealPoint nA(n.x,n.y, n.z);
    normV[*vb] = nA;
    
    ++vb;
    A =  geometry->vertexPositions[*vb];
    n =  geometry->vertexNormals[*vb];
    RealPoint pB(A.x,A.y, A.z);
    RealPoint nB(n.x,n.y, n.z);
    normV[*vb] = nB;
    ++vb;
    A =  geometry->vertexPositions[*vb];
    n =  geometry->vertexNormals[*vb];
    RealPoint pC(A.x,A.y, A.z);
    RealPoint nC(n.x,n.y, n.z);
    normV[*vb] = nC;
    
    m0[face] = CorrectedNormalCurrentFormula<RealPoint, RealPoint>::mu0InterpolatedU(pA, pB, pC, nA, nB, nC);
    m1[face] = CorrectedNormalCurrentFormula<RealPoint, RealPoint>::mu1InterpolatedU(pA, pB, pC, nA, nB, nC);
    m2[face] = CorrectedNormalCurrentFormula<RealPoint, RealPoint>::mu2InterpolatedU(pA, pB, pC, nA, nB, nC);
    cncMean[face] = clamp(m1[face]/m0[face]);
    cncGauss[face] = clamp(m2[face]/m0[face]);
    
    auto mxy =CorrectedNormalCurrentFormula<RealPoint, RealPoint>::muXYInterpolatedU(pA, pB, pC, nA, nB, nC);
    RealVector dd1,dd2;
    auto nf= geometry->faceNormal(face);
    RealVector nFace(nf.x,nf.y,nf.z);
    std::tie(dd1,dd2) = curvDirFromTensor(mxy, geometry->faceArea(face), nFace);
    d1[face] = dd1;
    d2[face] = dd2;
    auto vB = vBasisX[face];
    RealVector tan(vB.x,vB.y,vB.z);
    RealVector bitan(geometry->faceTangentBasis[face][1].x,
                     geometry->faceTangentBasis[face][1].y,
                     geometry->faceTangentBasis[face][1].z);

    double angle =std::acos(tan.dot(dd1));
    if (bitan.dot(dd1) < 0.0)
      angle  = 2*M_PI - angle;
    intd1[face] = std::polar(1.0, angle);
    
    angle =std::acos(tan.dot(dd2));
      if (bitan.dot(dd2) < 0.0)
        angle  = 2*M_PI - angle;
    intd2[face] = std::polar(1.0, angle);
    
    rusMean[face] = clamp(RusinkiewiczCurvatureFormula<RealPoint, RealPoint>::meanCurvature(pA, pB, pC, nA, nB, nC));
    rusGauss[face] = clamp(RusinkiewiczCurvatureFormula<RealPoint, RealPoint>::gaussianCurvature(pA, pB, pC, nA, nB, nC));
  }
  
  
  
  //NormalCycle per edge
  for(auto edge: mesh->edges())
  {
    auto vA = edge.halfedge().vertex();
    auto vB = edge.halfedge().twin().vertex();
    auto A = geometry->vertexPositions[vA];
    auto nA = geometry->faceNormal( edge.halfedge().face()  );
    auto B = geometry->vertexPositions[vB];
    auto nB = geometry->faceNormal(  edge.halfedge().twin().face() );
    RealPoint pA(A.x,A.y,A.z);
    RealPoint pB(B.x,B.y,B.z);
    RealPoint nnA(nA.x,nA.y,nA.z);
    RealPoint nnB(nB.x,nB.y,nB.z);
    ncMean[ edge ] = clamp(NormalCycleFormula<RealPoint, RealPoint>::meanCurvature(pA, pB, nnA, nnB));
  }
  
  psMesh->addVertexScalarQuantity("NC Gauss",ncGauss,
                                polyscope::DataType::SYMMETRIC);
  psMesh->addEdgeScalarQuantity("NC mean",ncMean,
                                polyscope::DataType::SYMMETRIC);
  
  psMesh->addFaceScalarQuantity("CNC mean",cncMean,
                                polyscope::DataType::SYMMETRIC);
  psMesh->addFaceScalarQuantity("CNC Gauss",cncGauss,
                                polyscope::DataType::SYMMETRIC);
  psMesh->addFaceScalarQuantity("Rusinkiewicz mean",rusMean,
                                polyscope::DataType::SYMMETRIC);
  psMesh->addFaceScalarQuantity("Rusinkiewicz Gauss",rusGauss,
                                polyscope::DataType::SYMMETRIC);
  psMesh->addFaceScalarQuantity("mu0",m0);
  psMesh->addFaceScalarQuantity("mu1",m1);
  psMesh->addFaceScalarQuantity("mu2",m2);
 
  psMesh->addFaceVectorQuantity("dir1 CNC",d1);
  psMesh->addFaceVectorQuantity("dir2 CNC",d2);

  psMesh->addFaceIntrinsicVectorQuantity("dir1 bis CNC",intd1);
  psMesh->addFaceIntrinsicVectorQuantity("dir2 bis CNC",intd2);

  
  psMesh->addVertexVectorQuantity("Normal at vertices", normV);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  ImGui::SliderFloat("Clamping value", &clampM, 0.0, 100);
  if (ImGui::Button("do work"))
  {
    doWork();
  }
}


int main(int argc, char **argv)
{
  // Configure the argument parser
   args::ArgumentParser parser("geometry-central & Polyscope example project");
   args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

   // Parse args
   try {
     parser.ParseCLI(argc, argv);
   } catch (args::Help) {
     std::cout << parser;
     return 0;
   } catch (args::ParseError e) {
     std::cerr << e.what() << std::endl;
     std::cerr << parser;
     return 1;
   }

   // Make sure a mesh name was given
   if (!inputFilename) {
     std::cerr << "Please specify a mesh file as argument" << std::endl;
     return EXIT_FAILURE;
   }

  
  polyscope::init();
  
  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  
  
  
  polyscope::show();
  return 0;
}
