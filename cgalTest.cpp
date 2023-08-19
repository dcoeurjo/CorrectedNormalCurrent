#include <iostream>
#include <thread>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/property_map.h>
#include <boost/graph/graph_traits.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> Surface_Mesh;
typedef boost::graph_traits<Surface_Mesh>::vertex_descriptor vertex_descriptor;


// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;


// Integration radius
float Radius;

void doWork()
{
  
}


void myCallback()
{
  ImGui::SliderFloat("Measuring ball radius", &Radius, 0.0, 1.0);
  
  if (ImGui::Button("do work"))
    doWork();
}


int main(int argc, char **argv)
{
  
  Surface_Mesh smesh;
  const std::string filename = argv[1] ;
  if (!CGAL::IO::read_polygon_mesh(filename, smesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  
  
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

 

  polyscope::show();
  return 0;
}
