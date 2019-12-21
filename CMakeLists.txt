project(polyscope-dgtal)


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_subdirectory(deps)

add_executable(simpleTest simpleTest)
target_link_libraries(simpleTest polyscope ${DGTAL_LIBRARIES})
