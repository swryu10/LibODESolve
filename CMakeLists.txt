cmake_minimum_required(VERSION 3.31)

project(ODESolve)
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
add_library(ODESolve STATIC)
target_sources(ODESolve PUBLIC ${PROJECT_SOURCE_DIR}/src/IntTrapezoid.cpp
                               ${PROJECT_SOURCE_DIR}/src/ODESolveBase.cpp
                               ${PROJECT_SOURCE_DIR}/src/nrutil.c
                               ${PROJECT_SOURCE_DIR}/src/ODESolveRK.cpp
                               ${PROJECT_SOURCE_DIR}/src/ODESolveRX.cpp)
target_include_directories(ODESolve PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_IntTrapezoid_pi)
if(OPENMP_FOUND)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(test_IntTrapezoid_pi ${PROJECT_SOURCE_DIR}/src/test_IntTrapezoid_pi.cpp)
set_target_properties(test_IntTrapezoid_pi PROPERTIES OUTPUT_NAME "test_IntTrapezoid_pi.exec")
target_link_libraries(test_IntTrapezoid_pi PUBLIC ODESolve)
target_include_directories(test_IntTrapezoid_pi PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_ODERK_cos)
add_executable(test_ODERK_cos ${PROJECT_SOURCE_DIR}/src/test_ODERK_cos.cpp)
set_target_properties(test_ODERK_cos PROPERTIES OUTPUT_NAME "test_ODERK_cos.exec")
target_link_libraries(test_ODERK_cos PUBLIC ODESolve)
target_include_directories(test_ODERK_cos PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_ODERX_cos)
add_executable(test_ODERX_cos ${PROJECT_SOURCE_DIR}/src/test_ODERX_cos.cpp)
set_target_properties(test_ODERX_cos PROPERTIES OUTPUT_NAME "test_ODERX_cos.exec")
target_link_libraries(test_ODERX_cos PUBLIC ODESolve)
target_include_directories(test_ODERX_cos PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(proj_OrbitKepler)
add_executable(proj_OrbitKepler ${PROJECT_SOURCE_DIR}/src/OrbitKepler.cpp
                                ${PROJECT_SOURCE_DIR}/src/main_proj_kepler.cpp)
set_target_properties(proj_OrbitKepler PROPERTIES OUTPUT_NAME "proj_OrbitKepler.exec")
target_link_libraries(proj_OrbitKepler PUBLIC ODESolve)
target_include_directories(proj_OrbitKepler PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(proj_GeodesicSphere)
find_package(Python COMPONENTS Interpreter Development)
add_executable(proj_GeodesicSphere ${PROJECT_SOURCE_DIR}/src/PythonHook.cpp
                                   ${PROJECT_SOURCE_DIR}/src/WrapGeoIATA.cpp
                                   ${PROJECT_SOURCE_DIR}/src/GeodesicSphere.cpp
                                   ${PROJECT_SOURCE_DIR}/src/main_proj_geodesic.cpp)
set_target_properties(proj_GeodesicSphere PROPERTIES OUTPUT_NAME "proj_GeodesicSphere.exec")
target_link_libraries(proj_GeodesicSphere PUBLIC ODESolve
                                                 ${Python_LIBRARIES})
target_include_directories(proj_GeodesicSphere PUBLIC "${PROJECT_SOURCE_DIR}/include"
                                                      "${Python_INCLUDE_DIRS}")
file(COPY "${PROJECT_SOURCE_DIR}/src/GeoIATA.py" DESTINATION "${CMAKE_BINARY_DIR}")
