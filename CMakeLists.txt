cmake_minimum_required(VERSION 3.31)

project(ODESolveRK)
add_library(ODESolveRK STATIC)
target_sources(ODESolveRK PUBLIC ${PROJECT_SOURCE_DIR}/src/ODESolveBase.cpp
                                 ${PROJECT_SOURCE_DIR}/src/ODESolveRK.cpp)
target_include_directories(ODESolveRK PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(ODESolveRX)
add_library(ODESolveRX STATIC)
target_sources(ODESolveRX PUBLIC ${PROJECT_SOURCE_DIR}/src/ODESolveBase.cpp
                                 ${PROJECT_SOURCE_DIR}/src/nrutil.c
                                 ${PROJECT_SOURCE_DIR}/src/ODESolveRX.cpp)
target_include_directories(ODESolveRX PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_ODERK_cos)
add_executable(test_ODERK_cos ${PROJECT_SOURCE_DIR}/src/test_ODERK_cos.cpp)
set_target_properties(test_ODERK_cos PROPERTIES OUTPUT_NAME "test_ODERK_cos.exec")
target_link_libraries(test_ODERK_cos PUBLIC ODESolveRK)
target_include_directories(test_ODERK_cos PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_ODERX_cos)
add_executable(test_ODERX_cos ${PROJECT_SOURCE_DIR}/src/test_ODERX_cos.cpp)
set_target_properties(test_ODERX_cos PROPERTIES OUTPUT_NAME "test_ODERX_cos.exec")
target_link_libraries(test_ODERX_cos PUBLIC ODESolveRX)
target_include_directories(test_ODERX_cos PUBLIC "${PROJECT_SOURCE_DIR}/include")
