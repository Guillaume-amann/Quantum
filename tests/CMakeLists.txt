# Include core headers
include_directories(${CMAKE_SOURCE_DIR})

# Find Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# Add each test executable
add_executable(test_density_matrix test_density_matrix.cpp)
target_link_libraries(test_density_matrix PRIVATE Eigen3::Eigen)

add_executable(test_gate test_gate.cpp)
target_link_libraries(test_gate PRIVATE Eigen3::Eigen)

add_executable(test_noisemodel test_noisemodel.cpp)
target_link_libraries(test_noisemodel PRIVATE Eigen3::Eigen)

add_executable(test_qbit test_qbit.cpp)
target_link_libraries(test_qbit PRIVATE Eigen3::Eigen)