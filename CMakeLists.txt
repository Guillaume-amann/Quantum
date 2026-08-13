cmake_minimum_required(VERSION 3.14)
project(Quantum)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directory for all executables to bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Enable CTest for automated test discovery
enable_testing()

# Include tests subdirectory
add_subdirectory(tests)

message(STATUS "")
message(STATUS "Build configuration:")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Runtime output: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "")
message(STATUS "To build everything:")
message(STATUS "  cd build && cmake --build . --target all")
message(STATUS "")
message(STATUS "To run unit tests:")
message(STATUS "  cd build && ctest --output-on-failure")
message(STATUS "  or from project root: ./run_tests.sh")
message(STATUS "")