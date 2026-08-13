# tests/CMakeLists.txt
# Unit tests for the Quantum simulator core infrastructure

# Include core headers
include_directories(${CMAKE_SOURCE_DIR})

# Find Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# Enable testing
enable_testing()

# Define test executables with shared configuration
set(TEST_SOURCES
    test_gate.cpp
    test_density_matrix.cpp
    test_qbit.cpp
    test_noisemodel.cpp
)

# Build each test
foreach(TEST_SOURCE ${TEST_SOURCES})
    # Extract test name from filename (e.g. test_gate.cpp -> test_gate)
    string(REPLACE ".cpp" "" TEST_NAME "${TEST_SOURCE}")
    
    # Create executable
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    
    # Link Eigen
    target_link_libraries(${TEST_NAME} PRIVATE Eigen3::Eigen)
    
    # Set compiler flags for this test (C++17, optimised)
    target_compile_features(${TEST_NAME} PRIVATE cxx_std_17)
    
    # Add to CTest registry (so `ctest` can discover and run it)
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${TEST_NAME}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    
    # Optional: set timeout per test (in seconds)
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 30)
endforeach()

# Print summary
message(STATUS "================================")
message(STATUS "  Unit Tests Configuration")
message(STATUS "================================")
message(STATUS "Test executables will be built to: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "Tests discovered by CTest:")
foreach(TEST_SOURCE ${TEST_SOURCES})
    string(REPLACE ".cpp" "" TEST_NAME "${TEST_SOURCE}")
    message(STATUS "  • ${TEST_NAME}")
endforeach()
message(STATUS "")
message(STATUS "To run all tests:")
message(STATUS "  cd build && ctest --output-on-failure")
message(STATUS "  or")
message(STATUS "  ../run_tests.sh")
message(STATUS "================================")