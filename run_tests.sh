#!/bin/bash
# =============================================================================
#  run_tests.sh — Automated test runner for Quantum simulator
#
#  Usage:
#    ./run_tests.sh              # Run all tests in sequence
#    ./run_tests.sh --verbose    # Print full output from each test
#    ./run_tests.sh --rebuild    # Rebuild before testing
#
#  Exit code:
#    0 = all tests passed
#    1 = one or more tests failed
# =============================================================================

set -e

VERBOSE=false
REBUILD=false
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build"
BIN_DIR="${PROJECT_ROOT}/bin"

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --verbose) VERBOSE=true ;;
        --rebuild) REBUILD=true ;;
        *)         echo "Unknown option: $arg"; exit 1 ;;
    esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Quantum Simulator Test Suite${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Create build directory if it doesn't exist
if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}→ Creating build directory...${NC}"
    mkdir -p "$BUILD_DIR"
    REBUILD=true
fi

# Rebuild if requested or if build files don't exist
if [ "$REBUILD" = true ]; then
    echo -e "${YELLOW}→ Rebuilding project...${NC}"
    cd "$BUILD_DIR"
    cmake ..
    cmake --build . --target all
    cd "$PROJECT_ROOT"
    echo ""
fi

# Define test executables (in dependency order)
TESTS=(
    "test_gate"
    "test_density_matrix"
    "test_qbit"
    "test_noisemodel"
)

# Track results
PASSED=0
FAILED=0
FAILED_TESTS=()

# Run each test
for test_name in "${TESTS[@]}"; do
    test_path="${BIN_DIR}/${test_name}"
    
    if [ ! -f "$test_path" ]; then
        echo -e "${RED}✗ FAIL${NC} : ${test_name} (executable not found)"
        FAILED=$((FAILED + 1))
        FAILED_TESTS+=("$test_name")
        continue
    fi
    
    echo -e "${YELLOW}Running${NC} ${test_name}..."
    
    if [ "$VERBOSE" = true ]; then
        # Print full output
        if "$test_path"; then
            echo -e "${GREEN}✓ PASS${NC} : ${test_name}\n"
            PASSED=$((PASSED + 1))
        else
            echo -e "${RED}✗ FAIL${NC} : ${test_name}\n"
            FAILED=$((FAILED + 1))
            FAILED_TESTS+=("$test_name")
        fi
    else
        # Suppress output, show only result
        if output=$("$test_path" 2>&1); then
            echo -e "${GREEN}✓ PASS${NC} : ${test_name}"
            PASSED=$((PASSED + 1))
        else
            echo -e "${RED}✗ FAIL${NC} : ${test_name}"
            echo "$output" | head -20  # Show first 20 lines of error
            FAILED=$((FAILED + 1))
            FAILED_TESTS+=("$test_name")
        fi
    fi
done

# Summary
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Test Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "Total:  $((PASSED + FAILED)) tests"
echo -e "${GREEN}Passed: ${PASSED}${NC}"

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}Failed: ${FAILED}${NC}\n"
    echo -e "${GREEN}✓ All tests passed!${NC}\n"
    exit 0
else
    echo -e "${RED}Failed: ${FAILED}${NC}\n"
    echo -e "${RED}Failed tests:${NC}"
    for test in "${FAILED_TESTS[@]}"; do
        echo -e "  ${RED}✗${NC} ${test}"
    done
    echo ""
    exit 1
fi