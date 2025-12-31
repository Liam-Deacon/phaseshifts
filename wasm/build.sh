#!/usr/bin/env bash
# Build script for WebAssembly version of phaseshifts
# 
# This script compiles libphsh.f to WebAssembly using Emscripten
# with MEMFS for virtual file I/O.
#
# Prerequisites:
#   - Emscripten SDK (emcc)
#   - f2c (Fortran to C converter)
#
# Usage:
#   ./build.sh [options]
#
# Options:
#   --clean     Remove build artifacts before building
#   --debug     Build with debug symbols
#   --help      Show this help message

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
FORTRAN_SRC="$PROJECT_ROOT/phaseshifts/lib/libphsh.f"
BUILD_DIR="$SCRIPT_DIR/build"
DIST_DIR="$SCRIPT_DIR/dist"
SRC_DIR="$SCRIPT_DIR/src"

# Build options
DEBUG=0
CLEAN=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --clean)
            CLEAN=1
            shift
            ;;
        --debug)
            DEBUG=1
            shift
            ;;
        --help)
            head -20 "$0" | tail -17
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check prerequisites
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed or not in PATH"
        echo "Please install $1 and try again."
        echo ""
        case $1 in
            emcc)
                echo "Install Emscripten:"
                echo "  brew install emscripten"
                echo "  # or"
                echo "  git clone https://github.com/emscripten-core/emsdk.git"
                echo "  cd emsdk && ./emsdk install latest && ./emsdk activate latest"
                echo "  source ./emsdk_env.sh"
                ;;
            f2c)
                echo "Install f2c:"
                echo "  brew install f2c  # macOS"
                echo "  apt-get install f2c  # Ubuntu/Debian"
                ;;
        esac
        exit 1
    fi
}

echo "=== phaseshifts WebAssembly Build ==="
echo ""

# Check for required tools
echo "Checking prerequisites..."
check_command emcc
check_command f2c
echo "  ✓ emcc found: $(emcc --version | head -1)"
echo "  ✓ f2c found"
echo ""

# Clean if requested
if [[ $CLEAN -eq 1 ]]; then
    echo "Cleaning build artifacts..."
    rm -rf "$BUILD_DIR" "$DIST_DIR"
fi

# Create directories
mkdir -p "$BUILD_DIR" "$DIST_DIR"

# Step 1: Convert Fortran to C using f2c
echo "Step 1: Converting Fortran to C..."
cd "$BUILD_DIR"

if [[ ! -f "$FORTRAN_SRC" ]]; then
    echo "Error: Fortran source not found: $FORTRAN_SRC"
    exit 1
fi

# Copy and convert Fortran source
cp "$FORTRAN_SRC" ./libphsh.f

# f2c conversion with automatic type sizing
f2c -a -A -C++ libphsh.f 2>&1 | head -20 || true

if [[ ! -f "libphsh.c" ]]; then
    echo "Error: f2c conversion failed"
    exit 1
fi
echo "  ✓ Generated libphsh.c"

# Step 2: Create wrapper functions for JavaScript interop
echo "Step 2: Creating JavaScript interop wrapper..."
cat > wrapper.c << 'EOF'
/*
 * JavaScript interop wrapper for phaseshifts WASM module
 * 
 * This file provides C-callable wrappers that can be invoked from JavaScript
 * via Emscripten's ccall/cwrap functions.
 */

#include <stdio.h>
#include <string.h>
#include <emscripten.h>

/* Forward declarations of Fortran functions (f2c naming convention) */
extern int hartfock_(char *input_file, int input_file_len);
extern int phsh_cav_(void);
extern int phsh_rel_(void);
extern int phsh_wil_(void);
extern int abinitio_(void);

/* 
 * Wrapper for hartfock subroutine
 * Calculates atomic radial charge density using Dirac-Fock method
 * 
 * @param input_file Path to input file in MEMFS
 * @return 0 on success, non-zero on error
 */
EMSCRIPTEN_KEEPALIVE
int run_hartfock(const char* input_file) {
    char fname[256];
    strncpy(fname, input_file, 255);
    fname[255] = '\0';
    int len = strlen(fname);
    return hartfock_(fname, len);
}

/*
 * Wrapper for phsh_cav (cavity LEED phase shifts)
 * Uses Loucks grid for potential-to-phase-shift calculation
 */
EMSCRIPTEN_KEEPALIVE
int run_phsh_cav(void) {
    return phsh_cav_();
}

/*
 * Wrapper for phsh_rel (relativistic phase shifts)
 * Calculates phase shifts with relativistic corrections
 */
EMSCRIPTEN_KEEPALIVE
int run_phsh_rel(void) {
    return phsh_rel_();
}

/*
 * Wrapper for phsh_wil (Williams method phase shifts)
 * A.R. Williams' phase shift calculation method
 */
EMSCRIPTEN_KEEPALIVE
int run_phsh_wil(void) {
    return phsh_wil_();
}

/*
 * Initialize the module and file system
 * Called automatically on module load
 */
EMSCRIPTEN_KEEPALIVE
void init_phaseshifts(void) {
    printf("phaseshifts WASM module initialized\n");
}

/*
 * Get version information
 */
EMSCRIPTEN_KEEPALIVE
const char* get_version(void) {
    return "phaseshifts-wasm 0.1.0";
}
EOF
echo "  ✓ Generated wrapper.c"

# Step 3: Compile with Emscripten
echo "Step 3: Compiling with Emscripten..."

# Set compiler flags based on debug mode
if [[ $DEBUG -eq 1 ]]; then
    EMCC_FLAGS="-O0 -g -s ASSERTIONS=2"
    echo "  (Debug build enabled)"
else
    EMCC_FLAGS="-O2"
fi

# Locate f2c library
F2C_LIB=""
if [[ -f "/usr/local/lib/libf2c.a" ]]; then
    F2C_LIB="/usr/local/lib/libf2c.a"
elif [[ -f "/opt/homebrew/lib/libf2c.a" ]]; then
    F2C_LIB="/opt/homebrew/lib/libf2c.a"
elif [[ -f "/usr/lib/libf2c.a" ]]; then
    F2C_LIB="/usr/lib/libf2c.a"
fi

# Check if we need to bundle f2c runtime
if [[ -z "$F2C_LIB" ]]; then
    echo "  Warning: libf2c.a not found, will try system linking"
    F2C_LINK="-lf2c"
else
    echo "  Found f2c library: $F2C_LIB"
    F2C_LINK="$F2C_LIB"
fi

# Main compilation
emcc $EMCC_FLAGS \
    libphsh.c \
    wrapper.c \
    $F2C_LINK \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME="createPhaseShiftsModule" \
    -s EXPORTED_FUNCTIONS="['_run_hartfock', '_run_phsh_cav', '_run_phsh_rel', '_run_phsh_wil', '_init_phaseshifts', '_get_version', '_main']" \
    -s EXPORTED_RUNTIME_METHODS="['FS', 'ccall', 'cwrap', 'UTF8ToString', 'stringToUTF8']" \
    -s FORCE_FILESYSTEM=1 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=67108864 \
    -s MAXIMUM_MEMORY=536870912 \
    -s NO_EXIT_RUNTIME=1 \
    -s ENVIRONMENT='web,worker' \
    -o "$DIST_DIR/phaseshifts.js"

echo "  ✓ Generated phaseshifts.js"
echo "  ✓ Generated phaseshifts.wasm"

# Step 4: Copy JavaScript wrapper API
echo "Step 4: Installing JavaScript API..."
cp "$SRC_DIR/phaseshifts.js" "$DIST_DIR/phaseshifts-api.js" 2>/dev/null || true
echo "  ✓ Installed phaseshifts-api.js"

# Step 5: Show build summary
echo ""
echo "=== Build Complete ==="
echo ""
echo "Output files:"
ls -lh "$DIST_DIR"/*.{js,wasm} 2>/dev/null || ls -lh "$DIST_DIR"
echo ""
echo "To test locally:"
echo "  cd $SCRIPT_DIR"
echo "  python3 -m http.server 8080"
echo "  # Open http://localhost:8080/web/"
echo ""
