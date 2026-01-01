**WebAssembly Build for phaseshifts**

This directory contains the build infrastructure for compiling the `libphsh.f`
Fortran phase shift code to WebAssembly, enabling browser-based calculations.

## Overview

The WASM build uses Emscripten to compile the Fortran code to WebAssembly with
MEMFS (in-memory filesystem) for virtualized file I/O. This allows the original
Fortran code to run in the browser with minimal modifications.

### Supported Algorithms

- **phsh_rel** - Relativistic phase shift calculation
- **phsh_cav** - Cavity LEED phase shift calculation
- **phsh_wil** - Williams' phase shift calculation

## Prerequisites

The build script supports three methods. Choose one:

### Option 1: LFortran (RECOMMENDED)

LFortran has native WASM support and is the easiest option:

```bash
# Via pip
pip install lfortran

# Or via conda
conda install -c conda-forge lfortran
```

### Option 2: Emscripten + f2c

```bash
# Install Emscripten
brew install emscripten
# Or from source:
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk && ./emsdk install latest && ./emsdk activate latest
source ./emsdk_env.sh

# f2c is NOT in Homebrew. Build from source automatically:
./build.sh --install-f2c

# Or on Ubuntu/Debian:
sudo apt-get install f2c libf2c2-dev
```

### Option 3: Docker (no local tools needed)

```bash
# Docker handles everything
brew install docker
# or: https://docs.docker.com/get-docker/

./build.sh --method=docker
```

## Building

### Quick Build (auto-detects method)

```bash
# From repository root
make wasm

# Or directly
cd wasm
./build.sh
```

### Build Options

```bash
./build.sh --help           # Show help
./build.sh --method=lfortran  # Force LFortran
./build.sh --method=f2c       # Force Emscripten+f2c
./build.sh --method=docker    # Force Docker
./build.sh --install-f2c      # Build f2c from source first
./build.sh --clean            # Clean before building
./build.sh --debug            # Debug build
```

### Manual Build Steps

1. Convert Fortran to C:

   ```bash
   f2c -a ../phaseshifts/lib/libphsh.f
   ```

2. Compile with Emscripten:

   ```bash
   emcc libphsh.c -lf2c \
       -s WASM=1 \
       -s MODULARIZE=1 \
       -s EXPORT_NAME="createPhaseShiftsModule" \
       -s EXPORTED_FUNCTIONS="['_hartfock_', '_phsh_cav_', \
       '_phsh_rel_', '_phsh_wil_']" \
       -s EXPORTED_RUNTIME_METHODS="['FS', 'ccall', 'cwrap']" \
       -s FORCE_FILESYSTEM=1 \
       -s ALLOW_MEMORY_GROWTH=1 \
       -o dist/phaseshifts.js
   ```

## Output Files

After building, the `dist/` directory will contain:

- `phaseshifts.js` - JavaScript loader and glue code
- `phaseshifts.wasm` - WebAssembly binary

## Usage

### In Browser

```html
<script src="dist/phaseshifts.js"></script>
<script>
  createPhaseShiftsModule().then((Module) => {
    // Write input file to virtual filesystem
    Module.FS.writeFile('/input.dat', inputData);

    // Run phase shift calculation
    Module.ccall('phsh_rel_', null, [], []);

    // Read output
    const output = Module.FS.readFile('/output.dat', { encoding: 'utf8' });
  });
</script>
```

### Using the Web Interface

Open `web/index.html` in a browser to use the graphical interface.

## Development

### Testing Locally

```bash
# Start a local server (required for WASM)
cd wasm
python3 -m http.server 8080

# Open http://localhost:8080/web/
```

### Running Tests

```bash
# From repository root
pytest tests/test_wasm_build.py -v
```

## Architecture

```text
wasm/
├── README.md           # This file
├── build.sh            # Main build script
├── CMakeLists.txt      # CMake configuration (optional)
├── src/
│   ├── phaseshifts.js  # JavaScript API wrapper
│   └── libf2c/         # f2c runtime library (if bundled)
├── web/
│   ├── index.html      # Browser interface
│   ├── app.js          # Application logic
│   └── style.css       # Styling
└── dist/               # Build output (generated)
    ├── phaseshifts.js
    └── phaseshifts.wasm
```

## Troubleshooting

### "emcc not found"

Ensure Emscripten is installed and the environment is activated:

```bash
source /path/to/emsdk/emsdk_env.sh
```

### "f2c not found"

f2c is not available in Homebrew on macOS. Use one of:

```bash
# Build from source (automated)
./build.sh --install-f2c

# Or use Docker instead
./build.sh --method=docker

# Or use LFortran (no f2c needed)
pip install lfortran
./build.sh --method=lfortran
```

### Memory errors in browser

Increase the memory limit in the build:

```bash
emcc ... -s INITIAL_MEMORY=64MB -s MAXIMUM_MEMORY=256MB ...
```

### File not found errors

Ensure you're writing files to the MEMFS before calling Fortran functions:

```javascript
Module.FS.writeFile('/input.dat', data);
```

## License

Same as the parent phaseshifts project (MIT License).

## References

- [Emscripten Documentation](https://emscripten.org/docs/)
- [f2c Manual](https://www.netlib.org/f2c/)
- [LEED Theory][leed-theory]

[leed-theory]: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/
