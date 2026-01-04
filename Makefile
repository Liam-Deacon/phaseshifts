# Makefile for assisting with some common development activities

SHELL = /bin/sh
PYTHON_VERSION ?= 3.13
PYTHON := $(shell command -v python$(PYTHON_VERSION) 2>/dev/null || command -v python3 2>/dev/null || echo python)
DOCKER ?= $(shell command -v docker 2>/dev/null || docker)
PLATFORM ?= $(shell uname -s | tr '[:upper:]' '[:lower:]' 2>/dev/null || echo $${OS:-unknown})
CIBW_PLATFORM ?= linux


PREFIX ?= /usr/local

.PHONY: binfmt build-deps cbuildwheel check check-cibuildwheel check-legacy-wheels check-types install install-deps libphsh sdist test wheel

#: Quickly generate binary wheel
wheel: install-deps
	@if $(PYTHON) -m build --help >/dev/null 2>&1; then \
		$(PYTHON) -m build --wheel --no-isolation; \
	else \
		$(PYTHON) setup.py build bdist_wheel; \
	fi

#: Meta target, will attempt to build all it can
build: sdist
	if command -v docker 1>/dev/null; then \
		$(MAKE) cibuildwheel; \
	else \
		$(MAKE) wheel; \
	fi

#: Install binfmt for cross-platform builds (only on macOS)
binfmt:
	@if [ "$(PLATFORM)" = "darwin" ]; then \
	  docker run --privileged --rm tonistiigi/binfmt --install amd64; \
	fi

#: Build a matrix of wheels for different OSs and CPU archs
cibuildwheel: build-deps $(DOCKER)
	$(PYTHON) -m cibuildwheel --platform=$(CIBW_PLATFORM) --output-dir=dist .

#: Build Linux x86_64 wheels for CPython 3.11/3.12 (+ host‑specific extras)
check-cibuildwheel: build-deps binfmt
	@if [ "$(PLATFORM)" = "linux" ]; then \
		CIBW_ARCHS_LINUX=x86_64 \
		CIBW_BUILD="cp311-* cp312-*" \
		$(PYTHON) -m cibuildwheel --platform linux --output-dir dist . ; \
	elif [ "$(PLATFORM)" = "darwin" ]; then \
		export TAR_OPTIONS="--warning=no-unknown-keyword"; \
		CIBW_ARCHS_LINUX=x86_64 \
		CIBW_BUILD="cp311-* cp312-*" \
		$(PYTHON) -m cibuildwheel --platform linux --output-dir dist . && \
		CIBW_ARCHS_MACOS=arm64 \
		CIBW_BUILD="cp311-*arm64 cp312-*arm64" \
		$(PYTHON) -m cibuildwheel --platform macos --output-dir dist . ; \
	elif [ "$(OS)" = "Windows_NT" ]; then \
		CIBW_ARCHS_LINUX=x86_64 \
		CIBW_BUILD="cp311-* cp312-*" \
		$(PYTHON) -m cibuildwheel --platform linux --output-dir dist . && \
		CIBW_BUILD="cp311-* cp312-*" \
		$(PYTHON) -m cibuildwheel --platform windows --output-dir dist . ; \
	else \
		echo "Unsupported host platform $(PLATFORM)"; \
	fi

###############################################################################
# Build legacy Linux Py‑2.7 wheels inside an isolated venv running
# cibuildwheel 1.12.0 so it never clashes with whatever cibuildwheel version
# the main environment uses.
#
# * Creates .venv-cibw112 if missing
# * Installs pip + cibuildwheel 1.12.0 (plus numpy for build)
# * Runs the legacy build matrix just like `legacy-wheel27`, but through the
#   venv’s python executable.
###############################################################################
check-legacy-wheels: binfmt
	@VENV_DIR=.venv-cibw112; \
	if [ "$(OS)" = "Windows_NT" ]; then \
		if [ ! -d "$$VENV_DIR" ]; then python -m venv "$$VENV_DIR"; fi; \
		PYEXEC="$$VENV_DIR\\Scripts\\python.exe"; \
	else \
		if [ ! -d "$$VENV_DIR" ]; then $(PYTHON) -m venv "$$VENV_DIR"; fi; \
		PYEXEC="$$VENV_DIR/bin/python"; \
	fi; \
	"$$PYEXEC" -m pip install --upgrade pip >/dev/null; \
	"$$PYEXEC" -m pip install --quiet 'cibuildwheel==1.12.0' 'numpy>=1.16.6'; \
	if [ "$(OS)" = "Windows_NT" ]; then \
		PIP_USE_PEP517=0 \
		PIP_NO_BUILD_ISOLATION=1 \
		SKBUILD_DISABLE=1 \
		CIBW_ARCHS_LINUX=x86_64 \
		CIBW_BEFORE_BUILD="rm pyproject.toml && python -m pip install 'pip<21' 'numpy>=1.16.6' 'setuptools<63'" \
		CIBW_BUILD="$(CIBW_BUILD_LEGACY)" \
		"$$PYEXEC" -m cibuildwheel --platform linux --output-dir dist . && \
		PIP_NO_BUILD_ISOLATION=1 SKBUILD_DISABLE=1 \
		CIBW_BUILD="$(CIBW_BUILD_LEGACY)" \
		"$$PYEXEC" -m cibuildwheel --platform windows --output-dir dist .; \
	else \
		export TAR_OPTIONS="--warning=no-unknown-keyword"; \
		PIP_USE_PEP517=0 \
		PIP_NO_BUILD_ISOLATION=1 \
		SKBUILD_DISABLE=1 \
		CIBW_ARCHS_LINUX=x86_64 \
		CIBW_BEFORE_BUILD="rm pyproject.toml; python -m pip install 'pip<21' 'numpy>=1.16.6' 'setuptools<=44' 'wheel<0.35'" \
		CIBW_BUILD="$(CIBW_BUILD_LEGACY)" \
		"$$PYEXEC" -m cibuildwheel --platform linux --output-dir dist .; \
	fi

#: Create source distributions
sdist: build-deps
	@if $(PYTHON) -m build --help >/dev/null 2>&1; then \
		$(PYTHON) -m build --sdist --no-isolation --formats=zip,tar; \
	else \
		$(PYTHON) setup.py sdist; \
	fi

#: Install build dependencies
build-deps:
	$(PYTHON) -m ensurepip && \
	$(PYTHON) -m pip install build cibuildwheel

#: Install setup_requires dependencies
install-deps:
	$(PYTHON) -m ensurepip && \
	$(PYTHON) -m pip install wheel numpy setuptools \
		'meson; python_version >= "3.5"' ninja pytest scikit-build-core

#: Install library into current virtualenv
pip-install:
	$(PYTHON) -m pip install .

libphsh.cmake:
	cmake -S . -B build \
		-DPYTHON_INCLUDE_DIR="$$($(PYTHON) -c "import sysconfig; print(sysconfig.get_path('include'))")"  \
		-DPYTHON_LIBRARY="$$($(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")"
	FFLAGS="-fallow-argument-mismatch -std=f95" \
		cmake --build build

#: Build the f2py wrapped libphsh shared library within source tree
libphsh:
	$(PYTHON) setup.py build_ext --inplace

#: Perform checks
check: libphsh
	$(PYTHON) -m pytest tests/ --verbose

test: check

#: Remove any artifacts
clean:
	rm -rf venv* .venv-cibw112
	rm -rf build dist _skbuild \
		phaseshifts/lib/libphshmodule.c phaseshifts/lib/libphsh-f2pywrappers.f \
		phaseshifts/lib/libphsh*.so phaseshifts/lib/libphsh*.pyd
	rm -rf phaseshifts/lib/_native_build*.so phaseshifts/lib/_native_build*.pyd
	rm -rf phaseshifts/lib/phshift2007 phaseshifts/lib/phshift2007.zip
	rm -rf .cmake/ CMakeInit.txt CMakeCache.txt .skbuild-info.json CMakeFiles/

#: Build docker image
docker:
	export DOCKER_TAG="$${DOCKER_TAG:-$$($(PYTHON) -c 'import phaseshifts; print(phaseshifts.__version__)')}" && \
	$(DOCKER) build \
		-t "ghcr.io/liam-deacon/phaseshifts:$$DOCKER_TAG" \
		-f dockerfiles/phaseshifts-phsh.dockerfile .

PHSHIFT2007_BUILD_ROOT ?= phaseshifts/lib/phshift2007

#: Obtain a copy of the phshift2007 package
phaseshifts/lib/phshift2007.zip:
	curl -sL -o "$@" "https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip" || (rm -f "$@" && false)

#: Extract the phshift2007 package
$(PHSHIFT2007_BUILD_ROOT): phaseshifts/lib/phshift2007.zip
	unzip -o -d "$@" "$<"

$(PHSHIFT2007_BUILD_ROOT)/psprog.ab3: $(PHSHIFT2007_BUILD_ROOT)

$(PHSHIFT2007_BUILD_ROOT)/psprog.ab4: $(PHSHIFT2007_BUILD_ROOT)

# NOTE: Using awk is crude way of splitting the file as it looses comment lines before the program statement

#: Obtain phsh0.for by splitting psprog.ab3 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh0.for: $(PHSHIFT2007_BUILD_ROOT)/psprog.ab3
	awk -v dir="$(@D)/" '/C  program / {filename = dir tolower($$3)}; filename && NF {if (prev) print prev >filename; print $$0 >filename}' "$<"

#: Obtain phsh1.for by splitting psprog.ab3 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh1.for: $(PHSHIFT2007_BUILD_ROOT)/phsh0.for

#: Obtain phsh2cav.for by splitting psprog.ab4 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh2cav.for: $(PHSHIFT2007_BUILD_ROOT)/psprog.ab4
	awk -v dir="$(@D)/" '/C  program / {filename = dir tolower($$3)}; filename && NF {if (prev) print prev >filename; print $$0 >filename}' "$<"

#: Obtain phsh2rel.for by splitting psprog.ab4 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh2rel.for: $(PHSHIFT2007_BUILD_ROOT)/phsh2cav.for

#: Obtain phsh2wil.for by splitting psprog.ab4 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh2wil.for: $(PHSHIFT2007_BUILD_ROOT)/phsh2cav.for

#: Obtain phsh3.for by splitting psprog.ab4 into separate files
$(PHSHIFT2007_BUILD_ROOT)/phsh3.for: $(PHSHIFT2007_BUILD_ROOT)/phsh2cav.for

# GFortran compiler flags
GFORTRAN_FLAGS ?= -Wall -Wno-unused-label -Wno-tabs -Wno-unused-variable -Wno-unused-dummy-argument -fcheck=bounds -frecursive -std=legacy -pie
GFORTRAN_FLAGS += -static-libgcc -static-libgfortran


GFORTRAN ?= $(shell command -v gfortran 2>/dev/null || echo gfortran)

#: Build the hartfock program
bin/phsh0: $(PHSHIFT2007_BUILD_ROOT)/phsh0.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

#: Build the cavpot program
bin/phsh1: $(PHSHIFT2007_BUILD_ROOT)/phsh1.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

#: Build the william's phase shift program
bin/phsh2wil: $(PHSHIFT2007_BUILD_ROOT)/phsh2wil.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

#: Build the CAVLEED phase shift program
bin/phsh2cav: $(PHSHIFT2007_BUILD_ROOT)/phsh2cav.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

#: Build the relativistic phase shift program
bin/phsh2rel: $(PHSHIFT2007_BUILD_ROOT)/phsh2rel.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

#: Build the conphas program
bin/phsh3: $(PHSHIFT2007_BUILD_ROOT)/phsh3.for
	@mkdir -p "$(@D)"
	$(GFORTRAN) $(GFORTRAN_FLAGS) -o "$@" "$<"

.PHONY: phshift2007
#: Build the phshift2007 package (meta target)
phshift2007: bin/phsh0 bin/phsh1 bin/phsh2wil bin/phsh2cav bin/phsh2rel bin/phsh3

#: Install the phshift2007 programs
install: phshift2007
	install -m 755 bin/phsh0 bin/phsh1 bin/phsh2wil bin/phsh2cav bin/phsh2rel bin/phsh3 "$(PREFIX)/bin"

#: Create a virtual environment and install the package in editable mode
venv:
	@echo "Creating virtual environment..."
	$(PYTHON) -m venv venv
	source venv/bin/activate && \
	python -m ensurepip && \
	python -m pip install uv && \
	uv pip install --upgrade pip setuptools wheel && \
	uv pip install -e '.[dev,test,doc]'
	@echo "Virtual environment created in 'venv/'"
	@echo "Activate it with: source venv/bin/activate"

#: Build WebAssembly version (requires Emscripten + f2c)
.PHONY: wasm wasm-clean wasm-serve
wasm:
	@echo "Building WebAssembly version..."
	cd wasm && ./build.sh

#: Clean WebAssembly build artifacts
wasm-clean:
	rm -rf wasm/build wasm/dist

#: Serve WASM web interface locally
wasm-serve:
	@echo "Starting local server at http://localhost:8080/web/"
	cd wasm && $(PYTHON) -m http.server 8080
