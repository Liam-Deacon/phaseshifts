# Makefile for assisting with some common development activities

PYTHON ?= $(shell pyenv which python 2>/dev/null || echo python)
DOCKER ?= $(shell command -v docker 2>/dev/null || docker)

.PHONY: build-deps cbuildwheel check install install-deps libphsh sdist test wheel

#: Quickly generate binary wheel
wheel: install-deps
	$(PYTHON) setup.py build bdist_wheel

#: Meta target, will attempt to build all it can
build: sdist
	if command -v docker 1>/dev/null; then $(MAKE) cibuildwheel; else $(MAKE) wheel; fi 

#: Build a matrix of wheels for different OSs and CPU archs
cibuildwheel: build-deps $(DOCKER)
	$(PYTHON) -m cibuildwheel --platform=auto --output-dir=dist .

#: Create source distributions
sdist: build-deps
	$(PYTHON) build --sdist --no-isolation --formats=zip,tar

#: Install build dependencies
build-deps:
	$(PYTHON) -m pip install build cibuildwheel

#: Install setup_requires dependencies
install-deps:
	$(PYTHON) -m pip install wheel numpy setuptools \
		'meson; python_version >= "3.5"' ninja pytest scikit-build

#: Install library into current virtualenv
install:
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
	rm -rf build dist _skbuild \
		phaseshifts/lib/libphshmodule.c phaseshifts/lib/libphsh-f2pywrappers.f \
		phaseshifts/lib/libphsh*.so phaseshifts/lib/libphsh*.pyd

#: Build docker image
docker:
	$(DOCKER) build . -t "ghcr.io/Liam-Deacon/phaseshifts:$${DOCKER_TAG:-$$($(PYTHON) -c 'import phaseshifts; print(phaseshifts.__version__)')}"
