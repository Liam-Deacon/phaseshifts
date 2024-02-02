# Makefile for assisting with some common development activities

PYTHON ?= $(shell pyenv which python 2>/dev/null || echo python)
DOCKER ?= $(shell command -v docker 2>/dev/null || docker)

PREFIX ?= /usr/local

.PHONY: build-deps cbuildwheel check clean install install-deps libphsh sdist test wheel

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
pip-install:
	$(PYTHON) -m pip install .

#: Build the project using cmake
cmake:
	(test -d build || $(MAKE) cmake.release) && cmake --build build

#: Build the project using cmake in release mode
cmake.release:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
	$(MAKE) cmake

#: Build the project using cmake in debug mode
cmake.debug:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
	$(MAKE) cmake

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
	rm -rf phaseshifts/lib/phshift2007 phaseshifts/lib/phshift2007.zip
	find phaseshifts -name '*.so' | xargs rm -f
	rm -rf htmlcov coverage.xml
	find phaseshifts -type d -name __pycache__ | xargs rm -rf
	find -name '*.pyc' -delete

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
