# Makefile for assisting with some common development activities

PYTHON ?= $(shell pyenv which python 2>/dev/null || echo python)

.PHONY: wheel install-deps libphsh check test clean

#: Quickly generate binary wheel
wheel: install-deps
	$(PYTHON) setup.py build bdist_wheel

#: Install dependencies
install-deps:
	$(PYTHON) -m pip install wheel numpy setuptools 'meson; python_version >= "3.5"' ninja pytest scikit-build

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
	docker build . -t "phaseshifts:$${DOCKER_TAG:-$$($(PYTHON) -c 'import phaseshifts; print(phaseshifts.__version__)')}"
