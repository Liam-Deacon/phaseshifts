ARG PYTHON_VERSION=3.11
FROM python:${PYTHON_VERSION}-slim as base

RUN apt-get update \
  && apt-get install --yes --no-install-recommends gfortran make \
  && apt-get clean --yes \
  && rm -rf /var/lib/apt/lists/*

# install key dependencies
RUN pip install --no-cache-dir --upgrade pip && \
pip install --no-cache-dir numpy==1.26.3 setuptools==69.0.3 wheel==0.42.0

FROM base as builder

WORKDIR /build

# TODO: Remove the manual install_requires step when migrating to pyproject.toml PEP-517 builds (see #8)
RUN pip install --no-cache-dir meson scikit-build

COPY . /build/

# Download dependencies and build phaseshifts wheel
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /build/wheels .

# Extras for phaseshifts must be provided in list syntax, e.g. "[atorb,dev,test]"
ARG PHASESHIFTS_EXTRAS=""
FROM base as runtime

WORKDIR /data

COPY --from=builder /build/wheels /wheels

RUN pip install --no-cache-dir /wheels/*
RUN pip install --no-cache-dir /wheels/phaseshifts*.whl${PHASESHIFTS_EXTRAS}

HEALTHCHECK NONE

ENTRYPOINT [ "phsh.py" ]
