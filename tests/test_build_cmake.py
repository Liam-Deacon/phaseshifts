"""Test building the extension using CMake."""
import pytest
from setuptools.dist import Distribution

from phaseshifts.build.cmake import CMakeBuild

@pytest.mark.build
@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.parametrize(
    ("extension",),
    [
        [None],
        ["phshift2007.zip"],
        ["libphsh"],
    ]
)
def test_cmake_build(extension):
    """Test use cmake build process for different cmake targets."""
    CMakeBuild(dist=Distribution()).build_extension(ext=extension)  # type: ignore


@pytest.mark.build
@pytest.mark.integration
def test_cmake_run():
    """Test cmake build run method, basically verifies that cmake is available."""
    builder = CMakeBuild(dist=Distribution())
    builder.extensions = []
    builder.run()
