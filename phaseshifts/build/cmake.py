"""Module for driving cmake builds from python as part of setup process."""

# pylint: disable=fixme

import os
import subprocess  # nosec: B404
import sys

import setuptools
import setuptools.command.build_ext

try:
    import cmake  # noqa
except ModuleNotFoundError:
    cmake = None  # pylint: disable=invalid-name


class CMakeBuild(setuptools.command.build_ext.build_ext):
    """Custom build_ext command class for driving CMake."""

    _CMAKE_EXE = "cmake" + (".exe" if sys.platform == "win32" else "")
    CMAKE = _CMAKE_EXE if not getattr(cmake, "CMAKE_BIN_DIR", None) else os.path.join(cmake.CMAKE_BIN_DIR, _CMAKE_EXE)

    def run(self):
        try:
            subprocess.check_output([self.CMAKE, "--version"])  # nosec: B603
        except subprocess.CalledProcessError:
            raise RuntimeError(  # noqa: B904  # TODO: Remove when we drop python2.7 support
                "CMake must be installed to build the following extensions: "
                + ", ".join(ext.name for ext in self.extensions)
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext=None):
        """Build extension via cmake.

        It assumes that any Extension.name has a target in CMakeLists.txt.
        Setting ext to None will build all targets (default).
        """
        target = getattr(ext, "name", None)
        cmake_args = list(filter(None, [] if not target else ["--target={}".format(target)]))
        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        env = os.environ.copy()

        build_dir = os.path.abspath("build")
        if not os.path.exists(build_dir):  # added for python 2.7 compatibility
            os.makedirs(build_dir, **({"exist_ok": True} if sys.version_info >= (3, 2) else {}))

        try:
            print("Running: " + " ".join([self.CMAKE, "-S", os.path.dirname(build_dir), "-B", build_dir] + cmake_args))
            subprocess.check_call(  # nosec: B603
                [self.CMAKE, "-S", os.path.dirname(build_dir), "-B", build_dir] + cmake_args,
                env=env,
            )
            subprocess.check_call([self.CMAKE, "--build", build_dir, "--verbose"] + build_args, env=env)  # nosec: B603
        except subprocess.CalledProcessError as err:
            raise RuntimeError(  # noqa: B904  # TODO: Remove when we drop python2.7 support
                "CMake error: {}".format(err)
            )


if __name__ == "__main__":
    from setuptools.dist import Distribution
    CMakeBuild(dist=Distribution()).build_extension()  # type: ignore
