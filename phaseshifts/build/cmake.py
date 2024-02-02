"""Module for driving cmake builds from python as part of setup process."""
import os
import platform
import subprocess
import sys

import setuptools
import setuptools.command.build_ext

try:
    import cmake  # noqa
except ModuleNotFoundError:
    cmake = None


class CMakeBuild(setuptools.command.build_ext.build_ext):
    """Custom build_ext command class for driving CMake."""

    _CMAKE_EXE = "cmake" + (".exe" if sys.platform == "win32" else "")
    CMAKE = _CMAKE_EXE if not cmake else os.path.join(cmake.CMAKE_BIN_DIR, _CMAKE_EXE)

    def run(self):
        try:
            subprocess.check_output([self.CMAKE, "--version"])
        except subprocess.CalledProcessError as err:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from err

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

        if platform.system() == "Windows":
            if sys.maxsize > 2**32:
                # We are on a Win 64-bit system
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            # These flags are passed as-is to the underlying build tool (make)
            # build_args += ["--", "-j2"]
            pass

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO="{}"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        build_dir = os.path.abspath("build")
        os.makedirs(build_dir, exist_ok=True)

        try:
            print("Running: " + " ".join([self.CMAKE, "-S", os.path.dirname(build_dir), "-B", build_dir] + cmake_args))
            subprocess.check_call(
                [self.CMAKE, "-S", os.path.dirname(build_dir), "-B", build_dir] + cmake_args,
                env=env,
            )
            subprocess.check_call([self.CMAKE, "--build", build_dir, "--verbose"] + build_args, env=env)
        except subprocess.CalledProcessError as err:
            raise RuntimeError(f"CMake error: {err}") from err


if __name__ == "__main__":
    from distutils.dist import Distribution
    CMakeBuild(dist=Distribution()).build_extension()