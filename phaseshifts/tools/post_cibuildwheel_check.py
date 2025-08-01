#!/usr/bin/env python3
"""
Comprehensive post-cibuildwheel wheel verification script.
Checks:
- Wheel(s) exist in wheelhouse
- For binary wheels: libphsh Fortran wrapper present
- Wheel is installable in a temp venv
- Main package and Fortran wrapper are importable

Usage:
    python post_cibuildwheel_check.py [--wheelhouse /path/to/wheelhouse]
"""
import argparse
import os
import sys
import zipfile
import tempfile
import shutil
import subprocess
import glob
import platform

# Configurable
FORTRAN_WRAPPER_BASENAME = "libphsh"
PACKAGE_IMPORT_NAME = "phaseshifts.lib.phsh"  # Update if needed


def is_binary_wheel(wheel_path):
    """Detect if wheel is binary (not pure Python) by filename or contents."""
    fname = os.path.basename(wheel_path)
    # If wheel name contains 'none-any', it's pure Python
    if "none-any.whl" in fname:
        return False
    # Otherwise, check for .so/.pyd/.dll inside wheel
    with zipfile.ZipFile(wheel_path, "r") as zf:
        for name in zf.namelist():
            if name.lower().endswith((".so", ".pyd", ".dll")):
                return True
    return False


def has_fortran_wrapper(wheel_path):
    """Check if wheel contains the compiled Fortran wrapper libphsh.*"""
    with zipfile.ZipFile(wheel_path, "r") as zf:
        for name in zf.namelist():
            base = os.path.basename(name).lower()
            if base.startswith(FORTRAN_WRAPPER_BASENAME) and base.endswith(
                (".so", ".pyd", ".dll")
            ):
                return True
    return False


def run_in_venv(wheel_path, import_names):
    """Create temp venv, install wheel, and try to import given names."""
    venv_dir = tempfile.mkdtemp(prefix="cibw_check_venv_")
    try:
        python_exe = shutil.which("python3") or shutil.which("python")
        if not python_exe:
            print("ERROR: No python executable found.")
            return False, "No python executable found"
        # Create venv
        subprocess.check_call([python_exe, "-m", "venv", venv_dir])
        vpy = os.path.join(
            venv_dir, "Scripts" if platform.system() == "Windows" else "bin", "python"
        )
        pip = os.path.join(
            venv_dir, "Scripts" if platform.system() == "Windows" else "bin", "pip"
        )
        # Upgrade pip
        subprocess.check_call(
            [vpy, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"]
        )
        # Install wheel
        subprocess.check_call([pip, "install", wheel_path])
        # Try imports
        for modname in import_names:
            code = f"import {modname}"
            try:
                subprocess.check_call([vpy, "-c", code])
            except subprocess.CalledProcessError:
                return False, f"Failed to import {modname}"
        return True, ""
    except Exception as e:
        return False, str(e)
    finally:
        shutil.rmtree(venv_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Post-cibuildwheel wheel verification script."
    )
    parser.add_argument(
        "--wheelhouse",
        type=str,
        default=None,
        help="Directory containing built wheels.",
    )
    args = parser.parse_args()

    search_dirs = []
    if args.wheelhouse:
        search_dirs = [args.wheelhouse]
    else:
        search_dirs = ["wheelhouse", "dist"]

    wheels = []
    for d in search_dirs:
        found = glob.glob(os.path.join(d, "*.whl"))
        wheels.extend([(w, d) for w in found])

    if not wheels:
        print(f"ERROR: No wheels found in {', '.join(search_dirs)}")
        sys.exit(1)
    print(f"Found {len(wheels)} wheel(s) in {', '.join(search_dirs)}.")
    all_ok = True
    for wheel, srcdir in wheels:
        print(f"\nChecking wheel: {os.path.basename(wheel)} (from {srcdir})")
        is_bin = is_binary_wheel(wheel)
        print(f"  Binary wheel: {'YES' if is_bin else 'NO'}")
        if is_bin:
            has_wrapper = has_fortran_wrapper(wheel)
            print(f"  Fortran wrapper present: {'YES' if has_wrapper else 'NO'}")
            if not has_wrapper:
                print(
                    f"ERROR: Fortran wrapper {FORTRAN_WRAPPER_BASENAME} not found in wheel."
                )
                all_ok = False
        # Install and import check
        import_names = [PACKAGE_IMPORT_NAME]
        if is_bin:
            # Try to import the wrapper directly if possible
            import_names.append(f"{PACKAGE_IMPORT_NAME}.{FORTRAN_WRAPPER_BASENAME}")
        ok, msg = run_in_venv(wheel, import_names)
        print(f"  Install & import test: {'PASS' if ok else 'FAIL'}")
        if not ok:
            print(f"ERROR: {msg}")
            all_ok = False
    print("\nSummary:")
    if all_ok:
        print("All checks passed.")
        sys.exit(0)
    else:
        print("Some checks failed.")
        sys.exit(2)


if __name__ == "__main__":
    main()
