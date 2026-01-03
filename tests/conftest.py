import glob
import os
import subprocess  # nosec
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
LIB_DIR = os.path.join(PROJECT_ROOT, "phaseshifts", "lib")
FORTRAN_COVERAGE_ENV = "PHASESHIFTS_FORTRAN_COVERAGE"
FORTRAN_COVERAGE_REPORT = os.path.join(PROJECT_ROOT, "fortran-coverage.xml")


def _has_local_compiled_lib():
    patterns = ("libphsh*.so", "libphsh*.pyd", "libphsh*.dll")
    return any(glob.glob(os.path.join(LIB_DIR, pattern)) for pattern in patterns)


use_source = os.environ.get("PHASESHIFTS_TEST_USE_SOURCE")

# Prefer the installed package (with compiled extension) unless explicitly
# told to use the source tree or a local compiled lib already exists.
if use_source or _has_local_compiled_lib():
    sys.path.insert(0, PROJECT_ROOT)
else:
    # Ensure the repo root does not shadow the installed wheel
    sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.curdir) != PROJECT_ROOT]


def pytest_addoption(parser):
    parser.addoption(
        "--fortran-coverage",
        action="store_true",
        help=(
            "Run gcovr to produce Fortran coverage (fortran-coverage.xml). "
            "Requires gfortran to be built with coverage flags."
        ),
    )


def _write_terminal(config, message):
    reporter = config.pluginmanager.get_plugin("terminalreporter")
    if reporter:
        reporter.write_line(message)


def _walk_gcda_files(paths):
    gcda_files = []
    for root in paths:
        if not os.path.isdir(root):
            continue
        for dirpath, _, filenames in os.walk(root):
            for name in filenames:
                if name.endswith(".gcda"):
                    gcda_files.append(os.path.join(dirpath, name))
    return gcda_files


def _iter_search_roots():
    for candidate in ("build", "_skbuild"):
        path = os.path.join(PROJECT_ROOT, candidate)
        if os.path.isdir(path):
            yield path


def _should_collect_fortran_coverage(config):
    return bool(config.getoption("--fortran-coverage") or os.environ.get(FORTRAN_COVERAGE_ENV))


def pytest_sessionfinish(session, exitstatus):
    if not _should_collect_fortran_coverage(session.config):
        return
    if exitstatus != 0:
        _write_terminal(
            session.config,
            "[fortran-coverage] Skipping gcovr because tests failed.",
        )
        return

    gcovr_exe = None
    try:
        import shutil

        gcovr_exe = shutil.which("gcovr") if hasattr(shutil, "which") else None
    except Exception:
        gcovr_exe = None

    if not gcovr_exe:
        _write_terminal(
            session.config,
            "[fortran-coverage] gcovr not installed; skipping Fortran coverage report.",
        )
        return

    search_roots = list(_iter_search_roots())
    if not search_roots:
        search_roots = [PROJECT_ROOT]

    gcda_files = _walk_gcda_files(search_roots)
    if not gcda_files:
        _write_terminal(
            session.config,
            "[fortran-coverage] No gcov data (.gcda) found. " "Ensure the build used PHASESHIFTS_FORTRAN_COVERAGE=1.",
        )
        return

    cmd = [
        gcovr_exe,
        "--root",
        PROJECT_ROOT,
        "--filter",
        os.path.join(PROJECT_ROOT, "phaseshifts", "lib"),
        "--xml",
        "--output",
        FORTRAN_COVERAGE_REPORT,
        "--gcov-executable",
        os.environ.get("GCOV", "gcov"),
        "--exclude-unreachable-branches",
        "--exclude-throw-branches",
    ]
    for root in search_roots:
        cmd.extend(["--search-path", root])

    result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # nosec
    stdout, stderr = result.communicate()

    if result.returncode != 0:
        _write_terminal(
            session.config,
            "[fortran-coverage] gcovr exited with code {code}. stdout:\n{stdout}\nstderr:\n{stderr}".format(
                code=result.returncode,
                stdout=stdout.decode("utf-8", "ignore"),
                stderr=stderr.decode("utf-8", "ignore"),
            ),
        )
        return

    summary = stdout.decode("utf-8", "ignore").strip()
    if summary:
        _write_terminal(
            session.config,
            "[fortran-coverage] gcovr summary:\n{summary}".format(summary=summary),
        )
    _write_terminal(
        session.config,
        "[fortran-coverage] Wrote {path}".format(path=FORTRAN_COVERAGE_REPORT),
    )
