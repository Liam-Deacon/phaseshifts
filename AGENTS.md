# AGENTS.md

## Build, Lint, and Test Commands

- **Install dependencies:** `pip install -r requirements.txt`
- **Best practices:**
  Always use a Python virtual environment (venv, conda, etc.) for development and installs.
  Prefer `uv` for managing virtual environments & packages.
- **Never:**
  Use `--break-system-packages` or install packages into the system Python unless absolutely necessary and explicitly confirmed by user. This prevents system breakage and ensures reproducible builds.
- **Build Fortran extension:** `make libphsh` or `python setup.py build_ext --inplace`
- **Run all tests:** `pytest tests/ test/ --verbose`
- **Run a single test:** `pytest tests/test_phshift2007.py::test_e2e`
- **Lint:** `flake8 . --max-line-length=127`
- **Build wheel:** `make wheel`
- **Build all:** `make build`
- **Clean artifacts:** `make clean`
- **Docker build:** `make docker`
- **Test coverage:** `pytest --cov=phaseshifts`

## Documentation Build Notes

- **Documentation is generated using Sphinx**. All source files are in the `docs/` directory.
- **Build HTML docs:**
  Run `make html` from the `docs/` directory. Output will be in `docs/_build/html`.
- **Build PDF docs:**
  Run `make latexpdf` from the `docs/` directory. Output will be in `docs/_build/latex`.
- **Clean build artifacts:**
  Run `make clean` in `docs/` to remove all generated files.
- **Other formats:**
  The `docs/Makefile` supports many targets (see `make help`), including `epub`, `man`, `text`, `json`, etc.
- **Prerequisites:**
  Sphinx must be installed (`pip install sphinx`). If `sphinx-build` is not found, follow the error message in the Makefile to install or set your PATH.
- **Troubleshooting:**
  If you see "The 'sphinx-build' command was not found", ensure Sphinx is installed and available in your PATH.
- **Advanced usage:**
  See `docs/Makefile` for all supported build targets and options.

## Documentation Tests

- **Test documentation build:**
  Run `pytest tests/test_docs_build.py` to verify that documentation builds successfully.
- The test will run `make html` in the `docs/` directory and check for the output in `docs/_build/html/index.html`.
- If Sphinx is not installed, the test will be skipped.
- If the build fails, the test will fail and print the error output.
- This helps ensure documentation is always buildable and up-to-date.

## Code Style Guidelines

- **Imports:** Use absolute imports; group stdlib, third-party, and local imports separately.
- **Formatting:** Follow PEP8; max line length 127 (flake8 config).
- **Types:** Use type hints where possible, especially for public APIs.
- **Naming:** Use snake_case for functions/variables, PascalCase for classes.
- **Error Handling:** Prefer explicit exceptions; use custom exceptions for CLI errors.
- **Tests:** Place in `tests/` or `test/` directories; use `pytest` or `unittest`.
- **Docstrings:** Use NumPy or Google style for functions/classes.
- **Compatibility:** Support Python 2.7+ and 3.5–3.11 (see README and requirements).
- **Fortran Extensions:** Ensure `libphsh` is built before running tests that require it. For fixed-format Fortran, use the `.f` extension and update all build references accordingly. For Fortran 90 free-format code, use `.f90`.
- **Pre-commit:** Use flake8 and pytest for linting and testing before commits.

## Fortran Code Notes

- The scientific core of phaseshifts is implemented in Fortran (see `phaseshifts/lib/libphsh.f`, `libphsh.f90`, and `.phsh.orig/phsh0.f`).
- These files provide atomic orbital, charge density, and pseudopotential calculations (Dirac-Fock, Hartree-Fock, relativistic corrections, etc.).
- Always ensure Fortran sources are present and compilable before running tests or builds. Use `make libphsh` or `python setup.py build_ext --inplace`.
- Python wrappers and tests depend on a working Fortran build; failures in Fortran compilation will break core functionality.
- For agentic workflows, automate Fortran build/test steps and document any scientific assumptions or changes in the Fortran code.
- Preserve and extend scientific comments in Fortran files—they are essential for future maintainers and agents.

## Agentic Coding Guidelines

- Always run lint and tests after code changes.
- Prefer small, incremental commits with clear messages.
- If modifying build/test infra, update this file.
- If Cursor or Copilot rules are added, include them here.
