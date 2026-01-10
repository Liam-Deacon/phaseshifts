"""
Tests for WebAssembly build infrastructure.

These tests verify:
1. Required build tools are properly detected
2. Build scripts are valid and executable
3. Web interface files are present and valid
4. (If WASM is built) The module loads correctly
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
WASM_DIR = PROJECT_ROOT / "wasm"
BUILD_SCRIPT = WASM_DIR / "build.sh"
DIST_DIR = WASM_DIR / "dist"
WEB_DIR = WASM_DIR / "web"
SRC_DIR = WASM_DIR / "src"


class TestWasmInfrastructure:
    """Tests for WASM build infrastructure files."""

    def test_wasm_directory_exists(self):
        """WASM directory should exist."""
        assert WASM_DIR.exists(), "wasm/ directory not found"

    def test_build_script_exists(self):
        """Build script should exist."""
        assert BUILD_SCRIPT.exists(), "wasm/build.sh not found"

    def test_build_script_is_executable(self):
        """Build script should be executable."""
        assert os.access(BUILD_SCRIPT, os.X_OK), "wasm/build.sh is not executable"

    def test_build_script_has_shebang(self):
        """Build script should have proper shebang."""
        with open(BUILD_SCRIPT, "r") as f:
            first_line = f.readline()
        assert first_line.startswith("#!/"), "build.sh should have a shebang line"
        assert "bash" in first_line, "build.sh should use bash"

    def test_readme_exists(self):
        """WASM README should exist."""
        readme = WASM_DIR / "README.md"
        assert readme.exists(), "wasm/README.md not found"

    def test_readme_has_content(self):
        """WASM README should have meaningful content."""
        readme = WASM_DIR / "README.md"
        content = readme.read_text(encoding="utf-8")
        assert len(content) > 500, "README.md seems too short"
        assert "WebAssembly" in content or "WASM" in content
        assert "Emscripten" in content


class TestWebInterface:
    """Tests for the browser interface files."""

    def test_web_directory_exists(self):
        """Web directory should exist."""
        assert WEB_DIR.exists(), "wasm/web/ directory not found"

    def test_index_html_exists(self):
        """index.html should exist."""
        index_html = WEB_DIR / "index.html"
        assert index_html.exists(), "wasm/web/index.html not found"

    def test_index_html_valid_structure(self):
        """index.html should have valid HTML structure."""
        index_html = WEB_DIR / "index.html"
        content = index_html.read_text(encoding="utf-8")
        content_lower = content.lower()

        # Check for essential HTML elements
        assert "<!doctype html>" in content_lower, "Missing DOCTYPE"
        assert "<html" in content_lower, "Missing html tag"
        assert "<head>" in content_lower, "Missing head tag"
        assert "<body>" in content_lower, "Missing body tag"
        assert "</html>" in content_lower, "Missing closing html tag"

    def test_index_html_has_required_elements(self):
        """index.html should have phase shift calculator elements."""
        index_html = WEB_DIR / "index.html"
        content = index_html.read_text(encoding="utf-8")

        # Check for key UI elements (updated for new tabbed structure builder UI)
        assert 'id="element-settings"' in content or 'id="element"' in content, "Missing element settings/selector"
        assert 'id="method"' in content, "Missing method selector"
        assert 'id="calculate-btn"' in content, "Missing calculate button"
        assert 'id="results-section"' in content, "Missing results section"

    def test_style_css_exists(self):
        """style.css should exist."""
        style_css = WEB_DIR / "style.css"
        assert style_css.exists(), "wasm/web/style.css not found"

    def test_style_css_has_content(self):
        """style.css should have CSS rules."""
        style_css = WEB_DIR / "style.css"
        content = style_css.read_text(encoding="utf-8")

        assert len(content) > 100, "style.css seems too short"
        assert "{" in content and "}" in content, "No CSS rules found"

    def test_app_js_exists(self):
        """app.js should exist."""
        app_js = WEB_DIR / "app.js"
        assert app_js.exists(), "wasm/web/app.js not found"

    def test_app_js_has_required_functions(self):
        """app.js should have required functions."""
        app_js = WEB_DIR / "app.js"
        content = app_js.read_text(encoding="utf-8")

        required_functions = [
            "runCalculation",
            "displayResults",
            "downloadResults",
            "loadPreset",
        ]

        for func in required_functions:
            assert func in content, f"Missing function: {func}"


class TestJavaScriptAPI:
    """Tests for the JavaScript API wrapper."""

    def test_src_directory_exists(self):
        """src directory should exist."""
        assert SRC_DIR.exists(), "wasm/src/ directory not found"

    def test_phaseshifts_js_exists(self):
        """phaseshifts.js API wrapper should exist."""
        api_js = SRC_DIR / "phaseshifts.js"
        assert api_js.exists(), "wasm/src/phaseshifts.js not found"

    def test_phaseshifts_js_exports(self):
        """phaseshifts.js should export required symbols."""
        api_js = SRC_DIR / "phaseshifts.js"
        content = api_js.read_text(encoding="utf-8")

        required_exports = [
            "PhaseShifts",
            "createPhaseShifts",
            "ELEMENTS",
        ]

        for export in required_exports:
            assert export in content, f"Missing export: {export}"

    def test_phaseshifts_js_has_methods(self):
        """PhaseShifts class should have required methods."""
        api_js = SRC_DIR / "phaseshifts.js"
        content = api_js.read_text(encoding="utf-8")

        required_methods = [
            "calculatePhaseShifts",
            "calculateRelativistic",
            "calculateCavity",
            "calculateWilliams",
            "init",
        ]

        for method in required_methods:
            assert method in content, f"Missing method: {method}"


class TestSharedModules:
    """Tests for shared WASM modules."""

    def test_shared_elements_exists(self):
        """Shared elements map should exist."""
        shared_elements = WASM_DIR / "shared" / "elements.js"
        assert shared_elements.exists(), "wasm/shared/elements.js not found"

    def test_web_app_imports_shared_elements(self):
        """web/app.js should import the shared elements map."""
        app_js = WEB_DIR / "app.js"
        content = app_js.read_text(encoding="utf-8")
        assert "./shared/elements.js" in content, "web/app.js should import ./shared/elements.js"

    def test_src_reexports_shared_elements(self):
        """src/elements.js should re-export the shared elements map."""
        src_elements = SRC_DIR / "elements.js"
        content = src_elements.read_text(encoding="utf-8")
        assert "../shared/elements.js" in content, "src/elements.js should re-export ../shared/elements.js"


class TestBuildTools:
    """Tests for build tool detection."""

    @pytest.mark.skipif(
        subprocess.run(["which", "emcc"], capture_output=True).returncode != 0,
        reason="Emscripten not installed",
    )
    def test_emcc_available(self):
        """emcc should be available when Emscripten is installed."""
        result = subprocess.run(["emcc", "--version"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "emcc" in result.stdout.lower()

    @pytest.mark.skipif(
        subprocess.run(["which", "f2c"], capture_output=True).returncode != 0,
        reason="f2c not installed",
    )
    def test_f2c_available(self):
        """f2c should be available when installed."""
        result = subprocess.run(["which", "f2c"], capture_output=True, text=True)
        assert result.returncode == 0

    def test_build_help_works(self):
        """Build script help should work."""
        cmd = [str(BUILD_SCRIPT), "--help"]
        if os.name == "nt":
            bash_path = shutil.which("bash")
            if not bash_path:
                pytest.skip("bash is not available on Windows")
            cmd = [bash_path, str(BUILD_SCRIPT), "--help"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=WASM_DIR)
        # Help should exit with 0
        assert result.returncode == 0
        assert "Usage" in result.stdout or "build" in result.stdout.lower()


class TestWasmBuild:
    """Tests for actual WASM build (only runs if tools are available)."""

    @pytest.fixture
    def has_build_tools(self):
        """Check if build tools are available."""
        emcc_available = subprocess.run(["which", "emcc"], capture_output=True).returncode == 0
        f2c_available = subprocess.run(["which", "f2c"], capture_output=True).returncode == 0
        return emcc_available and f2c_available

    @pytest.mark.skip(reason="Full WASM build requires Emscripten and f2c")
    def test_wasm_build_succeeds(self, has_build_tools):
        """WASM build should complete successfully."""
        if not has_build_tools:
            pytest.skip("Build tools not available")

        result = subprocess.run(
            [str(BUILD_SCRIPT), "--clean"],
            capture_output=True,
            text=True,
            cwd=WASM_DIR,
            timeout=300,  # 5 minute timeout
        )

        assert result.returncode == 0, f"Build failed: {result.stderr}"

        # Check output files exist
        assert (DIST_DIR / "phaseshifts.js").exists()
        assert (DIST_DIR / "phaseshifts.wasm").exists()

    @pytest.mark.skip(reason="Requires WASM build to be complete")
    def test_wasm_output_size_reasonable(self, has_build_tools):
        """WASM output should be a reasonable size."""
        if not has_build_tools:
            pytest.skip("Build tools not available")

        wasm_file = DIST_DIR / "phaseshifts.wasm"
        if not wasm_file.exists():
            pytest.skip("WASM not built")

        size = wasm_file.stat().st_size

        # Should be between 100KB and 50MB
        assert size > 100 * 1024, "WASM file too small"
        assert size < 50 * 1024 * 1024, "WASM file too large"


class TestFortranSource:
    """Tests that Fortran source is available for WASM build."""

    def test_libphsh_f_exists(self):
        """libphsh.f should exist."""
        fortran_src = PROJECT_ROOT / "phaseshifts" / "lib" / "libphsh.f"
        assert fortran_src.exists(), "phaseshifts/lib/libphsh.f not found"

    def test_libphsh_f90_exists(self):
        """libphsh.f90 should exist (free-form version)."""
        fortran_src = PROJECT_ROOT / "phaseshifts" / "lib" / "libphsh.f90"
        assert fortran_src.exists(), "phaseshifts/lib/libphsh.f90 not found"

    def test_fortran_has_required_subroutines(self):
        """Fortran source should have required subroutines."""
        fortran_src = PROJECT_ROOT / "phaseshifts" / "lib" / "libphsh.f"
        content = fortran_src.read_text()

        required_subroutines = [
            "hartfock",
            "phsh_cav",  # or phsh2cav
            "phsh_rel",  # or phsh_rel
            "phsh_wil",  # or phsh2wil
        ]

        content_lower = content.lower()
        for sub in required_subroutines:
            # Check for subroutine definition (may have different naming)
            assert sub.lower() in content_lower or sub.replace("_", "").lower() in content_lower, (
                f"Missing subroutine: {sub}"
            )
