"""Tests for WASM calculator web interface files."""

import os
import re

import pytest

WASM_DIR = os.path.join(os.path.dirname(__file__), "..", "wasm")
WEB_DIR = os.path.join(WASM_DIR, "web")
DIST_DIR = os.path.join(WASM_DIR, "dist")
SHARED_DIR = os.path.join(WASM_DIR, "shared")


class TestCalculatorWebFiles:
    """Tests for the calculator web interface files."""

    def test_index_html_exists(self):
        """Verify index.html exists."""
        index_path = os.path.join(WEB_DIR, "index.html")
        assert os.path.isfile(index_path), "wasm/web/index.html not found"

    def test_app_js_exists(self):
        """Verify app.js exists."""
        app_path = os.path.join(WEB_DIR, "app.js")
        assert os.path.isfile(app_path), "wasm/web/app.js not found"

    def test_style_css_exists(self):
        """Verify style.css exists."""
        style_path = os.path.join(WEB_DIR, "style.css")
        assert os.path.isfile(style_path), "wasm/web/style.css not found"

    def test_elements_js_exists(self):
        """Verify shared/elements.js exists."""
        elements_path = os.path.join(SHARED_DIR, "elements.js")
        assert os.path.isfile(elements_path), "wasm/shared/elements.js not found"


class TestCalculatorPaths:
    """Tests for correct path references in calculator files."""

    def test_index_html_script_paths(self):
        """Verify index.html references correct script paths for deployment."""
        index_path = os.path.join(WEB_DIR, "index.html")
        with open(index_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Should reference dist/phaseshifts.js (relative to calculator/)
        assert (
            'src="dist/phaseshifts.js"' in content or "src='dist/phaseshifts.js'" in content
        ), "index.html should reference dist/phaseshifts.js (not ../dist/)"

        # Should reference app.js (same directory)
        assert 'src="app.js"' in content or "src='app.js'" in content, "index.html should reference app.js"

    def test_app_js_import_paths(self):
        """Verify app.js uses correct import paths for deployment."""
        app_path = os.path.join(WEB_DIR, "app.js")
        with open(app_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Should import from ./shared/elements.js (relative to calculator/)
        # Match various import patterns
        import_pattern = r"from\s+['\"]\.\/shared\/elements\.js['\"]"
        assert re.search(import_pattern, content), "app.js should import from './shared/elements.js' (not '../shared/')"


class TestCalculatorSymlinks:
    """Tests for development symlinks (optional, for local dev convenience)."""

    @pytest.mark.skipif(
        not os.path.islink(os.path.join(WEB_DIR, "shared")), reason="Symlinks not set up (optional for local dev)"
    )
    def test_shared_symlink_target(self):
        """Verify shared symlink points to correct location."""
        shared_link = os.path.join(WEB_DIR, "shared")
        if os.path.islink(shared_link):
            target = os.readlink(shared_link)
            assert "shared" in target, f"shared symlink should point to ../shared, got: {target}"

    @pytest.mark.skipif(
        not os.path.islink(os.path.join(WEB_DIR, "dist")), reason="Symlinks not set up (optional for local dev)"
    )
    def test_dist_symlink_target(self):
        """Verify dist symlink points to correct location."""
        dist_link = os.path.join(WEB_DIR, "dist")
        if os.path.islink(dist_link):
            target = os.readlink(dist_link)
            assert "dist" in target, f"dist symlink should point to ../dist, got: {target}"


class TestWasmBuildDirectory:
    """Tests for WASM build output directory."""

    def test_dist_directory_exists_or_symlink(self):
        """Verify dist directory exists or is accessible via symlink."""
        # In the worktree, dist might not exist yet, but the symlink in web/ should work
        web_dist = os.path.join(WEB_DIR, "dist")
        # Either the actual dist dir exists, or the symlink in web/ exists
        assert os.path.isdir(DIST_DIR) or os.path.islink(
            web_dist
        ), "wasm/dist directory should exist or wasm/web/dist symlink should be present"

    @pytest.mark.skip(reason="WASM build is optional and may not be present locally")
    def test_wasm_binary_exists(self):
        """Verify WASM binary exists after build."""
        wasm_path = os.path.join(DIST_DIR, "phaseshifts.wasm")
        assert os.path.isfile(wasm_path), "wasm/dist/phaseshifts.wasm not found (run build.sh first)"
