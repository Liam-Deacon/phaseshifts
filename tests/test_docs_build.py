import os
import subprocess
import pytest

DOCS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs")
HTML_INDEX = os.path.join(DOCS_DIR, "_build", "html", "index.html")


def sphinx_installed():
    try:
        subprocess.run(
            ["sphinx-build", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


@pytest.mark.skipif(not sphinx_installed(), reason="Sphinx is not installed")
def test_docs_html_build():
    result = subprocess.run(
        ["make", "html"], cwd=DOCS_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    assert result.returncode == 0, f"make html failed: {result.stderr.decode()}"
    assert os.path.isfile(HTML_INDEX), "HTML index file not found after build"
