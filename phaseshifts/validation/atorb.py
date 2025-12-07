"""
Validation helpers for atorb input files.
"""

import os
from typing import List

try:
    from pydantic import BaseModel, ValidationError

    try:
        from pydantic import model_validator
    except ImportError:  # pragma: no cover - pydantic v1 fallback
        from pydantic import root_validator as model_validator  # type: ignore
except Exception:  # pragma: no cover - pydantic unavailable
    class ValidationError(ValueError):
        """Lightweight stand-in when pydantic is unavailable."""

    class BaseModel(object):
        """Minimal shim to keep type signatures consistent without pydantic."""

        def __init__(self, **data):
            for key, value in data.items():
                setattr(self, key, value)

        @classmethod
        def model_validate(cls, data):
            return cls(**data)

        def model_dump(self):
            return self.__dict__

    def model_validator(*args, **kwargs):  # type: ignore
        def decorator(func):
            return func

        return decorator


def coerce_model(model_cls, data):
    """Helper to support pydantic v1/v2 or fallback shim."""
    if hasattr(model_cls, "model_validate"):
        return model_cls.model_validate(data)  # pydantic v2 or shim
    if hasattr(model_cls, "parse_obj"):
        return model_cls.parse_obj(data)  # pragma: no cover - pydantic v1
    return model_cls(**data)  # pragma: no cover - shim safety net


class AtorbElectron(BaseModel):
    """Model describing a single orbital line in an atorb input file."""

    n: int
    l: int
    m: int
    j: float
    s: int
    occ: float

    def ensure_valid(self):
        """Run lightweight semantic checks that are independent of pydantic."""
        if self.n <= 0:
            raise ValidationError("n must be positive for an orbital entry")
        if self.l < 0 or self.l >= self.n:
            raise ValidationError("l must satisfy 0 <= l < n")
        if self.s not in (0, 1):
            raise ValidationError("s must be 0 or 1 for spin multiplicity")
        if self.occ < 0:
            raise ValidationError("occ must be non-negative")
        if abs(self.j * 2) % 1:
            raise ValidationError("j must be a multiple of 0.5")
        return self


class AtorbInputModel(BaseModel):
    """Structured representation of an atorb input file."""

    z: int
    nr: int
    rel: int
    method: str
    relic: float
    nlevels: int
    mixing_scf: float
    eigen_tol: float
    ech: int
    orbitals: List[AtorbElectron]
    output: str
    header: str = ""

    def ensure_valid(self):
        """Perform semantic validation independent of pydantic internals."""
        if self.z <= 0:
            raise ValidationError("Atomic number (Z) must be positive")
        if self.nr <= 0:
            raise ValidationError("Radial grid points (NR) must be positive")
        if self.rel not in (0, 1):
            raise ValidationError("rel must be 0 (non-relativistic) or 1 (relativistic)")
        if self.nlevels != len(self.orbitals):
            raise ValidationError(
                "nlevels ({0}) does not match number of orbital lines ({1})".format(
                    self.nlevels, len(self.orbitals)
                )
            )
        cleaned_orbitals = []
        for orbital in self.orbitals:
            if isinstance(orbital, dict):
                orbital = coerce_model(AtorbElectron, orbital)
            orbital.ensure_valid()
            cleaned_orbitals.append(orbital)
        self.orbitals = cleaned_orbitals
        if not str(self.output).strip():
            raise ValidationError("output filename must be non-empty")
        return self


def _validate_comment_spacing(raw_line, line_no):
    """Ensure inline comments are separated by whitespace for Fortran list reads."""
    if "!" not in raw_line:
        return
    bang_idx = raw_line.index("!")
    if bang_idx == 0:
        return
    if not raw_line[bang_idx - 1].isspace():
        raise ValueError(
            "Line {0} must contain whitespace before '!' to avoid breaking Fortran list input: {1}".format(
                line_no, raw_line.rstrip()
            )
        )


def _clean_atorb_lines(filename):
    """Yield non-empty, non-comment lines with comment spacing validated."""
    cleaned = []
    with open(filename, "r") as handle:
        for idx, raw_line in enumerate(handle.readlines()):
            _validate_comment_spacing(raw_line, idx + 1)
            stripped = raw_line.rstrip("\n").rstrip("\r")
            if not stripped.strip():
                continue
            if stripped.lstrip().startswith("C"):
                continue
            stripped = stripped.split("!")[0].split("#")[0].strip()
            if stripped:
                cleaned.append(stripped)
    return cleaned


def format_orbital_line(orbital):
    """Return a consistently spaced orbital entry string."""
    return "{0:>2d} {1:d} {2:d} {3:>5.1f} {4:d} {5:>14.8f}".format(
        int(orbital.n),
        int(orbital.l),
        int(orbital.m),
        float(orbital.j),
        int(orbital.s),
        float(orbital.occ),
    )


def _parse_atorb_file(filename):
    """Parse a user-provided atorb input file into a structured model."""
    lines = _clean_atorb_lines(filename)
    if not lines:
        raise ValueError("atorb input file {0} is empty".format(filename))

    cursor = 0
    if lines[cursor].lower() != "i":
        raise ValueError("atorb input must start with an 'i' line")
    cursor += 1

    try:
        z, nr = [int(tok) for tok in lines[cursor].split()[:2]]
    except (ValueError, IndexError):
        raise ValueError("Unable to parse Z and NR from line: {0}".format(lines[cursor]))
    cursor += 1

    if lines[cursor].lower() != "d":
        raise ValueError("Expected 'd' line specifying relativistic flag")
    cursor += 1

    try:
        rel = int(lines[cursor].split()[0])
    except (ValueError, IndexError):
        raise ValueError("Unable to parse relativistic flag from line: {0}".format(lines[cursor]))
    cursor += 1

    if lines[cursor].lower() != "x":
        raise ValueError("Expected 'x' line specifying exchange-correlation method")
    cursor += 1
    method = lines[cursor].split()[0]
    cursor += 1

    if lines[cursor].lower() != "a":
        raise ValueError("Expected 'a' line containing SCF parameters")
    cursor += 1
    try:
        relic, nlevels, mixing_scf, eigen_tol, ech = lines[cursor].split()[:5]
        relic = float(relic)
        nlevels = int(nlevels)
        mixing_scf = float(mixing_scf)
        eigen_tol = float(eigen_tol)
        ech = int(ech)
    except (ValueError, IndexError):
        raise ValueError("Unable to parse SCF parameters from line: {0}".format(lines[cursor]))
    cursor += 1

    orbitals = []
    for i in range(nlevels):
        try:
            entry = lines[cursor + i].split()
            if len(entry) < 6:
                raise ValueError
            (n, l, m, j, s, occ) = entry[:6]
            orbitals.append(
                {
                    "n": int(n),
                    "l": int(l),
                    "m": int(m),
                    "j": float(j),
                    "s": int(s),
                    "occ": float(occ),
                }
            )
        except Exception:
            raise ValueError(
                "Unable to parse orbital entry on line {0}: {1}".format(
                    i + cursor + 1, lines[cursor + i]
                )
            )
    cursor += nlevels

    if lines[cursor].lower() != "w":
        raise ValueError("Expected 'w' line specifying output filename")
    cursor += 1
    if cursor >= len(lines):
        raise ValueError("Output filename missing after 'w' line")
    output = lines[cursor].split()[0]

    data = {
        "z": z,
        "nr": nr,
        "rel": rel,
        "method": method,
        "relic": relic,
        "nlevels": nlevels,
        "mixing_scf": mixing_scf,
        "eigen_tol": eigen_tol,
        "ech": ech,
        "orbitals": orbitals,
        "output": output,
        "header": "",
    }

    model = coerce_model(AtorbInputModel, data)
    model.ensure_valid()
    return model


def validate_atorb_file(filename):
    """Validate a generated or user-supplied atorb file before invoking Fortran."""
    if not os.path.isfile(filename):
        raise IOError("atorb input file not found: {0}".format(filename))
    return _parse_atorb_file(filename)


def render_atorb_file(model, header, filename):
    """Render a validated AtorbInputModel to disk using safe spacing."""
    lines = []
    lines.append("C".ljust(70, "*"))
    lines.append("C {0}".format(header.strip() if header else "atorb input file"))
    lines.append("C".ljust(70, "*"))
    lines.append("i")
    lines.append(
        "{0} {1}".format(model.z, model.nr).ljust(30, " ")
        + " ! Z NR (number of points in radial grid)"
    )
    lines.append("d")
    lines.append("{0}".format(model.rel).ljust(30) + " ! 1=rel, 0=n.r.")
    lines.append("x")
    lines.append(
        "{0}".format(model.method).ljust(30)
        + " ! 0.d0=HF, 1.d0=LDA, -alfa = xalfa..."
    )
    lines.append("a")
    lines.append(
        "{0} {1} {2} {3} {4}".format(
            model.relic, model.nlevels, model.mixing_scf, model.eigen_tol, model.ech
        ).ljust(30)
        + " ! relic,levels,mixing SCF, eigen. tol,for ech."
    )
    for orbital in model.orbitals:
        lines.append(
            "{0}  ! n, l, l, -j, <1>, occupation".format(
                format_orbital_line(orbital)
            )
        )
    lines.append("w")
    lines.append("{0}".format(model.output))
    lines.append("q")

    with open(filename, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return filename
