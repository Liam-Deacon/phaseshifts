// eslint-disable-next-line -- Browser ESM uses relative paths for local modules.
import { elements } from './elements.js';

/**
 * Default parameters for phase shift input generation.
 * @type {{muffinTinRadius: number, energyMin: number, energyMax: number, energyStep: number, lmax: number}}
 */
const DEFAULT_PHSH_PARAMS = Object.freeze({
  muffinTinRadius: 2.5,
  energyMin: 1.0,
  energyMax: 12.0,
  energyStep: 0.25,
  lmax: 10,
});

/**
 * Format a float as a fixed-width field.
 * @param {number} value - Value to format.
 * @param {number} width - Field width.
 * @param {number} precision - Decimal precision.
 * @returns {string} Fixed-width string.
 */
function formatFloat(value, width, precision) {
  const fixed = Number(value).toFixed(precision);
  return fixed.padStart(width, ' ');
}

/**
 * Format an integer as a fixed-width field.
 * @param {number} value - Value to format.
 * @param {number} width - Field width.
 * @returns {string} Fixed-width string.
 */
function formatInt(value, width) {
  return String(Math.trunc(value)).padStart(width, ' ');
}

/**
 * Format a single orbital line for atorb input.
 * @param {Object} orbital - Orbital descriptor.
 * @param {number} orbital.n - Principal quantum number.
 * @param {number} orbital.l - Angular momentum.
 * @param {number} orbital.m - Magnetic quantum number.
 * @param {number} orbital.j - Total angular momentum.
 * @param {number} orbital.s - Spin.
 * @param {number} orbital.occ - Occupancy.
 * @returns {string} Fixed-width orbital line.
 */
function formatOrbitalLine(orbital) {
  const n = formatInt(orbital.n, 2);
  const l = String(Math.trunc(orbital.l));
  const m = String(Math.trunc(orbital.m));
  const j = formatFloat(orbital.j, 5, 1);
  const s = String(Math.trunc(orbital.s));
  const occ = formatFloat(orbital.occ, 14, 8);
  return `${n} ${l} ${m} ${j} ${s} ${occ}`;
}

/**
 * Create a default electron configuration up to 7p.
 * @param {number} atomicNumber - Atomic number (Z).
 * @returns {Array<{n: number, l: number, occupation: number}>} Orbital list.
 */
function buildDefaultOrbitals(atomicNumber) {
  const shells = [
    { n: 1, l: 0, capacity: 2 }, // 1s
    { n: 2, l: 0, capacity: 2 }, // 2s
    { n: 2, l: 1, capacity: 6 }, // 2p
    { n: 3, l: 0, capacity: 2 }, // 3s
    { n: 3, l: 1, capacity: 6 }, // 3p
    { n: 4, l: 0, capacity: 2 }, // 4s
    { n: 3, l: 2, capacity: 10 }, // 3d
    { n: 4, l: 1, capacity: 6 }, // 4p
    { n: 5, l: 0, capacity: 2 }, // 5s
    { n: 4, l: 2, capacity: 10 }, // 4d
    { n: 5, l: 1, capacity: 6 }, // 5p
    { n: 6, l: 0, capacity: 2 }, // 6s
    { n: 4, l: 3, capacity: 14 }, // 4f
    { n: 5, l: 2, capacity: 10 }, // 5d
    { n: 6, l: 1, capacity: 6 }, // 6p
    { n: 7, l: 0, capacity: 2 }, // 7s
    { n: 5, l: 3, capacity: 14 }, // 5f
    { n: 6, l: 2, capacity: 10 }, // 6d
    { n: 7, l: 1, capacity: 6 }, // 7p
  ];

  const orbitals = [];
  let remaining = Math.max(0, Math.trunc(atomicNumber));

  for (const shell of shells) {
    if (remaining <= 0) break;
    const occ = Math.min(remaining, shell.capacity);
    orbitals.push({
      n: shell.n,
      l: shell.l,
      occupation: occ,
    });
    remaining -= occ;
  }

  return orbitals;
}

/**
 * Expand compact orbitals into Dirac-split components.
 * @param {Array<{n: number, l: number, occupation: number}>} orbitals - Input list.
 * @returns {Array<Object>} Expanded orbital list.
 */
function expandOrbitals(orbitals) {
  const expanded = [];
  for (const orbital of orbitals) {
    const n = orbital.n;
    const l = orbital.l;
    const occ = orbital.occupation;
    const m = l;
    const s = 1;

    if (l === 0) {
      expanded.push({ n, l, m, j: -0.5, s, occ });
      continue;
    }

    const denom = 4 * l + 2;
    const occLower = (occ * (2 * l)) / denom;
    const occUpper = (occ * (2 * l + 2)) / denom;

    expanded.push({ n, l, m, j: -(l - 0.5), s, occ: occLower });
    expanded.push({ n, l, m, j: -(l + 0.5), s, occ: occUpper });
  }

  return expanded;
}

/**
 * Normalize orbital input or generate defaults.
 * @param {Object} params - Input parameters.
 * @returns {Array<{n: number, l: number, occupation: number}>} Normalized orbitals.
 */
function normalizeOrbitals(params) {
  if (Array.isArray(params.orbitals) && params.orbitals.length > 0) {
    return params.orbitals.map((orbital) => ({
      n: orbital.n,
      l: orbital.l,
      occupation: orbital.occupation ?? orbital.occ ?? orbital.occTotal ?? 0,
    }));
  }

  return buildDefaultOrbitals(params.atomicNumber);
}

/**
 * Resolve an element symbol from an atomic number.
 * @param {number} atomicNumber - Atomic number.
 * @returns {string|null} Element symbol or null if not found.
 */
function getElementSymbol(atomicNumber) {
  const target = Number(atomicNumber);
  if (!Number.isFinite(target)) {
    return null;
  }

  for (const [symbol, number] of Object.entries(elements)) {
    if (number === target) {
      return symbol;
    }
  }

  return null;
}

/**
 * Build an atorb input file payload.
 * @param {Object} params - Input parameters.
 * @returns {string} Formatted atorb input.
 */
function buildAtorbInput(params) {
  const atomicNumber = params.atomicNumber;
  const elementSymbol = getElementSymbol(atomicNumber);
  const ngrid = params.ngrid === undefined ? 1000 : params.ngrid;
  const rel = params.relativity === undefined ? 1 : params.relativity;
  const method =
    params.method === undefined ? '0.d0' : String(params.method).trim();
  const relic = params.relic === undefined ? 0 : params.relic;
  const mixingScf = params.mixingScf === undefined ? 0.5 : params.mixingScf;
  const eigenTol =
    params.eigenTol === undefined ? '0.0005' : String(params.eigenTol);
  const ech = params.ech === undefined ? 100 : params.ech;
  const output =
    params.output === undefined && elementSymbol
      ? `at_${elementSymbol}.i`
      : params.output || 'at_output.i';
  const header =
    params.header || `atorb input file: atorb_${elementSymbol || 'X'}.txt.`;

  const orbitals = normalizeOrbitals(params);
  const expanded = expandOrbitals(orbitals);
  const nlevels = expanded.length;

  const lines = [];
  lines.push('C'.padEnd(70, '*'));
  lines.push(`C ${header}`);
  lines.push('C'.padEnd(70, '*'));
  lines.push('i');
  lines.push(
    `${atomicNumber} ${ngrid}`.padEnd(30, ' ') +
      ' ! Z NR (number of points in radial grid)',
  );
  lines.push('d');
  lines.push(`${rel}`.padEnd(30, ' ') + ' ! 1=rel, 0=n.r.');
  lines.push('x');
  lines.push(
    `${method}`.padEnd(30, ' ') + ' ! 0.d0=HF, 1.d0=LDA, -alfa = xalfa...',
  );
  lines.push('a');
  lines.push(
    `${relic} ${nlevels} ${mixingScf} ${eigenTol} ${ech}`.padEnd(30, ' ') +
      ' ! relic,levels,mixing SCF, eigen. tol,for ech.',
  );

  for (const orbital of expanded) {
    lines.push(
      `${formatOrbitalLine(orbital)}`.padEnd(30, ' ') +
        ' ! n, l, l, -j, <1>, occupation',
    );
  }

  lines.push('w');
  lines.push(output);
  lines.push('q');

  return lines.join('\n') + '\n';
}

/**
 * Normalize phase shift parameters by applying default values.
 * @param {Object} params - Input parameters.
 * @returns {Object} Parameters with defaults applied.
 */
function normalizePhshParams(params) {
  const safe = params || {};
  return {
    atomicNumber: safe.atomicNumber,
    muffinTinRadius:
      safe.muffinTinRadius === undefined
        ? DEFAULT_PHSH_PARAMS.muffinTinRadius
        : safe.muffinTinRadius,
    energyMin:
      safe.energyMin === undefined
        ? DEFAULT_PHSH_PARAMS.energyMin
        : safe.energyMin,
    energyMax:
      safe.energyMax === undefined
        ? DEFAULT_PHSH_PARAMS.energyMax
        : safe.energyMax,
    energyStep:
      safe.energyStep === undefined
        ? DEFAULT_PHSH_PARAMS.energyStep
        : safe.energyStep,
    lmax: safe.lmax === undefined ? DEFAULT_PHSH_PARAMS.lmax : safe.lmax,
  };
}

/**
 * Build a phsh input file payload.
 * @param {Object} params - Input parameters.
 * @returns {string} Formatted phsh input.
 */
function buildPhshInput(params) {
  const {
    atomicNumber,
    muffinTinRadius,
    energyMin,
    energyMax,
    energyStep,
    lmax,
  } = normalizePhshParams(params);

  const nEnergies = Math.floor((energyMax - energyMin) / energyStep) + 1;

  const fields = [
    formatFloat(atomicNumber, 6, 1),
    formatFloat(muffinTinRadius, 10, 4),
    formatFloat(energyMin, 10, 4),
    formatFloat(energyMax, 10, 4),
    formatFloat(energyStep, 10, 4),
    formatInt(lmax, 4),
    formatInt(nEnergies, 6),
  ];

  return `${fields.join(' ')}\n`;
}

export { DEFAULT_PHSH_PARAMS, buildAtorbInput, buildPhshInput, normalizePhshParams };
