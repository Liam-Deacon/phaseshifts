import { elements } from './elements.js';

function formatFloat(value, width, precision) {
  const fixed = Number(value).toFixed(precision);
  return fixed.padStart(width, ' ');
}

function formatInt(value, width) {
  return String(Math.trunc(value)).padStart(width, ' ');
}

function formatOrbitalLine(orbital) {
  const n = formatInt(orbital.n, 2);
  const l = String(Math.trunc(orbital.l));
  const m = String(Math.trunc(orbital.m));
  const j = formatFloat(orbital.j, 5, 1);
  const s = String(Math.trunc(orbital.s));
  const occ = formatFloat(orbital.occ, 14, 8);
  return `${n} ${l} ${m} ${j} ${s} ${occ}`;
}

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

function buildAtorbInput(params) {
  const atomicNumber = params.atomicNumber;
  const elementSymbol = Object.keys(elements).find(
    (key) => elements[key] === atomicNumber,
  );
  const ngrid = params.ngrid === undefined ? 1000 : params.ngrid;
  const rel = params.relativity === undefined ? 1 : params.relativity;
  const method =
    params.method === undefined ? '0.d0' : String(params.method).trim();
  const relic = params.relic === undefined ? 0 : params.relic;
  const mixingScf =
    params.mixingScf === undefined ? 0.5 : params.mixingScf;
  const eigenTol =
    params.eigenTol === undefined ? 0.0005 : params.eigenTol;
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

function buildPhshInput(params) {
  const atomicNumber = params.atomicNumber;
  const muffinTinRadius =
    params.muffinTinRadius === undefined ? 2.5 : params.muffinTinRadius;
  const energyMin = params.energyMin === undefined ? 1.0 : params.energyMin;
  const energyMax = params.energyMax === undefined ? 12.0 : params.energyMax;
  const energyStep =
    params.energyStep === undefined ? 0.25 : params.energyStep;
  const lmax = params.lmax === undefined ? 10 : params.lmax;

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

export { buildAtorbInput, buildPhshInput };
