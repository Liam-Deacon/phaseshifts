/**
 * Parsing helpers for phase shift output.
 * @module phase_shift_parser
 */

/**
 * Parse phase shift data from text output.
 * @param {string} output - Raw text output containing phase shift data
 * @param {number} lmax - Maximum angular momentum quantum number
 * @returns {{energies: number[], data: number[][]}} Parsed energies and phase shifts indexed by L
 */
export function parsePhaseShiftData(output, lmax) {
  if (!Number.isInteger(lmax) || lmax < 0)
    throw new TypeError('lmax must be a non-negative integer');

  const lines = output.split('\n');
  const energies = [];
  const energySet = new Set();
  const data = [];

  for (let l = 0; l <= lmax; l++) data.push([]);

  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;

    const values = trimmed
      .split(/\s+/)
      .map(Number.parseFloat)
      .filter((v) => !Number.isNaN(v));
    if (values.length < 2) continue;

    const energy = values[0];
    if (!energySet.has(energy)) {
      energies.push(energy);
      energySet.add(energy);
    }

    const limit = Math.min(lmax + 1, values.length - 1);
    for (let l = 0; l < limit; l++)
      // eslint-disable-next-line security/detect-object-injection -- l is bounded by lmax.
      data[l].push(values[l + 1]);
  }

  return { energies, data };
}
