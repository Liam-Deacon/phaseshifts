/**
 * Parse phase shift data from text output.
 * @param {string} output - Raw text output containing phase shift data
 * @param {number} lmax - Maximum angular momentum quantum number
 * @returns {{energies: number[], data: number[][]}} Parsed energies and phase shifts indexed by L
 */
function parsePhaseShiftData(output, lmax) {
  if (!Number.isInteger(lmax) || lmax < 0) {
    throw new TypeError('lmax must be a non-negative integer');
  }

  const lines = output.split('\n');
  const energies = [];
  const energySet = new Set();
  const data = [];

  for (let l = 0; l <= lmax; l++) {
    data.push([]);
  }

  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;

    const values = trimmed
      .split(/\s+/)
      .map(parseFloat)
      .filter((v) => !isNaN(v));
    if (values.length >= 2) {
      const energy = values[0];
      if (!energySet.has(energy)) {
        energies.push(energy);
        energySet.add(energy);
      }

      for (let l = 0; l < values.length - 1 && l <= lmax; l++) {
        data.at(l).push(values[l + 1]);
      }
    }
  }

  return { energies, data };
}

export { parsePhaseShiftData };
