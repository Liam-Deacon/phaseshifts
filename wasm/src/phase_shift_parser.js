function parsePhaseShiftData(output, lmax) {
  const lines = output.split('\n');
  const energies = [];
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
      if (!energies.includes(energy)) {
        energies.push(energy);
      }

      for (let l = 0; l < values.length - 1 && l <= lmax; l++) {
        data[l].push(values[l + 1]);
      }
    }
  }

  return { energies, data };
}

export { parsePhaseShiftData };
