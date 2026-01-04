/**
 * Phase Shifts Calculator - Application Logic
 * Browser-based LEED/XPD phase shift calculations
 */
/* global Chart */

// Global state
let phaseShiftsModule = null;
let chart = null;
let currentResults = null;

let globalScope = null;
if (typeof globalThis === 'object') {
  globalScope = globalThis;
}

const listenerController = new AbortController();
const listenerSignal = listenerController.signal;

function addListener(element, eventName, handler, options = {}) {
  element.addEventListener(eventName, handler, {
    ...options,
    signal: listenerSignal,
  });
}

function handleBeforeUnload() {
  listenerController.abort();
}

// Element data
const elements = Object.freeze({
  H: 1,
  He: 2,
  Li: 3,
  Be: 4,
  B: 5,
  C: 6,
  N: 7,
  O: 8,
  F: 9,
  Ne: 10,
  Na: 11,
  Mg: 12,
  Al: 13,
  Si: 14,
  P: 15,
  S: 16,
  Cl: 17,
  Ar: 18,
  K: 19,
  Ca: 20,
  Sc: 21,
  Ti: 22,
  V: 23,
  Cr: 24,
  Mn: 25,
  Fe: 26,
  Co: 27,
  Ni: 28,
  Cu: 29,
  Zn: 30,
  Ga: 31,
  Ge: 32,
  As: 33,
  Se: 34,
  Br: 35,
  Kr: 36,
  Rb: 37,
  Sr: 38,
  Y: 39,
  Zr: 40,
  Nb: 41,
  Mo: 42,
  Tc: 43,
  Ru: 44,
  Rh: 45,
  Pd: 46,
  Ag: 47,
  Cd: 48,
  In: 49,
  Sn: 50,
  Sb: 51,
  Te: 52,
  I: 53,
  Xe: 54,
  Cs: 55,
  Ba: 56,
  La: 57,
  Ce: 58,
  Pr: 59,
  Nd: 60,
  Pm: 61,
  Sm: 62,
  Eu: 63,
  Gd: 64,
  Tb: 65,
  Dy: 66,
  Ho: 67,
  Er: 68,
  Tm: 69,
  Yb: 70,
  Lu: 71,
  Hf: 72,
  Ta: 73,
  W: 74,
  Re: 75,
  Os: 76,
  Ir: 77,
  Pt: 78,
  Au: 79,
  Hg: 80,
  Tl: 81,
  Pb: 82,
  Bi: 83,
  Po: 84,
  At: 85,
  Rn: 86,
  Fr: 87,
  Ra: 88,
  Ac: 89,
  Th: 90,
  Pa: 91,
  U: 92,
});

// Presets for common calculations
const presets = Object.freeze({
  'cu-fcc': {
    element: 29,
    muffinTinRadius: 2.41,
    lmax: 10,
    energyMin: 20,
    energyMax: 400,
    energyStep: 5,
    method: 'rel',
  },
  'ni-fcc': {
    element: 28,
    muffinTinRadius: 2.35,
    lmax: 10,
    energyMin: 20,
    energyMax: 400,
    energyStep: 5,
    method: 'rel',
  },
  'fe-bcc': {
    element: 26,
    muffinTinRadius: 2.38,
    lmax: 10,
    energyMin: 20,
    energyMax: 400,
    energyStep: 5,
    method: 'rel',
  },
  'si-diamond': {
    element: 14,
    muffinTinRadius: 2.22,
    lmax: 8,
    energyMin: 20,
    energyMax: 300,
    energyStep: 5,
    method: 'cav',
  },
});

// Method descriptions
const methodDescriptions = Object.freeze({
  rel: 'Relativistic: Full Dirac equation treatment, essential for heavy elements (Z > 30)',
  cav: 'Cavity LEED: Traditional cavity method using Loucks grid, suitable for most applications',
  wil: "Williams: A.R. Williams' method, good for comparison studies",
});

function getPreset(presetName) {
  if (!Object.hasOwn(presets, presetName)) {
    return null;
  }
  return presets[presetName];
}

function getMethodDescription(method) {
  if (!Object.hasOwn(methodDescriptions, method)) {
    return '';
  }
  return methodDescriptions[method];
}

function getPhaseShiftSeries(series, index) {
  if (!Number.isInteger(index) || index < 0 || index >= series.length) {
    return null;
  }
  return series[index];
}

function getPhaseShiftValue(series, lIndex, eIndex) {
  const target = getPhaseShiftSeries(series, lIndex);
  if (!target || eIndex < 0 || eIndex >= target.length) {
    return null;
  }
  return target[eIndex];
}

function getArrayItem(values, index) {
  if (
    !Array.isArray(values) ||
    !Number.isInteger(index) ||
    index < 0 ||
    index >= values.length
  ) {
    return null;
  }
  return values[index];
}

function pushPhaseShiftValue(series, index, value) {
  const target = getPhaseShiftSeries(series, index);
  if (target) {
    target.push(value);
  }
}

function handleCalculateClick() {
  runCalculation();
}

function handleClearClick() {
  clearResults();
}

function handleMethodChange(event) {
  const description = getMethodDescription(event.target.value);
  document.getElementById('method-help').textContent = description;
}

function handlePresetClick(event) {
  const presetName = event.currentTarget.dataset.preset;
  loadPreset(presetName);
}

function handleTabClick(event) {
  const tabName = event.currentTarget.dataset.tab;
  switchTab(tabName);
}

function handleDownloadCleed() {
  downloadResults('cleed');
}

function handleDownloadViperLeed() {
  downloadResults('viperleed');
}

function handleDownloadCsv() {
  downloadResults('csv');
}

function handleChartChange() {
  updateChart();
}

// Initialize application
async function handleDomContentLoaded() {
  setupEventListeners();
  populateElementSelect();
  await initializeWasmModule();
}

addListener(document, 'DOMContentLoaded', handleDomContentLoaded, {
  once: true,
});
if (globalScope && typeof globalScope.addEventListener === 'function') {
  addListener(globalScope, 'beforeunload', handleBeforeUnload, { once: true });
}

/**
 * Set up all event listeners
 */
function setupEventListeners() {
  // Calculate button
  addListener(
    document.getElementById('calculate-btn'),
    'click',
    handleCalculateClick,
  );

  // Clear button
  addListener(document.getElementById('clear-btn'), 'click', handleClearClick);

  // Method selection - update help text
  addListener(document.getElementById('method'), 'change', handleMethodChange);

  // Preset buttons
  document.querySelectorAll('.preset-btn').forEach((btn) => {
    addListener(btn, 'click', handlePresetClick);
  });

  // Tab switching
  document.querySelectorAll('.tab-btn').forEach((btn) => {
    addListener(btn, 'click', handleTabClick);
  });

  // Download buttons
  addListener(
    document.getElementById('download-cleed'),
    'click',
    handleDownloadCleed,
  );
  addListener(
    document.getElementById('download-viperleed'),
    'click',
    handleDownloadViperLeed,
  );
  addListener(
    document.getElementById('download-csv'),
    'click',
    handleDownloadCsv,
  );

  // Chart controls
  addListener(
    document.getElementById('show-all-l'),
    'change',
    handleChartChange,
  );
  addListener(document.getElementById('l-min'), 'change', handleChartChange);
  addListener(
    document.getElementById('l-max-display'),
    'change',
    handleChartChange,
  );
}

/**
 * Populate the element select with all elements
 */
function populateElementSelect() {
  const select = document.getElementById('element');
  const allGroup = select.querySelector('optgroup[label="All Elements"]');

  Object.entries(elements).forEach(([symbol, z]) => {
    const option = document.createElement('option');
    option.value = z;
    option.textContent = `${symbol} (Z=${z})`;
    allGroup.appendChild(option);
  });
}

/**
 * Initialize the WebAssembly module
 */
async function initializeWasmModule() {
  const statusBanner = document.getElementById('status-banner');
  const statusIcon = document.getElementById('status-icon');
  const statusText = document.getElementById('status-text');
  const calculateBtn = document.getElementById('calculate-btn');
  const versionInfo = document.getElementById('version-info');

  try {
    // Check if the WASM module loader is available
    const moduleFactory =
      globalScope && typeof globalScope.createPhaseShiftsModule === 'function'
        ? globalScope.createPhaseShiftsModule
        : null;
    if (typeof moduleFactory !== 'function') {
      throw new TypeError(
        'WASM module not found. Please build the WASM files first.',
      );
    }

    statusText.textContent = 'Initializing WebAssembly module...';

    // Create the module
    phaseShiftsModule = await moduleFactory();

    // Initialize filesystem
    phaseShiftsModule.FS.mkdir('/input');
    phaseShiftsModule.FS.mkdir('/output');
    phaseShiftsModule.FS.mkdir('/work');

    // Call init function
    phaseShiftsModule.ccall('init_phaseshifts', null, [], []);

    // Get version
    const version = phaseShiftsModule.ccall('get_version', 'string', [], []);
    versionInfo.textContent = version;

    // Update status
    statusBanner.className = 'status-banner ready';
    statusIcon.textContent = '✅';
    statusText.textContent = 'WebAssembly module ready!';
    calculateBtn.disabled = false;

    // Hide banner after 3 seconds
    scheduleStatusBannerHide(statusBanner, 3000);
  } catch (error) {
    console.error('Failed to initialize WASM module:', error);

    statusBanner.className = 'status-banner error';
    statusIcon.textContent = '❌';
    statusText.textContent = `Failed to load: ${error.message}`;
    versionInfo.textContent = 'Module not loaded';

    // Show demo mode option
    showDemoMode();
  }
}

function scheduleStatusBannerHide(element, delayMs) {
  const start = performance.now();

  function tick(now) {
    if (now - start >= delayMs) {
      element.classList.add('hidden');
      return;
    }
    // eslint-disable-next-line -- requestAnimationFrame used for UI animation timing
    requestAnimationFrame(tick);
  }

  // eslint-disable-next-line -- requestAnimationFrame used for UI animation timing
  requestAnimationFrame(tick);
}

/**
 * Show demo mode when WASM is not available
 */
function showDemoMode() {
  const calculateBtn = document.getElementById('calculate-btn');
  calculateBtn.disabled = false;
  calculateBtn.textContent = '🎭 Run Demo (No WASM)';
  calculateBtn.dataset.demoMode = 'true';
}

function setButtonLoading(button, label) {
  const spinner = document.createElement('span');
  spinner.className = 'spinner';
  spinner.setAttribute('aria-hidden', 'true');
  button.replaceChildren(spinner, document.createTextNode(` ${label}`));
}

/**
 * Load a preset configuration
 */
function loadPreset(presetName) {
  const preset = getPreset(presetName);
  if (!preset) return;

  document.getElementById('element').value = preset.element;
  document.getElementById('muffin-tin-radius').value = preset.muffinTinRadius;
  document.getElementById('lmax').value = preset.lmax;
  document.getElementById('energy-min').value = preset.energyMin;
  document.getElementById('energy-max').value = preset.energyMax;
  document.getElementById('energy-step').value = preset.energyStep;
  document.getElementById('method').value = preset.method;

  // Update method help text
  document.getElementById('method-help').textContent = getMethodDescription(
    preset.method,
  );
}

/**
 * Run the phase shift calculation
 */
async function runCalculation() {
  const calculateBtn = document.getElementById('calculate-btn');
  const originalText = calculateBtn.textContent || '';

  try {
    // Update button state
    calculateBtn.disabled = true;
    setButtonLoading(calculateBtn, 'Calculating...');

    // Get input parameters
    const params = {
      atomicNumber: Number.parseInt(
        document.getElementById('element').value,
        10,
      ),
      muffinTinRadius: Number.parseFloat(
        document.getElementById('muffin-tin-radius').value,
      ),
      lmax: Number.parseInt(document.getElementById('lmax').value, 10),
      energyMin: Number.parseFloat(document.getElementById('energy-min').value),
      energyMax: Number.parseFloat(document.getElementById('energy-max').value),
      energyStep: Number.parseFloat(
        document.getElementById('energy-step').value,
      ),
      method: document.getElementById('method').value,
    };

    // Validate inputs
    if (
      !params.atomicNumber ||
      params.atomicNumber < 1 ||
      params.atomicNumber > 92
    ) {
      throw new Error('Please select a valid element');
    }

    let results;

    // Check if we're in demo mode
    if (calculateBtn.dataset.demoMode === 'true') {
      results = generateDemoResults(params);
    } else {
      results = await calculatePhaseShifts(params);
    }

    // Store and display results
    currentResults = { params, results };
    displayResults(results, params);
  } catch (error) {
    alert(`Calculation error: ${error.message}`);
    console.error(error);
  } finally {
    calculateBtn.disabled = false;
    calculateBtn.textContent = originalText;
  }
}

/**
 * Calculate phase shifts using WASM module
 */
async function calculatePhaseShifts(params) {
  const method = params.method;

  // Generate input file for the calculation
  const inputData = generateInputFile(params);

  // Write input to MEMFS
  phaseShiftsModule.FS.writeFile('/input/phsh.i', inputData);

  // Run the appropriate calculation
  let result;
  switch (method) {
    case 'rel':
      result = phaseShiftsModule.ccall('run_phsh_rel', 'number', [], []);
      break;
    case 'cav':
      result = phaseShiftsModule.ccall('run_phsh_cav', 'number', [], []);
      break;
    case 'wil':
      result = phaseShiftsModule.ccall('run_phsh_wil', 'number', [], []);
      break;
    default:
      throw new Error(`Unknown method: ${method}`);
  }

  if (result !== 0) {
    throw new Error(`Calculation failed with error code ${result}`);
  }

  // Read output
  const output = phaseShiftsModule.FS.readFile('/output/phasout.o', {
    encoding: 'utf8',
  });

  return parsePhaseShiftOutput(output, params);
}

/**
 * Generate input file content
 */
function generateInputFile(params) {
  // This would generate the proper Fortran input format
  // Simplified for now
  return [
    params.atomicNumber,
    params.muffinTinRadius,
    params.energyMin,
    params.energyMax,
    params.energyStep,
    params.lmax,
  ].join(' ');
}

/**
 * Parse phase shift output
 */
function parsePhaseShiftOutput(output, params) {
  // Parse the output format from the Fortran code
  const lines = output.split('\n');
  const energies = [];
  const phaseShifts = [];

  // Initialize arrays for each L value
  for (let l = 0; l <= params.lmax; l++) {
    phaseShifts.push([]);
  }

  // Parse each line
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith('#')) continue;

    const values = trimmed
      .split(/\s+/)
      .map(Number.parseFloat)
      .filter((v) => !Number.isNaN(v));
    if (values.length >= 2) {
      const energy = values[0];
      energies.push(energy);

      for (let l = 0; l < Math.min(values.length - 1, params.lmax + 1); l++) {
        pushPhaseShiftValue(phaseShifts, l, values[l + 1]);
      }
    }
  }

  return {
    raw: output,
    energies,
    phaseShifts,
    success: true,
  };
}

/**
 * Generate demo results when WASM is not available
 */
function generateDemoResults(params) {
  const { atomicNumber, lmax, energyMin, energyMax, energyStep } = params;

  const energies = [];
  const phaseShifts = [];

  // Initialize arrays
  for (let l = 0; l <= lmax; l++) {
    phaseShifts.push([]);
  }

  // Generate synthetic phase shifts
  for (let e = energyMin; e <= energyMax; e += energyStep) {
    energies.push(e);

    for (let l = 0; l <= lmax; l++) {
      // Simplified model: phase shift depends on energy and L
      // Real phase shifts would come from the Fortran calculation
      const k = Math.sqrt(e / 13.6); // approximate wave vector
      const delta = Math.atan(Math.pow(atomicNumber / 10, 0.5) / (k * (l + 1)));

      // Add some realistic-looking energy dependence
      const phaseShift = delta * (1 + 0.1 * Math.sin(e / 50)) * (1 - l * 0.05);
      pushPhaseShiftValue(phaseShifts, l, phaseShift);
    }
  }

  // Generate raw output text
  let raw = '# Phase Shifts (DEMO MODE)\n';
  raw += `# Element: Z=${atomicNumber}, Lmax=${lmax}\n`;
  raw += `# Energy range: ${energyMin}-${energyMax} eV, step ${energyStep} eV\n`;
  raw += '#\n';
  raw +=
    '# E(eV)    ' +
    Array.from({ length: lmax + 1 }, (_, i) => `L=${i}`).join('      ') +
    '\n';

  for (let i = 0; i < energies.length; i++) {
    const energy = getArrayItem(energies, i) ?? 0;
    raw += energy.toFixed(2).padStart(8);
    for (let l = 0; l <= lmax; l++) {
      const phaseShift = getPhaseShiftValue(phaseShifts, l, i);
      raw += (phaseShift ?? 0).toFixed(4).padStart(10);
    }
    raw += '\n';
  }

  return {
    raw,
    energies,
    phaseShifts,
    success: true,
    isDemo: true,
  };
}

/**
 * Display calculation results
 */
function displayResults(results, params) {
  const resultsSection = document.getElementById('results-section');
  resultsSection.style.display = 'block';

  // Update L max display
  document.getElementById('l-max-display').value = params.lmax;

  // Update raw output
  document.getElementById('raw-output').textContent = results.raw;

  // Update table
  updateResultsTable(results, params);

  // Update chart
  createChart(results, params);

  // Scroll to results
  resultsSection.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Update the results table
 */
function updateResultsTable(results, params) {
  const headerRow = document.getElementById('table-header');
  const tableBody = document.getElementById('table-body');

  // Clear existing
  headerRow.replaceChildren();
  tableBody.replaceChildren();

  const energyHeader = document.createElement('th');
  energyHeader.textContent = 'Energy (eV)';
  headerRow.appendChild(energyHeader);

  // Add L headers
  for (let l = 0; l <= params.lmax; l++) {
    const th = document.createElement('th');
    th.textContent = `L=${l}`;
    headerRow.appendChild(th);
  }

  // Add data rows
  for (let i = 0; i < results.energies.length; i++) {
    const tr = document.createElement('tr');

    // Energy column
    const tdE = document.createElement('td');
    const energy = getArrayItem(results.energies, i);
    tdE.textContent = energy === null ? '-' : energy.toFixed(2);
    tr.appendChild(tdE);

    // Phase shift columns
    for (let l = 0; l <= params.lmax; l++) {
      const td = document.createElement('td');
      const phaseShift = getPhaseShiftValue(results.phaseShifts, l, i);
      td.textContent = phaseShift === null ? '-' : phaseShift.toFixed(4);
      tr.appendChild(td);
    }

    tableBody.appendChild(tr);
  }
}

/**
 * Create or update the phase shift chart
 */
function createChart(results, params) {
  const ctx = document.getElementById('phase-shift-chart').getContext('2d');

  // Destroy existing chart
  if (chart) {
    chart.destroy();
  }

  // Generate colors for each L value
  const colors = generateColors(params.lmax + 1);

  // Create datasets
  const datasets = [];
  for (let l = 0; l <= params.lmax; l++) {
    const series = getPhaseShiftSeries(results.phaseShifts, l) || [];
    const color = getArrayItem(colors, l) || '#888888';
    datasets.push({
      label: `L=${l}`,
      data: series,
      borderColor: color,
      backgroundColor: `${color}20`,
      borderWidth: 2,
      pointRadius: 0,
      tension: 0.3,
    });
  }

  chart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: results.energies.map((e) => e.toFixed(1)),
      datasets,
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        title: {
          display: true,
          text: `Phase Shifts vs Energy (${params.method.toUpperCase()} method)${
            results.isDemo ? ' [DEMO]' : ''
          }`,
        },
        legend: {
          position: 'right',
        },
      },
      scales: {
        x: {
          title: {
            display: true,
            text: 'Energy (eV)',
          },
        },
        y: {
          title: {
            display: true,
            text: 'Phase Shift (radians)',
          },
        },
      },
    },
  });
}

/**
 * Update chart based on controls
 */
function updateChart() {
  if (!chart || !currentResults) return;

  const showAll = document.getElementById('show-all-l').checked;
  const lMin = Number.parseInt(document.getElementById('l-min').value, 10) || 0;
  const lMax =
    Number.parseInt(document.getElementById('l-max-display').value, 10) || 0;

  chart.data.datasets.forEach((dataset, index) => {
    if (showAll) {
      dataset.hidden = false;
    } else {
      dataset.hidden = index < lMin || index > lMax;
    }
  });

  chart.update();
}

/**
 * Generate distinct colors for chart lines
 */
function generateColors(count) {
  const colors = [];
  for (let i = 0; i < count; i++) {
    const hue = ((i * 360) / count) % 360;
    colors.push(`hsl(${hue}, 70%, 50%)`);
  }
  return colors;
}

/**
 * Switch between tabs
 */
function switchTab(tabName) {
  // Update tab buttons
  document.querySelectorAll('.tab-btn').forEach((btn) => {
    btn.classList.toggle('active', btn.dataset.tab === tabName);
  });

  // Update tab content
  document.querySelectorAll('.tab-content').forEach((content) => {
    content.classList.toggle('active', content.id === `tab-${tabName}`);
  });
}

/**
 * Download results in specified format
 */
function downloadResults(format) {
  if (!currentResults) {
    alert('No results to download');
    return;
  }

  const { params, results } = currentResults;
  let content, filename, mimeType;

  switch (format) {
    case 'cleed':
      content = generateCleedFormat(results, params);
      filename = `phaseshifts_Z${params.atomicNumber}.phs`;
      mimeType = 'text/plain';
      break;

    case 'viperleed':
      content = generateViperLeedFormat(results, params);
      filename = `PHASESHIFTS_Z${params.atomicNumber}`;
      mimeType = 'text/plain';
      break;

    case 'csv':
      content = generateCsvFormat(results, params);
      filename = `phaseshifts_Z${params.atomicNumber}.csv`;
      mimeType = 'text/csv';
      break;

    default:
      content = results.raw;
      filename = 'phaseshifts.txt';
      mimeType = 'text/plain';
  }

  downloadFile(content, filename, mimeType);
}

/**
 * Generate CLEED format output
 */
function generateCleedFormat(results, params) {
  let output = `# Phase shifts for Z=${params.atomicNumber}\n`;
  output += `# Method: ${params.method}\n`;
  output += `# Muffin-tin radius: ${params.muffinTinRadius} Bohr\n`;
  output += `# Lmax: ${params.lmax}\n`;
  output += `# Energy range: ${params.energyMin}-${params.energyMax} eV\n`;
  output += '#\n';

  // Header
  output += `${params.energyMin.toFixed(4)} ${params.energyStep.toFixed(4)} ${
    results.energies.length
  } ${params.lmax + 1}\n`;

  // Phase shift data
  for (let i = 0; i < results.energies.length; i++) {
    const energy = getArrayItem(results.energies, i) ?? 0;
    let line = energy.toFixed(4);
    for (let l = 0; l <= params.lmax; l++) {
      const phaseShift = getPhaseShiftValue(results.phaseShifts, l, i) ?? 0;
      line += ` ${phaseShift.toFixed(6)}`;
    }
    output += line + '\n';
  }

  return output;
}

/**
 * Generate ViPErLEED format output
 */
function generateViperLeedFormat(results, params) {
  let output = `${params.lmax + 1} ${results.energies.length}\n`;

  for (let i = 0; i < results.energies.length; i++) {
    const energy = getArrayItem(results.energies, i) ?? 0;
    output += `${energy.toFixed(2)}`;
    for (let l = 0; l <= params.lmax; l++) {
      const phaseShift = getPhaseShiftValue(results.phaseShifts, l, i) ?? 0;
      output += ` ${phaseShift.toFixed(4)}`;
    }
    output += '\n';
  }

  return output;
}

/**
 * Generate CSV format output
 */
function generateCsvFormat(results, params) {
  let output = 'Energy(eV)';
  for (let l = 0; l <= params.lmax; l++) {
    output += `,L=${l}`;
  }
  output += '\n';

  for (let i = 0; i < results.energies.length; i++) {
    const energy = getArrayItem(results.energies, i) ?? 0;
    output += energy.toFixed(4);
    for (let l = 0; l <= params.lmax; l++) {
      const phaseShift = getPhaseShiftValue(results.phaseShifts, l, i) ?? 0;
      output += `,${phaseShift.toFixed(6)}`;
    }
    output += '\n';
  }

  return output;
}

/**
 * Download a file
 */
function downloadFile(content, filename, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

/**
 * Clear results
 */
function clearResults() {
  currentResults = null;
  document.getElementById('results-section').style.display = 'none';

  if (chart) {
    chart.destroy();
    chart = null;
  }
}
