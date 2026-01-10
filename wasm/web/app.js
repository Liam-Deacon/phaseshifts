/**
 * Phase Shifts Calculator - Application Logic
 * Browser-based LEED/XPD phase shift calculations
 */
/* global Chart */

// Note: shared/elements.js path works for both local dev and deployed structure
import { elements } from './shared/elements.js';
import { createStructureBuilder } from './shared/structure-builder.js';
import { createViewer } from './shared/viewer3d.js';

// Global state
let phaseShiftsModule = null;
let chart = null;
let currentResults = null;
let crystalStructure = null;
let structureBuilder = null;
let crystalViewer = null;

let globalScope = null;
if (typeof globalThis === 'object') {
  globalScope = globalThis;
}

const listenerController = new AbortController();
const listenerSignal = listenerController.signal;

function addListener(element, eventName, handler, options = {}) {
  if (!element) return;
  element.addEventListener(eventName, handler, {
    ...options,
    signal: listenerSignal,
  });
}

function handleBeforeUnload() {
  listenerController.abort();
  if (crystalViewer) {
    crystalViewer.dispose();
  }
}

// Method descriptions
const methodDescriptions = Object.freeze({
  rel: 'Relativistic: Full Dirac equation treatment, essential for heavy elements (Z > 30)',
  cav: 'Cavity LEED: Traditional cavity method using Loucks grid, suitable for most applications',
  wil: "Williams: A.R. Williams' method, good for comparison studies",
});

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

// Initialize application
async function handleDomContentLoaded() {
  setupMainTabs();
  setupEventListeners();
  initializeStructureBuilder();
  initializeViewer();
  await initializeWasmModule();
}

addListener(document, 'DOMContentLoaded', handleDomContentLoaded, {
  once: true,
});
if (globalScope && typeof globalScope.addEventListener === 'function') {
  addListener(globalScope, 'beforeunload', handleBeforeUnload, { once: true });
}

/**
 * Set up main tab navigation
 */
function setupMainTabs() {
  const tabs = document.querySelectorAll('.main-tab');

  tabs.forEach((tab) => {
    addListener(tab, 'click', () => {
      // Update active tab
      tabs.forEach((t) => t.classList.remove('active'));
      tab.classList.add('active');

      // Update active panel
      const panelId = `panel-${tab.dataset.tab}`;
      document.querySelectorAll('.tab-panel').forEach((panel) => {
        panel.classList.toggle('active', panel.id === panelId);
      });

      // Resize viewer when switching to structure tab
      if (tab.dataset.tab === 'structure' && crystalViewer) {
        crystalViewer._onResize();
      }
    });
  });
}

/**
 * Initialize the structure builder component
 */
function initializeStructureBuilder() {
  const container = document.getElementById('structure-builder');
  if (!container) return;

  structureBuilder = createStructureBuilder(container, {
    onStructureChange: handleStructureChange,
    showPresets: true,
    showImportExport: true,
  });

  // Load a default preset
  try {
    structureBuilder.loadPreset('Cu(111)');
  } catch (e) {
    console.warn('Could not load default preset:', e);
  }
}

/**
 * Initialize the 3D viewer
 */
function initializeViewer() {
  const container = document.getElementById('crystal-viewer');
  if (!container) return;

  try {
    crystalViewer = createViewer(container, {
      backgroundColor: 0x1a1a2e,
      showAxes: true,
      showBonds: true,
      showUnitCell: true,
      repeatX: 2,
      repeatY: 2,
    });
  } catch (error) {
    console.error('Failed to initialize 3D viewer:', error);
    container.innerHTML =
      '<p style="color: #888; padding: 2rem; text-align: center;">3D viewer requires WebGL support</p>';
  }

  // View control buttons
  addListener(document.getElementById('view-perspective'), 'click', () => {
    if (crystalViewer) crystalViewer.setView('perspective');
  });
  addListener(document.getElementById('view-top'), 'click', () => {
    if (crystalViewer) crystalViewer.setView('top');
  });
  addListener(document.getElementById('view-front'), 'click', () => {
    if (crystalViewer) crystalViewer.setView('front');
  });
  addListener(document.getElementById('view-side'), 'click', () => {
    if (crystalViewer) crystalViewer.setView('side');
  });
  addListener(document.getElementById('view-reset'), 'click', () => {
    if (crystalViewer) crystalViewer.resetView();
  });
}

/**
 * Handle structure changes from the builder
 */
function handleStructureChange(structure) {
  if (!structure) return;
  crystalStructure = structure;

  // Update 3D viewer
  if (crystalViewer) {
    crystalViewer.setStructure(structure);
  }

  // Update structure info display
  const nameEl = document.getElementById('structure-name');
  const elementsEl = document.getElementById('structure-elements');

  if (nameEl) {
    nameEl.textContent = structure.name;
  }

  if (elementsEl) {
    const uniqueElements = structure.getUniqueElements();
    elementsEl.textContent =
      uniqueElements.length > 0 ? `Elements: ${uniqueElements.join(', ')}` : '';
  }

  // Update structure summary in calculation tab
  updateStructureSummary(structure);

  // Update element-specific parameters
  updateElementParameters(structure);

  // Enable calculate button if structure has atoms
  const calculateBtn = document.getElementById('calculate-btn');
  const hasAtoms = structure.layers.some((l) => l.atoms.length > 0);
  if (calculateBtn) {
    calculateBtn.disabled = !hasAtoms;
  }
}

/**
 * Update the structure summary in the calculation tab
 */
function updateStructureSummary(structure) {
  const container = document.getElementById('structure-summary-content');
  if (!container) return;

  const uniqueElements = structure.getUniqueElements();
  const totalAtoms = structure.layers.reduce(
    (sum, l) => sum + l.atoms.length,
    0,
  );

  // Clear container
  container.replaceChildren();

  if (totalAtoms === 0) {
    const emptyMsg = document.createElement('p');
    emptyMsg.className = 'empty-message';
    emptyMsg.textContent =
      'No structure defined. Go to "Crystal Structure" tab to build one.';
    container.appendChild(emptyMsg);
    return;
  }

  const grid = document.createElement('div');
  grid.className = 'summary-grid';

  // Helper to create summary items
  const createSummaryItem = (label, value) => {
    const item = document.createElement('div');
    item.className = 'summary-item';

    const labelSpan = document.createElement('span');
    labelSpan.className = 'summary-label';
    labelSpan.textContent = label;
    item.appendChild(labelSpan);

    const valueSpan = document.createElement('span');
    valueSpan.className = 'summary-value';
    valueSpan.textContent = value;
    item.appendChild(valueSpan);

    return item;
  };

  grid.appendChild(createSummaryItem('Structure:', structure.name));
  grid.appendChild(
    createSummaryItem('Surface:', structure.millerIndices.toSimpleString()),
  );
  grid.appendChild(createSummaryItem('Elements:', uniqueElements.join(', ')));
  grid.appendChild(createSummaryItem('Layers:', structure.layers.length));
  grid.appendChild(createSummaryItem('Total Atoms:', totalAtoms));

  container.appendChild(grid);
}

/**
 * Update element-specific parameter inputs
 */
function updateElementParameters(structure) {
  const container = document.getElementById('element-params-list');
  if (!container) return;

  const uniqueElements = structure.getUniqueElements();

  // Clear container
  container.replaceChildren();

  if (uniqueElements.length === 0) {
    const emptyMsg = document.createElement('p');
    emptyMsg.className = 'empty-message';
    emptyMsg.textContent =
      'Add atoms to the structure to configure element parameters.';
    container.appendChild(emptyMsg);
    return;
  }

  for (const symbol of uniqueElements) {
    const z = elements[symbol] || 0;
    // Default muffin-tin radius based on atomic number
    const defaultMT = (1.5 + z * 0.02).toFixed(2);

    const row = document.createElement('div');
    row.className = 'element-param-row';

    // Element symbol span
    const symbolSpan = document.createElement('span');
    symbolSpan.className = 'element-symbol';
    symbolSpan.style.color = getElementColorForCSS(symbol);
    symbolSpan.textContent = symbol;
    row.appendChild(symbolSpan);

    // Muffin-tin radius input group
    const mtGroup = document.createElement('div');
    mtGroup.className = 'form-group';

    const mtLabel = document.createElement('label');
    mtLabel.textContent = 'Muffin-Tin Radius (Bohr)';
    mtGroup.appendChild(mtLabel);

    const mtInput = document.createElement('input');
    mtInput.type = 'number';
    mtInput.className = 'mt-radius-input';
    mtInput.dataset.element = symbol;
    mtInput.value = defaultMT;
    mtInput.step = '0.01';
    mtInput.min = '0.5';
    mtInput.max = '5.0';
    mtGroup.appendChild(mtInput);

    row.appendChild(mtGroup);

    // Inner potential input group
    const v0Group = document.createElement('div');
    v0Group.className = 'form-group';

    const v0Label = document.createElement('label');
    v0Label.textContent = 'Inner Potential V₀ (eV)';
    v0Group.appendChild(v0Label);

    const v0Input = document.createElement('input');
    v0Input.type = 'number';
    v0Input.className = 'v0-input';
    v0Input.dataset.element = symbol;
    v0Input.value = '10.0';
    v0Input.step = '0.5';
    v0Input.min = '0';
    v0Input.max = '30';
    v0Group.appendChild(v0Input);

    row.appendChild(v0Group);

    container.appendChild(row);
  }
}

/**
 * Get element color for CSS usage
 */
function getElementColorForCSS(symbol) {
  const colors = {
    Cu: '#C88033',
    Ni: '#50D050',
    Fe: '#E06633',
    Al: '#BFA6A6',
    Si: '#F0C8A0',
    Ag: '#C0C0C0',
    Au: '#FFD123',
    Pt: '#D0D0E0',
    Co: '#F090A0',
    Zn: '#7D80B0',
    C: '#909090',
    O: '#FF0D0D',
    N: '#3050F8',
    H: '#FFFFFF',
  };
  return colors[symbol] || '#888888';
}

/**
 * Set up all event listeners
 */
function setupEventListeners() {
  // Calculate button
  addListener(document.getElementById('calculate-btn'), 'click', () =>
    runCalculation(),
  );

  // Method selection - update help text
  addListener(document.getElementById('method'), 'change', (event) => {
    const description = getMethodDescription(event.target.value);
    document.getElementById('method-help').textContent = description;
  });

  // Tab switching for results
  document.querySelectorAll('.tabs .tab-btn').forEach((btn) => {
    addListener(btn, 'click', () => {
      const tabName = btn.dataset.tab;
      switchResultsTab(tabName);
    });
  });

  // Download buttons
  addListener(document.getElementById('download-cleed'), 'click', () =>
    downloadResults('cleed'),
  );
  addListener(document.getElementById('download-viperleed'), 'click', () =>
    downloadResults('viperleed'),
  );
  addListener(document.getElementById('download-csv'), 'click', () =>
    downloadResults('csv'),
  );

  // Chart controls
  addListener(document.getElementById('show-all-l'), 'change', () =>
    updateChart(),
  );
  addListener(document.getElementById('l-min'), 'change', () => updateChart());
  addListener(document.getElementById('l-max-display'), 'change', () =>
    updateChart(),
  );
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
    statusIcon.textContent = 'OK';
    statusText.textContent = 'WebAssembly module ready!';

    // Enable calculate if structure has atoms
    if (crystalStructure?.layers.some((l) => l.atoms.length > 0)) {
      calculateBtn.disabled = false;
    }

    // Hide banner after 3 seconds
    scheduleStatusBannerHide(statusBanner, 3000);
  } catch (error) {
    console.error('Failed to initialize WASM module:', error);

    statusBanner.className = 'status-banner error';
    statusIcon.textContent = 'X';
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
    requestAnimationFrame(tick);
  }

  requestAnimationFrame(tick);
}

/**
 * Show demo mode when WASM is not available
 */
function showDemoMode() {
  const calculateBtn = document.getElementById('calculate-btn');
  if (crystalStructure?.layers.some((l) => l.atoms.length > 0)) {
    calculateBtn.disabled = false;
  }
  calculateBtn.textContent = 'Run Demo (No WASM)';
  calculateBtn.dataset.demoMode = 'true';
}

function setButtonLoading(button, label) {
  const spinner = document.createElement('span');
  spinner.className = 'spinner';
  spinner.setAttribute('aria-hidden', 'true');
  button.replaceChildren(spinner, document.createTextNode(` ${label}`));
}

/**
 * Run the phase shift calculation
 */
async function runCalculation() {
  const calculateBtn = document.getElementById('calculate-btn');
  const originalText = calculateBtn.textContent || '';

  if (
    !crystalStructure ||
    crystalStructure.layers.every((l) => l.atoms.length === 0)
  ) {
    alert('Please define a crystal structure with atoms first.');
    return;
  }

  try {
    // Update button state
    calculateBtn.disabled = true;
    setButtonLoading(calculateBtn, 'Calculating...');

    // Get input parameters
    const params = {
      structure: crystalStructure,
      elements: crystalStructure.getUniqueElements(),
      lmax: Number.parseInt(document.getElementById('lmax').value, 10),
      energyMin: Number.parseFloat(document.getElementById('energy-min').value),
      energyMax: Number.parseFloat(document.getElementById('energy-max').value),
      energyStep: Number.parseFloat(
        document.getElementById('energy-step').value,
      ),
      method: document.getElementById('method').value,
      elementParams: getElementParams(),
    };

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

    // Switch to results tab
    document.querySelector('.main-tab[data-tab="results"]').click();
  } catch (error) {
    alert(`Calculation error: ${error.message}`);
    console.error(error);
  } finally {
    calculateBtn.disabled = false;
    calculateBtn.textContent = originalText;
  }
}

/**
 * Get element-specific parameters from UI
 */
function getElementParams() {
  const params = {};

  document.querySelectorAll('.mt-radius-input').forEach((input) => {
    const element = input.dataset.element;
    if (!params[element]) params[element] = {};
    params[element].muffinTinRadius = Number.parseFloat(input.value);
  });

  document.querySelectorAll('.v0-input').forEach((input) => {
    const element = input.dataset.element;
    if (!params[element]) params[element] = {};
    params[element].innerPotential = Number.parseFloat(input.value);
  });

  return params;
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
  const lines = [
    `# Phase shift calculation for ${params.structure.name}`,
    `# Method: ${params.method}`,
  ];

  // For each unique element
  for (const element of params.elements) {
    const z = elements[element];
    if (z === undefined) {
      throw new Error(
        `Unknown element symbol: ${element} in structure ${params.structure.name}`,
      );
    }
    const elementParams = params.elementParams[element] || {};
    const mt = Number(elementParams.muffinTinRadius) || 2.5;

    lines.push(
      `${z} ${mt} ${params.energyMin} ${params.energyMax} ${params.energyStep} ${params.lmax}`,
    );
  }

  return lines.join('\n');
}

/**
 * Parse phase shift output
 */
function parsePhaseShiftOutput(output, params) {
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
  const { lmax, energyMin, energyMax, energyStep } = params;

  const energies = [];
  const phaseShifts = [];

  // Initialize arrays
  for (let l = 0; l <= lmax; l++) {
    phaseShifts.push([]);
  }

  // Get representative atomic number from structure
  const primaryElement = params.elements[0] || 'Cu';
  const atomicNumber = elements[primaryElement] || 29;

  // Generate synthetic phase shifts
  for (let e = energyMin; e <= energyMax; e += energyStep) {
    energies.push(e);

    for (let l = 0; l <= lmax; l++) {
      // Simplified model: phase shift depends on energy and L
      const k = Math.sqrt(e / 13.6);
      const delta = Math.atan(Math.pow(atomicNumber / 10, 0.5) / (k * (l + 1)));
      const phaseShift = delta * (1 + 0.1 * Math.sin(e / 50)) * (1 - l * 0.05);
      pushPhaseShiftValue(phaseShifts, l, phaseShift);
    }
  }

  // Generate raw output text
  let raw = '# Phase Shifts (DEMO MODE)\n';
  raw += `# Structure: ${params.structure.name}\n`;
  raw += `# Elements: ${params.elements.join(', ')}\n`;
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
  const noResults = document.getElementById('no-results');

  resultsSection.style.display = 'block';
  if (noResults) noResults.style.display = 'none';

  // Update L max display
  document.getElementById('l-max-display').value = params.lmax;

  // Update raw output
  document.getElementById('raw-output').textContent = results.raw;

  // Update table
  updateResultsTable(results, params);

  // Update chart
  createChart(results, params);
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
          text: `Phase Shifts - ${params.structure.name} (${params.method.toUpperCase()})${
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
 * Switch between results tabs
 */
function switchResultsTab(tabName) {
  // Update tab buttons
  document.querySelectorAll('.tabs .tab-btn').forEach((btn) => {
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
      filename = `phaseshifts_${params.structure.name.replaceAll(/[^a-z0-9]/gi, '_')}.phs`;
      mimeType = 'text/plain';
      break;

    case 'viperleed':
      content = generateViperLeedFormat(results, params);
      filename = `PHASESHIFTS_${params.structure.name.replaceAll(/[^a-z0-9]/gi, '_')}`;
      mimeType = 'text/plain';
      break;

    case 'csv':
      content = generateCsvFormat(results, params);
      filename = `phaseshifts_${params.structure.name.replaceAll(/[^a-z0-9]/gi, '_')}.csv`;
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
  let output = `# Phase shifts for ${params.structure.name}\n`;
  output += `# Method: ${params.method}\n`;
  output += `# Elements: ${params.elements.join(', ')}\n`;
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
