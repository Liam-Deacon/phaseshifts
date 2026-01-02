/**
 * phaseshifts.js - JavaScript API for browser-based phase shift calculations
 *
 * This module provides a high-level JavaScript API for the WebAssembly-compiled
 * Fortran phase shift calculation library.
 *
 * @module phaseshifts
 * @version 0.1.0
 * @license MIT
 */

/* global createPhaseShiftsModule */

import { elements } from './elements.js';
import { parsePhaseShiftData } from './phase_shift_parser.js';

/**
 * PhaseShifts calculator class
 * Wraps the WASM module and provides a clean JavaScript API
 */
class PhaseShifts {
  /**
   * Create a PhaseShifts calculator instance
   * @param {Object} wasmModule - The initialized Emscripten module
   */
  constructor(wasmModule) {
    this.Module = wasmModule;
    this.FS = wasmModule.FS;
    this._initialized = false;
  }

  /**
   * Initialize the PhaseShifts module
   * Must be called before any calculations
   * @returns {Promise<PhaseShifts>} This instance for chaining
   */
  async init() {
    if (this._initialized) return this;

    // Call the WASM init function
    this.Module.ccall('init_phaseshifts', null, [], []);

    // Create working directories in MEMFS
    this._createDir('/input');
    this._createDir('/output');
    this._createDir('/work');

    this._initialized = true;
    return this;
  }

  /**
   * Create a directory in MEMFS if it doesn't exist
   * @private
   */
  _createDir(path) {
    const info = this.FS.analyzePath(path);
    if (!info.exists) {
      this.FS.mkdir(path);
    }
  }

  /**
   * Get the library version
   * @returns {string} Version string
   */
  getVersion() {
    const ptr = this.Module.ccall('get_version', 'string', [], []);
    return ptr;
  }

  /**
   * Calculate atomic radial charge density using Hartree-Fock/Dirac-Fock method
   *
   * @param {Object} params - Calculation parameters
   * @param {number} params.atomicNumber - Atomic number (Z)
   * @param {number} [params.relativity=1] - Relativity flag (0=non-rel, 1=rel)
   * @param {number} [params.exchangeAlpha=0.7] - Exchange correlation alpha
   * @param {Array<Object>} [params.orbitals] - Orbital configuration
   * @returns {Object} Calculation results including charge density
   */
  calculateChargeDensity(params) {
    this._ensureInitialized();

    const inputFile = this._generateAtorbInput(params);
    this.FS.writeFile('/input/atorb.i', inputFile);

    // Run hartfock calculation
    const result = this.Module.ccall(
      'run_hartfock',
      'number',
      ['string'],
      ['/input/atorb.i'],
    );

    if (result !== 0) {
      throw new Error(`Hartree-Fock calculation failed with code ${result}`);
    }

    // Read output files
    return this._parseChargeDensityOutput();
  }

  /**
   * Calculate phase shifts using the specified method
   *
   * @param {Object} params - Calculation parameters
   * @param {string} [params.method='rel'] - Method: 'rel', 'cav', or 'wil'
   * @param {number} params.atomicNumber - Atomic number
   * @param {number} params.muffinTinRadius - Muffin-tin radius in Bohr
   * @param {number} [params.energyMin=1.0] - Minimum energy in Hartrees
   * @param {number} [params.energyMax=12.0] - Maximum energy in Hartrees
   * @param {number} [params.energyStep=0.25] - Energy step in Hartrees
   * @param {number} [params.lmax=10] - Maximum angular momentum quantum number
   * @param {string} [params.potentialFile] - Pre-calculated potential file content
   * @returns {Object} Phase shift results
   */
  calculatePhaseShifts(params) {
    this._ensureInitialized();

    const method = (params.method || 'rel').toLowerCase();

    // Generate or use provided potential
    if (params.potentialFile) {
      this.FS.writeFile('/input/mufftin.o', params.potentialFile);
    } else {
      // Need to run hartfock first to generate potential
      this.calculateChargeDensity(params);
    }

    // Generate phase shift input file
    const phshInput = this._generatePhshInput(params);
    this.FS.writeFile('/input/phsh.i', phshInput);

    // Run the appropriate phase shift calculation
    let result;
    switch (method) {
      case 'rel':
        result = this.Module.ccall('run_phsh_rel', 'number', [], []);
        break;
      case 'cav':
        result = this.Module.ccall('run_phsh_cav', 'number', [], []);
        break;
      case 'wil':
        result = this.Module.ccall('run_phsh_wil', 'number', [], []);
        break;
      default:
        throw new Error(
          `Unknown method: ${method}. Use 'rel', 'cav', or 'wil'.`,
        );
    }

    if (result !== 0) {
      throw new Error(
        `Phase shift calculation (${method}) failed with code ${result}`,
      );
    }

    return this._parsePhaseShiftOutput(params);
  }

  /**
   * Calculate phase shifts using relativistic method
   * @param {Object} params - Same as calculatePhaseShifts
   * @returns {Object} Phase shift results
   */
  calculateRelativistic(params) {
    return this.calculatePhaseShifts({ ...params, method: 'rel' });
  }

  /**
   * Calculate phase shifts using cavity method
   * @param {Object} params - Same as calculatePhaseShifts
   * @returns {Object} Phase shift results
   */
  calculateCavity(params) {
    return this.calculatePhaseShifts({ ...params, method: 'cav' });
  }

  /**
   * Calculate phase shifts using Williams method
   * @param {Object} params - Same as calculatePhaseShifts
   * @returns {Object} Phase shift results
   */
  calculateWilliams(params) {
    return this.calculatePhaseShifts({ ...params, method: 'wil' });
  }

  /**
   * Write a file to the virtual filesystem
   * @param {string} path - File path in MEMFS
   * @param {string|Uint8Array} content - File content
   */
  writeFile(path, content) {
    this._ensureInitialized();
    this.FS.writeFile(path, content);
  }

  /**
   * Read a file from the virtual filesystem
   * @param {string} path - File path in MEMFS
   * @param {Object} [options] - Read options
   * @param {string} [options.encoding='utf8'] - File encoding
   * @returns {string|Uint8Array} File content
   */
  readFile(path, options = { encoding: 'utf8' }) {
    this._ensureInitialized();
    return this.FS.readFile(path, options);
  }

  /**
   * List files in a directory
   * @param {string} path - Directory path
   * @returns {string[]} List of filenames
   */
  listDir(path) {
    this._ensureInitialized();
    return this.FS.readdir(path).filter((f) => f !== '.' && f !== '..');
  }

  /**
   * Ensure the module is initialized
   * @private
   */
  _ensureInitialized() {
    if (!this._initialized) {
      throw new Error('PhaseShifts not initialized. Call init() first.');
    }
  }

  /**
   * Generate atorb input file content
   * @private
   */
  _generateAtorbInput(params) {
    const atomicNumber = params.atomicNumber;
    const relativity = params.relativity === undefined ? 1 : params.relativity;
    const exchangeAlpha =
      params.exchangeAlpha === undefined ? 0.7 : params.exchangeAlpha;

    let input = `d                           ! Dirac calculation\n`;
    input += `${relativity}                           ! relativity flag\n`;
    input += `x                           ! exchange correlation\n`;
    input += `${exchangeAlpha.toFixed(4)}                     ! alpha\n`;
    input += `i                           ! initialize\n`;
    input += `${atomicNumber.toFixed(1)}                       ! atomic number\n`;
    input += `a                           ! ab initio calculation\n`;
    input += `w                           ! write output\n`;
    input += `q                           ! quit\n`;

    return input;
  }

  /**
   * Generate phase shift input file content
   * @private
   */
  _generatePhshInput(params) {
    const atomicNumber = params.atomicNumber;
    const muffinTinRadius =
      params.muffinTinRadius === undefined ? 2.5 : params.muffinTinRadius;
    const energyMin = params.energyMin === undefined ? 1.0 : params.energyMin;
    const energyMax = params.energyMax === undefined ? 12.0 : params.energyMax;
    const energyStep =
      params.energyStep === undefined ? 0.25 : params.energyStep;
    const lmax = params.lmax === undefined ? 10 : params.lmax;

    // Number of energy points
    const nEnergies = Math.floor((energyMax - energyMin) / energyStep) + 1;

    return `${atomicNumber.toFixed(1)} ${muffinTinRadius.toFixed(
      4,
    )} ${energyMin.toFixed(4)} ${energyMax.toFixed(4)} ${energyStep.toFixed(
      4,
    )} ${lmax} ${nEnergies}\n`;
  }

  /**
   * Parse charge density output
   * @private
   */
  _parseChargeDensityOutput() {
    const outputPath = '/output/atorb.o';
    if (!this.FS.analyzePath(outputPath).exists) {
      const parsed = {
        raw: '',
        success: false,
        error: 'Charge density output not found.',
      };
      return parsed;
    }

    const output = this.FS.readFile(outputPath, { encoding: 'utf8' });
    const parsed = {
      raw: output,
      success: true,
      // Additional parsing can be added here
    };
    return parsed;
  }

  /**
   * Parse phase shift output
   * @private
   */
  _parsePhaseShiftOutput(params) {
    const lmax = params.lmax === undefined ? 10 : params.lmax;
    const energyMin = params.energyMin === undefined ? 1.0 : params.energyMin;
    const energyMax = params.energyMax === undefined ? 12.0 : params.energyMax;
    const energyStep =
      params.energyStep === undefined ? 0.25 : params.energyStep;

    const outputPath = '/output/phasout.o';
    if (!this.FS.analyzePath(outputPath).exists) {
      const parsed = {
        raw: '',
        success: false,
        error: 'Phase shift output not found.',
        energies: [],
        phaseShifts: [],
      };
      return parsed;
    }

    const output = this.FS.readFile(outputPath, {
      encoding: 'utf8',
    });

    // Parse phase shifts from output
    const phaseShifts = parsePhaseShiftData(output, lmax);

    const parsed = {
      raw: output,
      success: true,
      energies: phaseShifts.energies,
      phaseShifts: phaseShifts.data,
      lmax: lmax,
      energyRange: { min: energyMin, max: energyMax, step: energyStep },
    };
    return parsed;
  }

  /**
   * Get default orbital configuration for an element
   * @private
   */
  _getDefaultOrbitals(Z) {
    // Simplified orbital filling based on atomic number
    // Full implementation would use periodic table data
    const orbitals = [];
    let remaining = Z;

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
    ];

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
}

/**
 * Create and initialize a PhaseShifts calculator
 *
 * @param {Object} [options] - Configuration options
 * @param {string} [options.wasmPath='./phaseshifts.js'] - Path to WASM loader
 * @returns {Promise<PhaseShifts>} Initialized PhaseShifts instance
 *
 * @example
 * const phsh = await createPhaseShifts();
 * const result = await phsh.calculatePhaseShifts({
 *     atomicNumber: 29,  // Copper
 *     muffinTinRadius: 2.5,
 *     method: 'rel'
 * });
 */
async function createPhaseShifts(options = {}) {
  const wasmPath = options.wasmPath || './phaseshifts.js';

  // Dynamically load the WASM module
  // In browser, this relies on the Emscripten-generated loader
  let createModule;

  if (typeof createPhaseShiftsModule !== 'undefined') {
    // Module already loaded (via script tag)
    createModule = createPhaseShiftsModule;
  } else {
    // Dynamic import
    const script = await import(/* webpackChunkName: 'phaseshifts' */ wasmPath);
    createModule = script.default || script.createPhaseShiftsModule;
  }

  const wasmModule = await createModule();
  const phsh = new PhaseShifts(wasmModule);
  await phsh.init();

  return phsh;
}

// Export for different module systems
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    PhaseShifts,
    createPhaseShifts,
    elements,
    ELEMENTS: elements,
  };
}

if (typeof window !== 'undefined') {
  window.PhaseShifts = PhaseShifts;
  window.createPhaseShifts = createPhaseShifts;
  window.elements = elements;
  window.ELEMENTS = elements;
}

export { PhaseShifts, createPhaseShifts, elements };
export { elements as ELEMENTS };
