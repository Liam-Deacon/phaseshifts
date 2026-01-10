/**
 * Crystal Structure Builder UI Component
 * Provides interactive UI for building and editing crystal structures
 */

import {
  CrystalStructure,
  Layer,
  Atom,
  Vector3,
  UnitCell3D,
  MillerIndices,
  getPresetNames,
  createFromPreset,
} from './crystal.js';
import { elements } from './elements.js';

/**
 * Structure Builder UI Component
 */
export class StructureBuilder {
  constructor(container, options = {}) {
    this.container =
      typeof container === 'string'
        ? document.getElementById(container)
        : container;

    if (!this.container) {
      throw new Error('Container element not found');
    }

    this.options = {
      onStructureChange: options.onStructureChange || (() => {}),
      showPresets: options.showPresets ?? true,
      showImportExport: options.showImportExport ?? true,
      compact: options.compact ?? false,
      ...options,
    };

    this.structure = new CrystalStructure({ name: 'New Structure' });
    this._selectedLayerIndex = -1;
    this._selectedAtomIndex = -1;

    this._render();
    this._attachEventListeners();
  }

  /**
   * Get the current crystal structure
   * @returns {CrystalStructure}
   */
  getStructure() {
    return this.structure;
  }

  /**
   * Set the crystal structure
   * @param {CrystalStructure} structure
   */
  setStructure(structure) {
    this.structure = structure;
    this._render();
    this._notifyChange();
  }

  /**
   * Load a preset structure
   * @param {string} presetName
   */
  loadPreset(presetName) {
    try {
      this.structure = createFromPreset(presetName);
      this._render();
      this._notifyChange();
    } catch (error) {
      console.error('Failed to load preset:', error);
      alert(`Failed to load preset: ${error.message}`);
    }
  }

  _notifyChange() {
    this.options.onStructureChange(this.structure);
  }

  _render() {
    this.container.innerHTML = '';
    this.container.className =
      'structure-builder' + (this.options.compact ? ' compact' : '');

    // Build main sections
    if (this.options.showPresets) {
      this.container.appendChild(this._createPresetsSection());
    }

    this.container.appendChild(this._createUnitCellSection());
    this.container.appendChild(this._createLayersSection());

    if (this.options.showImportExport) {
      this.container.appendChild(this._createImportExportSection());
    }
  }

  _createPresetsSection() {
    const section = document.createElement('div');
    section.className = 'builder-section presets-section';

    const header = document.createElement('h3');
    header.textContent = 'Load Preset Structure';
    section.appendChild(header);

    const grid = document.createElement('div');
    grid.className = 'presets-grid';

    const presets = getPresetNames();
    for (const name of presets) {
      const btn = document.createElement('button');
      btn.className = 'preset-btn';
      btn.textContent = name;
      btn.dataset.preset = name;
      btn.addEventListener('click', () => this.loadPreset(name));
      grid.appendChild(btn);
    }

    section.appendChild(grid);
    return section;
  }

  _createUnitCellSection() {
    const section = document.createElement('div');
    section.className = 'builder-section unit-cell-section';

    const header = document.createElement('h3');
    header.textContent = 'Unit Cell Parameters';
    section.appendChild(header);

    // Bulk unit cell
    const bulkGroup = document.createElement('div');
    bulkGroup.className = 'form-group';

    const bulkLabel = document.createElement('label');
    bulkLabel.textContent = 'Bulk Lattice Constant (Å)';
    bulkGroup.appendChild(bulkLabel);

    const bulkInputs = document.createElement('div');
    bulkInputs.className = 'input-row';

    const createInput = (label, value, onChange) => {
      const wrapper = document.createElement('div');
      wrapper.className = 'input-wrapper';

      const lbl = document.createElement('span');
      lbl.textContent = label;
      wrapper.appendChild(lbl);

      const input = document.createElement('input');
      input.type = 'number';
      input.step = '0.001';
      input.value = value;
      input.addEventListener('change', onChange);
      wrapper.appendChild(input);

      return wrapper;
    };

    bulkInputs.appendChild(
      createInput('a:', this.structure.bulkUnitCell.aLength.toFixed(3), (e) => {
        const value = Number.parseFloat(e.target.value);
        if (Number.isNaN(value)) return;
        this.structure.bulkUnitCell = new UnitCell3D(
          value,
          this.structure.bulkUnitCell.bLength,
          this.structure.bulkUnitCell.cLength,
        );
        this._notifyChange();
      }),
    );

    bulkInputs.appendChild(
      createInput('b:', this.structure.bulkUnitCell.bLength.toFixed(3), (e) => {
        const value = Number.parseFloat(e.target.value);
        if (Number.isNaN(value)) return;
        this.structure.bulkUnitCell = new UnitCell3D(
          this.structure.bulkUnitCell.aLength,
          value,
          this.structure.bulkUnitCell.cLength,
        );
        this._notifyChange();
      }),
    );

    bulkInputs.appendChild(
      createInput('c:', this.structure.bulkUnitCell.cLength.toFixed(3), (e) => {
        const value = Number.parseFloat(e.target.value);
        if (Number.isNaN(value)) return;
        this.structure.bulkUnitCell = new UnitCell3D(
          this.structure.bulkUnitCell.aLength,
          this.structure.bulkUnitCell.bLength,
          value,
        );
        this._notifyChange();
      }),
    );

    bulkGroup.appendChild(bulkInputs);
    section.appendChild(bulkGroup);

    // Miller indices
    const millerGroup = document.createElement('div');
    millerGroup.className = 'form-group';

    const millerLabel = document.createElement('label');
    millerLabel.textContent = 'Surface Orientation (Miller Indices)';
    millerGroup.appendChild(millerLabel);

    const millerInputs = document.createElement('div');
    millerInputs.className = 'input-row';

    const mi = this.structure.millerIndices;
    millerInputs.appendChild(
      createInput('h:', mi.h, (e) => {
        const value = Number.parseInt(e.target.value, 10);
        if (Number.isNaN(value)) return;
        this.structure.millerIndices = new MillerIndices(
          value,
          this.structure.millerIndices.k,
          this.structure.millerIndices.l,
        );
        this._notifyChange();
      }),
    );

    millerInputs.appendChild(
      createInput('k:', mi.k, (e) => {
        const value = Number.parseInt(e.target.value, 10);
        if (Number.isNaN(value)) return;
        this.structure.millerIndices = new MillerIndices(
          this.structure.millerIndices.h,
          value,
          this.structure.millerIndices.l,
        );
        this._notifyChange();
      }),
    );

    millerInputs.appendChild(
      createInput('l:', mi.l, (e) => {
        const value = Number.parseInt(e.target.value, 10);
        if (Number.isNaN(value)) return;
        this.structure.millerIndices = new MillerIndices(
          this.structure.millerIndices.h,
          this.structure.millerIndices.k,
          value,
        );
        this._notifyChange();
      }),
    );

    millerGroup.appendChild(millerInputs);
    section.appendChild(millerGroup);

    // Structure name
    const nameGroup = document.createElement('div');
    nameGroup.className = 'form-group';

    const nameLabel = document.createElement('label');
    nameLabel.textContent = 'Structure Name';
    nameGroup.appendChild(nameLabel);

    const nameInput = document.createElement('input');
    nameInput.type = 'text';
    nameInput.value = this.structure.name;
    nameInput.addEventListener('change', (e) => {
      this.structure.name = e.target.value;
      this._notifyChange();
    });
    nameGroup.appendChild(nameInput);

    section.appendChild(nameGroup);

    return section;
  }

  _createLayersSection() {
    const section = document.createElement('div');
    section.className = 'builder-section layers-section';

    const header = document.createElement('div');
    header.className = 'section-header';

    const title = document.createElement('h3');
    title.textContent = 'Layers';
    header.appendChild(title);

    const addBtn = document.createElement('button');
    addBtn.className = 'btn btn-small btn-primary';
    addBtn.textContent = '+ Add Layer';
    addBtn.addEventListener('click', () => this._addLayer());
    header.appendChild(addBtn);

    section.appendChild(header);

    // Layers list
    const layersList = document.createElement('div');
    layersList.className = 'layers-list';

    if (this.structure.layers.length === 0) {
      const emptyMsg = document.createElement('p');
      emptyMsg.className = 'empty-message';
      emptyMsg.textContent =
        'No layers defined. Click "Add Layer" to start building.';
      layersList.appendChild(emptyMsg);
    } else {
      this.structure.layers.forEach((layer, index) => {
        layersList.appendChild(this._createLayerCard(layer, index));
      });
    }

    section.appendChild(layersList);
    return section;
  }

  _createLayerCard(layer, index) {
    const card = document.createElement('div');
    card.className =
      'layer-card' + (index === this._selectedLayerIndex ? ' selected' : '');
    card.dataset.layerIndex = index;

    // Layer header
    const header = document.createElement('div');
    header.className = 'layer-header';

    const nameInput = document.createElement('input');
    nameInput.type = 'text';
    nameInput.value = layer.name;
    nameInput.className = 'layer-name-input';
    nameInput.addEventListener('change', (e) => {
      layer.name = e.target.value;
      this._notifyChange();
    });
    header.appendChild(nameInput);

    const zInput = document.createElement('input');
    zInput.type = 'number';
    zInput.step = '0.01';
    zInput.value = layer.zPosition.toFixed(3);
    zInput.className = 'layer-z-input';
    zInput.title = 'Z position (Å)';
    zInput.addEventListener('change', (e) => {
      const value = Number.parseFloat(e.target.value);
      if (Number.isNaN(value)) return;
      layer.zPosition = value;
      this._notifyChange();
    });
    header.appendChild(zInput);

    const deleteBtn = document.createElement('button');
    deleteBtn.className = 'btn btn-small btn-danger';
    deleteBtn.textContent = '×';
    deleteBtn.title = 'Delete layer';
    deleteBtn.addEventListener('click', (e) => {
      e.stopPropagation();
      this._removeLayer(index);
    });
    header.appendChild(deleteBtn);

    card.appendChild(header);

    // Atoms in layer
    const atomsList = document.createElement('div');
    atomsList.className = 'atoms-list';

    if (layer.atoms.length === 0) {
      const emptyMsg = document.createElement('span');
      emptyMsg.className = 'empty-atoms';
      emptyMsg.textContent = 'No atoms';
      atomsList.appendChild(emptyMsg);
    } else {
      layer.atoms.forEach((atom, atomIndex) => {
        atomsList.appendChild(this._createAtomTag(atom, index, atomIndex));
      });
    }

    const addAtomBtn = document.createElement('button');
    addAtomBtn.className = 'btn btn-tiny';
    addAtomBtn.textContent = '+';
    addAtomBtn.title = 'Add atom';
    addAtomBtn.addEventListener('click', () => this._showAddAtomDialog(index));
    atomsList.appendChild(addAtomBtn);

    card.appendChild(atomsList);

    // Click to select layer
    card.addEventListener('click', () => {
      this._selectedLayerIndex = index;
      this._render();
    });

    return card;
  }

  _createAtomTag(atom, layerIndex, atomIndex) {
    const tag = document.createElement('span');
    tag.className = 'atom-tag';
    tag.style.backgroundColor = atom.getColor();
    tag.style.color = this._getContrastColor(atom.getColor());
    tag.textContent = `${atom.symbol} (${atom.position.x.toFixed(2)}, ${atom.position.y.toFixed(2)})`;
    tag.title = `${atom.symbol} at (${atom.position.x.toFixed(3)}, ${atom.position.y.toFixed(3)}, ${atom.position.z.toFixed(3)})`;

    const removeBtn = document.createElement('span');
    removeBtn.className = 'atom-remove';
    removeBtn.textContent = '×';
    removeBtn.addEventListener('click', (e) => {
      e.stopPropagation();
      this._removeAtom(layerIndex, atomIndex);
    });
    tag.appendChild(removeBtn);

    return tag;
  }

  _getContrastColor(hexColor) {
    // Convert hex to RGB
    const hex = hexColor.replace('#', '');
    const r = Number.parseInt(hex.substr(0, 2), 16);
    const g = Number.parseInt(hex.substr(2, 2), 16);
    const b = Number.parseInt(hex.substr(4, 2), 16);

    // Calculate luminance
    const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;

    return luminance > 0.5 ? '#000000' : '#FFFFFF';
  }

  _addLayer() {
    const zPos =
      this.structure.layers.length > 0
        ? this.structure.layers[this.structure.layers.length - 1].zPosition -
          2
        : 0;

    const layer = new Layer([], {
      name: `Layer ${this.structure.layers.length + 1}`,
      zPosition: zPos,
    });

    this.structure.addLayer(layer);
    this._render();
    this._notifyChange();
  }

  _removeLayer(index) {
    if (confirm(`Delete "${this.structure.layers[index].name}"?`)) {
      this.structure.removeLayer(index);
      this._selectedLayerIndex = -1;
      this._render();
      this._notifyChange();
    }
  }

  _showAddAtomDialog(layerIndex) {
    // Create modal dialog
    const modal = document.createElement('div');
    modal.className = 'modal-overlay';

    const dialog = document.createElement('div');
    dialog.className = 'modal-dialog';

    const title = document.createElement('h4');
    title.textContent = 'Add Atom';
    dialog.appendChild(title);

    // Element select
    const elementGroup = document.createElement('div');
    elementGroup.className = 'form-group';

    const elementLabel = document.createElement('label');
    elementLabel.textContent = 'Element';
    elementGroup.appendChild(elementLabel);

    const elementSelect = document.createElement('select');
    Object.entries(elements).forEach(([symbol, z]) => {
      const opt = document.createElement('option');
      opt.value = symbol;
      opt.textContent = `${symbol} (Z=${z})`;
      elementSelect.appendChild(opt);
    });
    elementGroup.appendChild(elementSelect);
    dialog.appendChild(elementGroup);

    // Position inputs
    const posGroup = document.createElement('div');
    posGroup.className = 'form-group';

    const posLabel = document.createElement('label');
    posLabel.textContent = 'Position (x, y in fractional or Å)';
    posGroup.appendChild(posLabel);

    const posInputs = document.createElement('div');
    posInputs.className = 'input-row';

    const xInput = document.createElement('input');
    xInput.type = 'number';
    xInput.step = '0.001';
    xInput.value = '0';
    xInput.placeholder = 'x';
    posInputs.appendChild(xInput);

    const yInput = document.createElement('input');
    yInput.type = 'number';
    yInput.step = '0.001';
    yInput.value = '0';
    yInput.placeholder = 'y';
    posInputs.appendChild(yInput);

    posGroup.appendChild(posInputs);
    dialog.appendChild(posGroup);

    // Buttons
    const buttons = document.createElement('div');
    buttons.className = 'dialog-buttons';

    const cancelBtn = document.createElement('button');
    cancelBtn.className = 'btn btn-secondary';
    cancelBtn.textContent = 'Cancel';
    cancelBtn.addEventListener('click', () => modal.remove());
    buttons.appendChild(cancelBtn);

    const addBtn = document.createElement('button');
    addBtn.className = 'btn btn-primary';
    addBtn.textContent = 'Add';
    addBtn.addEventListener('click', () => {
      const symbol = elementSelect.value;
      const x = Number.parseFloat(xInput.value) || 0;
      const y = Number.parseFloat(yInput.value) || 0;
      const layer = this.structure.layers[layerIndex];
      const z = layer.zPosition;

      const atom = new Atom(symbol, new Vector3(x, y, z));
      layer.addAtom(atom);

      modal.remove();
      this._render();
      this._notifyChange();
    });
    buttons.appendChild(addBtn);

    dialog.appendChild(buttons);
    modal.appendChild(dialog);
    document.body.appendChild(modal);

    // Close on background click
    modal.addEventListener('click', (e) => {
      if (e.target === modal) modal.remove();
    });
  }

  _removeAtom(layerIndex, atomIndex) {
    this.structure.layers[layerIndex].removeAtom(atomIndex);
    this._render();
    this._notifyChange();
  }

  _createImportExportSection() {
    const section = document.createElement('div');
    section.className = 'builder-section import-export-section';

    const header = document.createElement('h3');
    header.textContent = 'Import / Export';
    section.appendChild(header);

    const buttons = document.createElement('div');
    buttons.className = 'button-row';

    const exportBtn = document.createElement('button');
    exportBtn.className = 'btn btn-secondary';
    exportBtn.textContent = 'Export JSON';
    exportBtn.addEventListener('click', () => this._exportJSON());
    buttons.appendChild(exportBtn);

    const importBtn = document.createElement('button');
    importBtn.className = 'btn btn-secondary';
    importBtn.textContent = 'Import JSON';
    importBtn.addEventListener('click', () => this._importJSON());
    buttons.appendChild(importBtn);

    section.appendChild(buttons);
    return section;
  }

  _exportJSON() {
    const json = JSON.stringify(this.structure.toJSON(), null, 2);
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${this.structure.name.replaceAll(/[^a-z0-9]/gi, '_')}.json`;
    a.click();
    URL.revokeObjectURL(url);
  }

  _importJSON() {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.json';
    input.addEventListener('change', async (e) => {
      const file = e.target.files[0];
      if (!file) return;

      try {
        const text = await file.text();
        const json = JSON.parse(text);
        this.structure = CrystalStructure.fromJSON(json);
        this._render();
        this._notifyChange();
      } catch (error) {
        console.error('Failed to import structure:', error);
        alert(`Failed to import: ${error.message}`);
      }
    });
    input.click();
  }

  _attachEventListeners() {
    // Could add keyboard shortcuts, drag-drop, etc.
  }

  /**
   * Dispose and clean up
   */
  dispose() {
    this.container.innerHTML = '';
  }
}

/**
 * Create a structure builder instance
 * @param {string|HTMLElement} container
 * @param {Object} options
 * @returns {StructureBuilder}
 */
export function createStructureBuilder(container, options = {}) {
  return new StructureBuilder(container, options);
}
