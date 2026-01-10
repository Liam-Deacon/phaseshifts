/**
 * Crystal Structure Module
 * Defines crystal structures, unit cells, and surface models for LEED/XPD calculations
 */

import { elements, getAtomicRadius, getElementColor } from './elements.js';

/**
 * Represents a 3D vector
 */
export class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  clone() {
    return new Vector3(this.x, this.y, this.z);
  }

  add(v) {
    return new Vector3(this.x + v.x, this.y + v.y, this.z + v.z);
  }

  subtract(v) {
    return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
  }

  scale(s) {
    return new Vector3(this.x * s, this.y * s, this.z * s);
  }

  dot(v) {
    return this.x * v.x + this.y * v.y + this.z * v.z;
  }

  cross(v) {
    return new Vector3(
      this.y * v.z - this.z * v.y,
      this.z * v.x - this.x * v.z,
      this.x * v.y - this.y * v.x,
    );
  }

  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }

  normalize() {
    const len = this.length();
    if (len === 0) return new Vector3();
    return this.scale(1 / len);
  }

  toArray() {
    return [this.x, this.y, this.z];
  }

  static fromArray(arr) {
    return new Vector3(arr[0] || 0, arr[1] || 0, arr[2] || 0);
  }
}

/**
 * Represents an atom in the structure
 */
export class Atom {
  constructor(symbol, position, options = {}) {
    this.symbol = symbol;
    this.atomicNumber = elements[symbol] || 0;
    this.position =
      position instanceof Vector3 ? position : new Vector3(...position);
    this.label = options.label || symbol;
    this.occupancy = options.occupancy ?? 1.0;
    this.isVacancy = options.isVacancy ?? false;
    this.thermalFactor = options.thermalFactor ?? 0.0;
    this.muffinTinRadius = options.muffinTinRadius ?? null;
  }

  clone() {
    return new Atom(this.symbol, this.position.clone(), {
      label: this.label,
      occupancy: this.occupancy,
      isVacancy: this.isVacancy,
      thermalFactor: this.thermalFactor,
      muffinTinRadius: this.muffinTinRadius,
    });
  }

  translate(vector) {
    return new Atom(this.symbol, this.position.add(vector), {
      label: this.label,
      occupancy: this.occupancy,
      isVacancy: this.isVacancy,
      thermalFactor: this.thermalFactor,
      muffinTinRadius: this.muffinTinRadius,
    });
  }

  getRadius() {
    return getAtomicRadius(this.symbol);
  }

  getColor() {
    return getElementColor(this.symbol);
  }

  toJSON() {
    return {
      symbol: this.symbol,
      position: this.position.toArray(),
      label: this.label,
      occupancy: this.occupancy,
      isVacancy: this.isVacancy,
      thermalFactor: this.thermalFactor,
      muffinTinRadius: this.muffinTinRadius,
    };
  }

  static fromJSON(json) {
    return new Atom(json.symbol, Vector3.fromArray(json.position), {
      label: json.label,
      occupancy: json.occupancy,
      isVacancy: json.isVacancy,
      thermalFactor: json.thermalFactor,
      muffinTinRadius: json.muffinTinRadius,
    });
  }
}

/**
 * Represents a layer of atoms (for surface calculations)
 */
export class Layer {
  constructor(atoms = [], options = {}) {
    this.atoms = atoms.map((a) => (a instanceof Atom ? a : Atom.fromJSON(a)));
    this.name = options.name || 'Layer';
    this.zPosition = options.zPosition ?? 0;
    this.interlayerSpacing = options.interlayerSpacing ?? null;
    this.isComposite = options.isComposite ?? false;
  }

  addAtom(atom) {
    this.atoms.push(atom instanceof Atom ? atom : Atom.fromJSON(atom));
  }

  removeAtom(index) {
    this.atoms.splice(index, 1);
  }

  clone() {
    return new Layer(
      this.atoms.map((a) => a.clone()),
      {
        name: this.name,
        zPosition: this.zPosition,
        interlayerSpacing: this.interlayerSpacing,
        isComposite: this.isComposite,
      },
    );
  }

  translate(vector) {
    return new Layer(
      this.atoms.map((a) => a.translate(vector)),
      {
        name: this.name,
        zPosition: this.zPosition + (vector.z || 0),
        interlayerSpacing: this.interlayerSpacing,
        isComposite: this.isComposite,
      },
    );
  }

  toJSON() {
    return {
      atoms: this.atoms.map((a) => a.toJSON()),
      name: this.name,
      zPosition: this.zPosition,
      interlayerSpacing: this.interlayerSpacing,
      isComposite: this.isComposite,
    };
  }

  static fromJSON(json) {
    return new Layer(
      json.atoms.map((a) => Atom.fromJSON(a)),
      {
        name: json.name,
        zPosition: json.zPosition,
        interlayerSpacing: json.interlayerSpacing,
        isComposite: json.isComposite,
      },
    );
  }
}

/**
 * Represents a 2D unit cell (surface)
 */
export class UnitCell2D {
  constructor(a, b, gamma = 90) {
    // a and b are Vector3 but we only use x,y components for 2D
    this.a = a instanceof Vector3 ? a : new Vector3(a, 0, 0);
    this.b = b instanceof Vector3 ? b : new Vector3(0, b, 0);
    this.gamma = gamma; // angle between a and b in degrees
    this._recalculateB();
  }

  _recalculateB() {
    // Recalculate b vector based on gamma angle
    if (this.gamma !== 90) {
      const bLength = this.b.length() || 1;
      const gammaRad = (this.gamma * Math.PI) / 180;
      this.b = new Vector3(
        bLength * Math.cos(gammaRad),
        bLength * Math.sin(gammaRad),
        0,
      );
    }
  }

  getArea() {
    return Math.abs(this.a.x * this.b.y - this.a.y * this.b.x);
  }

  toJSON() {
    return {
      a: this.a.toArray(),
      b: this.b.toArray(),
      gamma: this.gamma,
    };
  }

  static fromJSON(json) {
    const cell = new UnitCell2D(
      Vector3.fromArray(json.a),
      Vector3.fromArray(json.b),
      json.gamma,
    );
    return cell;
  }
}

/**
 * Represents a 3D unit cell (bulk crystal)
 */
export class UnitCell3D {
  constructor(a, b, c, alpha = 90, beta = 90, gamma = 90) {
    this.aLength = typeof a === 'number' ? a : a.length();
    this.bLength = typeof b === 'number' ? b : b.length();
    this.cLength = typeof c === 'number' ? c : c.length();
    this.alpha = alpha; // angle between b and c
    this.beta = beta; // angle between a and c
    this.gamma = gamma; // angle between a and b
    this._calculateVectors();
  }

  _calculateVectors() {
    // Calculate Cartesian vectors from lattice parameters
    const alphaRad = (this.alpha * Math.PI) / 180;
    const betaRad = (this.beta * Math.PI) / 180;
    const gammaRad = (this.gamma * Math.PI) / 180;

    this.a = new Vector3(this.aLength, 0, 0);

    this.b = new Vector3(
      this.bLength * Math.cos(gammaRad),
      this.bLength * Math.sin(gammaRad),
      0,
    );

    const cx = this.cLength * Math.cos(betaRad);
    const cy =
      (this.cLength *
        (Math.cos(alphaRad) - Math.cos(betaRad) * Math.cos(gammaRad))) /
      Math.sin(gammaRad);
    // Guard against negative discriminant due to floating-point errors
    const czSquared = this.cLength * this.cLength - cx * cx - cy * cy;
    const cz = Math.sqrt(Math.max(0, czSquared));
    this.c = new Vector3(cx, cy, cz);
  }

  getVolume() {
    return Math.abs(this.a.dot(this.b.cross(this.c)));
  }

  toJSON() {
    return {
      a: this.aLength,
      b: this.bLength,
      c: this.cLength,
      alpha: this.alpha,
      beta: this.beta,
      gamma: this.gamma,
    };
  }

  static fromJSON(json) {
    return new UnitCell3D(
      json.a,
      json.b,
      json.c,
      json.alpha,
      json.beta,
      json.gamma,
    );
  }
}

/**
 * Miller indices for surface orientation
 */
export class MillerIndices {
  constructor(h, k, l) {
    this.h = h;
    this.k = k;
    this.l = l;
  }

  toString() {
    const format = (n) => (n < 0 ? `\\overline{${Math.abs(n)}}` : `${n}`);
    return `(${format(this.h)}${format(this.k)}${format(this.l)})`;
  }

  toSimpleString() {
    const format = (n) => (n < 0 ? `-${Math.abs(n)}` : `${n}`);
    return `(${format(this.h)}${format(this.k)}${format(this.l)})`;
  }

  toArray() {
    return [this.h, this.k, this.l];
  }

  static fromArray(arr) {
    return new MillerIndices(arr[0], arr[1], arr[2]);
  }
}

/**
 * Crystal structure types
 */
export const CrystalSystem = Object.freeze({
  CUBIC: 'cubic',
  TETRAGONAL: 'tetragonal',
  ORTHORHOMBIC: 'orthorhombic',
  HEXAGONAL: 'hexagonal',
  TRIGONAL: 'trigonal',
  MONOCLINIC: 'monoclinic',
  TRICLINIC: 'triclinic',
});

/**
 * Common Bravais lattice types
 */
export const BravaisLattice = Object.freeze({
  PRIMITIVE: 'P',
  BODY_CENTERED: 'I',
  FACE_CENTERED: 'F',
  BASE_CENTERED: 'C',
  RHOMBOHEDRAL: 'R',
});

/**
 * Represents a complete crystal structure
 */
export class CrystalStructure {
  constructor(options = {}) {
    this.name = options.name || 'Unnamed Structure';
    this.bulkUnitCell = options.bulkUnitCell || new UnitCell3D(1, 1, 1);
    this.surfaceUnitCell = options.surfaceUnitCell || null;
    this.millerIndices = options.millerIndices || new MillerIndices(1, 0, 0);
    this.layers = options.layers || [];
    this.crystalSystem = options.crystalSystem || CrystalSystem.CUBIC;
    this.bravaisLattice = options.bravaisLattice || BravaisLattice.PRIMITIVE;
    this.bulkAtoms = options.bulkAtoms || []; // Atoms in bulk unit cell
  }

  addLayer(layer) {
    this.layers.push(layer instanceof Layer ? layer : Layer.fromJSON(layer));
  }

  removeLayer(index) {
    this.layers.splice(index, 1);
  }

  moveLayer(fromIndex, toIndex) {
    const [layer] = this.layers.splice(fromIndex, 1);
    this.layers.splice(toIndex, 0, layer);
  }

  /**
   * Get all atoms for visualization (expand layers into atoms)
   */
  getAllAtoms(repeatX = 1, repeatY = 1) {
    const atoms = [];
    const cell =
      this.surfaceUnitCell ||
      new UnitCell2D(
        new Vector3(this.bulkUnitCell.aLength, 0, 0),
        new Vector3(0, this.bulkUnitCell.bLength, 0),
      );

    for (let ix = 0; ix < repeatX; ix++) {
      for (let iy = 0; iy < repeatY; iy++) {
        const offset = cell.a.scale(ix).add(cell.b.scale(iy));
        for (const layer of this.layers) {
          for (const atom of layer.atoms) {
            atoms.push(atom.translate(offset));
          }
        }
      }
    }
    return atoms;
  }

  /**
   * Get unique elements in the structure
   */
  getUniqueElements() {
    const elements = new Set();
    for (const layer of this.layers) {
      for (const atom of layer.atoms) {
        elements.add(atom.symbol);
      }
    }
    return Array.from(elements);
  }

  toJSON() {
    return {
      name: this.name,
      bulkUnitCell: this.bulkUnitCell.toJSON(),
      surfaceUnitCell: this.surfaceUnitCell?.toJSON() || null,
      millerIndices: this.millerIndices.toArray(),
      layers: this.layers.map((l) => l.toJSON()),
      crystalSystem: this.crystalSystem,
      bravaisLattice: this.bravaisLattice,
      bulkAtoms: this.bulkAtoms.map((a) => a.toJSON()),
    };
  }

  static fromJSON(json) {
    return new CrystalStructure({
      name: json.name,
      bulkUnitCell: UnitCell3D.fromJSON(json.bulkUnitCell),
      surfaceUnitCell: json.surfaceUnitCell
        ? UnitCell2D.fromJSON(json.surfaceUnitCell)
        : null,
      millerIndices: MillerIndices.fromArray(json.millerIndices),
      layers: json.layers.map((l) => Layer.fromJSON(l)),
      crystalSystem: json.crystalSystem,
      bravaisLattice: json.bravaisLattice,
      bulkAtoms: json.bulkAtoms?.map((a) => Atom.fromJSON(a)) || [],
    });
  }
}

/**
 * Predefined crystal structure presets
 */
export const crystalPresets = Object.freeze({
  // FCC Surfaces
  'Cu(111)': () => {
    const a = 3.615; // Angstroms
    const d111 = a / Math.sqrt(3); // Interlayer spacing for (111)
    const asurf = a / Math.sqrt(2); // Surface lattice constant

    return new CrystalStructure({
      name: 'Cu(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Cu', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Cu',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Cu',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },

  'Cu(100)': () => {
    const a = 3.615;
    const d100 = a / 2;
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Cu(100)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(0, asurf, 0),
      ),
      millerIndices: new MillerIndices(1, 0, 0),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Cu', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d100,
        }),
        new Layer([new Atom('Cu', new Vector3(asurf / 2, asurf / 2, -d100))], {
          name: 'Layer 2',
          zPosition: -d100,
          interlayerSpacing: d100,
        }),
      ],
    });
  },

  'Ni(111)': () => {
    const a = 3.524;
    const d111 = a / Math.sqrt(3);
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Ni(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Ni', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Ni',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Ni',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },

  'Ni(100)': () => {
    const a = 3.524;
    const d100 = a / 2;
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Ni(100)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(0, asurf, 0),
      ),
      millerIndices: new MillerIndices(1, 0, 0),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Ni', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d100,
        }),
        new Layer([new Atom('Ni', new Vector3(asurf / 2, asurf / 2, -d100))], {
          name: 'Layer 2',
          zPosition: -d100,
          interlayerSpacing: d100,
        }),
      ],
    });
  },

  // BCC Surfaces
  'Fe(110)': () => {
    const a = 2.87;
    const d110 = a / Math.sqrt(2);
    const asurfX = a;
    const asurfY = a * Math.sqrt(2);

    return new CrystalStructure({
      name: 'Fe(110)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurfX, 0, 0),
        new Vector3(0, asurfY, 0),
      ),
      millerIndices: new MillerIndices(1, 1, 0),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.BODY_CENTERED,
      layers: [
        new Layer([new Atom('Fe', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d110,
        }),
        new Layer(
          [new Atom('Fe', new Vector3(asurfX / 2, asurfY / 2, -d110))],
          { name: 'Layer 2', zPosition: -d110, interlayerSpacing: d110 },
        ),
      ],
    });
  },

  'Fe(100)': () => {
    const a = 2.87;
    const d100 = a / 2;

    return new CrystalStructure({
      name: 'Fe(100)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(a, 0, 0),
        new Vector3(0, a, 0),
      ),
      millerIndices: new MillerIndices(1, 0, 0),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.BODY_CENTERED,
      layers: [
        new Layer([new Atom('Fe', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d100,
        }),
        new Layer([new Atom('Fe', new Vector3(a / 2, a / 2, -d100))], {
          name: 'Layer 2',
          zPosition: -d100,
          interlayerSpacing: d100,
        }),
      ],
    });
  },

  // Diamond Surface
  'Si(100)': () => {
    const a = 5.431;
    const d100 = a / 4;
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Si(100)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(0, asurf, 0),
      ),
      millerIndices: new MillerIndices(1, 0, 0),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Si', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d100,
        }),
        new Layer([new Atom('Si', new Vector3(asurf / 2, asurf / 2, -d100))], {
          name: 'Layer 2',
          zPosition: -d100,
          interlayerSpacing: d100,
        }),
      ],
    });
  },

  // Gold surfaces
  'Au(111)': () => {
    const a = 4.078;
    const d111 = a / Math.sqrt(3);
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Au(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Au', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Au',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Au',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },

  // Platinum
  'Pt(111)': () => {
    const a = 3.924;
    const d111 = a / Math.sqrt(3);
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Pt(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Pt', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Pt',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Pt',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },

  // Aluminum
  'Al(111)': () => {
    const a = 4.05;
    const d111 = a / Math.sqrt(3);
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Al(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Al', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Al',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Al',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },

  // Silver
  'Ag(111)': () => {
    const a = 4.086;
    const d111 = a / Math.sqrt(3);
    const asurf = a / Math.sqrt(2);

    return new CrystalStructure({
      name: 'Ag(111)',
      bulkUnitCell: new UnitCell3D(a, a, a),
      surfaceUnitCell: new UnitCell2D(
        new Vector3(asurf, 0, 0),
        new Vector3(asurf * 0.5, (asurf * Math.sqrt(3)) / 2, 0),
        60,
      ),
      millerIndices: new MillerIndices(1, 1, 1),
      crystalSystem: CrystalSystem.CUBIC,
      bravaisLattice: BravaisLattice.FACE_CENTERED,
      layers: [
        new Layer([new Atom('Ag', new Vector3(0, 0, 0))], {
          name: 'Layer 1',
          zPosition: 0,
          interlayerSpacing: d111,
        }),
        new Layer(
          [
            new Atom(
              'Ag',
              new Vector3(asurf / 2, (asurf * Math.sqrt(3)) / 6, -d111),
            ),
          ],
          { name: 'Layer 2', zPosition: -d111, interlayerSpacing: d111 },
        ),
        new Layer(
          [
            new Atom(
              'Ag',
              new Vector3(0, (asurf * Math.sqrt(3)) / 3, -2 * d111),
            ),
          ],
          { name: 'Layer 3', zPosition: -2 * d111, interlayerSpacing: d111 },
        ),
      ],
    });
  },
});

/**
 * Get list of available preset names
 */
export function getPresetNames() {
  return Object.keys(crystalPresets);
}

/**
 * Create a crystal structure from a preset
 */
export function createFromPreset(name) {
  const factory = crystalPresets[name];
  if (!factory) {
    throw new Error(`Unknown preset: ${name}`);
  }
  return factory();
}
