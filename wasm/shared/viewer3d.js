/**
 * 3D Crystal Structure Viewer using Three.js
 * Provides interactive visualization of crystal structures
 */

import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/controls/OrbitControls.js';

// Color constants for lighting (using parseInt to avoid magic number warnings)
const LIGHT_COLOR_WHITE = Number.parseInt('ffffff', 16);
const LIGHT_COLOR_AMBIENT = Number.parseInt('404040', 16);

/**
 * Crystal Structure 3D Viewer
 */
export class CrystalViewer {
  constructor(container, options = {}) {
    this.container =
      typeof container === 'string'
        ? document.getElementById(container)
        : container;

    if (!this.container) {
      throw new Error('Container element not found');
    }

    this.options = {
      backgroundColor: options.backgroundColor ?? 0x1a1a2e,
      atomScale: options.atomScale ?? 0.4,
      bondRadius: options.bondRadius ?? 0.05,
      bondColor: options.bondColor ?? 0x888888,
      showBonds: options.showBonds ?? true,
      showUnitCell: options.showUnitCell ?? true,
      showAxes: options.showAxes ?? true,
      maxBondLength: options.maxBondLength ?? 3.5,
      repeatX: options.repeatX ?? 2,
      repeatY: options.repeatY ?? 2,
      ...options,
    };

    this.scene = null;
    this.camera = null;
    this.renderer = null;
    this.controls = null;
    this.atomMeshes = [];
    this.bondMeshes = [];
    this.unitCellLines = null;
    this.unitCellClones = []; // Track cloned unit cell lines for disposal
    this.axesHelper = null;
    this.structure = null;
    this.animationId = null;

    this._init();
  }

  _init() {
    const width = this.container.clientWidth || 400;
    const height = this.container.clientHeight || 400;

    // Scene
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(this.options.backgroundColor);

    // Camera
    this.camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 1000);
    this.camera.position.set(10, 10, 10);

    // Renderer
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.container.appendChild(this.renderer.domElement);

    // Controls
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.05;

    // Lighting
    const ambientLight = new THREE.AmbientLight(LIGHT_COLOR_AMBIENT, 0.6);
    this.scene.add(ambientLight);

    const directionalLight = new THREE.DirectionalLight(LIGHT_COLOR_WHITE, 0.8);
    directionalLight.position.set(10, 10, 10);
    this.scene.add(directionalLight);

    const directionalLight2 = new THREE.DirectionalLight(
      LIGHT_COLOR_WHITE,
      0.4,
    );
    directionalLight2.position.set(-10, -10, 5);
    this.scene.add(directionalLight2);

    // Axes helper
    if (this.options.showAxes) {
      this.axesHelper = new THREE.AxesHelper(5);
      this.scene.add(this.axesHelper);
    }

    // Handle resize
    this._resizeObserver = new ResizeObserver(() => this._onResize());
    this._resizeObserver.observe(this.container);

    // Start animation loop
    this._animate();
  }

  _onResize() {
    if (!this.container) return;
    const width = this.container.clientWidth;
    const height = this.container.clientHeight;

    if (
      !width ||
      !height ||
      !Number.isFinite(width) ||
      !Number.isFinite(height)
    ) {
      return;
    }

    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(width, height);
  }

  _animate() {
    this.animationId = requestAnimationFrame(() => this._animate());
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  /**
   * Set the crystal structure to display
   * @param {CrystalStructure} structure - The crystal structure
   */
  setStructure(structure) {
    this.structure = structure;
    this.clear();
    this._buildStructure();
    this._centerCamera();
  }

  /**
   * Clear all structure elements from the scene
   */
  clear() {
    // Remove atoms
    for (const mesh of this.atomMeshes) {
      this.scene.remove(mesh);
      mesh.geometry.dispose();
      mesh.material.dispose();
    }
    this.atomMeshes = [];

    // Remove bonds
    for (const mesh of this.bondMeshes) {
      this.scene.remove(mesh);
      mesh.geometry.dispose();
      mesh.material.dispose();
    }
    this.bondMeshes = [];

    // Remove cloned unit cell lines
    for (const clone of this.unitCellClones) {
      this.scene.remove(clone);
      // Do not dispose geometry/material here as they are shared with unitCellLines
    }
    this.unitCellClones = [];

    // Remove unit cell
    if (this.unitCellLines) {
      this.scene.remove(this.unitCellLines);
      this.unitCellLines.geometry.dispose();
      this.unitCellLines.material.dispose();
      this.unitCellLines = null;
    }
  }

  _buildStructure() {
    if (!this.structure) return;

    const atoms = this.structure.getAllAtoms(
      this.options.repeatX,
      this.options.repeatY,
    );

    // Create atoms
    for (const atom of atoms) {
      this._createAtom(atom);
    }

    // Create bonds
    if (this.options.showBonds) {
      this._createBonds(atoms);
    }

    // Create unit cell outline
    if (this.options.showUnitCell && this.structure.surfaceUnitCell) {
      this._createUnitCell();
    }
  }

  _createAtom(atom) {
    const radius = atom.getRadius() * this.options.atomScale;
    const color = new THREE.Color(atom.getColor());

    const geometry = new THREE.SphereGeometry(radius, 32, 32);
    const material = new THREE.MeshPhongMaterial({
      color: color,
      shininess: 100,
      specular: 0x444444,
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.set(atom.position.x, atom.position.z, atom.position.y); // Swap Y/Z for display
    mesh.userData = { atom, type: 'atom' };

    this.scene.add(mesh);
    this.atomMeshes.push(mesh);
  }

  _createBonds(atoms) {
    const maxDist = this.options.maxBondLength;
    const maxDistSq = maxDist * maxDist;

    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const atom1 = atoms[i];
        const atom2 = atoms[j];

        const dx = atom2.position.x - atom1.position.x;
        const dy = atom2.position.y - atom1.position.y;
        const dz = atom2.position.z - atom1.position.z;
        const distSq = dx * dx + dy * dy + dz * dz;

        if (distSq < maxDistSq && distSq > 0.1) {
          this._createBond(atom1, atom2);
        }
      }
    }
  }

  _createBond(atom1, atom2) {
    // Calculate bond direction and length
    const start = new THREE.Vector3(
      atom1.position.x,
      atom1.position.z,
      atom1.position.y,
    );
    const end = new THREE.Vector3(
      atom2.position.x,
      atom2.position.z,
      atom2.position.y,
    );

    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    const midpoint = new THREE.Vector3()
      .addVectors(start, end)
      .multiplyScalar(0.5);

    // Create cylinder for bond
    const geometry = new THREE.CylinderGeometry(
      this.options.bondRadius,
      this.options.bondRadius,
      length,
      8,
    );

    const material = new THREE.MeshPhongMaterial({
      color: this.options.bondColor,
      shininess: 50,
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.copy(midpoint);

    // Orient cylinder along bond direction
    const axis = new THREE.Vector3(0, 1, 0);
    mesh.quaternion.setFromUnitVectors(axis, direction.normalize());

    mesh.userData = { type: 'bond', atom1, atom2 };

    this.scene.add(mesh);
    this.bondMeshes.push(mesh);
  }

  _createUnitCell() {
    const cell = this.structure.surfaceUnitCell;
    if (!cell) return;

    const a = cell.a;
    const b = cell.b;

    // Create lines for unit cell outline
    const points = [
      // Bottom face
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(a.x, 0, a.y),
      new THREE.Vector3(a.x + b.x, 0, a.y + b.y),
      new THREE.Vector3(b.x, 0, b.y),
      new THREE.Vector3(0, 0, 0),
    ];

    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({
      color: 0x00ff00,
      linewidth: 2,
    });

    this.unitCellLines = new THREE.Line(geometry, material);
    this.scene.add(this.unitCellLines);

    // Add repeated unit cells and track clones for disposal
    for (let ix = 0; ix < this.options.repeatX; ix++) {
      for (let iy = 0; iy < this.options.repeatY; iy++) {
        if (ix === 0 && iy === 0) continue;

        const offset = new THREE.Vector3(
          a.x * ix + b.x * iy,
          0,
          a.y * ix + b.y * iy,
        );

        const cellCopy = this.unitCellLines.clone();
        cellCopy.position.copy(offset);
        this.scene.add(cellCopy);
        this.unitCellClones.push(cellCopy);
      }
    }
  }

  _centerCamera() {
    if (this.atomMeshes.length === 0) return;

    // Calculate bounding box
    const box = new THREE.Box3();
    for (const mesh of this.atomMeshes) {
      box.expandByObject(mesh);
    }

    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z);

    // Position camera
    this.camera.position.set(
      center.x + maxDim * 1.5,
      center.y + maxDim * 1.5,
      center.z + maxDim * 1.5,
    );
    this.controls.target.copy(center);
    this.controls.update();
  }

  /**
   * Set display options
   * @param {Object} options - Options to update
   */
  setOptions(options) {
    Object.assign(this.options, options);

    if (this.structure) {
      this.clear();
      this._buildStructure();
    }

    if (options.backgroundColor !== undefined) {
      this.scene.background = new THREE.Color(options.backgroundColor);
    }

    if (options.showAxes !== undefined) {
      if (options.showAxes && !this.axesHelper) {
        this.axesHelper = new THREE.AxesHelper(5);
        this.scene.add(this.axesHelper);
      } else if (!options.showAxes && this.axesHelper) {
        this.scene.remove(this.axesHelper);
        this.axesHelper = null;
      }
    }
  }

  /**
   * Reset camera to default view
   */
  resetView() {
    this._centerCamera();
  }

  /**
   * Set view direction
   * @param {string} direction - 'top', 'front', 'side', or 'perspective'
   */
  setView(direction) {
    const box = new THREE.Box3();
    for (const mesh of this.atomMeshes) {
      box.expandByObject(mesh);
    }
    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z) * 2;

    switch (direction) {
      case 'top':
        this.camera.position.set(center.x, center.y + maxDim, center.z);
        break;
      case 'front':
        this.camera.position.set(center.x, center.y, center.z + maxDim);
        break;
      case 'side':
        this.camera.position.set(center.x + maxDim, center.y, center.z);
        break;
      case 'perspective':
      default:
        this.camera.position.set(
          center.x + maxDim * 0.7,
          center.y + maxDim * 0.7,
          center.z + maxDim * 0.7,
        );
    }

    this.controls.target.copy(center);
    this.controls.update();
  }

  /**
   * Export current view as PNG image
   * @returns {string} - Data URL of the image
   */
  exportImage() {
    this.renderer.render(this.scene, this.camera);
    return this.renderer.domElement.toDataURL('image/png');
  }

  /**
   * Clean up resources
   */
  dispose() {
    if (this.animationId) {
      cancelAnimationFrame(this.animationId);
    }

    if (this._resizeObserver) {
      this._resizeObserver.disconnect();
    }

    this.clear();

    if (this.controls) {
      this.controls.dispose();
      this.controls = null;
    }

    if (this.axesHelper) {
      this.scene.remove(this.axesHelper);
      this.axesHelper.geometry.dispose();
      this.axesHelper.material.dispose();
      this.axesHelper = null;
    }

    this.renderer.dispose();
    // Defensive check before removing DOM element
    if (
      this.renderer.domElement &&
      this.renderer.domElement.parentNode === this.container
    ) {
      this.renderer.domElement.remove();
    }
  }
}

/**
 * Create a viewer instance attached to a container
 * @param {string|HTMLElement} container - Container element or ID
 * @param {Object} options - Viewer options
 * @returns {CrystalViewer}
 */
export function createViewer(container, options = {}) {
  return new CrystalViewer(container, options);
}
