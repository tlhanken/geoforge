// The geoforge viewer's page: it asks the server for the scene at an address and draws it.
// Everything shown is worked out by the server (geoforge-viewer); this only draws and picks.
import * as THREE from '/three.js';

const canvas = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.outputColorSpace = THREE.SRGBColorSpace;
const camera = new THREE.PerspectiveCamera(55, 1, 0.001, 1e6);
const controls = new THREE.OrbitControls(camera, canvas);
controls.enableDamping = true;
let scene = new THREE.Scene();
let pickables = [];   // { object, index?, go, label }
let labels = [];      // { position: Vector3, text, dim }
let spinning = null;  // a mesh that turns slowly
const state = { address: null, stars: false, scene: null, layer: 'appearance', flat: false, start: null };

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const rgb = (c) => new THREE.Color(c[0] / 255, c[1] / 255, c[2] / 255);

async function api(path, method = 'GET') {
  const r = await fetch(path, { method });
  const j = await r.json();
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}

function resize() {
  const w = window.innerWidth, h = window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);

// ---------------------------------------------------------------- navigation

async function go(address, opts = {}) {
  const stars = !!opts.stars;
  $('busy').style.display = 'block';
  try {
    const s = await api(`/api/scene?at=${encodeURIComponent(address)}${stars ? '&stars=1' : ''}`);
    state.address = address;
    state.stars = stars;
    state.scene = s;
    if (!opts.keepLayer) { state.layer = 'appearance'; state.flat = false; }
    const hash = `#at=${encodeURIComponent(address)}${stars ? '&stars=1' : ''}`;
    if (location.hash !== hash) history.pushState(null, '', hash);
    draw(s, !opts.keepView);
    crumbs();
  } catch (e) {
    tip(`${e.message}`);
  } finally {
    $('busy').style.display = 'none';
  }
}

window.addEventListener('popstate', () => fromHash());
function fromHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  const at = p.get('at');
  if (at) go(at, { stars: p.has('stars') });
  else go(state.start.home);
}

function crumbs() {
  const out = (state.scene.crumbs || []).map((c) => ({ a: c.address, label: c.name, stars: c.stars }));
  $('crumbs').innerHTML = out.map((c, k) => {
    const here = c.a === state.address && c.stars === state.stars;
    return `${k ? '<span>›</span>' : ''}<a href="#" data-a="${esc(c.a)}" data-s="${c.stars ? 1 : ''}" class="${here ? 'here' : ''}" title="${esc(c.a)}">${esc(c.label)}</a>`;
  }).join('');
  for (const el of $('crumbs').querySelectorAll('a')) el.onclick = (ev) => { ev.preventDefault(); go(el.dataset.a, { stars: !!el.dataset.s }); };
}

function tip(t) { $('tip').textContent = t; }

// ---------------------------------------------------------------- drawing

function fresh(far) {
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05070d);
  pickables = [];
  labels = [];
  spinning = null;
  camera.near = far / 1e5;
  camera.far = far * 20;
  camera.updateProjectionMatrix();
}

function aim(radius, from = [0.55, 0.5, 0.8]) {
  const d = new THREE.Vector3(...from).normalize().multiplyScalar(radius * 2.2);
  camera.position.copy(d);
  controls.target.set(0, 0, 0);
  controls.update();
}

function dotTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const g = c.getContext('2d');
  const grad = g.createRadialGradient(32, 32, 0, 32, 32, 32);
  grad.addColorStop(0, 'rgba(255,255,255,1)');
  grad.addColorStop(0.25, 'rgba(255,255,255,0.8)');
  grad.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grad;
  g.fillRect(0, 0, 64, 64);
  return new THREE.CanvasTexture(c);
}
const DOT = dotTexture();

// A thin ring, to circle the stars whose systems hold a world with open water.
function ringTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const g = c.getContext('2d');
  g.strokeStyle = 'rgba(255,255,255,1)';
  g.lineWidth = 4;
  g.beginPath();
  g.arc(32, 32, 26, 0, 2 * Math.PI);
  g.stroke();
  return new THREE.CanvasTexture(c);
}
const RING = ringTexture();

// A diamond outline, to mark the stars whose systems hold a fictional material.
function diamondTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const g = c.getContext('2d');
  g.strokeStyle = 'rgba(255,255,255,1)';
  g.lineWidth = 4;
  g.beginPath();
  g.moveTo(32, 4); g.lineTo(60, 32); g.lineTo(32, 60); g.lineTo(4, 32); g.closePath();
  g.stroke();
  return new THREE.CanvasTexture(c);
}
const DIAMOND = diamondTexture();
const HABITABLE = { temperate: [110, 230, 140], habitable: [90, 170, 220] };

function points(positions, colours, size, opts = {}) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  geo.setAttribute('color', new THREE.Float32BufferAttribute(colours, 3));
  if (opts.sizes) geo.setAttribute('size', new THREE.Float32BufferAttribute(opts.sizes, 1));
  const mat = new THREE.PointsMaterial({ size, vertexColors: true, map: DOT, transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, sizeAttenuation: opts.attenuate ?? true, opacity: opts.opacity ?? 1 });
  const p = new THREE.Points(geo, mat);
  scene.add(p);
  return p;
}

function marker(position, colour, size, text, goTo, opts = {}) {
  const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: DOT, color: rgb(colour), depthWrite: false, blending: THREE.AdditiveBlending }));
  m.position.copy(position);
  m.scale.setScalar(size);
  scene.add(m);
  if (goTo) pickables.push({ object: m, go: goTo, label: text });
  if (text) labels.push({ position: position.clone(), text, dim: opts.dim });
  return m;
}

function draw(s, reframe) {
  ({ universe: drawUniverse, galaxy: drawGalaxy, stars: drawStars, system: drawSystem, body: drawBody })[s.level](s, reframe);
  panel(s);
}

function drawUniverse(s, reframe) {
  fresh(s.radius_mpc);
  const pos = [], col = [];
  for (const g of s.galaxies) {
    pos.push(...g.position_mpc);
    const c = rgb(g.colour);
    col.push(c.r, c.g, c.b);
  }
  const p = points(pos, col, s.radius_mpc * 0.012, { opacity: 0.9 });
  s.galaxies.forEach((g, i) => {
    const big = g.stellar_msun > 2e10;
    if (g.home || big) marker(new THREE.Vector3(...g.position_mpc), g.colour, s.radius_mpc * (g.home ? 0.05 : 0.03), g.home ? `${g.name} (home)` : g.name, g.address ? () => go(g.address) : null, { dim: !g.home });
  });
  pickables.push({ object: p, pick: (i) => s.galaxies[i].address && go(s.galaxies[i].address), label: (i) => `${s.galaxies[i].name} · ${s.galaxies[i].kind} · ${s.galaxies[i].stellar_msun.toExponential(1)} M☉${s.galaxies[i].active ? ' · ' + s.galaxies[i].active : ''}` });
  // Active nuclei ringed: violet for Seyferts and quasars, dim for LINERs.
  for (const g of s.galaxies.filter((x) => x.active && !/LINER/.test(x.active))) {
    const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: RING, color: rgb([200, 140, 255]), depthWrite: false, transparent: true, sizeAttenuation: false }));
    m.position.set(...g.position_mpc);
    m.scale.setScalar(/quasar/.test(g.active) ? 0.04 : 0.025);
    scene.add(m);
  }
  if (reframe) aim(s.radius_mpc * 0.6);
}

function drawGalaxy(s, reframe) {
  fresh(s.radius_kpc);
  const cloud = (flat, bright, size) => {
    const pos = [], col = [];
    for (let i = 0; i < flat.length; i += 6) {
      pos.push(flat[i], flat[i + 2], -flat[i + 1]);
      col.push(flat[i + 3] * bright, flat[i + 4] * bright, flat[i + 5] * bright);
    }
    points(pos, col, s.radius_kpc * size, { opacity: 0.85 });
  };
  cloud(s.cloud, 0.45, 0.006);
  cloud(s.young, 0.9, 0.009);
  // Systems made close together share one mark: within a tenth of a kiloparsec.
  const groups = [];
  for (const sys of s.systems) {
    const p = new THREE.Vector3(sys.position_kpc[0], sys.position_kpc[2], -sys.position_kpc[1]);
    const g = groups.find((x) => x.p.distanceTo(p) < 0.1);
    if (g) g.names.push(sys.name); else groups.push({ p, names: [sys.name], address: sys.address });
  }
  for (const g of groups) marker(g.p, [255, 230, 120], s.radius_kpc * 0.04, g.names.join(', '), () => go(g.address, { stars: true }));
  // Globular clusters: warm points, each telling about itself.
  for (const c of s.globulars) {
    const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: DOT, color: rgb([255, 215, 150]), depthWrite: false, blending: THREE.AdditiveBlending }));
    m.position.set(c.position_kpc[0], c.position_kpc[2], -c.position_kpc[1]);
    m.scale.setScalar(s.radius_kpc * 0.012);
    scene.add(m);
    pickables.push({ object: m, go: () => showNebula(s, c), label: `${c.name} · globular cluster` });
  }
  // The centre's black hole.
  const hole = new THREE.Sprite(new THREE.SpriteMaterial({ map: RING, color: rgb([255, 170, 90]), depthWrite: false, transparent: true, sizeAttenuation: false }));
  hole.scale.setScalar(0.03);
  scene.add(hole);
  pickables.push({ object: hole, go: () => {}, label: `Centre: ${s.centre}` });
  if (reframe) aim(s.radius_kpc * 0.6, [0.15, 1.0, 0.35]);
  markNebulae(s);
}

// ---------------------------------------------------------------- nebulae

const NEBULA_COLOUR = { emission: [255, 80, 120], reflection: [100, 150, 255], dark: [120, 95, 70], planetary: [80, 230, 200], supernova_remnant: [255, 150, 170] };
const NEBULA_KIND = { emission: 'emission nebula', reflection: 'reflection nebula', dark: 'dark cloud', planetary: 'planetary nebula', supernova_remnant: 'supernova remnant', globular: 'globular cluster' };

// Fetch the nebulae round the scene (a galaxy's home star, or a system) and draw them: in a
// galaxy as glows by kind, round a system as thin spheres where they reach into the view.
async function markNebulae(s) {
  let found;
  try {
    found = await api(`/api/nebulae?at=${encodeURIComponent(s.address)}`);
  } catch (e) {
    return;
  }
  if (state.scene !== s) return;
  s.nebulae = found;
  for (const n of found) {
    const c = NEBULA_COLOUR[n.kind];
    const label = `${n.name} · ${NEBULA_KIND[n.kind]}${n.materials.length ? ' · ' + n.materials.join(', ') + ' (fictional)' : ''}`;
    if (s.level === 'galaxy') {
      const p = new THREE.Vector3(n.position_kpc[0], n.position_kpc[2], -n.position_kpc[1]);
      const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: DOT, color: rgb(c), depthWrite: false, transparent: true, opacity: n.kind === 'dark' ? 0.5 : 0.9, blending: n.kind === 'dark' ? THREE.NormalBlending : THREE.AdditiveBlending }));
      m.position.copy(p);
      m.scale.setScalar(s.radius_kpc * 0.008 + 4 * n.radius_pc / 1000);
      scene.add(m);
      pickables.push({ object: m, go: () => showNebula(s, n), label });
    } else if (n.distance_pc - n.radius_pc <= s.radius_pc) {
      const p = new THREE.Vector3(n.offset_pc[0], n.offset_pc[2], -n.offset_pc[1]);
      const m = new THREE.Mesh(new THREE.SphereGeometry(n.radius_pc, 32, 16), new THREE.MeshBasicMaterial({ color: rgb(c), transparent: true, opacity: 0.07, depthWrite: false, side: THREE.DoubleSide }));
      m.position.copy(p);
      scene.add(m);
      labels.push({ position: p, text: n.name, dim: true });
    }
  }
  panel(s);
}

// Tell about nebula `n` in the side panel, until the scene changes or Back is pressed.
function showNebula(s, n) {
  s.picked = n;
  panel(s);
}

// The side panel's rows for the nebulae listed round a stars scene: the nearest few, each a
// link that tells about it.
function nebulaList(s) {
  if (!s.nebulae || !s.nebulae.length) return '';
  const near = s.nebulae.filter((n) => n.kind !== 'dark' || n.distance_pc < 200).slice(0, 6);
  return `<h3>Nebulae within ${Math.round(500 * 3.2616)} ly</h3><ul class="finds">${near.map((n, i) => `<li><a href="#" data-neb="${i}">${esc(n.name)}</a>: ${esc(NEBULA_KIND[n.kind])}, ${Math.round(n.distance_pc * 3.2616)} ly</li>`).join('')}</ul>`;
}

const KIND_NAME = { star: 'star', white_dwarf: 'white dwarf', brown_dwarf: 'brown dwarf', neutron_star: 'neutron star', black_hole: 'black hole' };

function drawStars(s, reframe) {
  fresh(s.radius_pc);
  const pos = [], col = [];
  for (const st of s.stars) {
    pos.push(st.position_pc[0], st.position_pc[2], -st.position_pc[1]);
    const f = st.absolute_v === null ? 1 : Math.max(0.25, Math.min(1, 1.2 - (st.absolute_v - 2) / 14));
    const c = rgb(st.colour);
    col.push(c.r * f, c.g * f, c.b * f);
  }
  const label = (st) => `${st.name} · ${st.kind === 'star' ? st.spectral : KIND_NAME[st.kind]} · ${(Math.hypot(...st.position_pc) * 3.2616).toFixed(1)} ly`;
  const layer = (keep, size) => {
    const idx = s.stars.map((st, i) => i).filter((i) => keep(s.stars[i]));
    if (!idx.length) return;
    const p = points(idx.flatMap((i) => pos.slice(3 * i, 3 * i + 3)), idx.flatMap((i) => col.slice(3 * i, 3 * i + 3)), s.radius_pc * size);
    pickables.push({ object: p, pick: (k) => go(s.stars[idx[k]].address), label: (k) => label(s.stars[idx[k]]) });
  };
  // Shining stars in three sizes by brightness: the bright few, the Sun-like, the dim many.
  const lit = (st) => st.kind === 'star' || st.kind === 'white_dwarf';
  layer((st) => lit(st) && st.absolute_v < 2, 0.05);
  layer((st) => lit(st) && st.absolute_v >= 2 && st.absolute_v < 8, 0.025);
  layer((st) => lit(st) && st.absolute_v >= 8, 0.012);
  // Brown dwarfs and neutron stars small, in their own colours; black holes as a violet ring.
  layer((st) => st.kind === 'brown_dwarf', 0.01);
  layer((st) => st.kind === 'neutron_star', 0.014);
  for (const st of s.stars.filter((x) => x.kind === 'black_hole')) {
    const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: RING, color: rgb(st.colour), depthWrite: false, transparent: true, sizeAttenuation: false }));
    m.position.set(st.position_pc[0], st.position_pc[2], -st.position_pc[1]);
    m.scale.setScalar(0.018);
    scene.add(m);
    pickables.push({ object: m, go: () => go(st.address), label: label(st) });
  }
  for (const st of s.stars) {
    if (st.centre || st.made || (st.absolute_v !== null && st.absolute_v < 3)) marker(new THREE.Vector3(st.position_pc[0], st.position_pc[2], -st.position_pc[1]), st.centre ? [255, 240, 160] : st.colour, s.radius_pc * (st.centre ? 0.08 : 0.04), st.name, () => go(st.address), { dim: !st.centre && !st.made });
  }
  // A ring at the centre's plane, for scale: 10 light-years.
  const ring = new THREE.LineLoop(circle(10 / 3.2616, 96), new THREE.LineBasicMaterial({ color: 0x26304a }));
  scene.add(ring);
  if (reframe) aim(s.radius_pc * 0.35);
  markHabitable(s);
  markExotics(s);
  markNebulae(s);
}

// Diamond the systems holding a fictional material (in its colour), once the server has
// worked them out; nothing if the universe lists none.
async function markExotics(s) {
  let found;
  try {
    found = await api(`/api/exotics?at=${encodeURIComponent(s.address)}`);
  } catch (e) {
    return;
  }
  if (state.scene !== s) return;
  s.exotics = found;
  const colour = new Map((found.materials || []).map((m) => [m.name, m.colour]));
  for (const h of found.systems || []) {
    const st = s.stars.find((x) => x.address === h.address);
    if (!st) continue;
    const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: DIAMOND, color: rgb(colour.get(h.materials[0]) || [190, 110, 255]), depthWrite: false, transparent: true, sizeAttenuation: false }));
    m.position.set(st.position_pc[0], st.position_pc[2], -st.position_pc[1]);
    m.scale.setScalar(0.036);
    scene.add(m);
    pickables.push({ object: m, go: () => go(st.address), label: `${st.name} · ${h.materials.join(', ')} (fictional)` });
  }
  panel(s);
}

// Circle the systems with a world of open water (green: temperate, like Earth; blue: water
// open somewhere), once the server has worked them out.
async function markHabitable(s) {
  let found;
  try {
    found = await api(`/api/habitable?at=${encodeURIComponent(s.address)}`);
  } catch (e) {
    return;
  }
  if (state.scene !== s) return;
  s.habitable = found;
  for (const h of found) {
    const st = s.stars.find((x) => x.address === h.address);
    if (!st) continue;
    const p = new THREE.Vector3(st.position_pc[0], st.position_pc[2], -st.position_pc[1]);
    // A fixed size on screen, near or far.
    const m = new THREE.Sprite(new THREE.SpriteMaterial({ map: RING, color: rgb(HABITABLE[h.habitability]), depthWrite: false, transparent: true, sizeAttenuation: false }));
    m.position.copy(p);
    m.scale.setScalar(h.habitability === 'temperate' ? 0.045 : 0.028);
    scene.add(m);
    pickables.push({ object: m, go: () => go(st.address), label: `${st.name} · ${h.habitability} world` });
    if (h.habitability === 'temperate' && !st.centre) labels.push({ position: p, text: `${st.name} ✓` });
  }
  panel(s);
}

function circle(r, n) {
  const g = new THREE.BufferGeometry();
  const pts = [];
  for (let k = 0; k < n; k++) { const a = (k / n) * Math.PI * 2; pts.push(r * Math.cos(a), 0, r * Math.sin(a)); }
  g.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
  return g;
}

// Distances in a system are squeezed so inner and outer worlds both show: a planet at a AU
// is drawn at 12·√a; moons on rings just outside their planet.
const squeeze = (v) => { const d = Math.hypot(...v); return d > 0 ? v.map((x) => (x / d) * 12 * Math.sqrt(d)) : [0, 0, 0]; };
const toScene = (v) => new THREE.Vector3(v[0], v[2], -v[1]);
let systemScale = 1;
const bodySize = (km) => systemScale * Math.max(0.25, 0.6 * Math.cbrt(km / 6371));

function texture(address, layer, width) {
  const t = new THREE.TextureLoader().load(`/api/texture?at=${encodeURIComponent(address)}&layer=${encodeURIComponent(layer)}&width=${width}`);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

function drawSystem(s, reframe) {
  // Bodies are drawn far larger than they are, in proportion to the system's size.
  let extent = 4;
  for (const b of s.bodies) if (!s.bodies.some((x) => x.key === b.parent)) extent = Math.max(extent, Math.hypot(...squeeze(b.position_au)));
  systemScale = Math.max(1, extent / 40);
  fresh(60);
  scene.add(new THREE.AmbientLight(0xffffff, 0.08));
  const where = new Map();
  for (const st of s.stars) {
    const p = toScene(squeeze(st.position_au));
    where.set(st.key, p);
    const r = systemScale * Math.max(0.6, 0.6 * Math.cbrt(st.radius_solar));
    // A black hole is a dark disc in a ring of bent light; it lights nothing.
    const dark = st.kind === 'black_hole';
    const m = new THREE.Mesh(new THREE.SphereGeometry(r, 32, 16), new THREE.MeshBasicMaterial({ color: dark ? 0x000000 : rgb(st.colour) }));
    m.position.copy(p);
    scene.add(m);
    marker(p, st.colour, r * 7, st.name, null);
    if (dark) {
      const halo = new THREE.Sprite(new THREE.SpriteMaterial({ map: RING, color: rgb(st.colour), depthWrite: false, transparent: true }));
      halo.position.copy(p);
      halo.scale.setScalar(r * 3.2);
      scene.add(halo);
    } else {
      const sun = new THREE.PointLight(rgb(st.colour), st.kind === 'neutron_star' || st.kind === 'brown_dwarf' ? 0.6 : 2.5, 0, 0);
      sun.position.copy(p);
      scene.add(sun);
    }
    // A pulsar's two beams sweep out from its magnetic poles.
    if (st.kind === 'neutron_star' && st.notes.some((n) => /pulsar|magnetar/.test(n))) {
      for (const sgn of [1, -1]) {
        const beam = new THREE.Mesh(new THREE.ConeGeometry(r * 1.2, r * 14, 16, 1, true), new THREE.MeshBasicMaterial({ color: 0x9fc4ff, transparent: true, opacity: 0.25, side: THREE.DoubleSide, depthWrite: false }));
        beam.position.copy(p).add(new THREE.Vector3(0, sgn * r * 7, 0).applyAxisAngle(new THREE.Vector3(0, 0, 1), 0.4));
        beam.rotation.z = 0.4 + (sgn > 0 ? Math.PI : 0);
        scene.add(beam);
      }
    }
  }
  // Bodies around stars first (their places), then moons around them.
  const around = (b) => s.bodies.find((x) => x.key === b.parent);
  const moonIndex = new Map();
  for (const b of s.bodies) {
    const parent = around(b);
    let ringPts;
    if (!parent) {
      const pk = where.get(b.parent) || new THREE.Vector3();
      ringPts = b.orbit_au.map((v) => toScene(squeeze(v)).add(pk));
      const p = toScene(squeeze(b.position_au.map((v, k) => v - (s.stars.find((x) => x.key === b.parent)?.position_au[k] ?? 0)))).add(pk);
      where.set(b.key, p);
      extent = Math.max(extent, p.length());
    }
    if (ringPts && ringPts.length) {
      const g = new THREE.BufferGeometry().setFromPoints(ringPts);
      scene.add(new THREE.LineLoop(g, new THREE.LineBasicMaterial({ color: b.focus ? 0x5a7ad0 : 0x26304a })));
    }
  }
  for (const b of s.bodies) {
    const parent = around(b);
    if (!parent) continue;
    const pp = where.get(parent.key) || new THREE.Vector3();
    const k = (moonIndex.get(parent.key) || 0) + 1;
    moonIndex.set(parent.key, k);
    const r = bodySize(parent.radius_km) * (1.8 + 0.7 * k);
    const a = Math.atan2(b.position_au[1] - parent.position_au[1], b.position_au[0] - parent.position_au[0]);
    const p = pp.clone().add(new THREE.Vector3(r * Math.cos(a), 0, -r * Math.sin(a)));
    where.set(b.key, p);
    const ring = new THREE.LineLoop(circle(r, 64), new THREE.LineBasicMaterial({ color: 0x1c2438 }));
    ring.position.copy(pp);
    scene.add(ring);
  }
  for (const b of s.bodies) {
    const p = where.get(b.key);
    const size = bodySize(b.radius_km) * (around(b) ? 0.7 : 1);
    const m = new THREE.Mesh(new THREE.SphereGeometry(size, 32, 16), new THREE.MeshStandardMaterial({ map: texture(b.address, 'appearance', 256), roughness: 1, metalness: 0 }));
    m.position.copy(p);
    scene.add(m);
    pickables.push({ object: m, go: () => go(b.address), label: `${b.name} · ${b.class}${b.made ? ' · made' : ''}` });
    const r = (s.rings || []).find((x) => x.parent === b.key);
    if (r) addRing(r, size, p);
    if (!around(b) || b.made || b.focus || b.habitable) labels.push({ position: p.clone().add(new THREE.Vector3(0, size, 0)), text: b.name + (b.focus ? ' ★' : '') + (b.habitable === 'temperate' ? ' ✓' : b.habitable ? ' ~' : ''), dim: !b.made && !b.focus && !b.habitable });
  }
  for (const belt of s.belts) {
    const pk = where.get(belt.parent) || new THREE.Vector3();
    const inner = 12 * Math.sqrt(belt.inner_au), outer = 12 * Math.sqrt(belt.outer_au);
    if (where.has(belt.parent) && s.bodies.some((b) => b.key === belt.parent)) continue;
    const m = new THREE.Mesh(new THREE.RingGeometry(inner, Math.max(outer, inner + 0.05), 128), new THREE.MeshBasicMaterial({ color: 0x8a7f6a, transparent: true, opacity: 0.06, side: THREE.DoubleSide, depthWrite: false }));
    m.rotation.x = -Math.PI / 2;
    m.position.copy(pk);
    scene.add(m);
  }
  if (reframe) aim(Math.min(extent, 60) * 0.75, [0.2, 0.7, 0.9]);
}

// A planet's rings around a body drawn `size` across at `at`: in its equator, tilted with it;
// bright ice as opaque as its optical depth makes it, in bands with a gap two thirds out (as
// moons' resonances clear Saturn's Cassini division); faint dust as a whisper.
function addRing(r, size, at) {
  const group = new THREE.Group();
  const bands = r.icy ? 28 : 1;
  const base = r.icy ? Math.min(0.8, 1 - Math.exp(-2 * r.depth)) : 0.07;
  for (let k = 0; k < bands; k++) {
    const t0 = k / bands, t1 = (k + 1) / bands;
    const wave = r.icy ? 0.55 + 0.45 * Math.sin(k * 2.7 + r.depth * 10) ** 2 : 1;
    const gap = r.icy && Math.abs((t0 + t1) / 2 - 0.66) < 0.025 ? 0.12 : 1;
    const ri = (r.inner_r + (r.outer_r - r.inner_r) * t0) * size, ro = (r.inner_r + (r.outer_r - r.inner_r) * t1) * size;
    group.add(new THREE.Mesh(new THREE.RingGeometry(ri, ro, 160), new THREE.MeshBasicMaterial({ color: r.icy ? 0xd2c3a5 : 0x96785f, transparent: true, opacity: base * wave * gap, side: THREE.DoubleSide, depthWrite: false })));
  }
  group.rotation.x = -Math.PI / 2 + (r.tilt_deg * Math.PI) / 180;
  group.position.copy(at);
  scene.add(group);
  return group;
}

function drawBody(s, reframe) {
  fresh(10);
  scene.add(new THREE.AmbientLight(0xffffff, state.layer === 'appearance' ? 0.25 : 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, state.layer === 'appearance' ? 2.4 : 0.6);
  sun.position.set(5, 2, 4);
  scene.add(sun);
  const width = Math.min(4096, Math.max(512, Math.round(window.innerWidth * renderer.getPixelRatio())));
  const map = texture(s.address, state.layer, s.made ? width : 1024);
  if (state.flat) {
    const m = new THREE.Mesh(new THREE.PlaneGeometry(4, 2), new THREE.MeshBasicMaterial({ map }));
    scene.add(m);
    for (const d of s.deposits) {
      const pin = new THREE.Sprite(new THREE.SpriteMaterial({ map: DIAMOND, color: rgb(d.colour), depthWrite: false, transparent: true, sizeAttenuation: false }));
      pin.position.set(d.lon / 90, d.lat / 90, 0.01);
      pin.scale.setScalar(0.03);
      scene.add(pin);
    }
    controls.enableRotate = false;
    if (reframe) { camera.position.set(0, 0, 3.2); controls.target.set(0, 0, 0); controls.update(); }
  } else {
    controls.enableRotate = true;
    const m = new THREE.Mesh(new THREE.SphereGeometry(1, 128, 64), new THREE.MeshStandardMaterial({ map, roughness: 1, metalness: 0 }));
    scene.add(m);
    spinning = m;
    // Fictional deposits, pinned where they are (they turn with the globe).
    for (const d of s.deposits) {
      const lat = d.lat * Math.PI / 180, lon = d.lon * Math.PI / 180;
      const pin = new THREE.Sprite(new THREE.SpriteMaterial({ map: DIAMOND, color: rgb(d.colour), depthWrite: false, transparent: true, sizeAttenuation: false }));
      pin.position.set(1.01 * Math.cos(lat) * Math.cos(lon), 1.01 * Math.sin(lat), -1.01 * Math.cos(lat) * Math.sin(lon));
      pin.scale.setScalar(0.03);
      m.add(pin);
    }
    if (s.ring) addRing(s.ring, 1, new THREE.Vector3());
    if (reframe) aim(s.ring && s.ring.icy ? 2.2 : 1.1, [0.0, 0.25, 1.0]);
  }
}

// ---------------------------------------------------------------- the side panel

// The stars panel's rows for fictional materials: how many systems hold each, or the fault
// in exotics.toml.
function exoticRows(s) {
  const x = s.exotics;
  if (!x) return [];
  if (x.error) return [['Fictional', x.error]];
  return (x.materials || []).map((m) => [m.name, x.systems.filter((h) => h.materials.includes(m.name)).length]).filter(([, n]) => n > 0).map(([name, n]) => [`${name} (fictional)`, `${n} systems ◇`]);
}

function panel(s) {
  const p = $('panel');
  const facts = (rows) => `<dl>${rows.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('')}</dl>`;
  const head = (title, sub) => `<h2>${esc(title)}</h2><div class="sub">${esc(sub)}</div>`;
  if (s.picked) {
    const n = s.picked;
    const note = n.kind === 'globular' ? 'Its light and place are drawn; its stars are not yet generated one by one.' : 'Nebulae are not places of their own: each is named by the stars or cloud that make it.';
    p.innerHTML = head(n.name, NEBULA_KIND[n.kind]) + facts([...n.facts, ['Distance', `${Math.round(n.distance_pc * 3.2616).toLocaleString()} ly from ${s.level === 'galaxy' ? 'the home star' : 'here'}`]]) +
      `<div class="row">${n.star ? '<button id="nebstar">Its star</button>' : ''}<button id="back">Back</button></div><div class="note">${note}</div>`;
    if ($('nebstar')) $('nebstar').onclick = () => go(n.star, { stars: true });
    $('back').onclick = () => { s.picked = null; panel(s); };
    return;
  }
  if (s.level === 'universe') {
    const active = (re) => s.galaxies.filter((g) => g.active && re.test(g.active)).length;
    p.innerHTML = head(s.name, s.address) + facts([['Galaxies shown', `${s.galaxies.length} within ${s.radius_mpc} Mpc`], ['Active nuclei', `${active(/quasar/)} quasars, ${active(/Seyfert/)} Seyferts, ${active(/LINER/)} LINERs`], ['Pick', 'a galaxy to go in']]) +
      `<div class="note">Ringed in <span style="color:#c88cff">violet</span>: galaxies with a Seyfert nucleus or a quasar.</div>`;
  } else if (s.level === 'galaxy') {
    p.innerHTML = head(s.name, s.address) + facts([['Type', s.kind], ['Stars', `${s.stellar_msun.toExponential(2)} M☉`], ['Disk to', `${s.radius_kpc.toFixed(1)} kpc`], ['Centre', s.centre], ['Globular clusters', String(s.globulars.length)]]) +
      (s.nebulae ? facts([['Nebulae', `${s.nebulae.length} within ${Math.round(3000 * 3.2616).toLocaleString()} ly of the home star`]]) : '') +
      `<div class="note">Gold marks: systems made here (or the default star). Pick one to see the stars around it. Round the home star, nebulae glow by kind: <span style="color:#ff5078">emission</span>, <span style="color:#6496ff">reflection</span>, <span style="color:#50e6c8">planetary</span>, <span style="color:#ff96aa">supernova remnants</span>, <span style="color:#a08060">dark clouds</span>.</div>`;
  } else if (s.level === 'stars') {
    const n = (k) => (s.habitable ? String(s.habitable.filter((h) => h.habitability === k).length) : '…');
    const count = (k) => s.stars.filter((st) => st.kind === k).length;
    p.innerHTML = head('Neighbourhood', s.address) + facts([['Systems', `${s.stars.length} within ${(s.radius_pc * 3.2616).toFixed(0)} ly`], ['Of them', `${count('brown_dwarf')} brown dwarfs, ${count('white_dwarf')} white dwarfs, ${count('neutron_star')} neutron stars, ${count('black_hole')} black holes`], ['Temperate', n('temperate')], ['Open water', n('habitable')], ...exoticRows(s), ['Ring', '10 light-years']]) +
      nebulaList(s) +
      `<div class="row"><button id="tosys">Open this system</button></div><div class="note">Pick any star to open its system: it is generated as you look. Circled: systems with a world of open water, <span style="color:#6ee68c">green</span> temperate like Earth, <span style="color:#5aaadc">blue</span> colder or hotter.</div>`;
    $('tosys').onclick = () => go(s.address);
    for (const a of p.querySelectorAll('[data-neb]')) a.onclick = (e) => { e.preventDefault(); showNebula(s, s.nebulae.filter((n) => n.kind !== 'dark' || n.distance_pc < 200)[Number(a.dataset.neb)]); };
  } else if (s.level === 'system') {
    const star = s.stars[0];
    const what = (st) => (st.kind === 'black_hole' ? 'black hole' : `${st.spectral}, ${Math.round(st.temperature_k)} K`);
    const notes = s.stars.flatMap((st) => st.notes.map((n) => [st.name, n]));
    p.innerHTML = head(s.name, s.address) + facts([['Star', star ? what(star) : '—'], ...s.stars.slice(1).map((st) => [st.name, what(st)]), ...notes, ...s.belts.filter((b) => b.mix).map((b) => [b.name, `${b.inner_au.toFixed(1)}–${b.outer_au.toFixed(1)} AU: ${b.mix}`]), ['Age', `${s.age_gyr.toFixed(2)} Gyr`], ['Bodies', String(s.bodies.length)]]) +
      (s.exotics.length ? `<h3>Fictional materials</h3><ul class="finds">${s.exotics.map((f, i) => `<li>${f.address ? `<a href="#" data-find="${i}">` : ''}${esc(f.material)}: ${esc(f.place)}${f.address ? '</a>' : ''}</li>`).join('')}</ul>` : '') +
      `<div class="row"><button id="near">Stars around</button></div><div class="note">Distances are squeezed (√ of AU); moons on rings by their planet. ★ is the system's most Earth-like body; ✓ a temperate world, ~ one with open water but colder or hotter.</div>`;
    $('near').onclick = () => go(s.address, { stars: true });
    for (const a of p.querySelectorAll('[data-find]')) a.onclick = (e) => { e.preventDefault(); go(s.exotics[Number(a.dataset.find)].address); };
  } else if (s.level === 'body') {
    const layerOpts = s.layers.map((l) => `<option ${l === state.layer ? 'selected' : ''}>${esc(l)}</option>`).join('');
    const make = s.made
      ? (s.layers.length > 1 ? '' : `<button id="gen">Generate its layers</button>`)
      : `<button id="make">Make this world</button><button id="gen">Make and generate</button>`;
    p.innerHTML = head(s.name, s.address) + facts(s.facts) +
      `<div class="row"><select id="layer">${layerOpts}</select><button id="sphere" class="${state.flat ? '' : 'on'}">Globe</button><button id="flat" class="${state.flat ? 'on' : ''}">Map</button></div>` +
      `<div class="row">${make}</div>` +
      `<div class="note">${s.job ? `Generating: ${esc(s.job)}` : s.made ? (s.layers.length > 1 ? 'Generated here: every layer can be shown.' : 'Made, not generated yet: true colour is drawn from what it is.') : 'Not made: true colour is drawn from what it is.'}</div>`;
    $('layer').onchange = (e) => { state.layer = e.target.value; draw(s, false); };
    $('sphere').onclick = () => { state.flat = false; draw(s, true); };
    $('flat').onclick = () => { state.flat = true; draw(s, true); };
    if ($('make')) $('make').onclick = async () => { await act('make', s.address); };
    if ($('gen')) $('gen').onclick = async () => { await act('generate', s.address); poll(s.address); };
  }
}

async function act(what, address) {
  try { await api(`/api/${what}?at=${encodeURIComponent(address)}`, 'POST'); } catch (e) { tip(e.message); }
  go(address, { keepLayer: true, keepView: true });
}

function poll(address) {
  const t = setInterval(async () => {
    if (state.address !== address) { clearInterval(t); return; }
    const s = await api(`/api/scene?at=${encodeURIComponent(address)}`);
    state.scene = s;
    panel(s);
    if (!s.job || s.job === 'done' || s.job.startsWith('failed')) { clearInterval(t); }
  }, 3000);
}

// ---------------------------------------------------------------- search

let searchTimer = null, hits = [], chosen = -1;
$('search').addEventListener('input', () => {
  clearTimeout(searchTimer);
  searchTimer = setTimeout(async () => {
    const q = $('search').value.trim();
    if (!q) { $('hits').style.display = 'none'; return; }
    try { hits = await api(`/api/search?q=${encodeURIComponent(q)}`); } catch (e) { hits = []; tip(e.message); }
    chosen = hits.length ? 0 : -1;
    showHits();
  }, 220);
});
$('search').addEventListener('keydown', (e) => {
  if (e.key === 'ArrowDown') { chosen = Math.min(hits.length - 1, chosen + 1); showHits(); e.preventDefault(); }
  if (e.key === 'ArrowUp') { chosen = Math.max(0, chosen - 1); showHits(); e.preventDefault(); }
  if (e.key === 'Enter' && hits[chosen]) pick(hits[chosen]);
  if (e.key === 'Escape') $('hits').style.display = 'none';
});
function showHits() {
  const el = $('hits');
  el.innerHTML = hits.length ? hits.map((h, i) => `<div data-i="${i}" class="${i === chosen ? 'on' : ''}">${esc(h.name)}<small>${esc(h.kind)} · ${esc(h.note)}</small></div>`).join('') : '<div>No match among what is known nearby. Type an address to go anywhere.</div>';
  el.style.display = 'block';
  for (const d of el.querySelectorAll('[data-i]')) d.onclick = () => pick(hits[+d.dataset.i]);
}
function pick(h) {
  $('hits').style.display = 'none';
  $('search').value = '';
  go(h.address);
}

// ---------------------------------------------------------------- picking and the loop

const ray = new THREE.Raycaster();
const pointer = new THREE.Vector2();
let downAt = null;
canvas.addEventListener('pointerdown', (e) => { downAt = [e.clientX, e.clientY]; });
canvas.addEventListener('pointerup', (e) => {
  if (!downAt || Math.hypot(e.clientX - downAt[0], e.clientY - downAt[1]) > 6) return;
  const h = hit(e);
  if (h) h.go();
});
canvas.addEventListener('pointermove', (e) => { const h = hit(e); tip(h ? h.label : ''); canvas.style.cursor = h ? 'pointer' : ''; });

function hit(e) {
  pointer.set((e.clientX / window.innerWidth) * 2 - 1, -(e.clientY / window.innerHeight) * 2 + 1);
  ray.setFromCamera(pointer, camera);
  ray.params.Points.threshold = camera.position.distanceTo(controls.target) * 0.012;
  let best = null;
  for (const p of pickables) {
    const found = ray.intersectObject(p.object, false)[0];
    if (!found || (best && found.distance >= best.distance)) continue;
    if (p.pick) best = { distance: found.distance, go: () => p.pick(found.index), label: p.label(found.index) };
    else best = { distance: found.distance, go: p.go, label: p.label };
  }
  return best;
}

function placeLabels() {
  const w = window.innerWidth, h = window.innerHeight;
  const out = [];
  const v = new THREE.Vector3();
  for (const l of labels) {
    v.copy(l.position).project(camera);
    if (v.z > 1 || v.x < -1.1 || v.x > 1.1 || v.y < -1.1 || v.y > 1.1) continue;
    out.push(`<div class="${l.dim ? 'dim' : ''}" style="left:${((v.x + 1) / 2) * w}px;top:${((1 - v.y) / 2) * h}px">${esc(l.text)}</div>`);
  }
  $('labels').innerHTML = out.join('');
}

function loop() {
  controls.update();
  if (spinning && !state.flat) spinning.rotation.y += 0.0008;
  renderer.render(scene, camera);
  placeLabels();
  requestAnimationFrame(loop);
}

resize();
loop();
api('/api/start').then((s) => { state.start = s; fromHash(); }).catch((e) => tip(e.message));
