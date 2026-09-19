/* ============================================================================
   Drop-3D desktop — application core
   ----------------------------------------------------------------------------
   No framework and no build step.  A UI this size does not need one, and the
   packaged application has to run from a frozen directory where a bundler's
   assumptions about the filesystem stop holding.

   Layout of this file:
     i18n + DOM helpers -> bridge -> presentation pieces -> charts -> router
   Everything the pages use is re-exported on `window.Drop3D`, so pages.js stays
   declarative.
   ========================================================================== */

(function () {
  'use strict';

  /* ------------------------------------------------------------------- state */

  const state = {
    lang: 'zh',
    page: null,
    info: null,
    ready: false,
    busy: false,
  };

  const pages = [];

  /* -------------------------------------------------------------------- i18n */

  function t(key, vars) {
    const table = I18N[state.lang] || I18N.zh;
    let value = table[key];
    if (value === undefined) value = I18N.zh[key];
    if (value === undefined) return key;
    if (vars) {
      value = String(value).replace(/\{(\w+)\}/g, (m, name) =>
        Object.prototype.hasOwnProperty.call(vars, name) ? vars[name] : m);
    }
    return value;
  }

  function setLang(lang) {
    if (!I18N[lang]) return;
    state.lang = lang;
    document.documentElement.lang = lang === 'zh' ? 'zh' : 'en';
    document.querySelectorAll('#lang-switch button').forEach((b) => {
      b.classList.toggle('active', b.dataset.lang === lang);
    });
    try { localStorage.setItem('drop3d.lang', lang); } catch (e) { /* not fatal */ }
    renderSidebar();
    render();
  }

  /* --------------------------------------------------------------- DOM helper */

  function h(tag, attrs, ...kids) {
    const el = document.createElement(tag);
    if (attrs) {
      for (const [k, v] of Object.entries(attrs)) {
        if (v === null || v === undefined || v === false) continue;
        if (k === 'class') el.className = v;
        else if (k === 'text') el.textContent = v;
        else if (k === 'html') el.innerHTML = v;
        else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
        else if (k.startsWith('on') && typeof v === 'function') {
          el.addEventListener(k.slice(2).toLowerCase(), v);
        } else if (k === 'dataset') Object.assign(el.dataset, v);
        else el.setAttribute(k, v);
      }
    }
    for (const kid of kids.flat(Infinity)) {
      if (kid === null || kid === undefined || kid === false) continue;
      el.appendChild(kid instanceof Node ? kid : document.createTextNode(String(kid)));
    }
    return el;
  }

  const NS = 'http://www.w3.org/2000/svg';
  function s(tag, attrs, ...kids) {
    const el = document.createElementNS(NS, tag);
    if (attrs) {
      for (const [k, v] of Object.entries(attrs)) {
        if (v === null || v === undefined || v === false) continue;
        if (k === 'text') el.textContent = v;
        else el.setAttribute(k, v);
      }
    }
    for (const kid of kids.flat(Infinity)) {
      if (kid === null || kid === undefined || kid === false) continue;
      el.appendChild(kid instanceof Node ? kid : document.createTextNode(String(kid)));
    }
    return el;
  }

  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }

  /* --------------------------------------------------------------- formatting */

  function num(value, digits) {
    if (value === null || value === undefined || value === '') return t('common.none');
    if (typeof value !== 'number' || !isFinite(value)) return String(value);
    if (value !== 0 && (Math.abs(value) < 1e-4 || Math.abs(value) >= 1e6)) {
      return value.toExponential(2);
    }
    const d = digits === undefined ? (Math.abs(value) >= 100 ? 2 : 4) : digits;
    return Number(value.toFixed(d)).toString();
  }

  function pct(value, digits) {
    if (typeof value !== 'number' || !isFinite(value)) return t('common.none');
    return (value * 100).toFixed(digits === undefined ? 1 : digits) + '%';
  }

  /* --------------------------------------------------------------- the bridge */

  async function bridge(method, payload) {
    const api = window.pywebview && window.pywebview.api;
    if (!api || typeof api[method] !== 'function') {
      throw new Error(t('err.noBridge') + ' (' + method + ')');
    }
    // Methods that take no arguments must be called with none.  Passing an empty
    // object looks harmless and is not: it lands as a positional argument and
    // raises TypeError, which is how app_info and probe_liquids silently returned
    // nothing while every payload-taking call worked.
    const result = payload === undefined ? await api[method]() : await api[method](payload);
    if (result && result.ok === false) {
      const err = new Error(result.error || t('err.title'));
      err.payload = result;
      throw err;
    }
    return result;
  }

  /* ----------------------------------------------------------------- feedback */

  function toast(message, kind) {
    const host = document.getElementById('toasts');
    const el = h('div', { class: 'toast ' + (kind || ''), text: message });
    host.appendChild(el);
    setTimeout(() => {
      el.style.transition = 'opacity 240ms, transform 240ms';
      el.style.opacity = '0';
      el.style.transform = 'translateY(10px)';
      setTimeout(() => el.remove(), 260);
    }, kind === 'error' ? 5200 : 2600);
  }

  function busyRow(label) {
    return h('div', { class: 'busy' },
      h('div', { class: 'spinner' }),
      h('span', { text: label || t('common.calculating') }));
  }

  /* ------------------------------------------------------------ presentation */

  function card(title, hint, ...body) {
    const kids = [];
    if (title) {
      kids.push(h('div', { class: 'card-title' },
        h('span', { text: title }),
        hint ? h('span', { class: 'hint', text: hint }) : null));
    }
    kids.push(...body);
    return h('div', { class: 'card' }, kids);
  }

  function stat(label, value, unit, extra) {
    return h('div', { class: 'stat' + (extra ? ' ' + extra : '') },
      h('div', { class: 'k', text: label }),
      h('div', { class: 'v', text: value }),
      unit ? h('div', { class: 'u', text: unit }) : null);
  }

  const VERDICT_CLASS = {
    accept: 'accept',
    accept_with_warning: 'warn',
    reject: 'reject',
  };

  function verdictPill(verdict) {
    if (!verdict) return h('span', { class: 'verdict idle', text: t('verdict.idle') });
    const cls = VERDICT_CLASS[verdict] || 'idle';
    const icon = cls === 'accept' ? '✓' : (cls === 'warn' ? '!' : '×');
    return h('span', { class: 'verdict ' + cls }, icon + ' ' + t('verdict.' + verdict));
  }

  function severityPill(severity) {
    const s = String(severity || '').toLowerCase();
    const cls = s === 'hard' || s === 'blocking' ? 'hard'
      : s === 'soft' || s === 'serious' ? 'soft' : 'skip';
    return h('span', { class: 'pill ' + cls, text: t('sev.' + s) !== 'sev.' + s ? t('sev.' + s) : s });
  }

  function provenancePill(provenance) {
    const p = String(provenance || '').toLowerCase();
    if (!['literature', 'verified', 'engineering'].includes(p)) return null;
    return h('span', { class: 'prov ' + p, text: p });
  }

  function callout(tone, icon, ...body) {
    return h('div', { class: 'callout ' + tone },
      icon ? h('div', { class: 'ico', text: icon }) : null,
      h('div', { class: 'body' }, body));
  }

  /* Text that comes from the library verbatim is English-only: its messages are
     written for a developer reading a traceback, and rewriting them here would
     put words in the library's mouth.  So they are shown, but labelled -- the
     interface chrome is translated, the library's own sentences are quoted. */
  function libraryNote(...body) {
    return h('div', { class: 'callout' },
      h('div', { class: 'ico', text: '⌘' }),
      h('div', { class: 'body' },
        h('div', { class: 'lib-tag', text: t('common.libraryOutput') }),
        ...body));
  }

  /* Translate the library's error codes into an explanation plus a remedy.  The
     remedy is the point: "rejected" on its own tells the user nothing they can
     act on. */
  function remedyList(codes) {
    if (!codes || !codes.length) return null;
    return h('ul', { class: 'tight' }, codes.map((code) => {
      const label = t('code.' + code) !== 'code.' + code ? t('code.' + code) : code;
      const fix = t('code.' + code + '.fix');
      return h('li', {},
        h('b', { text: label }),
        fix !== 'code.' + code + '.fix' ? h('span', { text: ' — ' + fix }) : null);
    }));
  }

  /* The gate table is the centrepiece: every gate, whether it ran, what it
     measured, the threshold it was compared against, and where that threshold
     came from.  A failure is highlighted rather than hidden in a log. */
  function gatesTable(validity) {
    const checks = (validity && validity.checks) || [];
    if (!checks.length) return h('div', { class: 'faint small', text: t('common.none') });
    const rows = checks.map((c) => {
      const passed = c.passed === true
        ? h('span', { class: 'pill ok', text: t('common.passed') })
        : c.passed === false
          ? h('span', { class: 'pill fail', text: t('common.no') })
          : h('span', { class: 'pill skip', text: t('common.notEvaluated') });
      const name = t('gate.' + c.name) !== 'gate.' + c.name ? t('gate.' + c.name) : c.name;
      return h('tr', { class: c.passed === false ? 'row-fail' : '' },
        h('td', {}, h('div', { text: name }), h('div', { class: 'faint mono', style: { fontSize: '10.5px' }, text: c.name })),
        h('td', {}, passed),
        h('td', {}, severityPill(c.severity)),
        h('td', { class: 'num' }, num(c.value)),
        h('td', { class: 'num' }, num(c.threshold)),
        h('td', {}, provenancePill(c.provenance) || h('span', { class: 'faint', text: '—' })));
    });
    return h('div', { class: 'scroll-x' }, h('table', { class: 'data' },
      h('thead', {}, h('tr', {},
        h('th', { text: t('common.gate') }),
        h('th', { text: t('common.passed') }),
        h('th', { text: t('common.severity') }),
        h('th', { class: 'num', text: t('common.value') }),
        h('th', { class: 'num', text: t('common.threshold') }),
        h('th', { text: t('common.provenance') }))),
      h('tbody', {}, rows)));
  }

  /* A refusal is only useful if it says what to change, so the error code and
     its remedy travel together. */
  function errorCodesBlock(codes) {
    if (!codes || !codes.length) return null;
    const rows = codes.map((code) => h('tr', {},
      h('td', { class: 'mono' }, code),
      h('td', { text: t('code.' + code) !== 'code.' + code ? t('code.' + code) : '' }),
      h('td', { class: 'muted', text: t('code.' + code + '.fix') !== 'code.' + code + '.fix' ? t('code.' + code + '.fix') : '' })));
    return card(t('pendant.errorCodes'), null,
      h('div', { class: 'scroll-x' }, h('table', { class: 'data' },
        h('thead', {}, h('tr', {},
          h('th', { text: t('pendant.code') }),
          h('th', { text: '' }),
          h('th', { text: t('pendant.fix') }))),
        h('tbody', {}, rows))));
  }

  function field(label, control, unit) {
    return h('div', { class: 'field' },
      h('label', {}, label, unit ? h('span', { class: 'unit', text: '· ' + unit }) : null),
      control);
  }

  function numberField(value, onchange, opts) {
    const o = opts || {};
    const input = h('input', {
      type: 'number',
      value: value === null || value === undefined ? '' : value,
      step: o.step === undefined ? 'any' : o.step,
      min: o.min, max: o.max,
    });
    input.addEventListener('change', () => {
      const v = input.value === '' ? null : Number(input.value);
      onchange(v);
    });
    return input;
  }

  function selectField(options, selected, onchange) {
    const sel = h('select', {});
    options.forEach((opt) => {
      const value = typeof opt === 'string' ? opt : opt.value;
      const label = typeof opt === 'string' ? opt : opt.label;
      sel.appendChild(h('option', { value, selected: String(value) === String(selected) }, label));
    });
    sel.addEventListener('change', () => onchange(sel.value));
    return sel;
  }

  function checkField(label, checked, onchange) {
    const input = h('input', { type: 'checkbox', checked: checked ? 'checked' : null });
    input.addEventListener('change', () => onchange(input.checked));
    return h('label', { class: 'check' }, input, h('span', { text: label }));
  }

  function dataTable(columns, rows) {
    return h('div', { class: 'scroll-x' }, h('table', { class: 'data' },
      h('thead', {}, h('tr', {}, columns.map((c) =>
        h('th', { class: c.num ? 'num' : '', text: c.label })))),
      h('tbody', {}, rows)));
  }

  function rawJson(value) {
    return h('details', {},
      h('summary', { class: 'faint small', style: { cursor: 'pointer', padding: '4px 0' }, text: t('common.rawJson') }),
      h('pre', {
        class: 'mono small',
        style: { maxHeight: '260px', overflow: 'auto', background: 'rgba(0,0,0,0.25)', padding: '10px', borderRadius: '10px', margin: '8px 0 0' },
        text: JSON.stringify(value, null, 1),
      }));
  }

  /* ------------------------------------------------------------------ charts */

  function plotFrame(width, height) {
    return { width, height, left: 46, right: 12, top: 12, bottom: 30 };
  }

  function lineChart(series, opts) {
    const o = Object.assign({ width: 520, height: 220, xLabel: '', yLabel: '' }, opts);
    const f = plotFrame(o.width, o.height);
    const all = series.flatMap((sr) => sr.points);
    if (!all.length) return h('svg', { class: 'chart' });
    const xs = all.map((p) => p[0]);
    const ys = all.map((p) => p[1]);
    let x0 = Math.min(...xs), x1 = Math.max(...xs);
    let y0 = Math.min(...ys), y1 = Math.max(...ys);
    if (x0 === x1) x1 = x0 + 1;
    if (y0 === y1) { y0 -= 1; y1 += 1; }
    const pad = (y1 - y0) * 0.08;
    y0 -= pad; y1 += pad;
    // Padding must not invent negative values: P_s and uncertainties are
    // non-negative, and an axis running below zero suggests otherwise.
    if (Math.min(...ys) >= 0 && y0 < 0) y0 = 0;

    const px = (x) => f.left + (x - x0) / (x1 - x0) * (o.width - f.left - f.right);
    const py = (y) => o.height - f.bottom - (y - y0) / (y1 - y0) * (o.height - f.top - f.bottom);

    const g = s('svg', { class: 'chart', viewBox: `0 0 ${o.width} ${o.height}`, preserveAspectRatio: 'xMidYMid meet' });

    // gridlines + y labels
    const ticks = 4;
    for (let i = 0; i <= ticks; i++) {
      const v = y0 + (y1 - y0) * (i / ticks);
      const y = py(v);
      g.appendChild(s('line', { class: 'gridline', x1: f.left, x2: o.width - f.right, y1: y, y2: y }));
      g.appendChild(s('text', { class: 'label', x: f.left - 7, y: y + 3.5, 'text-anchor': 'end', text: num(v, 3) }));
    }
    // x labels
    for (let i = 0; i <= ticks; i++) {
      const v = x0 + (x1 - x0) * (i / ticks);
      g.appendChild(s('text', {
        class: 'label', x: px(v), y: o.height - 10, 'text-anchor': 'middle', text: num(v, 3),
      }));
    }
    g.appendChild(s('line', { class: 'axis', x1: f.left, x2: o.width - f.right, y1: py(y0), y2: py(y0) }));

    series.forEach((sr, idx) => {
      const d = sr.points.map((p, i) => (i ? 'L' : 'M') + px(p[0]).toFixed(2) + ' ' + py(p[1]).toFixed(2)).join(' ');
      const stroke = sr.color || (idx === 0 ? 'var(--accent)' : 'var(--teal)');
      g.appendChild(s('path', {
        d,
        fill: 'none',
        stroke,
        'stroke-width': sr.width || 2,
        'stroke-linejoin': 'round',
        'stroke-dasharray': sr.dashed ? '5 4' : null,
      }));
      (sr.markers || []).forEach((m) => {
        // Markers are {x, y, color}; reading them as a pair would silently
        // produce NaN coordinates, and SVG draws a NaN circle at the origin --
        // which looks like a deliberate annotation in the top-left corner.
        g.appendChild(s('circle', {
          cx: px(m.x), cy: py(m.y), r: m.r || 4,
          fill: m.color || stroke, stroke: 'rgba(255,255,255,0.85)', 'stroke-width': 1.5,
        }));
      });
    });

    const wrap = h('div', {});
    wrap.appendChild(g);
    if (o.xLabel || o.yLabel) {
      wrap.appendChild(h('div', { class: 'legend' },
        o.xLabel ? h('span', { text: o.xLabel }) : null,
        o.yLabel ? h('span', { text: o.yLabel }) : null));
    }
    if (series.some((sr) => sr.label)) {
      wrap.appendChild(h('div', { class: 'legend' }, series.filter((sr) => sr.label).map((sr) =>
        h('span', {}, h('i', { style: { background: sr.color || 'var(--accent)' } }), sr.label))));
    }
    return wrap;
  }

  function barChart(items, opts) {
    const o = Object.assign({ width: 520, rowHeight: 30, labelWidth: 150 }, opts);
    if (!items.length) return h('div', { class: 'faint small', text: t('common.none') });
    const height = items.length * o.rowHeight + 10;
    const width = o.width;
    const max = Math.max(...items.map((i) => Math.abs(i.value) || 0), 1e-12);
    const trackX = o.labelWidth + 8;
    const trackW = width - trackX - 62;

    const g = s('svg', { class: 'chart', viewBox: `0 0 ${width} ${height}`, preserveAspectRatio: 'xMidYMid meet' });
    const defs = s('defs', {});
    const grad = s('linearGradient', { id: 'barGrad', x1: '0', y1: '0', x2: '1', y2: '0' });
    grad.appendChild(s('stop', { offset: '0%', 'stop-color': 'rgba(10,132,255,0.95)' }));
    grad.appendChild(s('stop', { offset: '100%', 'stop-color': 'rgba(64,200,224,0.85)' }));
    defs.appendChild(grad);
    g.appendChild(defs);

    items.forEach((item, i) => {
      const y = i * o.rowHeight + 6;
      g.appendChild(s('text', {
        class: 'label', x: o.labelWidth, y: y + 14, 'text-anchor': 'end', text: item.label,
      }));
      g.appendChild(s('rect', {
        x: trackX, y: y + 4, width: trackW, height: 14, rx: 7, fill: 'rgba(255,255,255,0.07)',
      }));
      const w = Math.max(2, (Math.abs(item.value) / max) * trackW);
      g.appendChild(s('rect', {
        x: trackX, y: y + 4, width: w, height: 14, rx: 7,
        fill: item.color || 'url(#barGrad)',
      }));
      g.appendChild(s('text', {
        class: 'label', x: trackX + trackW + 8, y: y + 15, text: item.display || num(item.value, 3),
      }));
    });
    return h('div', {}, g);
  }

  /* ---------------------------------------------------------- image overlay */

  /* Draws the measurement image with the extracted profile and the fitted curve
     on top, in the image's own pixel coordinates.  `step` is the display
     downsample stride reported by the backend; ignoring it would misplace the
     overlay on any image larger than the display limit. */
  function overlayCanvas(payload, opts) {
    const o = Object.assign({ showProfile: true, showFit: true, showAxis: true }, opts);
    const wrap = h('div', { class: 'canvas-wrap' });
    const display = payload.display || {};
    const step = display.step || 1;
    const w = display.width || 320;
    const hgt = display.height || 240;
    const dpr = window.devicePixelRatio || 1;

    const canvas = h('canvas', {});
    canvas.width = w * dpr;
    canvas.height = hgt * dpr;
    // Display size is left to CSS: a 240x320 test image shown at 1:1 would sit
    // marooned in the middle of a wide card.  The backing store stays at the
    // image's own resolution, so scaling up costs nothing in accuracy and the
    // overlay coordinates never change.
    wrap.appendChild(canvas);

    const ctx = canvas.getContext('2d');
    const img = new Image();
    img.onload = () => {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, w, hgt);
      ctx.drawImage(img, 0, 0, w, hgt);

      const overlay = payload.overlay || {};
      const scale = 1 / step;
      const drawPath = (points, color, width, dash) => {
        if (!points || points.length < 2) return;
        ctx.save();
        ctx.strokeStyle = color;
        ctx.lineWidth = width;
        ctx.setLineDash(dash || []);
        ctx.beginPath();
        points.forEach((p, i) => {
          const x = p[0] * scale;
          const y = p[1] * scale;
          if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        });
        ctx.stroke();
        ctx.restore();
      };

      if (o.showAxis && overlay.axis_x) {
        ctx.save();
        ctx.strokeStyle = 'rgba(255,214,10,0.55)';
        ctx.setLineDash([6, 5]);
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(overlay.axis_x * scale, 0);
        ctx.lineTo(overlay.axis_x * scale, hgt);
        ctx.stroke();
        ctx.restore();
      }
      if (o.showProfile) drawPath(overlay.profile, 'rgba(64,200,224,0.95)', 1.6);
      if (o.showFit) drawPath(overlay.fit_curve, 'rgba(255,105,180,0.95)', 2, [7, 4]);

      ctx.font = '600 10px "Segoe UI", sans-serif';
      ctx.fillStyle = 'rgba(64,200,224,0.98)';
      ctx.fillText('— profile', 8, 14);
      ctx.fillStyle = 'rgba(255,105,180,0.98)';
      ctx.fillText('— fit', 8, 27);
    };
    img.src = payload.image;
    return wrap;
  }

  /* ---------------------------------------------------------------- dropzone */

  function dropzone(onImage, onPick) {
    const zone = h('div', { class: 'dropzone' },
      h('div', { class: 'big', text: '⇩' }),
      h('div', { text: t('common.dropHint') }),
      h('button', { class: 'btn small ghost', onclick: (e) => { e.stopPropagation(); onPick(); } },
        t('common.loadImage')));

    zone.addEventListener('dragover', (e) => { e.preventDefault(); zone.classList.add('over'); });
    zone.addEventListener('dragleave', () => zone.classList.remove('over'));
    zone.addEventListener('drop', (e) => {
      e.preventDefault();
      zone.classList.remove('over');
      const file = e.dataTransfer.files && e.dataTransfer.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = () => onImage(reader.result, file.name);
      reader.readAsDataURL(file);
    });
    return zone;
  }

  /* ------------------------------------------------------------------ router */

  function registerPage(page) { pages.push(page); }

  function renderSidebar() {
    const nav = document.getElementById('sidebar');
    clear(nav);
    nav.appendChild(h('div', { class: 'brand' },
      h('div', { class: 'name' }, h('span', { class: 'dot' }), 'Drop-3D'),
      h('div', { class: 'ver', text: 'v' + ((state.info && state.info.drop3d_version) || '0.1.0') })));

    let lastGroup = null;
    pages.forEach((page) => {
      if (page.group !== lastGroup) {
        lastGroup = page.group;
        nav.appendChild(h('div', { class: 'nav-group', text: t(page.group) }));
      }
      nav.appendChild(h('div', {
        class: 'nav-item' + (state.page === page.id ? ' active' : ''),
        onclick: () => navigate(page.id),
      }, h('span', { class: 'ico', text: page.icon }), h('span', { text: t(page.title) })));
    });

    nav.appendChild(h('div', { class: 'foot' }, t('app.tagline')));
  }

  function currentPage() { return pages.find((p) => p.id === state.page) || pages[0]; }

  function navigate(id) {
    state.page = id;
    try { localStorage.setItem('drop3d.page', id); } catch (e) { /* not fatal */ }
    renderSidebar();
    render();
  }

  function render() {
    const view = document.getElementById('view');
    const page = currentPage();
    if (!page) return;
    clear(view);
    const content = page.render();
    content.classList.add('fade-in');
    view.appendChild(content);
    view.scrollTop = 0;
  }

  /* -------------------------------------------------------------------- boot */

  async function loadInfo(attempts) {
    // The bridge is injected asynchronously, and on a loaded machine it can take
    // a moment after the DOM is ready.  Retrying here (rather than blocking the
    // first paint on it) means the interface appears immediately and simply fills
    // in the version numbers when they arrive.
    const tries = attempts || 14;
    for (let i = 0; i < tries; i++) {
      try {
        state.info = await bridge('app_info');
        return true;
      } catch (e) {
        await new Promise((r) => setTimeout(r, 150));
      }
    }
    state.info = null;
    return false;
  }

  async function boot() {
    try {
      const saved = localStorage.getItem('drop3d.lang');
      if (saved && I18N[saved]) state.lang = saved;
    } catch (e) { /* not fatal */ }

    document.querySelectorAll('#lang-switch button').forEach((b) => {
      b.classList.toggle('active', b.dataset.lang === state.lang);
      b.addEventListener('click', () => setLang(b.dataset.lang));
    });

    const winClose = document.getElementById('win-close');
    const winMin = document.getElementById('win-min');
    const winMax = document.getElementById('win-max');
    if (winClose) winClose.addEventListener('click', () => bridge('window_close').catch(() => {}));
    if (winMin) winMin.addEventListener('click', () => bridge('window_minimise').catch(() => {}));
    if (winMax) winMax.addEventListener('click', () => bridge('window_toggle_maximise').catch(() => {}));

    let start = pages.length ? pages[0].id : null;
    try {
      const saved = localStorage.getItem('drop3d.page');
      if (saved && pages.some((p) => p.id === saved)) start = saved;
    } catch (e) { /* not fatal */ }
    state.page = start;

    // First paint does not wait on the backend.
    renderSidebar();
    render();
    state.ready = true;

    if (await loadInfo()) {
      renderSidebar();
      render();
    }

    // Pages that need backend data before first use declare an init(); it runs
    // here, once, rather than on every render.  Kept generic so app.js never has
    // to know a page's name -- the page modules are closed over their own state.
    pages.forEach((page) => {
      if (typeof page.init !== 'function') return;
      Promise.resolve()
        .then(() => page.init())
        .catch(() => {});
    });
  }

  window.Drop3D = {
    state, t, h, s, clear, num, pct,
    bridge, toast, busyRow,
    card, stat, verdictPill, severityPill, provenancePill, callout, libraryNote, remedyList,
    gatesTable, errorCodesBlock, field, numberField, selectField, checkField,
    dataTable, rawJson,
    lineChart, barChart, overlayCanvas, dropzone,
    registerPage, navigate, render, renderSidebar, boot, setLang,
    VERDICT_CLASS,
  };
})();
