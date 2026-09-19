/* ============================================================================
   Drop-3D desktop — the pages
   ----------------------------------------------------------------------------
   Seven working pages plus a system-information page.  Each page owns its own
   local state, talks to the bridge, and renders itself; the router in app.js
   knows nothing about any of them.

   A note on how refusals are presented, because it is the one design decision
   that runs through every page: when a hard gate fails, the headline number is
   *replaced* rather than shown with a warning next to it.  A number on screen
   gets quoted; a refusal with a reason gets fixed.  That asymmetry is the whole
   point of building the interface this way.
   ========================================================================== */

(function () {
  'use strict';

  const D = window.Drop3D;
  const { h, t, num, pct, bridge, card, stat, verdictPill, callout, libraryNote,
    remedyList, field, numberField, selectField, checkField, dataTable, rawJson,
    toast, lineChart, barChart, overlayCanvas, dropzone, gatesTable,
    errorCodesBlock } = D;

  /* ------------------------------------------------------------------ helpers */

  function pageHead(titleKey, leadKey, right) {
    return h('div', { class: 'page-head' },
      h('div', { class: 'grow' },
        h('h1', { text: t(titleKey) }),
        h('p', { text: t(leadKey) })),
      right || null);
  }

  function runButton(label, onRun) {
    const btn = h('button', { class: 'btn primary' }, label || t('common.run'));
    btn.addEventListener('click', async () => {
      if (btn.disabled) return;
      btn.disabled = true;
      const original = btn.textContent;
      btn.textContent = t('common.running');
      try {
        await onRun();
      } catch (err) {
        toast(err.message || String(err), 'error');
      } finally {
        btn.disabled = false;
        btn.textContent = original;
      }
    });
    return btn;
  }

  function errorBlock(err) {
    const payload = err.payload || {};
    return card(t('err.title'), null,
      callout('error', '×',
        h('div', {}, h('b', { text: err.message || String(err) })),
        payload.traceback
          ? h('details', {}, h('summary', { class: 'small', style: { cursor: 'pointer' }, text: 'traceback' }),
            h('pre', { class: 'mono small', style: { overflow: 'auto', maxHeight: '200px' }, text: payload.traceback }))
          : null));
  }

  function synthBadge() {
    return h('span', { class: 'pill skip', text: t('common.syntheticBadge') });
  }

  function parseSeries(text) {
    const rows = String(text || '').split(/\r?\n/)
      .map((line) => line.trim())
      .filter((line) => line && !line.startsWith('#'))
      .map((line) => line.split(/[,\s;]+/).map(Number).filter((v) => !Number.isNaN(v)))
      .filter((parts) => parts.length >= 3);
    return {
      t: rows.map((r) => r[0]),
      area: rows.map((r) => r[1]),
      gamma: rows.map((r) => r[2]),
    };
  }

  /* ==================================================================== 1. overview */

  const overview = {
    id: 'overview', group: 'nav.group.start', icon: '◈', title: 'nav.overview',
    render() {
      const info = D.state.info || {};
      const workflows = [
        ['pendant', '①', 'ov.wf.pendant'],
        ['surface-energy', '②', 'ov.wf.surface_energy'],
        ['uncertainty', '③', 'ov.wf.uncertainty'],
        ['oscillation', '④', 'ov.wf.oscillation'],
        ['sliding', '⑤', 'ov.wf.sliding'],
        ['design', '⑥', 'ov.wf.design'],
      ];
      return h('div', {},
        pageHead('ov.title', 'ov.lead'),
        h('div', { class: 'grid c3 mb' }, workflows.map(([id, n, key]) =>
          h('div', {
            class: 'card',
            style: { cursor: 'pointer' },
            onclick: () => D.navigate(id),
          },
          h('div', { class: 'card-title' }, h('span', { style: { opacity: 0.7 }, text: n }), t(key)),
          h('div', { class: 'small muted', text: t(key + '.d') })))),
        card(t('ov.dataNote'), null,
          callout('warn', '!', h('div', { text: t('ov.dataNoteBody') }))),
        card(t('ov.modules'), t('ov.modulesNote'),
          h('div', { class: 'row' }, (info.modules || []).map((m) =>
            h('span', { class: 'pill skip mono', text: m })))),
        h('div', { class: 'row mt' },
          h('button', { class: 'btn primary', onclick: () => D.navigate('pendant') },
            '▶ ' + t('ov.startHere')),
          h('button', { class: 'btn ghost', onclick: () => D.navigate('about') },
            t('nav.about'))));
    },
  };

  /* ===================================================================== 2. pendant */

  const pendant = (() => {
    const st = {
      pending: null,        // {kind:'synthetic'|'file'|'dataurl', ...}
      preview: null,        // last load_image / make_synthetic_drop payload
      result: null,
      error: null,
      params: {
        bond: 0.3, radius_px: 60.0, noise: 3.0, seed: 0,
        needle_mm: 1.5, needle_px: 60.6,
        liquid_density: 998.0, air_density: 1.184,
        volume_m3: '',
        u_px_frac: 0.005, u_bond_frac: 0.005,
      },
    };

    async function generate() {
      const p = st.params;
      st.pending = { kind: 'synthetic', bond: p.bond, radius_px: p.radius_px, noise: p.noise, seed: p.seed };
      st.preview = await bridge('make_synthetic_drop', {
        bond: p.bond, radius_px: p.radius_px, noise: p.noise, seed: p.seed,
      });
      st.result = null; st.error = null;
      D.render();
    }

    async function pick() {
      const picked = await bridge('pick_image_file');
      if (picked.cancelled || !picked.path) return;
      st.pending = { kind: 'file', path: picked.path };
      st.preview = await bridge('load_image', { source: st.pending });
      st.result = null; st.error = null;
      D.render();
    }

    async function dropped(dataUrl) {
      st.pending = { kind: 'dataurl', data: dataUrl };
      st.preview = await bridge('load_image', { source: st.pending });
      st.result = null; st.error = null;
      D.render();
    }

    async function run() {
      const p = st.params;
      if (!st.pending) { toast(t('pendant.needImage'), 'error'); return; }
      st.error = null;
      try {
        const payload = {
          source: st.pending,
          needle_diameter_mm: p.needle_mm,
          needle_diameter_px: p.needle_px,
          liquid_density: p.liquid_density,
          air_density: p.air_density,
          u_px_size_frac: p.u_px_frac,
          u_bond_frac: p.u_bond_frac,
        };
        if (p.volume_m3 !== '' && p.volume_m3 !== null) payload.volume_m3 = p.volume_m3;
        st.result = await bridge('analyse_pendant_drop', payload);
        D.state.lastPendant = st.result;
      } catch (err) {
        st.error = err;
      }
      D.render();
    }

    async function exportReport() {
      if (!st.result) return;
      const rep = await bridge('export_report', { measurement: st.result });
      const saved = await bridge('save_text_file', { filename: rep.filename, content: rep.content });
      if (!saved.cancelled) toast(t('common.export') + ' → ' + saved.path, 'ok');
    }

    function render() {
      const p = st.params;
      const r = st.result;
      const v = r && r.validity;

      /* ---- image column */
      const imageBody = st.preview
        ? (r ? overlayCanvas(r) : overlayCanvas({
          image: st.preview.image,
          display: st.preview.display,
          overlay: { profile: [], fit_curve: [], axis_x: 0 },
        }, { showProfile: false, showFit: false, showAxis: false }))
        : dropzone((dataUrl) => dropped(dataUrl).catch((e) => toast(e.message, 'error')),
          () => pick().catch((e) => toast(e.message, 'error')));

      const imageCard = card(t('pendant.image'),
        st.result ? t('pendant.overlay') : null,
        imageBody,
        h('div', { class: 'row mt' },
          h('button', {
            class: 'btn small',
            onclick: () => generate().catch((e) => toast(e.message, 'error')),
          }, '⚗ ' + t('common.synthetic')),
          h('button', {
            class: 'btn small ghost',
            onclick: () => pick().catch((e) => toast(e.message, 'error')),
          }, t('common.loadImage')),
          st.preview ? synthBadge() : null,
          r && r.source_label ? h('span', { class: 'faint small', text: r.source_label }) : null));

      /* ---- parameter column */
      const paramCard = card(t('common.parameters'),
        r ? (r.seconds + ' s') : null,
        h('div', { class: 'grid c2' },
          field(t('pendant.bondParam'), numberField(p.bond, (val) => { p.bond = val; }, { step: 0.01, min: 0.001, max: 0.6 })),
          field(t('pendant.needleMm'), numberField(p.needle_mm, (val) => { p.needle_mm = val; }, { step: 0.1 }), 'mm'),
          field(t('pendant.needlePx'), numberField(p.needle_px, (val) => { p.needle_px = val; }, { step: 0.1 }), 'px'),
          field(t('pendant.liquidDensity'), numberField(p.liquid_density, (val) => { p.liquid_density = val; }, { step: 1 }), 'kg/m³'),
          field(t('pendant.airDensity'), numberField(p.air_density, (val) => { p.air_density = val; }, { step: 0.001 }), 'kg/m³'),
          field(t('pendant.volume'), numberField(p.volume_m3, (val) => { p.volume_m3 = val; }, {}), 'm³')),
        callout('info', 'i', h('div', { class: 'small', text: t('pendant.scaleNote') })),
        h('div', { class: 'row mt' },
          runButton(null, run),
          r ? h('button', { class: 'btn ghost small', onclick: () => exportReport().catch((e) => toast(e.message, 'error')) }, t('common.export')) : null));

      /* ---- result column */
      let resultCard = null;
      if (st.error) {
        resultCard = errorBlock(st.error);
      } else if (r) {
        const refused = v && v.verdict === 'reject';
        const unc = r.uncertainty || {};
        const heroStat = refused
          ? h('div', { class: 'stat hero refused' },
            h('div', { class: 'k', text: t('hero.gamma') }),
            h('div', { class: 'v', text: t('hero.refused') }),
            h('div', { class: 'u', text: t('hero.refusedHint') }))
          : stat(t('hero.gamma'), r.gamma_mN_m === null ? t('common.none') : num(r.gamma_mN_m, 3), 'mN/m', 'hero');

        const stats = h('div', { class: 'grid c4' },
          heroStat,
          stat(t('hero.uncertainty'), unc.std === undefined ? t('common.none') : '± ' + num(unc.std, 3), 'mN/m · k=' + num(unc.k, 2)),
          stat(t('hero.bond'), num(r.fit && r.fit.bond, 4), ''),
          stat(t('hero.radius'), num(r.fit && r.fit.radius_px, 2), 'px'));

        const stats2 = h('div', { class: 'grid c3 mt' },
          stat(t('hero.shapeParam'), num(r.fit && r.fit.shape_parameter, 4), ''),
          stat(t('hero.rms'), num(r.fit && r.fit.rms_px, 4), 'px'),
          stat(t('hero.reliability'), r.reliability_class || t('common.none'), ''));

        /* The verdict is explained in the interface's own language first; the
           library's English sentence is quoted underneath, labelled as such. */
        const verdictBody = refused
          ? callout('error', '×', h('div', {},
            h('b', { text: t('hero.refused') }), remedyList(r.error_codes)))
          : (r.error_codes && r.error_codes.length
            ? callout('warn', '!', h('div', {}, remedyList(r.error_codes)))
            : null);

        resultCard = card(t('common.result'), null,
          h('div', { class: 'row mb' }, verdictPill(v && v.verdict), synthBadge(),
            h('span', { class: 'faint small', text: 'Δρ = ' + num(r.calibration && r.calibration.delta_rho, 4) + ' kg/m³ · ' +
              'px = ' + num(r.calibration && r.calibration.px_size_mm, 5) + ' mm/px' })),
          verdictBody,
          v && v.reason ? h('div', { class: 'mt' }, libraryNote(h('div', { class: 'small', text: v.reason }))) : null,
          stats, stats2);
      }

      /* ---- gates + codes + residuals */
      const body = [
        pageHead('pendant.title', 'pendant.lead'),
        h('div', { class: 'grid side' }, h('div', {}, imageCard), h('div', {}, paramCard)),
        resultCard,
        r && v ? card(t('common.gates'), t('pendant.gatesNote'), gatesTable(v)) : null,
        r ? errorCodesBlock(r.error_codes) : null,
        r && r.overlay && r.overlay.residual_px && r.overlay.residual_px.length
          ? card(t('pendant.residual'),
            'RMS = ' + num(r.overlay.residual_rms, 4) + ' px',
            lineChart([{
              points: r.overlay.residual_px.map((val, i) => [i, val]),
              label: t('pendant.residual') + ' (px)',
            }, {
              points: [[0, r.overlay.residual_rms], [r.overlay.residual_px.length - 1, r.overlay.residual_rms]],
              label: 'RMS', color: 'rgba(255,214,10,0.85)', dashed: true, width: 1.5,
            }], { xLabel: 'point index', yLabel: 'px' }),
            h('div', { class: 'small faint mt', text: t('pendant.residualNote') }))
          : null,
        r ? rawJson(r) : null,
      ];
      return h('div', {}, body);
    }

    return { id: 'pendant', group: 'nav.group.measure', icon: '◐', title: 'nav.pendant', render };
  })();

  /* ============================================================ 3. surface energy */

  const surfaceEnergy = (() => {
    const st = {
      rows: [
        { liquid: 'water', angle: 48 },
        { liquid: 'diiodomethane', angle: 72 },
        { liquid: 'ethylene_glycol', angle: 65 },
      ],
      liquids: [],
      result: null,
      error: null,
    };

    async function loadLiquids() {
      const res = await bridge('probe_liquids');
      st.liquids = Object.keys(res.liquids || {});
    }

    async function run() {
      st.error = null;
      try {
        st.result = await bridge('analyse_surface_energy', {
          angles: st.rows.map((r) => Number(r.angle)),
          liquids: st.rows.map((r) => r.liquid),
        });
      } catch (err) { st.error = err; }
      D.render();
    }

    function setPreset(kind) {
      if (kind === 'two') {
        st.rows = [{ liquid: 'water', angle: 48 }, { liquid: 'diiodomethane', angle: 72 }];
      } else {
        st.rows = [
          { liquid: 'water', angle: 48 },
          { liquid: 'diiodomethane', angle: 72 },
          { liquid: 'ethylene_glycol', angle: 65 },
        ];
      }
      D.render();
    }

    function render() {
      const options = st.liquids.length
        ? st.liquids
        : ['water', 'diiodomethane', 'ethylene_glycol'];

      const rows = st.rows.map((row, idx) => h('tr', {},
        h('td', {}, selectField(options, row.liquid, (val) => { row.liquid = val; })),
        h('td', { style: { width: '130px' } }, numberField(row.angle, (val) => { row.angle = val; }, { step: 0.1 }), ' °'),
        h('td', { style: { width: '40px' } }, h('button', {
          class: 'btn small ghost',
          onclick: () => { st.rows.splice(idx, 1); D.render(); },
        }, '−'))));

      const liquidCard = card(t('se.liquids'), t('se.angle') + ' °',
        dataTable([{ label: t('se.liquid') }, { label: t('se.angle') }, { label: '' }], rows),
        h('div', { class: 'row mt' },
          h('button', { class: 'btn small', onclick: () => { st.rows.push({ liquid: options[0], angle: 60 }); D.render(); } }, '+ ' + t('se.addRow')),
          h('button', { class: 'btn small ghost', onclick: () => setPreset('two') }, t('se.preset.waterDii')),
          h('button', { class: 'btn small ghost', onclick: () => setPreset('three') }, t('se.preset.threeLiq')),
          h('div', { class: 'spacer' }),
          runButton(null, run)));

      const models = (st.result && st.result.models) || {};
      const spread = (st.result && st.result.spread) || null;

      let resultCards = null;
      if (st.error) {
        resultCards = errorBlock(st.error);
      } else if (st.result) {
        const keys = Object.keys(models).filter((k) => !k.startsWith('_'));
        const modelRows = keys.map((k) => {
          const m = models[k];
          return h('tr', { class: (m.ok === false) ? 'row-fail' : '' },
            h('td', {}, h('b', { text: k })),
            h('td', { class: 'num' }, num(m.total, 4)),
            h('td', { class: 'num' }, m.dispersive === null ? t('common.none') : num(m.dispersive, 4)),
            h('td', { class: 'num' }, m.polar === null ? t('common.none') : num(m.polar, 4)),
            h('td', { class: 'num' }, num(m.rms_deg, 3)));
        });

        resultCards = h('div', {},
          card(t('common.result'), null,
            h('div', { class: 'grid c3' },
              stat(t('se.spread'), spread ? num(spread.range, 3) : t('common.none'), 'mN/m'),
              stat(t('common.result') + ' min', spread ? num(spread.min, 3) : t('common.none'), 'mN/m'),
              stat(t('common.result') + ' max', spread ? num(spread.max, 3) : t('common.none'), 'mN/m')),
            h('div', { class: 'mt' }, callout('warn', '!', h('div', { text: t('se.spreadNote') }))),
            st.result.recommended
              ? libraryNote(h('div', { text: String(st.result.recommended) }))
              : null),
          card(t('se.model'), null,
            dataTable([
              { label: t('se.model') },
              { label: t('se.total'), num: true },
              { label: t('se.dispersive'), num: true },
              { label: t('se.polar'), num: true },
              { label: t('se.rms'), num: true },
            ], modelRows),
            h('div', { class: 'small faint mt', text: t('se.zismanNote') }),
            h('div', { class: 'small faint', text: t('se.probeNote') })));
      }

      return h('div', {},
        pageHead('se.title', 'se.lead'),
        liquidCard,
        resultCards,
        st.result ? rawJson(st.result) : null);
    }

    return {
      id: 'surface-energy', group: 'nav.group.measure', icon: '◑', title: 'nav.surface-energy',
      render,
      // Fetched once at boot so the liquid dropdowns are populated before the
      // page is first opened; re-rendered only if it is the visible page.
      init: async () => {
        await loadLiquids();
        if (D.state.page === 'surface-energy') D.render();
      },
    };
  })();

  /* ============================================================== 4. uncertainty */

  const uncertainty = (() => {
    const st = {
      inputs: {
        delta_rho: { nominal: 997.0, rel: 0.001 },
        radius_px: { nominal: 60.94, rel: 0.005 },
        px_size_mm: { nominal: 0.02475, rel: 0.005 },
        bond: { nominal: 0.2934, rel: 0.005 },
      },
      k: 2.0,
      useMc: true,
      nSamples: 20000,
      result: null,
      error: null,
      conformal: {
        scores: '', alpha: 0.05, maxHalfWidth: 1.0, value: '',
        result: null, error: null,
      },
    };

    function pullFromPendant() {
      const last = D.state.lastPendant;
      if (!last || !last.fit) { toast(t('pendant.needImage'), 'error'); return; }
      st.inputs.delta_rho.nominal = last.calibration.delta_rho;
      st.inputs.radius_px.nominal = last.fit.radius_px;
      st.inputs.px_size_mm.nominal = last.calibration.px_size_mm;
      st.inputs.bond.nominal = last.fit.bond;
      D.render();
      toast(t('common.result') + ' ← ' + t('nav.pendant'), 'ok');
    }

    async function run() {
      st.error = null;
      try {
        st.result = await bridge('analyse_uncertainty', {
          delta_rho: st.inputs.delta_rho.nominal,
          radius_px: st.inputs.radius_px.nominal,
          px_size_mm: st.inputs.px_size_mm.nominal,
          bond: st.inputs.bond.nominal,
          u_delta_rho_frac: st.inputs.delta_rho.rel,
          u_radius_frac: st.inputs.radius_px.rel,
          u_px_size_frac: st.inputs.px_size_mm.rel,
          u_bond_frac: st.inputs.bond.rel,
          k: st.k,
          monte_carlo: st.useMc,
          n_samples: st.nSamples,
        });
      } catch (err) { st.error = err; }
      D.render();
    }

    async function runConformal() {
      st.conformal.error = null;
      try {
        const scores = st.conformal.scores.split(/[\s,;]+/).map(Number).filter((v) => !Number.isNaN(v));
        if (!scores.length) throw new Error(t('unc.scores'));
        if (!st.conformal.value) {
          // No measurement supplied: report what the calibration set supports.
          const cal = await bridge('conformal_calibrate', { scores, alpha: st.conformal.alpha });
          st.conformal.result = { calibrate_only: true, calibration: cal.calibration, required_size: cal.required_size };
        } else {
          st.conformal.result = await bridge('conformal_predict', {
            scores,
            alpha: st.conformal.alpha,
            value: Number(st.conformal.value),
            max_half_width: st.conformal.maxHalfWidth,
          });
        }
      } catch (err) { st.conformal.error = err; }
      D.render();
    }

    function inputRow(key, labelKey, unit) {
      const item = st.inputs[key];
      return h('tr', {},
        h('td', {}, h('div', { text: t(labelKey) }), h('div', { class: 'faint mono', style: { fontSize: '10.5px' }, text: key })),
        h('td', { style: { width: '160px' } }, numberField(item.nominal, (val) => { item.nominal = val; }, {})),
        h('td', { class: 'faint small', style: { width: '60px' }, text: unit || '' }),
        h('td', { style: { width: '150px' } }, numberField((item.rel * 100).toFixed(3), (val) => { item.rel = (val || 0) / 100; }, { step: 0.01 })),
        h('td', { class: 'faint small', style: { width: '40px' }, text: '%' }));
    }

    function render() {
      const gum = st.result && st.result.gum;
      const mc = st.result && st.result.monte_carlo;

      const contributionItems = gum && gum.contributions
        ? Object.entries(gum.contributions)
          .filter(([, val]) => val > 0)
          .sort((a, b) => b[1] - a[1])
          .map(([key, val]) => ({
            label: t('input.' + key) !== 'input.' + key ? t('input.' + key) : key,
            value: val * 100,
            display: pct(val, 1),
          }))
        : [];

      let resultCards = null;
      if (st.error) resultCards = errorBlock(st.error);
      else if (st.result) {
        resultCards = card(t('common.result'), null,
          h('div', { class: 'grid c4' },
            stat(t('hero.gamma'), num(gum.value, 4), 'mN/m', 'hero'),
            stat('u (k=' + num(gum.k, 2) + ')', '± ' + num(gum.std, 4), 'mN/m'),
            stat('relative', pct(gum.relative, 3), ''),
            stat(t('unc.mc'), mc ? '± ' + num(mc.std, 4) : t('common.none'), 'mN/m')),
          h('div', { class: 'small faint mt', text: gum.method || '' }),
          contributionItems.length
            ? h('div', { class: 'mt' },
              h('div', { class: 'card-title', text: t('unc.contrib') }),
              barChart(contributionItems),
              h('div', { class: 'small faint mt', text: t('unc.contribNote') }))
            : null,
          mc ? h('div', { class: 'mt' }, callout('info', 'i', h('div', {},
            h('b', { text: t('unc.mc') + ': ' }),
            'μ = ' + num(mc.value, 4) + ', σ = ' + num(mc.std, 4) + ' mN/m · ' + (mc.method || ''),
            h('div', { class: 'small faint', text: 'GUM σ = ' + num(gum.std, 4) + ' mN/m' })))) : null);
      }

      const conf = st.conformal;
      const confResult = conf.result;
      const confCard = card(t('unc.conformal'), null,
        h('div', { class: 'grid c3' },
          field(t('unc.scores'), h('textarea', {
            placeholder: t('unc.scoresHint'),
            oninput: (e) => { conf.scores = e.target.value; },
          }, conf.scores)),
          field(t('unc.alpha'), numberField(conf.alpha, (val) => { conf.alpha = val; }, { step: 0.01, min: 0.001, max: 0.5 })),
          field(t('unc.maxHalfWidth'), numberField(conf.maxHalfWidth, (val) => { conf.maxHalfWidth = val; }, { step: 0.1 }), 'mN/m')),
        h('div', { class: 'grid c2 mt' },
          field(t('unc.predictValue'), numberField(conf.value, (val) => { conf.value = val; }, {})),
          h('div', { class: 'field' }, h('label', { text: ' ' }), runButton(null, runConformal))),
        h('div', { class: 'mt' }, callout('info', 'i', h('div', { class: 'small', text: t('unc.conformalNote') }))),
        conf.error ? h('div', { class: 'mt' }, errorBlock(conf.error)) : null,
        confResult ? h('div', { class: 'mt' },
          confResult.calibrate_only
            ? dataTable([{ label: '' }, { label: '' }], [
              h('tr', {}, h('td', { text: t('unc.requiredSize') }), h('td', { class: 'num', text: String(confResult.required_size) })),
              h('tr', {}, h('td', { text: 'half width' }), h('td', { class: 'num', text: num(confResult.calibration.half_width, 4) })),
              h('tr', {}, h('td', { text: 'level achieved' }), h('td', { class: 'num', text: num(confResult.calibration.level_achieved, 4) })),
            ])
            : h('div', {},
              h('div', { class: 'grid c3' },
                stat(t('unc.interval'), '[' + num(confResult.interval.lo, 3) + ', ' + num(confResult.interval.hi, 3) + ']', 'mN/m'),
                stat('half width', num((confResult.interval.hi - confResult.interval.lo) / 2, 4), 'mN/m'),
                stat(t('unc.widthDecision'), confResult.decision.reportable ? '✓' : '×', '')),
              confResult.decision.reportable
                ? callout('ok', '✓', h('div', { text: 'reportable' }))
                : callout('warn', '!', h('div', {}, h('b', { text: 'not reportable' }), ' — ' + (confResult.decision.reason || ''))),
              h('div', { class: 'small faint mt', text: t('unc.widthNote') })))
          : null);

      return h('div', {},
        pageHead('unc.title', 'unc.lead',
          h('button', { class: 'btn ghost small', onclick: pullFromPendant }, '← ' + t('nav.pendant'))),
        card(t('unc.inputs'), t('common.unit') + ': SI',
          dataTable([
            { label: t('unc.quantity') },
            { label: t('unc.nominal') },
            { label: '' },
            { label: t('unc.relUnc') },
            { label: '' },
          ], [
            inputRow('delta_rho', 'pendant.deltaRho', 'kg/m³'),
            inputRow('radius_px', 'hero.radius', 'px'),
            inputRow('px_size_mm', 'pendant.pxMm', 'mm/px'),
            inputRow('bond', 'hero.bond', ''),
          ]),
          h('div', { class: 'row mt' },
            field(t('unc.k'), numberField(st.k, (val) => { st.k = val; }, { step: 0.1, min: 1 })),
            field(t('unc.nSamples'), numberField(st.nSamples, (val) => { st.nSamples = val; }, { step: 1000, min: 1000 })),
            h('div', { class: 'field' }, h('label', { text: ' ' }),
              checkField(t('unc.runMc'), st.useMc, (val) => { st.useMc = val; })),
            h('div', { class: 'field' }, h('label', { text: ' ' }), runButton(null, run)))),
        resultCards,
        confCard,
        st.result ? rawJson(st.result) : null);
    }

    return { id: 'uncertainty', group: 'nav.group.measure', icon: '◒', title: 'nav.uncertainty', render };
  })();

  /* ============================================================== 5. oscillation */

  const oscillation = (() => {
    const st = {
      text: '', frequency: 1.0, fitLag: true, fitFreq: false,
      lag: 15.0, deltaDeg: 20.0, nCycles: 3.0, samplesPerCycle: 120,
      result: null, error: null,
    };

    async function generate() {
      const res = await bridge('make_synthetic_oscillation', {
        frequency_hz: st.frequency,
        n_cycles: st.nCycles,
        samples_per_cycle: st.samplesPerCycle,
        instrument_lag_deg: st.lag,
        delta_deg: st.deltaDeg,
      });
      st.text = res.t.map((tv, i) => [tv.toFixed(5), res.area[i].toFixed(5), res.gamma[i].toFixed(5)].join(', ')).join('\n');
      st.result = null; st.error = null;
      D.render();
    }

    async function run() {
      st.error = null;
      try {
        const series = parseSeries(st.text);
        if (series.t.length < 8) throw new Error(t('osc.csvHint'));
        st.result = await bridge('analyse_oscillation', {
          t: series.t, area: series.area, gamma: series.gamma,
          frequency_hz: st.frequency,
          fit_instrument_lag: st.fitLag,
          fit_frequency: st.fitFreq,
        });
      } catch (err) { st.error = err; }
      D.render();
    }

    function render() {
      const r = st.result;
      let resultCards = null;
      if (st.error) resultCards = errorBlock(st.error);
      else if (r) {
        const series = parseSeries(st.text);
        resultCards = h('div', {},
          card(t('common.result'), null,
            h('div', { class: 'grid c4' },
              stat(t('osc.storage'), num(r.storage, 4), 'mN/m', 'hero'),
              stat(t('osc.loss'), num(r.loss, 4), 'mN/m'),
              stat(t('osc.delta'), num(r.delta_deg, 3), '°'),
              stat(t('osc.lag'), num(r.instrument_lag_deg, 3), '°')),
            h('div', { class: 'grid c3 mt' },
              stat(t('osc.modulus'), num(r.modulus, 4), 'mN/m'),
              stat(t('osc.fpc'), num(r.frames_per_cycle, 1), ''),
              stat('n', String(r.n_points))),
            r.warnings && r.warnings.length
              ? h('div', { class: 'mt' }, libraryNote(h('ul', { class: 'tight' },
                r.warnings.map((w) => h('li', { text: w })))))
              : null,
            h('div', { class: 'mt' }, callout('info', 'i', h('div', { class: 'small', text: t('osc.lagWarning') }))),
            h('div', { class: 'small faint mt', text: t('osc.fpcNote') })),
          card(t('osc.series'), null,
            h('div', { class: 'grid c2' },
              h('div', {}, h('div', { class: 'small muted mb', text: 'γ(t) mN/m' }),
                lineChart([{ points: series.t.map((tv, i) => [tv, series.gamma[i]]) }], { height: 170 })),
              h('div', {}, h('div', { class: 'small muted mb', text: 'A(t)' }),
                lineChart([{ points: series.t.map((tv, i) => [tv, series.area[i]]), color: 'var(--teal)' }], { height: 170 })))));
      }

      return h('div', {},
        pageHead('osc.title', 'osc.lead'),
        card(t('osc.data'), t('osc.csvHint'),
          h('textarea', {
            style: { minHeight: '150px' },
            placeholder: 't, area, gamma\n0.00000, 20.12941, 72.22943\n…',
            oninput: (e) => { st.text = e.target.value; },
          }, st.text),
          h('div', { class: 'grid c4 mt' },
            field(t('osc.freq'), numberField(st.frequency, (val) => { st.frequency = val; }, { step: 0.1 }), 'Hz'),
            field(t('osc.lag') + ' (synthetic)', numberField(st.lag, (val) => { st.lag = val; }, { step: 1 }), '°'),
            field('δ (synthetic)', numberField(st.deltaDeg, (val) => { st.deltaDeg = val; }, { step: 1 }), '°'),
            field('cycles', numberField(st.nCycles, (val) => { st.nCycles = val; }, { step: 0.5 }))),
          h('div', { class: 'row mt' },
            h('div', { class: 'field' }, h('label', { text: ' ' }),
              checkField(t('osc.fitLag'), st.fitLag, (val) => { st.fitLag = val; })),
            h('div', { class: 'field' }, h('label', { text: ' ' }),
              checkField(t('osc.fitFreq'), st.fitFreq, (val) => { st.fitFreq = val; })),
            h('div', { class: 'spacer' }),
            h('button', { class: 'btn small', onclick: () => generate().catch((e) => toast(e.message, 'error')) },
              '⚗ ' + t('common.synthetic')),
            runButton(null, run))),
        resultCards,
        r ? rawJson(r) : null);
    }

    return { id: 'oscillation', group: 'nav.group.measure', icon: '◓', title: 'nav.oscillation', render };
  })();

  /* ================================================================= 6. sliding */

  const sliding = (() => {
    const st = {
      acq: {
        scenario: 'contact_line_dynamics', fps: 100.0, um_per_px: 0.7,
        contact_line_resolved: true, temperature_c: 22.0, relative_humidity_pct: 45.0,
        drop_volume_ul: 45.0, tilt_rate_deg_s: 1.0, duration_s: 2.0,
        refilled: false, axisymmetric_fit_used: false,
      },
      acqResult: null, acqError: null,
      fur: { width_m: 0.002, gamma_mN_m: 72.0, theta_a_deg: 100.0, theta_r_deg: 80.0 },
      furResult: null,
      trk: { vx: 50, noise: 0.2, n_frames: 50, duration_s: 1.0, max_speed: 200, position_sigma: 0.5 },
      trkResult: null,
    };

    async function runAcq() {
      st.acqError = null;
      try {
        st.acqResult = await bridge('assess_acquisition', st.acq);
      } catch (err) { st.acqError = err; }
      D.render();
    }

    async function runFur() {
      st.furResult = await bridge('analyse_furmidge', st.fur);
      D.render();
    }

    async function runTrk() {
      const synth = await bridge('make_synthetic_trajectory', {
        vx: st.trk.vx, noise: st.trk.noise, n_frames: st.trk.n_frames, duration_s: st.trk.duration_s,
      });
      st.trkResult = await bridge('analyse_tracking', {
        times: synth.data.times,
        detections: synth.data.detections,
        max_speed: st.trk.max_speed,
        position_sigma: st.trk.position_sigma,
      });
      st.trkResult.true_velocity = synth.data.true_velocity;
      D.render();
    }

    function render() {
      const a = st.acqResult;
      const hazardRows = a && a.hazards
        ? a.hazards.map((hz) => h('tr', { class: hz.severity === 'blocking' ? 'row-fail' : '' },
          h('td', { class: 'mono' }, hz.name),
          h('td', {}, D.severityPill(hz.severity)),
          h('td', { class: 'small muted' }, hz.message)))
        : [];

      let acqCard = card(t('sl.acq'), t('sl.verdict'),
        h('div', { class: 'grid c4' },
          field(t('sl.scenario'), selectField([
            { value: 'contact_line_dynamics', label: t('sl.scenario.contact_line_dynamics') },
            { value: 'sliding', label: t('sl.scenario.sliding') },
            { value: 'oscillating', label: t('sl.scenario.oscillating') },
            { value: 'impact', label: t('sl.scenario.impact') },
          ], st.acq.scenario, (val) => { st.acq.scenario = val; })),
          field(t('sl.fps'), numberField(st.acq.fps, (val) => { st.acq.fps = val; }, { step: 10 }), 'fps'),
          field(t('sl.umPerPx'), numberField(st.acq.um_per_px, (val) => { st.acq.um_per_px = val; }, { step: 0.1 }), 'µm/px'),
          field(t('sl.dropVolume'), numberField(st.acq.drop_volume_ul, (val) => { st.acq.drop_volume_ul = val; }, { step: 1 }), 'µL')),
        h('div', { class: 'grid c4 mt' },
          field(t('sl.temperature'), numberField(st.acq.temperature_c, (val) => { st.acq.temperature_c = val; }, { step: 0.5 }), '°C'),
          field(t('sl.humidity'), numberField(st.acq.relative_humidity_pct, (val) => { st.acq.relative_humidity_pct = val; }, { step: 1 }), '%'),
          field(t('sl.tiltRate'), numberField(st.acq.tilt_rate_deg_s, (val) => { st.acq.tilt_rate_deg_s = val; }, { step: 0.1 }), '°/s'),
          field(t('sl.duration'), numberField(st.acq.duration_s, (val) => { st.acq.duration_s = val; }, { step: 0.5 }), 's')),
        h('div', { class: 'row mt' },
          checkField(t('sl.contactLineResolved'), st.acq.contact_line_resolved, (val) => { st.acq.contact_line_resolved = val; }),
          checkField(t('sl.axisymUsed'), st.acq.axisymmetric_fit_used, (val) => { st.acq.axisymmetric_fit_used = val; }),
          h('div', { class: 'spacer' }),
          runButton(null, runAcq)),
        st.acqError ? h('div', { class: 'mt' }, errorBlock(st.acqError)) : null,
        a ? h('div', { class: 'mt' },
          h('div', { class: 'row mb' },
            h('span', { class: 'verdict ' + (a.verdict && a.verdict.indexOf('unusable') === 0 ? 'reject' : a.verdict && a.verdict.indexOf('marginal') === 0 ? 'warn' : 'accept') },
              a.verdict || ''),
            h('span', { class: 'faint small', text: t('sl.blocking') + ': ' + (a.blocking ? a.blocking.length : 0) + ' · ' + t('sl.serious') + ': ' + (a.serious ? a.serious.length : 0) })),
          hazardRows.length
            ? h('div', {},
              dataTable([{ label: t('sl.hazards') }, { label: t('common.severity') }, { label: '' }], hazardRows),
              h('div', { class: 'lib-tag mt', text: t('common.libraryOutput') }))
            : callout('ok', '✓', h('div', { text: 'no hazards raised' }))) : null);

      const f = st.furResult;
      const furCard = card(t('sl.furmidge'), t('sl.furmidgeNote'),
        h('div', { class: 'grid c4' },
          field(t('sl.width'), numberField(st.fur.width_m, (val) => { st.fur.width_m = val; }, { step: 0.0001 }), 'm'),
          field(t('sl.gamma'), numberField(st.fur.gamma_mN_m, (val) => { st.fur.gamma_mN_m = val; }, { step: 1 }), 'mN/m'),
          field(t('sl.thetaA'), numberField(st.fur.theta_a_deg, (val) => { st.fur.theta_a_deg = val; }, { step: 1 }), '°'),
          field(t('sl.thetaR'), numberField(st.fur.theta_r_deg, (val) => { st.fur.theta_r_deg = val; }, { step: 1 }), '°')),
        h('div', { class: 'row mt' },
          h('div', { class: 'spacer' }),
          h('button', { class: 'btn small ghost', onclick: () => runFur().catch((e) => toast(e.message, 'error')) }, t('common.run'))),
        f ? h('div', { class: 'grid c4 mt' },
          stat(t('sl.force'), num(f.force_mN, 4), 'mN', 'hero'),
          stat(t('sl.forceRange'), num(f.force_min_mN, 4) + ' … ' + num(f.force_max_mN, 4), 'mN'),
          stat(t('sl.k'), num(f.k_used, 4), f.k_is_physical ? '✓ physical' : '⚠ outside [0.5, 0.884]'),
          stat('hysteresis', num(f.hysteresis_deg, 2), '°'),
          f.warnings && f.warnings.length
            ? h('div', { style: { gridColumn: '1 / -1' } }, libraryNote(h('ul', { class: 'tight' },
              f.warnings.map((w) => h('li', { text: w })))))
            : null) : null);

      const tr = st.trkResult;
      const trkCard = card(t('sl.tracking'), t('sl.trackingNote'),
        h('div', { class: 'grid c4' },
          field('vx (synthetic)', numberField(st.trk.vx, (val) => { st.trk.vx = val; }, { step: 5 }), 'px/s'),
          field('noise', numberField(st.trk.noise, (val) => { st.trk.noise = val; }, { step: 0.1 }), 'px'),
          field(t('sl.maxSpeed'), numberField(st.trk.max_speed, (val) => { st.trk.max_speed = val; }, { step: 10 }), 'px/s'),
          field(t('sl.posSigma'), numberField(st.trk.position_sigma, (val) => { st.trk.position_sigma = val; }, { step: 0.1 }), 'px')),
        h('div', { class: 'row mt' },
          h('div', { class: 'spacer' }),
          h('button', { class: 'btn small', onclick: () => runTrk().catch((e) => toast(e.message, 'error')) },
            '⚗ ' + t('common.run'))),
        tr ? h('div', { class: 'mt' },
          h('div', { class: 'grid c3' },
            stat(t('sl.tracks'), String(tr.n_tracks), ''),
            tr.tracks.length ? stat(t('sl.speed'), num(tr.tracks[0].velocity.speed, 3), '± ' + num(tr.tracks[0].velocity.std_speed, 3) + ' px/s') : null,
            tr.true_velocity ? stat('true vx', num(tr.true_velocity[0], 3), 'px/s') : null),
          tr.tracks.length
            ? dataTable([
              { label: '#' }, { label: 'n detected' }, { label: 'gaps' },
              { label: 'vx', num: true }, { label: 'vy', num: true }, { label: t('sl.trackAssessment') },
            ], tr.tracks.slice(0, 12).map((item, i) => h('tr', {},
              h('td', { text: String(i) }),
              h('td', { class: 'num', text: String(item.track.n_detected) }),
              h('td', { class: 'num', text: String(item.velocity.n_missing_frames) }),
              h('td', { class: 'num', text: num(item.velocity.vx, 3) }),
              h('td', { class: 'num', text: num(item.velocity.vy, 3) }),
              h('td', { text: (item.assessment && item.assessment.verdict) || '' }))))
            : null) : null);

      return h('div', {},
        pageHead('sl.title', 'sl.lead'),
        acqCard, furCard, trkCard,
        (a || f) ? rawJson({ acquisition: a, furmidge: f, tracking: st.trkResult }) : null);
    }

    return { id: 'sliding', group: 'nav.group.measure', icon: '◔', title: 'nav.sliding', render };
  })();

  /* ================================================================== 7. design */

  const design = (() => {
    const st = {
      n: 25, target: 0.15, gamma: 72.0, deltaRho: 997.0,
      curve: null, decision: null, error: null,
    };

    async function scan() {
      st.error = null;
      try {
        const res = await bridge('scan_shape_parameter', { n: st.n, bo_min: 0.01, bo_max: 0.45 });
        st.curve = res;
        st.decision = await bridge('design_drop', {
          target_ps: st.target, gamma_mN_m: st.gamma, delta_rho: st.deltaRho,
        });
      } catch (err) { st.error = err; }
      D.render();
    }

    function render() {
      const c = st.curve;
      // Mark the P_s peak on the curve: past it, a larger drop *lowers* P_s, and
      // that reversal is the single most counter-intuitive fact on this page.
      // Taken as the largest sampled P_s rather than the nearest to the nominal
      // peak, so the marker is always a real point on the drawn curve.
      let markers = [];
      if (c && c.curve.length) {
        const peakRow = c.curve.reduce((best, row) => (row.ps > best.ps ? row : best), c.curve[0]);
        markers = [{ x: peakRow.bond, y: peakRow.ps, color: 'rgba(255,214,10,0.95)' }];
      }
      const chart = c
        ? lineChart([{
          points: c.curve.map((row) => [row.bond, row.ps]),
          label: 'P_s',
          markers,
        }], { xLabel: 'Bond number', yLabel: 'P_s', height: 260 })
        : null;

      const d = st.decision;
      return h('div', {},
        pageHead('ds.title', 'ds.lead'),
        card(t('ds.scan'), t('ds.samples') + ': ' + st.n,
          h('div', { class: 'grid c4' },
            field(t('ds.samples'), numberField(st.n, (val) => { st.n = val; }, { step: 1, min: 5, max: 200 })),
            field(t('ds.target'), numberField(st.target, (val) => { st.target = val; }, { step: 0.01, min: 0.01 })),
            field(t('ds.gammaInput'), numberField(st.gamma, (val) => { st.gamma = val; }, { step: 1 }), 'mN/m'),
            field(t('ds.rhoInput'), numberField(st.deltaRho, (val) => { st.deltaRho = val; }, { step: 1 }), 'kg/m³')),
          h('div', { class: 'row mt' },
            h('div', { class: 'spacer' }),
            runButton(null, scan)),
          st.error ? h('div', { class: 'mt' }, errorBlock(st.error)) : null,
          c ? h('div', { class: 'mt' },
            chart,
            h('div', { class: 'small faint', text: t('ds.peak') + ' Bo ≈ ' + num(c.peak_bond, 3) },
              h('br'), t('ds.maxValidated') + ' Bo = ' + num(c.max_validated_bond, 3))) : null),
        d ? card(t('common.conclusion'), null,
          d.reachable
            ? h('div', {},
              h('div', { class: 'grid c3' },
                stat(t('ds.bondNeeded'), num(d.bond, 4), ''),
                stat(t('ds.radiusNeeded'), d.required_radius_mm ? num(d.required_radius_mm, 3) : t('common.none'), 'mm'),
                stat(t('ds.target'), num(d.target_ps, 3), '')),
              h('div', { class: 'mt' }, callout('info', 'i', h('div', { class: 'small', text: t('ds.lead') }))))
            : callout('warn', '!', h('div', {},
              h('b', { text: t('ds.notReachable') }),
              h('div', { class: 'small', text: t('ds.notReachableNote') })))) : null,
        c ? rawJson({ curve: c, decision: d }) : null);
    }

    return { id: 'design', group: 'nav.group.design', icon: '◆', title: 'nav.design', render };
  })();

  /* =================================================================== 8. about */

  const about = {
    id: 'about', group: 'nav.group.about', icon: 'ⓘ', title: 'nav.about',
    render() {
      const info = D.state.info || {};
      const pairs = [
        ['drop3d', info.drop3d_version],
        ['Python', info.python],
        ['numpy', info.numpy],
        ['scipy', info.scipy],
        ['Pillow', info.pillow],
        ['pywebview', info.webview],
        ['air density', info.air_density + ' kg/m³'],
        ['gravity', info.gravity + ' m/s²'],
        ['platform', info.python_platform],
      ];
      return h('div', {},
        pageHead('about.title', null),
        card(t('about.env'), null,
          dataTable([{ label: '' }, { label: '' }], pairs.map(([k, v]) =>
            h('tr', {}, h('td', { class: 'mono small', text: k }), h('td', { class: 'mono small', text: String(v) }))))),
        card(t('about.gatesExplained'), null,
          dataTable([{ label: t('common.verdict') }, { label: '' }], [
            h('tr', {}, h('td', {}, verdictPill('accept')), h('td', { class: 'small muted', text: 'every hard gate ran and passed' })),
            h('tr', {}, h('td', {}, verdictPill('accept_with_warning')), h('td', { class: 'small muted', text: 'usable, but something is marginal or a hard gate could not be evaluated' })),
            h('tr', {}, h('td', {}, verdictPill('reject')), h('td', { class: 'small muted', text: 'at least one hard gate failed; no number is reported, only the reason' })),
          ])),
        card(t('about.license'), null,
          callout('warn', '!', h('div', { text: t('about.licenseBody') }))),
        card(t('about.repo'), null,
          h('div', { class: 'small muted', text: t('about.repoBody') })));
    },
  };

  /* ------------------------------------------------------------------ register */

  D.registerPage(overview);
  D.registerPage(pendant);
  D.registerPage(surfaceEnergy);
  D.registerPage(uncertainty);
  D.registerPage(oscillation);
  D.registerPage(sliding);
  D.registerPage(design);
  D.registerPage(about);
  // surfaceEnergy.init() is called by boot() once the bridge answers, so the
  // probe-liquid list is filled before the page is first opened.
})();
