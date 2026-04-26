    const allGenes = [];
    const geneInfo = [];
    const allAntibiotics = [];
    const allImportances = [];

    let aiDetected = {};
    let activeGeneScope = [...allGenes];

    /* â”€â”€ Tabs â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    function switchTab(t) {
      ['manual', 'auto'].forEach(id => {
        const pane = document.getElementById('pane-' + id);
        const tab = document.getElementById('tab-' + id);
        if (pane) pane.classList.toggle('active', id === t);
        if (tab) tab.classList.toggle('active', id === t);
      });
    }

    function normalizeGeneSearchValue(value) {
      return String(value || '').trim().toLowerCase().replace(/[\s_-]+/g, ' ');
    }

    function updateGeneSearchCount(visible, total, queryActive) {
      const count = document.getElementById('gene-search-count');
      if (!count) return;
      count.textContent = queryActive
        ? `${visible} result${visible === 1 ? '' : 's'}`
        : `${total} genes total`;
    }

    function updateGeneSearchResults(query, matches) {
      const results = document.getElementById('gene-search-results');
      if (!results) return;
      if (!query) {
        results.innerHTML = '';
        return;
      }
      if (!matches.length) {
        results.innerHTML = `<strong>No direct matches</strong> for "${query}".`;
        return;
      }
      results.innerHTML = `<strong>Showing matches</strong> for "${query}".`;
    }

    function filterGeneOptions(value = '') {
      const query = normalizeGeneSearchValue(value);
      const searchWrap = document.getElementById('gene-search-wrap');
      const empty = document.getElementById('genes-empty');
      const matches = [];
      let visible = 0;
      let defaultShown = 0;

      allGenes.forEach(gene => {
        const wrap = document.getElementById('wrap_' + gene);
        if (!wrap) return;
        const cb = document.getElementById('cb_' + gene);
        const isChecked = cb && cb.checked;
        const haystack = normalizeGeneSearchValue(
          wrap.dataset.search || `${gene} ${(geneInfo[gene] || {}).antibiotic_class || ''}`
        );
        // Show if: matches search query, OR is selected
        const matchesSearch = query && haystack.includes(query);
        let show = matchesSearch || isChecked;

        // Show all genes by default if there's no search query
        if (!query && !show) {
          show = true;
        }

        wrap.style.display = show ? 'block' : 'none';
        wrap.setAttribute('aria-hidden', show ? 'false' : 'true');
        wrap.classList.toggle('search-hit', Boolean(query) && matchesSearch);
        if (show) visible++;
        if (matchesSearch) matches.push(gene);
      });

      if (searchWrap) searchWrap.classList.toggle('has-value', Boolean(query));
      if (empty) empty.classList.toggle('visible', Boolean(query) && visible === 0);
      updateGeneSearchCount(visible, allGenes.length, Boolean(query));
      updateGeneSearchResults(query, matches);
    }

    function clearGeneSearch() {
      const input = document.getElementById('gene-search');
      if (!input) return;
      input.value = '';
      filterGeneOptions('');
      input.focus();
    }

    /* â”€â”€ Gene sync â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    function syncGene(g, val, skipFilter) {
      const cb = document.getElementById('cb_' + g);
      const hid = document.getElementById('hid_' + g);
      const st = document.getElementById('st_' + g);
      const v = (val !== undefined) ? val : cb.checked;
      if (cb) cb.checked = v;
      if (hid) hid.value = v ? '1' : '0';
      if (st) st.textContent = v ? 'Present' : 'Absent';
      if (!skipFilter) {
        const searchInput = document.getElementById('gene-search');
        filterGeneOptions(searchInput ? searchInput.value : '');
      }
    }
    function selectAll(s) {
      allGenes.forEach(g => syncGene(g, s, true));
      const searchInput = document.getElementById('gene-search');
      filterGeneOptions(searchInput ? searchInput.value : '');
    }

    /* â”€â”€ File load â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    function loadFile(inp) {
      const f = inp.files[0]; if (!f) return;
      const r = new FileReader();
      r.onload = e => {
        document.getElementById('fasta-input').value = e.target.result;
        detectGenes(allGenes);
      };
      r.readAsText(f);
    }
    async function loadExample() {
      const r = await fetch('/example-fasta?genes=blaNDM,gyrA_mut,tetA,mcr1');
      const d = await r.json();
      if (d.fasta) {
        document.getElementById('fasta-input').value = d.fasta;
        detectGenes(allGenes);
      }
    }
    function clearFasta() {
      document.getElementById('fasta-input').value = '';
      document.getElementById('fasta-file').value = '';
      document.getElementById('detect-banner').style.display = 'none';
      document.getElementById('predict-banner').style.display = 'none';
      document.getElementById('seq-pills').innerHTML = '';
      allGenes.forEach(g => { const l = document.getElementById('lbl_' + g); if (l) l.classList.remove('auto-hit'); });
      document.getElementById('ai-badge').style.display = 'none';
      aiDetected = {};
      activeGeneScope = [...allGenes];
      allGenes.forEach(g => syncGene(g, false, true));
      const searchInput = document.getElementById('gene-search');
      filterGeneOptions(searchInput ? searchInput.value : '');
      document.getElementById('result-panel').style.display = 'none';
    }

    /* â”€â”€ Detect genes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    async function detectGenes(geneScope = allGenes) {
      const fasta = document.getElementById('fasta-input').value.trim();
      if (!fasta) return;
      const btn = document.getElementById('detect-btn');
      if (btn) btn.classList.add('loading');
      const scopedGenes = Array.isArray(geneScope) && geneScope.length ? geneScope : allGenes;

      const fd = new FormData();
      fd.append('fasta', fasta);
      fd.append('feature_subset', JSON.stringify(scopedGenes));
      const fileInp = document.getElementById('fasta-file');
      if (fileInp.files[0]) fd.append('fasta_file', fileInp.files[0]);

      try {
        const resp = await fetch('/analyze-genome', { method: 'POST', body: fd });
        const data = await resp.json();
        if (data.error) {
          clearDetectedMarkers();
          showBanner('err', '', 'Analysis Failed', data.error);
          return;
        }

        clearDetectedMarkers(false);
        aiDetected = data.detected;
        activeGeneScope = Array.isArray(data.screened_genes) && data.screened_genes.length
          ? data.screened_genes
          : [...scopedGenes];

        let found = 0;
        allGenes.forEach(g => {
          const present = data.detected[g] === 1;
          syncGene(g, present, true);
          if (present && activeGeneScope.includes(g)) {
            found++;
            const l = document.getElementById('lbl_' + g);
            if (l) l.classList.add('auto-hit');
          }
        });

        const searchInput = document.getElementById('gene-search');
        filterGeneOptions(searchInput ? searchInput.value : '');

        const badge = document.getElementById('ai-badge');
        if (badge) {
          badge.textContent = 'Screened';
          badge.style.display = 'inline-flex';
        }

        const pills = document.getElementById('seq-pills');
        if (pills) {
          pills.innerHTML = `<span class="seq-pill">${(data.sequence_length || 0).toLocaleString()} bp</span>
          <span class="seq-pill">${found}/${activeGeneScope.length} markers flagged</span>
          <span class="seq-pill">Heuristic Screen</span>`;
        }

        const foundNames = activeGeneScope.filter(g => data.detected[g] === 1)
          .map(g => (data.meta || {})[g]?.family?.split(' ')[0] || g).join(' &middot; ');
        showBanner('ok', '',
          `${found} AMR marker${found !== 1 ? 's' : ''} flagged`,
          foundNames || 'No curated AMR markers passed the heuristic screen. Review genes, then click Predict.');
        showPredictBanner('ok', 'Screen Complete', 'Review the flagged markers, adjust them if needed, then click Predict Antibiotic Resistance.');
      } catch (e) {
        showBanner('err', '', 'Network Error', e.message);
      } finally {
        if (btn) btn.classList.remove('loading');
      }
    }

    function clearDetectedMarkers(resetGenes = true) {
      allGenes.forEach(g => {
        const l = document.getElementById('lbl_' + g);
        if (l) l.classList.remove('auto-hit');
        if (resetGenes) syncGene(g, false, true);
      });
      if (resetGenes) {
        const searchInput = document.getElementById('gene-search');
        filterGeneOptions(searchInput ? searchInput.value : '');
      }
      aiDetected = {};
      document.getElementById('ai-badge').style.display = 'none';
      document.getElementById('seq-pills').innerHTML = '';
    }

    function showBanner(type, icon, title, sub) {
      const b = document.getElementById('detect-banner');
      if (!b) return;
      b.className = 'detect-banner' + (type === 'err' ? ' err' : '');
      document.getElementById('b-icon').textContent = icon;
      document.getElementById('b-title').textContent = title;
      document.getElementById('b-sub').textContent = sub;
      b.style.display = 'flex';
    }

    function showPredictBanner(type, title, sub) {
      const banner = document.getElementById('predict-banner');
      if (!banner) return;
      banner.className = 'detect-banner' + (type === 'err' ? ' err' : '');
      document.getElementById('predict-title').textContent = title;
      document.getElementById('predict-sub').textContent = sub;
      banner.style.display = 'flex';
    }

    /* â”€â”€ ML Predict â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    async function runPredict(e, options = {}) {
      if (e) e.preventDefault();
      const saveHistory = options.saveHistory ?? true;

      const btn = document.getElementById('predict-btn');
      if (btn) btn.classList.add('loading');

      try {
        const formData = new FormData();

        // Use aiDetected as base if it exists, then override with current form values
        const currentData = { ...aiDetected };
        allGenes.forEach(g => {
          const hid = document.getElementById('hid_' + g);
          if (hid) currentData[g] = hid.value;
        });

        Object.entries(currentData).forEach(([key, val]) => {
          formData.append(key, val);
        });

        formData.append('feature_subset', JSON.stringify(activeGeneScope));
        formData.append('save_history', saveHistory ? '1' : '0');

        const resp = await fetch('/predict', { method: 'POST', body: formData });
        const data = await resp.json();
        if (!resp.ok || data.error) {
          showPredictBanner('err', 'Prediction Failed', data.error || 'The server could not complete the prediction.');
          return;
        }
        showPredictBanner('ok', 'Prediction Complete', data.saved_to_history ? 'Prediction saved to history.' : 'Prediction generated without saving to history.');
        renderResults(data);
      } catch (err) {
        showPredictBanner('err', 'Network Error', err.message);
      } finally {
        if (btn) btn.classList.remove('loading');
      }
    }

    function getPredictionMeta(prediction) {
      const raw = String(prediction || '').trim();
      const normalized = raw.toLowerCase().replace(/[\s_-]+/g, '');

      if (normalized === 'resistant') {
        return { tone: 'resistant', label: 'Resistant', normalized };
      }
      if (normalized === 'intermediate') {
        return { tone: 'intermediate', label: 'Intermediate', normalized };
      }
      if (['susceptible', 'notresistant', 'nonresistant', 'sensitive'].includes(normalized)) {
        const label = normalized === 'susceptible' ? 'Susceptible' : 'Non-Resistant';
        return { tone: 'susceptible', label, normalized };
      }

      return { tone: 'intermediate', label: raw || 'Unknown', normalized };
    }

    function renderResults(data) {
      const panel = document.getElementById('result-panel');
      const resultGeneScope = Array.isArray(data.screened_genes) && data.screened_genes.length
        ? data.screened_genes
        : activeGeneScope;
      activeGeneScope = [...resultGeneScope];
      const results = data.results || [];
      const normalizedResults = results.map(r => ({ ...r, meta: getPredictionMeta(r.prediction) }));
      const nRes = normalizedResults.filter(r => r.meta.tone === 'resistant').length;
      const nInt = normalizedResults.filter(r => r.meta.tone === 'intermediate').length;
      const nSus = normalizedResults.filter(r => r.meta.tone === 'susceptible').length;
      const avgConf = results.length ? Math.round(results.reduce((a, r) => a + r.confidence, 0) / results.length) : 0;
      const nonResLabel = normalizedResults.some(r => r.meta.normalized === 'notresistant') ? 'Non-Resistant' : 'Susceptible';

      /* Summary boxes */
      document.getElementById('result-summary').innerHTML = `
    <div class="sum-box ${nRes > 0 ? 'danger' : 'ok'}">
      <div class="sum-val">${nRes}</div>
      <div class="sum-lbl">Resistant</div>
    </div>
    <div class="sum-box ${nInt > 0 ? 'warning' : 'ok'}">
      <div class="sum-val">${nInt}</div>
      <div class="sum-lbl">Intermediate</div>
    </div>
    <div class="sum-box ok">
      <div class="sum-val">${nSus}</div>
      <div class="sum-lbl">${nonResLabel}</div>
    </div>
    <div class="sum-box info">
      <div class="sum-val">${avgConf}%</div>
      <div class="sum-lbl">Avg Confidence</div>
    </div>`;

      /* Antibiotic table */
      const tbody = document.getElementById('ab-tbody');
      tbody.innerHTML = '';
      normalizedResults.forEach(r => {
        const key = r.meta.tone;
        const tr = document.createElement('tr');
        tr.className = 'ab-row';
        tr.innerHTML = `
      <td><span class="ab-name">${r.antibiotic}</span></td>
      <td><span class="ab-badge ${key}">${r.meta.label}</span></td>
      <td><span class="ab-conf">${r.confidence}%</span></td>
      <td class="ab-bar-wrap">
        <div class="ab-bar-track">
          <div class="ab-bar-fill fill-${key}" data-pct="${r.confidence}"></div>
        </div>
      </td>`;
        tbody.appendChild(tr);
      });
      requestAnimationFrame(() => {
        tbody.querySelectorAll('.ab-bar-fill').forEach(el => { el.style.width = el.dataset.pct + '%'; });
      });

      /* Feature importance */
      const scopedGenes = new Set(resultGeneScope);
      const scopedImportances = Object.entries(allImportances).filter(([gene]) => scopedGenes.has(gene));
      const visibleImportances = (scopedImportances.length ? scopedImportances : Object.entries(allImportances)).slice(0, 8);
      const maxImp = Math.max(...visibleImportances.map(([, imp]) => imp), 1);
      const fi = document.getElementById('fi-section');
      fi.innerHTML = '';
      visibleImportances.forEach(([gene, imp]) => {
        const pct = Math.round(imp * 100);
        const pctOfMax = Math.round((imp / maxImp) * 100);
        const row = document.createElement('div');
        row.className = 'fi-row';
        row.innerHTML = `
      <span class="fi-name">${gene}</span>
      <div class="fi-track"><div class="fi-fill" data-pct="${pctOfMax}"></div></div>
      <span class="fi-pct">${pct}%</span>`;
        fi.appendChild(row);
      });
      requestAnimationFrame(() => {
        fi.querySelectorAll('.fi-fill').forEach(el => { el.style.width = el.dataset.pct + '%'; });
      });

      /* Gene chips â€“ only show genes that are present/detected */
      const chips = document.getElementById('gene-chips');
      chips.innerHTML = '';
      resultGeneScope.forEach(g => {
        const cb = document.getElementById('cb_' + g);
        const present = cb && cb.checked;
        if (!present) return;
        const isAI = aiDetected[g] === 1;
        const chip = document.createElement('span');
        chip.className = `gene-chip present${isAI ? ' ai' : ''}`;
        chip.innerHTML = (isAI ? '' : ' &#10003; ') + g;
        chip.title = geneInfo[g]?.antibiotic_class || '';
        chips.appendChild(chip);
      });

      panel.style.display = 'block';
      panel.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    /* â”€â”€ Smart Input Module Logic â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
    let currentSmartSequence = { bacteria: '', disease: '' };
    let currentSmartScreenedGenes = { bacteria: [], disease: [] };
    let smartDropdownOptions = { bacteria: [], disease: [] };
    let smartDropdownState = {
      bacteria: { open: false, highlighted: -1, filtered: [], closeTimer: null },
      disease: { open: false, highlighted: -1, filtered: [], closeTimer: null }
    };

    function normalizeSmartOption(value = '') {
      return String(value).trim().toLowerCase();
    }

    function getSmartDropdownElements(mode) {
      return {
        combo: document.getElementById(`${mode}-combobox`),
        input: document.getElementById(`${mode}-input`),
        list: document.getElementById(`${mode}-dropdown-list`),
        meta: document.getElementById(`${mode}-dropdown-meta`),
      };
    }

    function filterSmartDropdownOptions(mode, query = '') {
      const options = smartDropdownOptions[mode] || [];
      const normalizedQuery = normalizeSmartOption(query);
      if (!normalizedQuery) return [...options];

      const starts = options.filter(option => normalizeSmartOption(option).startsWith(normalizedQuery));
      const contains = options.filter(option =>
        !normalizeSmartOption(option).startsWith(normalizedQuery) &&
        normalizeSmartOption(option).includes(normalizedQuery)
      );
      return [...starts, ...contains];
    }

    function renderSmartDropdown(mode, query = '') {
      const state = smartDropdownState[mode];
      const { list, meta } = getSmartDropdownElements(mode);
      if (!list || !meta) return;

      const filtered = filterSmartDropdownOptions(mode, query);
      state.filtered = filtered;
      if (filtered.length === 0) {
        state.highlighted = -1;
      } else if (state.highlighted >= filtered.length || state.highlighted < 0) {
        state.highlighted = 0;
      }

      meta.textContent = `${filtered.length} option${filtered.length === 1 ? '' : 's'}`;
      list.innerHTML = '';

      if (!filtered.length) {
        list.innerHTML = `
          <div class="smart-dropdown-empty">
            <strong>No mapped matches</strong>
            Keep typing, or press Enter to try the exact value you entered.
          </div>`;
        return;
      }

      filtered.forEach((option, index) => {
        const item = document.createElement('button');
        item.type = 'button';
        item.className = `smart-dropdown-item${index === state.highlighted ? ' active' : ''}`;
        item.textContent = option;
        item.onmousedown = evt => {
          evt.preventDefault();
          selectSmartDropdownOption(mode, option);
        };
        list.appendChild(item);
      });
    }

    function openSmartDropdown(mode) {
      const state = smartDropdownState[mode];
      const { combo, input } = getSmartDropdownElements(mode);
      if (!combo || !input) return;
      clearTimeout(state.closeTimer);
      state.open = true;
      combo.classList.add('open');
      renderSmartDropdown(mode, input.value);
    }

    function closeSmartDropdown(mode) {
      const state = smartDropdownState[mode];
      const { combo } = getSmartDropdownElements(mode);
      if (!combo) return;
      state.open = false;
      combo.classList.remove('open');
    }

    function scheduleSmartDropdownClose(mode) {
      const state = smartDropdownState[mode];
      clearTimeout(state.closeTimer);
      state.closeTimer = setTimeout(() => closeSmartDropdown(mode), 140);
    }

    function toggleSmartDropdown(mode) {
      const state = smartDropdownState[mode];
      const { input } = getSmartDropdownElements(mode);
      if (state.open) {
        closeSmartDropdown(mode);
      } else {
        if (input) input.focus();
        openSmartDropdown(mode);
      }
    }

    async function commitSmartInput(mode, value) {
      const trimmed = String(value || '').trim();
      if (!trimmed) return;
      if (mode === 'bacteria') {
        await fetchSequenceForBacteria(trimmed, mode);
      } else {
        await fetchBacteriaForDisease(trimmed);
      }
    }

    function handleSmartInput(mode, value) {
      const state = smartDropdownState[mode];
      state.highlighted = 0;
      openSmartDropdown(mode);
      renderSmartDropdown(mode, value);
    }

    async function selectSmartDropdownOption(mode, value) {
      const { input } = getSmartDropdownElements(mode);
      if (input) input.value = value;
      closeSmartDropdown(mode);
      await commitSmartInput(mode, value);
    }

    async function handleSmartInputKeydown(mode, event) {
      const state = smartDropdownState[mode];
      if (!state.open && ['ArrowDown', 'ArrowUp'].includes(event.key)) {
        openSmartDropdown(mode);
      }

      if (event.key === 'ArrowDown') {
        event.preventDefault();
        if (state.filtered.length) {
          state.highlighted = Math.min(state.highlighted + 1, state.filtered.length - 1);
          renderSmartDropdown(mode, event.target.value);
        }
        return;
      }

      if (event.key === 'ArrowUp') {
        event.preventDefault();
        if (state.filtered.length) {
          state.highlighted = Math.max(state.highlighted - 1, 0);
          renderSmartDropdown(mode, event.target.value);
        }
        return;
      }

      if (event.key === 'Enter') {
        event.preventDefault();
        const selected = state.filtered[state.highlighted];
        await selectSmartDropdownOption(mode, selected || event.target.value);
        return;
      }

      if (event.key === 'Escape') {
        closeSmartDropdown(mode);
      }
    }

    function resetSmartSequence(mode, message = 'Genomic sequence will appear here...') {
      currentSmartSequence[mode] = '';
      currentSmartScreenedGenes[mode] = [];
      const display = document.getElementById(`${mode}-seq-display`);
      if (display) display.textContent = message;
    }

    function setSmartMode(mode) {
      document.querySelectorAll('.smart-tab-btn').forEach(b => b.classList.remove('active'));
      document.querySelectorAll('.smart-panel').forEach(p => p.classList.remove('active'));

      const btn = document.querySelector(`.smart-tab-btn.${mode}-tab`);
      const panel = document.getElementById(`panel-${mode}`);
      if (btn) btn.classList.add('active');
      if (panel) panel.classList.add('active');
    }

    async function initSmartInput() {
      try {
        const [diseasesRes, bacteriaRes] = await Promise.all([
          fetch('/api/get-disease-suggestions'),
          fetch('/api/get-bacteria-suggestions')
        ]);
        const diseases = await diseasesRes.json();
        const bacteria = await bacteriaRes.json();

        smartDropdownOptions.bacteria = Array.isArray(bacteria) ? bacteria : [];
        smartDropdownOptions.disease = Array.isArray(diseases) ? diseases : [];
        renderSmartDropdown('bacteria');
        renderSmartDropdown('disease');
      } catch (e) { console.error("Error init smart input:", e); }
    }

    async function fetchSequenceForBacteria(bacteria, mode = 'bacteria') {
      if (!bacteria) return;
      const overlayId = `${mode}-loading`;
      const displayId = `${mode}-seq-display`;

      document.getElementById(overlayId).classList.add('visible');
      try {
        const res = await fetch(`/api/get-sequence?bacteria=${encodeURIComponent(bacteria)}`);
        const data = await res.json();
        if (res.ok && data.sequence) {
          document.getElementById(displayId).textContent = data.sequence;
          currentSmartSequence[mode] = data.sequence;
          currentSmartScreenedGenes[mode] = Array.isArray(data.screened_genes) ? data.screened_genes : [];
        } else {
          document.getElementById(displayId).textContent = `Error: ${data.error || 'Sequence not found'}`;
          currentSmartSequence[mode] = '';
          currentSmartScreenedGenes[mode] = [];
        }
      } catch (e) {
        document.getElementById(displayId).textContent = `Network error trying to fetch sequence.`;
        currentSmartSequence[mode] = '';
        currentSmartScreenedGenes[mode] = [];
      } finally {
        document.getElementById(overlayId).classList.remove('visible');
      }
    }

    async function fetchBacteriaForDisease(disease) {
      if (!disease) return;
      const suggContainer = document.getElementById('disease-bacteria-suggestions');
      const pillsDiv = document.getElementById('disease-suggestion-pills');
      pillsDiv.innerHTML = '';
      resetSmartSequence('disease', 'Select a suggested bacterium to load its genomic sequence.');

      try {
        const res = await fetch(`/api/get-bacteria-from-disease?disease=${encodeURIComponent(disease)}`);
        const bacteriaList = await res.json();
        if (res.ok && bacteriaList.length > 0) {
          suggContainer.classList.add('visible');
          let firstPill = null;
          bacteriaList.forEach(b => {
            const pill = document.createElement('div');
            pill.className = 'suggestion-pill';
            pill.textContent = b;
            pill.onclick = async () => {
              pillsDiv.querySelectorAll('.suggestion-pill').forEach(node => node.classList.remove('active'));
              pill.classList.add('active');
              fetchSequenceForBacteria(b, 'disease');
            };
            pillsDiv.appendChild(pill);
            if (!firstPill) firstPill = pill;
          });
          if (firstPill) firstPill.click();
        } else {
          suggContainer.classList.remove('visible');
          resetSmartSequence('disease', 'No mapped bacteria found for this disease.');
        }
      } catch (e) { console.error(e); }
    }

    function copySmartSequence(mode) {
      const seq = currentSmartSequence[mode];
      if (seq) {
        navigator.clipboard.writeText(seq);
        alert('Sequence copied to clipboard!');
      } else {
        alert('No sequence to copy.');
      }
    }

    function downloadSmartSequence(mode) {
      const seq = currentSmartSequence[mode];
      if (!seq) return alert('No sequence to download.');
      const blob = new Blob([seq], { type: 'text/plain' });
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `genome_sequence_${mode}.txt`;
      a.click();
    }

    async function sendToPredictor(mode) {
      const seq = currentSmartSequence[mode];
      if (!seq) return alert('Retrieve a sequence first.');
      clearFasta();
      document.getElementById('fasta-input').value = seq;
      switchTab('auto');
      document.getElementById('pane-auto').scrollIntoView({ behavior: 'smooth', block: 'start' });
      await detectGenes(currentSmartScreenedGenes[mode]);
    }

    // Initialize on load (Main parts)
    initSmartInput();
    filterGeneOptions('');

