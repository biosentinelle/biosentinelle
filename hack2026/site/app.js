const INSERTIONS_URLS = [
  "/output/insertions.json",
  "../output/insertions.json",
];
const CHR_IV = "ref|NC_001136|";
const CHR_IV_LENGTH = 1531933;
const COLORS = {
  Kan_gene: "#2F9E73",
  Promotor_VCF: "#4E79A7",
};

const state = {
  insertions: [],
  binSize: 500,
  viewWidth: 200000,
  visible: {
    Kan_gene: true,
    Promotor_VCF: true,
  },
  isLoaded: false,
  loadPromise: null,
};

const statusEl = document.getElementById("status");
const modal = document.getElementById("modal");
const openChrIV = document.getElementById("openChrIV");
const closeModal = document.getElementById("closeModal");
const closeBackdrop = document.getElementById("closeBackdrop");
const binSlider = document.getElementById("binSlider");
const binValue = document.getElementById("binValue");
const viewSlider = document.getElementById("viewSlider");
const viewValue = document.getElementById("viewValue");
const plotTarget = document.getElementById("chrIVPlot");
const modCenter = document.getElementById("modCenter");

function formatInt(value) {
  return Math.round(value).toLocaleString("en-US");
}

function updateSummary() {
  const totalSignals = document.getElementById("totalSignals");
  const kanReads = document.getElementById("kanReads");
  const vcfReads = document.getElementById("vcfReads");
  if (!totalSignals || !kanReads || !vcfReads) return;

  const totals = state.insertions.reduce(
    (acc, row) => {
      const count = row.junction_positions?.length ?? 0;
      acc.signals += 1;
      if (row.element === "Kan_gene") acc.kan += count;
      if (row.element === "Promotor_VCF") acc.vcf += count;
      return acc;
    },
    { signals: 0, kan: 0, vcf: 0 }
  );

  totalSignals.textContent = formatInt(totals.signals);
  kanReads.textContent = formatInt(totals.kan);
  vcfReads.textContent = formatInt(totals.vcf);
}

function positionsFor(element) {
  return state.insertions
    .filter((row) => row.chromosome === CHR_IV && row.element === element)
    .flatMap((row) => row.junction_positions || []);
}

function estimateGapCenter() {
  const kanRows = state.insertions
    .filter((row) => row.chromosome === CHR_IV && row.element === "Kan_gene")
    .sort((a, b) => (b.junction_positions?.length || 0) - (a.junction_positions?.length || 0));

  const best = kanRows[0];
  if (!best) return Math.floor(CHR_IV_LENGTH / 2);

  if (best.left_breakpoint && best.right_breakpoint) {
    return Math.round((best.left_breakpoint + best.right_breakpoint) / 2);
  }

  if (best.insertion_position) return best.insertion_position;
  return Math.floor(CHR_IV_LENGTH / 2);
}

function centeredRegion(center, width) {
  let start = Math.round(center - width / 2);
  let end = start + width - 1;

  if (start < 1) {
    start = 1;
    end = width;
  }

  if (end > CHR_IV_LENGTH) {
    end = CHR_IV_LENGTH;
    start = Math.max(1, end - width + 1);
  }

  return { start, end };
}

function binPositions(positions, binSize, regionStart, regionEnd) {
  const regionLength = regionEnd - regionStart + 1;
  const binCount = Math.ceil(regionLength / binSize);
  const bins = new Array(binCount).fill(0);

  for (const pos of positions) {
    if (pos >= regionStart && pos <= regionEnd) {
      bins[Math.floor((pos - regionStart) / binSize)] += 1;
    }
  }

  return bins;
}

function drawChrIVPlot() {
  if (!state.insertions.length) {
    plotTarget.innerHTML = '<div class="empty-plot">Chargement des donnees JSON...</div>';
    return;
  }

  const binSize = state.binSize;
  const gapCenter = estimateGapCenter();
  if (modCenter) {
    modCenter.textContent = `Centre modification: chrIV:${formatInt(gapCenter)} bp`;
  }
  const { start: regionStart, end: regionEnd } = centeredRegion(gapCenter, state.viewWidth);
  const regionLength = regionEnd - regionStart + 1;
  const kan = binPositions(positionsFor("Kan_gene"), binSize, regionStart, regionEnd);
  const vcf = binPositions(positionsFor("Promotor_VCF"), binSize, regionStart, regionEnd);
  const visibleBins = [];
  if (state.visible.Kan_gene) visibleBins.push(...kan);
  if (state.visible.Promotor_VCF) visibleBins.push(...vcf);
  const maxBin = Math.max(1, ...visibleBins);

  const width = 1260;
  const height = 430;
  const left = 76;
  const right = 1228;
  const top = 52;
  const bottom = 360;
  const plotW = right - left;
  const plotH = bottom - top;
  const barW = plotW / kan.length;

  const bars = [];
  const drawSeries = (values, color, opacity) => {
    values.forEach((count, idx) => {
      if (!count) return;
      const x = left + idx * barW;
      const h = Math.max(2.4, (count / maxBin) * plotH);
      const visibleBarW = Math.max(1.6, barW - 0.25);
      bars.push(
        `<rect x="${x.toFixed(2)}" y="${(bottom - h).toFixed(2)}" width="${visibleBarW.toFixed(2)}" height="${h.toFixed(2)}" fill="${color}" opacity="${opacity}"/>`
      );
    });
  };

  if (state.visible.Promotor_VCF) drawSeries(vcf, COLORS.Promotor_VCF, 0.66);
  if (state.visible.Kan_gene) drawSeries(kan, COLORS.Kan_gene, 0.72);

  const legendOpacity = (element) => state.visible[element] ? 1 : 0.32;
  const legendLabel = (element) => state.visible[element] ? "masquer" : "afficher";

  const ticks = Array.from({ length: 7 }, (_, idx) => {
    const x = left + (idx / 6) * plotW;
    const pos = regionStart + (idx / 6) * regionLength;
    return `
      <line x1="${x}" y1="${bottom}" x2="${x}" y2="${bottom + 6}" stroke="#536176"/>
      <text x="${x}" y="${bottom + 24}" text-anchor="middle" font-size="12" fill="#536176">${formatInt(pos)}</text>
    `;
  }).join("");

  const yTicks = Array.from({ length: 5 }, (_, idx) => {
    const value = Math.round((idx / 4) * maxBin);
    const y = bottom - (idx / 4) * plotH;
    return `
      <line x1="${left}" y1="${y}" x2="${right}" y2="${y}" stroke="#DDE5EF"/>
      <text x="${left - 12}" y="${y + 4}" text-anchor="end" font-size="12" fill="#6B778C">${value}</text>
    `;
  }).join("");

  plotTarget.innerHTML = `
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Histogramme chrIV">
      <rect width="100%" height="100%" fill="#ffffff"/>
      <text x="0" y="21" font-size="18" font-weight="700" fill="#172033">chrIV - zoom du pic Kan_gene + Promotor_VCF</text>
      <text x="0" y="42" font-size="13" fill="#6B778C">centre du vide=${formatInt(gapCenter)} bp | region=${formatInt(regionStart)}-${formatInt(regionEnd)} bp | bin=${formatInt(binSize)} bp | max bin=${formatInt(maxBin)}</text>
      <line x1="${left + ((gapCenter - regionStart) / regionLength) * plotW}" y1="${top}" x2="${left + ((gapCenter - regionStart) / regionLength) * plotW}" y2="${bottom}" stroke="#172033" stroke-width="2" stroke-dasharray="5 5" opacity="0.68"/>
      <rect x="${left}" y="${top}" width="${plotW}" height="${plotH}" fill="#F6F8FB" stroke="#DDE5EF"/>
      ${yTicks}
      ${bars.join("")}
      <line x1="${left}" y1="${bottom}" x2="${right}" y2="${bottom}" stroke="#9EADBF"/>
      ${ticks}
      <text x="${left}" y="${height - 14}" font-size="13" fill="#6B778C">${formatInt(regionStart)}</text>
      <text x="${right}" y="${height - 14}" text-anchor="end" font-size="13" fill="#6B778C">position sur chrIV</text>
      <g data-toggle="Kan_gene" style="cursor:pointer" opacity="${legendOpacity("Kan_gene")}">
        <title>Cliquer pour ${legendLabel("Kan_gene")} Kan_gene</title>
        <rect x="700" y="8" width="104" height="27" rx="6" fill="#F6F8FB" stroke="#DDE5EF"/>
        <rect x="713" y="16" width="14" height="14" rx="3" fill="${COLORS.Kan_gene}" opacity="0.72"/>
        <text x="735" y="28" font-size="13" fill="#334155">Kan_gene</text>
      </g>
      <g data-toggle="Promotor_VCF" style="cursor:pointer" opacity="${legendOpacity("Promotor_VCF")}">
        <title>Cliquer pour ${legendLabel("Promotor_VCF")} Promotor_VCF</title>
        <rect x="816" y="8" width="138" height="27" rx="6" fill="#F6F8FB" stroke="#DDE5EF"/>
        <rect x="829" y="16" width="14" height="14" rx="3" fill="${COLORS.Promotor_VCF}" opacity="0.66"/>
        <text x="851" y="28" font-size="13" fill="#334155">Promotor_VCF</text>
      </g>
    </svg>
  `;
}

function loadInsertions() {
  if (state.isLoaded) return Promise.resolve(state.insertions);
  if (state.loadPromise) return state.loadPromise;

  if (Array.isArray(window.BIOSENTINEL_INSERTIONS)) {
    state.insertions = window.BIOSENTINEL_INSERTIONS;
    state.isLoaded = true;
    if (statusEl) statusEl.textContent = "Prêt";
    updateSummary();
    return Promise.resolve(state.insertions);
  }

  state.loadPromise = INSERTIONS_URLS.reduce(
    (promise, url) => promise.catch(() => fetch(url)),
    Promise.reject()
  )
    .then((response) => {
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      return response.json();
    })
    .then((data) => {
      state.insertions = data;
      state.isLoaded = true;
      if (statusEl) statusEl.textContent = "Prêt";
      updateSummary();
      return data;
    })
    .catch((error) => {
      state.loadPromise = null;
      if (statusEl) statusEl.textContent = "Erreur JSON";
      throw error;
    });

  return state.loadPromise;
}

function openModal() {
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  if (state.isLoaded) {
    drawChrIVPlot();
    return;
  }

  plotTarget.innerHTML = '<div class="empty-plot">Chargement des donnees JSON...</div>';
  loadInsertions().then(drawChrIVPlot).catch((error) => {
    plotTarget.innerHTML = `
      <div class="empty-plot">
        Impossible de charger <code>output/insertions.json</code>.<br>
        Ouvre la page via <code>http://127.0.0.1:8765/site/index.html</code>, pas en <code>file://</code>.
      </div>
    `;
    console.error(error);
  });
}

function hideModal() {
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
}

binSlider.addEventListener("input", () => {
  state.binSize = Number(binSlider.value);
  binValue.textContent = `${formatInt(state.binSize)} bp`;
  if (state.isLoaded) drawChrIVPlot();
});

viewSlider.addEventListener("input", () => {
  state.viewWidth = Number(viewSlider.value);
  viewValue.textContent = `${formatInt(state.viewWidth)} bp`;
  if (state.isLoaded) drawChrIVPlot();
});

plotTarget.addEventListener("click", (event) => {
  const toggle = event.target.closest("[data-toggle]");
  if (!toggle) return;

  const element = toggle.dataset.toggle;
  state.visible[element] = !state.visible[element];
  if (!state.visible.Kan_gene && !state.visible.Promotor_VCF) {
    state.visible[element] = true;
  }
  drawChrIVPlot();
});

openChrIV.addEventListener("click", openModal);
closeModal.addEventListener("click", hideModal);
closeBackdrop.addEventListener("click", hideModal);

document.addEventListener("keydown", (event) => {
  if (event.key === "Escape") hideModal();
});

loadInsertions()
  .then(() => {
    if (modal.classList.contains("is-open")) drawChrIVPlot();
  })
  .catch((error) => {
    console.error(error);
  });
