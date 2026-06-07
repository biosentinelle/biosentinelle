const DEFAULT_INSERTIONS_URLS = ["insertions.json", "./insertions.json"];
const COLORS = {
  Kan_gene: "#2F9E73",
  Promotor_VCF: "#4E79A7",
};
const CHROMS = [
  ["ref|NC_001133|", "chrI", 230218],
  ["ref|NC_001134|", "chrII", 813184],
  ["ref|NC_001135|", "chrIII", 316620],
  ["ref|NC_001136|", "chrIV", 1531933],
  ["ref|NC_001137|", "chrV", 576874],
  ["ref|NC_001138|", "chrVI", 270161],
  ["ref|NC_001139|", "chrVII", 1090940],
  ["ref|NC_001140|", "chrVIII", 562643],
  ["ref|NC_001141|", "chrIX", 439888],
  ["ref|NC_001142|", "chrX", 745751],
  ["ref|NC_001143|", "chrXI", 666816],
  ["ref|NC_001144|", "chrXII", 1078177],
  ["ref|NC_001145|", "chrXIII", 924431],
  ["ref|NC_001146|", "chrXIV", 784333],
  ["ref|NC_001147|", "chrXV", 1091291],
  ["ref|NC_001148|", "chrXVI", 948066],
];

const state = {
  insertions: [],
  selectedChrom: "ref|NC_001136|",
  isLoaded: false,
};

const overview = document.getElementById("genomeOverview");
const metadata = document.getElementById("metadata");
const modal = document.getElementById("modal");
const modalChrom = document.getElementById("modalChrom");
const closeModal = document.getElementById("closeModal");
const closeBackdrop = document.getElementById("closeBackdrop");
const plotTarget = document.getElementById("chromPlot");
const modCenter = document.getElementById("modCenter");

function formatInt(value) {
  return Math.round(value).toLocaleString("en-US");
}

function chromInfo(chrom) {
  const known = CHROMS.find((row) => row[0] === chrom);
  if (known) return { id: known[0], label: known[1], length: known[2] };
  const positions = state.insertions
    .filter((row) => row.chromosome === chrom)
    .flatMap((row) => row.junction_positions || []);
  return { id: chrom, label: chrom, length: Math.max(1, ...positions) };
}

function insertionUrls() {
  const params = new URLSearchParams(window.location.search);
  const dataUrl = params.get("data");
  return dataUrl ? [dataUrl, ...DEFAULT_INSERTIONS_URLS] : DEFAULT_INSERTIONS_URLS;
}

function fetchJson(url) {
  return fetch(url, { cache: "no-store" })
    .then((response) => {
      if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
      return response.json();
    })
    .then((data) => {
      if (!Array.isArray(data)) throw new Error(`${url}: JSON root is not an array`);
      return data;
    });
}

function loadInsertions() {
  return insertionUrls().reduce(
    (promise, url) => promise.catch(() => fetchJson(url)),
    Promise.reject()
  ).catch((error) => {
    if (Array.isArray(window.BIOSENTINEL_INSERTIONS)) return window.BIOSENTINEL_INSERTIONS;
    throw error;
  });
}

function rowsForChrom(chrom) {
  return state.insertions.filter((row) => row.chromosome === chrom);
}

function positionsFor(chrom, element) {
  return rowsForChrom(chrom)
    .filter((row) => row.element === element)
    .flatMap((row) => row.junction_positions || []);
}

function bestCenter(chrom) {
  const rows = rowsForChrom(chrom)
    .slice()
    .sort((a, b) => (b.junction_positions?.length || 0) - (a.junction_positions?.length || 0));
  const best = rows[0];
  if (!best) return Math.floor(chromInfo(chrom).length / 2);
  if (best.left_breakpoint && best.right_breakpoint) {
    return Math.round((best.left_breakpoint + best.right_breakpoint) / 2);
  }
  return best.insertion_position || Math.floor(chromInfo(chrom).length / 2);
}

function renderOverview() {
  if (!state.insertions.length) {
    overview.innerHTML = '<div class="empty-plot">Aucun signal insertion detecte.</div>';
    metadata.textContent = "0 insertion signals";
    return;
  }

  const chroms = Array.from(new Set([
    ...CHROMS.map((row) => row[0]),
    ...state.insertions.map((row) => row.chromosome),
  ]));
  const rowH = 34;
  const top = 42;
  const left = 90;
  const right = 1180;
  const height = top + chroms.length * rowH + 48;
  const maxLen = Math.max(...chroms.map((chrom) => chromInfo(chrom).length));
  const totalReads = state.insertions.reduce((sum, row) => sum + (row.junction_positions?.length || 0), 0);
  metadata.textContent = `${state.insertions.length} insertion signal(s) | ${formatInt(totalReads)} junction read(s)`;

  const rows = chroms.map((chrom, idx) => {
    const info = chromInfo(chrom);
    const y = top + idx * rowH;
    const w = (info.length / maxLen) * (right - left);
    const markers = rowsForChrom(chrom).map((row) => {
      const pos = row.insertion_position || bestCenter(chrom);
      const x = left + (pos / info.length) * w;
      const color = COLORS[row.element] || "#172033";
      const size = Math.min(14, 4 + Math.sqrt(row.junction_positions?.length || 1));
      return `<circle cx="${x.toFixed(1)}" cy="${y}" r="${size.toFixed(1)}" fill="${color}" opacity="0.76"><title>${info.label} ${row.element} ${formatInt(pos)} bp</title></circle>`;
    }).join("");
    return `
      <g data-chrom="${chrom}" tabindex="0" role="button" style="cursor:pointer">
        <text x="16" y="${y + 5}" font-size="13" font-weight="700" fill="#172033">${info.label}</text>
        <rect x="${left}" y="${y - 8}" width="${w.toFixed(1)}" height="16" rx="4" fill="#F6F8FB" stroke="#DDE5EF"/>
        ${markers}
      </g>
    `;
  }).join("");

  overview.innerHTML = `
    <svg viewBox="0 0 1240 ${height}" role="img" aria-label="Insertion overview">
      <rect width="100%" height="100%" fill="#ffffff"/>
      <text x="16" y="22" font-size="18" font-weight="700" fill="#172033">Insertion overview</text>
      <g transform="translate(930 7)">
        <rect x="0" y="0" width="16" height="16" rx="3" fill="${COLORS.Kan_gene}"/>
        <text x="24" y="13" font-size="13" fill="#334155">Kan_gene</text>
        <rect x="118" y="0" width="16" height="16" rx="3" fill="${COLORS.Promotor_VCF}"/>
        <text x="142" y="13" font-size="13" fill="#334155">Promotor_VCF</text>
      </g>
      ${rows}
      <text x="${left}" y="${height - 18}" font-size="12" fill="#6B778C">0 bp</text>
      <text x="${right}" y="${height - 18}" text-anchor="end" font-size="12" fill="#6B778C">scaled chromosome length</text>
    </svg>
  `;
}

function centeredRegion(center, width, chromLength) {
  let start = Math.round(center - width / 2);
  let end = start + width - 1;
  if (start < 1) {
    start = 1;
    end = Math.min(width, chromLength);
  }
  if (end > chromLength) {
    end = chromLength;
    start = Math.max(1, end - width + 1);
  }
  return { start, end };
}

function binPositions(positions, binSize, regionStart, regionEnd) {
  const binCount = Math.ceil((regionEnd - regionStart + 1) / binSize);
  const bins = new Array(binCount).fill(0);
  for (const pos of positions) {
    if (pos >= regionStart && pos <= regionEnd) {
      bins[Math.floor((pos - regionStart) / binSize)] += 1;
    }
  }
  return bins;
}

function drawChromPlot() {
  const info = chromInfo(state.selectedChrom);
  const center = bestCenter(info.id);
  const region = centeredRegion(center, Math.min(200000, info.length), info.length);
  const binSize = 500;
  const kan = binPositions(positionsFor(info.id, "Kan_gene"), binSize, region.start, region.end);
  const vcf = binPositions(positionsFor(info.id, "Promotor_VCF"), binSize, region.start, region.end);
  const maxBin = Math.max(1, ...kan, ...vcf);
  const width = 1260;
  const height = 430;
  const left = 76;
  const right = 1228;
  const top = 52;
  const bottom = 360;
  const plotW = right - left;
  const plotH = bottom - top;
  const barW = plotW / kan.length;

  modalChrom.textContent = info.label;
  modCenter.textContent = `Centre modification: ${info.label}:${formatInt(center)} bp`;

  const bars = [];
  const drawSeries = (values, color, opacity) => {
    values.forEach((count, idx) => {
      if (!count) return;
      const x = left + idx * barW;
      const h = Math.max(2.4, (count / maxBin) * plotH);
      bars.push(`<rect x="${x.toFixed(2)}" y="${(bottom - h).toFixed(2)}" width="${Math.max(1.6, barW - 0.25).toFixed(2)}" height="${h.toFixed(2)}" fill="${color}" opacity="${opacity}"/>`);
    });
  };
  drawSeries(vcf, COLORS.Promotor_VCF, 0.66);
  drawSeries(kan, COLORS.Kan_gene, 0.72);

  const ticks = Array.from({ length: 7 }, (_, idx) => {
    const x = left + (idx / 6) * plotW;
    const pos = region.start + (idx / 6) * (region.end - region.start + 1);
    return `<line x1="${x}" y1="${bottom}" x2="${x}" y2="${bottom + 6}" stroke="#536176"/><text x="${x}" y="${bottom + 24}" text-anchor="middle" font-size="12" fill="#536176">${formatInt(pos)}</text>`;
  }).join("");

  plotTarget.innerHTML = `
    <svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Histogramme ${info.label}">
      <rect width="100%" height="100%" fill="#ffffff"/>
      <text x="0" y="21" font-size="18" font-weight="700" fill="#172033">${info.label} - zoom insertion signals</text>
      <text x="0" y="42" font-size="13" fill="#6B778C">region=${formatInt(region.start)}-${formatInt(region.end)} bp | bin=${formatInt(binSize)} bp | max bin=${formatInt(maxBin)}</text>
      <rect x="${left}" y="${top}" width="${plotW}" height="${plotH}" fill="#F6F8FB" stroke="#DDE5EF"/>
      <line x1="${left + ((center - region.start) / (region.end - region.start + 1)) * plotW}" y1="${top}" x2="${left + ((center - region.start) / (region.end - region.start + 1)) * plotW}" y2="${bottom}" stroke="#172033" stroke-width="2" stroke-dasharray="5 5" opacity="0.68"/>
      ${bars.join("")}
      <line x1="${left}" y1="${bottom}" x2="${right}" y2="${bottom}" stroke="#9EADBF"/>
      ${ticks}
    </svg>
  `;
}

function openChrom(chrom) {
  state.selectedChrom = chrom;
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  drawChromPlot();
}

function hideModal() {
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
}

overview.addEventListener("click", (event) => {
  const target = event.target.closest("[data-chrom]");
  if (target) openChrom(target.dataset.chrom);
});
closeModal.addEventListener("click", hideModal);
closeBackdrop.addEventListener("click", hideModal);
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape") hideModal();
});

loadInsertions()
  .then((data) => {
    state.insertions = data;
    state.isLoaded = true;
    renderOverview();
  })
  .catch((error) => {
    overview.innerHTML = '<div class="empty-plot">Impossible de charger les insertions.</div>';
    metadata.textContent = "JSON load error";
    console.error(error);
  });
