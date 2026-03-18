window.__DOCS_DATA__ = window.__DOCS_DATA__ || {};
window.__DOCS_STATE__ = window.__DOCS_STATE__ || {};

function formatNumber(value, digits = 3) {
  return Number(value).toFixed(digits);
}

function createGalleryCard(title, description, imagePath) {
  return `
    <div class="gallery-card">
      <img src="${imagePath}" alt="${title}">
      <h3>${title}</h3>
      <p>${description}</p>
      <p><a class="button-link" href="${imagePath}" target="_blank" rel="noopener">Open image</a></p>
    </div>
  `;
}

function renderTable(container, rows) {
  if (!rows.length) {
    container.textContent = "No rows available.";
    return;
  }
  const headers = Object.keys(rows[0]);
  const thead = `<thead><tr>${headers.map((header) => `<th>${header}</th>`).join("")}</tr></thead>`;
  const tbody = `<tbody>${rows
    .map(
      (row) =>
        `<tr>${headers
          .map((header) => `<td>${typeof row[header] === "number" ? formatNumber(row[header]) : row[header]}</td>`)
          .join("")}</tr>`,
    )
    .join("")}</tbody>`;
  container.innerHTML = `<table>${thead}${tbody}</table>`;
}

function renderJsonViewer(container, text) {
  try {
    const data = JSON.parse(text);
    container.textContent = JSON.stringify(data, null, 2);
  } catch {
    container.textContent = text;
  }
}

function byId(id) {
  return document.getElementById(id);
}

function loadScript(path) {
  return new Promise((resolve, reject) => {
    const script = document.createElement("script");
    script.src = path;
    script.onload = resolve;
    script.onerror = () => reject(new Error(`Failed to load script: ${path}`));
    document.head.appendChild(script);
  });
}

async function loadDataEntry(entry) {
  if (window.__DOCS_DATA__[entry.data_key] !== undefined) {
    return window.__DOCS_DATA__[entry.data_key];
  }
  await loadScript(entry.script);
  return window.__DOCS_DATA__[entry.data_key];
}

function ensurePlotly() {
  if (!window.Plotly) {
    throw new Error("Plotly is not available. Run scripts/prepare_docs_site.py before building the site.");
  }
  return window.Plotly;
}

function loadManifest() {
  const manifest = window.__DOCS_DATA__.site_manifest;
  if (!manifest) {
    throw new Error("Site manifest is not available. Run scripts/prepare_docs_site.py first.");
  }
  return manifest;
}

function preferredSystemTheme() {
  return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
}

function getThemeMode() {
  const saved = window.localStorage.getItem("docs-theme-mode");
  if (saved === "light" || saved === "dark") {
    return saved;
  }
  const attr = document.documentElement.getAttribute("data-bs-theme") || document.body.getAttribute("data-bs-theme");
  if (attr === "light" || attr === "dark") {
    return attr;
  }
  return preferredSystemTheme();
}

function applyThemeMode(mode, persist = true) {
  document.documentElement.setAttribute("data-bs-theme", mode);
  document.body.setAttribute("data-bs-theme", mode);
  if (persist) {
    window.localStorage.setItem("docs-theme-mode", mode);
  }
  const toggleButton = byId("docs-theme-toggle-button");
  const toggleIcon = byId("docs-theme-toggle-icon");
  if (toggleButton && toggleIcon) {
    const nextMode = mode === "dark" ? "light" : "dark";
    toggleButton.setAttribute("aria-label", `Switch to ${nextMode} mode`);
    toggleButton.setAttribute("title", `Switch to ${nextMode} mode`);
    toggleButton.setAttribute("data-mode", mode);
    toggleIcon.className = mode === "dark" ? "fa-solid fa-sun" : "fa-solid fa-moon";
  }
}

function getThemePalette() {
  const style = getComputedStyle(document.documentElement);
  const mode = getThemeMode();
  return {
    mode,
    paper: style.getPropertyValue("--surface").trim() || "#ffffff",
    plot: style.getPropertyValue("--surface").trim() || "#ffffff",
    grid: style.getPropertyValue("--border").trim() || "#d9e2ec",
    font: style.getPropertyValue("--heading").trim() || "#102a43",
    accent: style.getPropertyValue("--accent").trim() || "#0f4c81",
    textSoft: style.getPropertyValue("--text-soft").trim() || "#5b6b79",
    colorway:
      mode === "dark"
        ? ["#7cc4ff", "#ffb259", "#ffd166", "#ff6b6b", "#7bd389", "#c792ea", "#f4a261"]
        : ["#0f4c81", "#d97706", "#c18b00", "#c92a2a", "#2f855a", "#6b46c1", "#8c4f1d"],
  };
}

function plotConfig() {
  return {
    responsive: true,
    displaylogo: false,
  };
}

function baseAxis(title, palette) {
  return {
    title: { text: title, font: { color: palette.font, size: 15 } },
    tickfont: { color: palette.font, size: 12 },
    gridcolor: palette.grid,
    zerolinecolor: palette.grid,
    automargin: true,
  };
}

function baseLayout(title, xTitle, yTitle, extra = {}) {
  const palette = getThemePalette();
  return Object.assign(
    {
      title: { text: title, font: { color: palette.font, size: 20 } },
      autosize: true,
      paper_bgcolor: palette.paper,
      plot_bgcolor: palette.plot,
      font: { color: palette.font, size: 13 },
      colorway: palette.colorway,
      xaxis: baseAxis(xTitle, palette),
      yaxis: baseAxis(yTitle, palette),
      legend: {
        orientation: "h",
        x: 0,
        y: -0.2,
        bgcolor: "rgba(0,0,0,0)",
        font: { color: palette.font, size: 11 },
      },
      margin: { t: 72, r: 28, b: 112, l: 88 },
    },
    extra,
  );
}

function rightLegendLayout(title, xTitle, yTitle, extra = {}) {
  const palette = getThemePalette();
  return baseLayout(title, xTitle, yTitle, Object.assign(
    {
      legend: {
        orientation: "v",
        x: 1.02,
        y: 1,
        xanchor: "left",
        yanchor: "top",
        bgcolor: "rgba(0,0,0,0)",
        font: { color: palette.font, size: 11 },
      },
      margin: { t: 72, r: 168, b: 88, l: 88 },
    },
    extra,
  ));
}

function withPlotHeight(container, layout) {
  const height = Math.max(Math.round(container.clientHeight || 0) - 4, 360);
  return Object.assign({}, layout, { height });
}

function customHeaderLinks() {
  return [
    { href: "index.html", label: "Home" },
    { href: "validation-explorer.html", label: "Validation Explorer" },
    { href: "publication-figures.html", label: "Publication Figures" },
    { href: "supporting-runs.html", label: "Supporting Runs" },
    { href: "citation.html", label: "Citation" },
    { href: "fair-data.html", label: "FAIR Data" },
  ];
}

function installCustomHeader() {
  if (byId("docs-header")) {
    return;
  }

  const current = window.location.pathname.split("/").pop() || "index.html";
  const links = customHeaderLinks()
    .map((entry) => {
      const active = entry.href === current ? ' class="active" aria-current="page"' : "";
      return `<a href="${entry.href}"${active}>${entry.label}</a>`;
    })
    .join("");

  const header = document.createElement("header");
  header.id = "docs-header";
  header.className = "docs-header";
  header.innerHTML = `
    <div class="docs-header-inner">
      <a class="docs-header-brand" href="index.html">Schmitt 2026 Battery Data Explorer</a>
      <nav class="docs-header-nav" aria-label="Primary">
        ${links}
      </nav>
      <div class="docs-header-controls" id="docs-header-controls"></div>
    </div>
  `;

  document.body.insertBefore(header, document.body.firstChild);
}

function installThemeToggle() {
  if (byId("docs-theme-toggle")) {
    return;
  }

  const host = byId("docs-header-controls") || document.body;

  const wrapper = document.createElement("div");
  wrapper.className = "docs-theme-toggle";
  wrapper.id = "docs-theme-toggle";
  wrapper.innerHTML = `
    <button id="docs-theme-toggle-button" type="button" aria-label="Toggle theme" title="Toggle theme">
      <i id="docs-theme-toggle-icon" class="fa-solid fa-moon"></i>
    </button>
  `;
  host.appendChild(wrapper);

  byId("docs-theme-toggle-button").addEventListener("click", async () => {
    const nextMode = getThemeMode() === "dark" ? "light" : "dark";
    applyThemeMode(nextMode);
    await rerenderInteractiveContent();
  });

  applyThemeMode(getThemeMode(), false);
}

async function initOverview(manifest) {
  const gallery = byId("overview-gallery");
  if (gallery) {
    gallery.innerHTML = manifest.publication_gallery
      .map((item) => createGalleryCard(item.title, item.description, item.image))
      .join("");
  }
}

async function initValidationExplorer(manifest) {
  const select = byId("validation-case-select");
  const plotDiv = byId("validation-plot");
  const tableDiv = byId("validation-summary-table");
  if (!select || !plotDiv || !tableDiv) {
    return;
  }

  const Plotly = ensurePlotly();
  const [reference, summary] = await Promise.all([
    loadDataEntry(manifest.data.validation_reference),
    loadDataEntry(manifest.data.validation_summary),
  ]);

  const options = [{ label: "All discharge rates", value: "all" }].concat(
    reference.cases.map((entry, index) => ({
      label: `${entry.case_name} (${formatNumber(entry.current_a, 3)} A)`,
      value: String(index),
    })),
  );
  select.innerHTML = options.map((option) => `<option value="${option.value}">${option.label}</option>`).join("");

  const makeTraces = (caseIndices) =>
    caseIndices.flatMap((caseIndex) => {
      const entry = reference.cases[caseIndex];
      return [
        {
          x: entry.experimental.time_s.map((time) => (time * entry.current_a) / 3600.0),
          y: entry.experimental.voltage_v,
          mode: "lines",
          name: `${entry.case_name} experiment`,
          line: { dash: "dash", width: 2 },
        },
        {
          x: entry.battmo.time_s.map((time) => (time * entry.current_a) / 3600.0),
          y: entry.battmo.voltage_v,
          mode: "lines",
          name: `${entry.case_name} BattMo`,
          line: { width: 2.5 },
        },
      ];
    });

  const render = async () => {
    const selected = select.value;
    const caseIndices = selected === "all" ? reference.cases.map((_, index) => index) : [Number(selected)];
    const title = selected === "all" ? "BattMo validation curves across all rates" : reference.cases[Number(selected)].case_name;
    const layout =
      selected === "all"
        ? rightLegendLayout(title, "Capacity / Ah", "Voltage / V", {
            legend: {
              orientation: "v",
              x: 1.02,
              y: 1,
              xanchor: "left",
              yanchor: "top",
              bgcolor: "rgba(0,0,0,0)",
              font: { color: getThemePalette().font, size: 10 },
            },
            margin: { t: 72, r: 210, b: 88, l: 88 },
          })
        : baseLayout(title, "Capacity / Ah", "Voltage / V", {
            legend: {
              orientation: "h",
              x: 0,
              y: -0.17,
              bgcolor: "rgba(0,0,0,0)",
              font: { color: getThemePalette().font, size: 12 },
            },
            margin: { t: 72, r: 24, b: 94, l: 88 },
          });
    await Plotly.react(
      plotDiv,
      makeTraces(caseIndices),
      withPlotHeight(plotDiv, layout),
      plotConfig(),
    );
  };

  if (!select.dataset.bound) {
    select.addEventListener("change", render);
    select.dataset.bound = "true";
  }
  await render();

  renderTable(
    tableDiv,
    summary.summary_metrics.map((row) => ({
      case: row.case_name,
      current_A: row.current_a,
      rmse_mV: row.rmse_v * 1000.0,
      mae_mV: row.mae_v * 1000.0,
      max_abs_mV: row.max_abs_v * 1000.0,
      final_diff_mV: row.final_voltage_diff_v * 1000.0,
    })),
  );
}

async function initPublicationFigures(manifest) {
  const fig12Discharge = byId("figure12-discharge");
  const fig12Ocv = byId("figure12-ocv");
  const fig13Plot = byId("figure13-plot");
  const fig14Plot = byId("figure14-plot");
  const gallery = byId("publication-gallery");

  if (!fig12Discharge && !fig12Ocv && !fig13Plot && !fig14Plot && !gallery) {
    return;
  }

  const Plotly = ensurePlotly();
  const [fig12, fig13, fig14] = await Promise.all([
    loadDataEntry(manifest.data.figure12),
    loadDataEntry(manifest.data.figure13),
    loadDataEntry(manifest.data.figure14),
  ]);

  if (fig12Discharge) {
    await Plotly.react(
      fig12Discharge,
      [
        { x: fig12.graphite_init.capacity_mAh_cm2, y: fig12.graphite_init.voltage_v, mode: "lines", name: "Graphite init", line: { dash: "dash", width: 2 } },
        { x: fig12.lnmo_init.capacity_mAh_cm2, y: fig12.lnmo_init.voltage_v, mode: "lines", name: "LNMO init", line: { dash: "dash", width: 2 } },
        { x: fig12.full_cell_init.capacity_mAh_cm2, y: fig12.full_cell_init.voltage_v, mode: "lines", name: "Full Cell init", line: { dash: "dash", width: 2 } },
        { x: fig12.graphite_opt.capacity_mAh_cm2, y: fig12.graphite_opt.voltage_v, mode: "lines", name: "Graphite opt", line: { width: 2.5 } },
        { x: fig12.lnmo_opt.capacity_mAh_cm2, y: fig12.lnmo_opt.voltage_v, mode: "lines", name: "LNMO opt", line: { width: 2.5 } },
        { x: fig12.full_cell_opt.capacity_mAh_cm2, y: fig12.full_cell_opt.voltage_v, mode: "lines", name: "Full Cell opt", line: { width: 2.5 } },
        { x: fig12.experiment.capacity_mAh_cm2, y: fig12.experiment.voltage_v, mode: "lines", name: "Experiment 0.05 C", line: { dash: "dot", width: 2 } },
      ],
      withPlotHeight(fig12Discharge, rightLegendLayout(
        "Half-cell and full-cell voltages",
        "Capacity / mAh cm^-2",
        "Voltage / V",
      )),
      plotConfig(),
    );
  }

  if (fig12Ocv) {
    await Plotly.react(
      fig12Ocv,
      [
        { x: fig12.full_cell_init.capacity_mAh_cm2, y: fig12.full_cell_init.voltage_v, mode: "lines", name: "Full Cell init", line: { dash: "dash", width: 2 } },
        { x: fig12.full_cell_opt.capacity_mAh_cm2, y: fig12.full_cell_opt.voltage_v, mode: "lines", name: "Full Cell opt", line: { width: 2.5 } },
        { x: fig12.experiment.capacity_mAh_cm2, y: fig12.experiment.voltage_v, mode: "lines", name: "Experiment 0.05 C", line: { dash: "dot", width: 2 } },
      ],
      withPlotHeight(fig12Ocv, baseLayout("Full-cell zoom", "Capacity / mAh cm^-2", "Voltage / V", {
        yaxis: Object.assign(baseAxis("Voltage / V", getThemePalette()), { range: [3.35, 4.95] }),
      })),
      plotConfig(),
    );
  }

  if (fig13Plot) {
    await Plotly.react(
      fig13Plot,
      [
        { x: fig13.experiment.time_h, y: fig13.experiment.voltage_v, mode: "lines", name: "Experiment 2C", line: { dash: "dash", width: 2 } },
        { x: fig13.initial_guess_dne_1e_14.time_h, y: fig13.initial_guess_dne_1e_14.voltage_v, mode: "lines", name: "Initial D_NE=1e-14", line: { width: 2 } },
        { x: fig13.initial_guess_dne_1e_13.time_h, y: fig13.initial_guess_dne_1e_13.voltage_v, mode: "lines", name: "Initial D_NE=1e-13", line: { width: 2 } },
        { x: fig13.calibrated_from_dne_1e_14.time_h, y: fig13.calibrated_from_dne_1e_14.voltage_v, mode: "lines", name: "Calibrated from 1e-14", line: { width: 2.5 } },
        { x: fig13.calibrated_from_dne_1e_13.time_h, y: fig13.calibrated_from_dne_1e_13.voltage_v, mode: "lines", name: "Calibrated from 1e-13", line: { dash: "dash", width: 2.5 } },
      ],
      withPlotHeight(fig13Plot, rightLegendLayout("High-rate calibration at 2C", "Time / h", "Voltage / V")),
      plotConfig(),
    );
  }

  if (fig14Plot) {
    const traces = fig14.cases.flatMap((entry) => [
      { x: entry.exp_capacity_ah, y: entry.experimental_voltage_v, mode: "lines", name: `${entry.case_name} experiment`, line: { dash: "dash", width: 2 } },
      { x: entry.sim_capacity_ah, y: entry.sim_voltage_v, mode: "lines", name: `${entry.case_name} BattMo`, line: { width: 2.5 } },
    ]);
    await Plotly.react(
      fig14Plot,
      traces,
      withPlotHeight(fig14Plot, rightLegendLayout("Discharge-rate validation", "Capacity / Ah", "Voltage / V")),
      plotConfig(),
    );
  }

  if (gallery) {
    gallery.innerHTML = manifest.publication_gallery
      .map((item) => createGalleryCard(item.title, item.description, item.image))
      .join("");
  }
}

async function initSupportingRuns(manifest) {
  const caseSelect = byId("supporting-case-select");
  const variableSelect = byId("supporting-variable-select");
  const metaDiv = byId("supporting-case-meta");
  const voltagePlot = byId("supporting-voltage-plot");
  const heatmapPlot = byId("supporting-state-heatmap");
  const gallery = byId("supporting-gallery");

  if (!caseSelect || !variableSelect || !metaDiv || !voltagePlot || !heatmapPlot || !gallery) {
    return;
  }

  const Plotly = ensurePlotly();
  const stateData = await loadDataEntry(manifest.data.supporting_states);

  const variables = [
    { key: "elyte_c", label: "Electrolyte concentration", unit: "mol m^-3", xKey: "x_elyte_um" },
    { key: "elyte_phi", label: "Electrolyte potential", unit: "V", xKey: "x_elyte_um" },
    { key: "ne_phi", label: "Negative electrode potential", unit: "V", xKey: "x_ne_um" },
    { key: "pe_phi", label: "Positive electrode potential", unit: "V", xKey: "x_pe_um" },
    { key: "ne_theta", label: "Negative particle surface stoichiometry", unit: "-", xKey: "x_ne_um" },
    { key: "pe_theta", label: "Positive particle surface stoichiometry", unit: "-", xKey: "x_pe_um" },
  ];

  caseSelect.innerHTML = stateData.cases
    .map((entry, index) => `<option value="${index}">${entry.case_name}</option>`)
    .join("");
  variableSelect.innerHTML = variables
    .map((entry) => `<option value="${entry.key}">${entry.label}</option>`)
    .join("");

  const render = async () => {
    const caseEntry = stateData.cases[Number(caseSelect.value)];
    const variable = variables.find((entry) => entry.key === variableSelect.value);
    const palette = getThemePalette();

    metaDiv.innerHTML = `
      <strong>${caseEntry.case_name}</strong><br>
      Current / A: ${formatNumber(caseEntry.current_a, 4)}<br>
      Discharge rate / C: ${formatNumber(caseEntry.drate, 3)}
    `;

    await Plotly.react(
      voltagePlot,
      [
        { x: caseEntry.exp_capacity_ah, y: caseEntry.experimental_voltage_v, mode: "lines", name: "Experiment", line: { dash: "dash", width: 2 } },
        { x: caseEntry.sim_capacity_ah, y: caseEntry.sim_voltage_v, mode: "lines", name: "BattMo", line: { width: 2.5 } },
      ],
      withPlotHeight(voltagePlot, baseLayout(`${caseEntry.case_name} voltage curve`, "Capacity / Ah", "Voltage / V", {
        legend: {
          orientation: "h",
          x: 0,
          y: -0.17,
          bgcolor: "rgba(0,0,0,0)",
          font: { color: getThemePalette().font, size: 12 },
        },
        margin: { t: 72, r: 24, b: 94, l: 88 },
      })),
      plotConfig(),
    );

    await Plotly.react(
      heatmapPlot,
      [
        {
          x: caseEntry.time_h,
          y: caseEntry[variable.xKey],
          z: caseEntry[variable.key],
          type: "heatmap",
          colorscale: palette.mode === "dark" ? "Cividis" : "Viridis",
          colorbar: {
            title: { text: variable.unit, side: "right" },
            tickfont: { color: palette.font },
            thickness: 14,
            len: 0.8,
          },
        },
      ],
      withPlotHeight(heatmapPlot, baseLayout(`${caseEntry.case_name}: ${variable.label}`, "Time / h", "x / um", {
        xaxis: Object.assign(baseAxis("Time / h", palette), { nticks: 6, tickangle: 0 }),
        yaxis: Object.assign(baseAxis("x / um", palette), { nticks: 6 }),
        margin: { t: 72, r: 84, b: 88, l: 84 },
      })),
      plotConfig(),
    );
  };

  if (!caseSelect.dataset.bound) {
    caseSelect.addEventListener("change", render);
    caseSelect.dataset.bound = "true";
  }
  if (!variableSelect.dataset.bound) {
    variableSelect.addEventListener("change", render);
    variableSelect.dataset.bound = "true";
  }
  await render();

  gallery.innerHTML = manifest.supporting_gallery
    .map((item) => createGalleryCard(item.title, item.description, item.image))
    .join("");
}

async function initFairData(manifest) {
  const select = byId("fair-data-select");
  const meta = byId("fair-data-meta");
  const viewer = byId("fair-data-viewer");
  if (!select || !meta || !viewer) {
    return;
  }

  select.innerHTML = manifest.fair_documents
    .map((entry, index) => `<option value="${index}">${entry.label}</option>`)
    .join("");

  const render = async () => {
    const entry = manifest.fair_documents[Number(select.value)];
    meta.innerHTML = `<strong>${entry.label}</strong><br>${entry.description}`;
    const text = await loadDataEntry(entry);
    renderJsonViewer(viewer, text);
  };

  if (!select.dataset.bound) {
    select.addEventListener("change", render);
    select.dataset.bound = "true";
  }
  await render();
}

async function rerenderInteractiveContent() {
  const manifest = window.__DOCS_STATE__.manifest || loadManifest();
  window.__DOCS_STATE__.manifest = manifest;
  await Promise.all([
    initOverview(manifest),
    initValidationExplorer(manifest),
    initPublicationFigures(manifest),
    initSupportingRuns(manifest),
    initFairData(manifest),
  ]);
}

document.addEventListener("DOMContentLoaded", async () => {
  try {
    installCustomHeader();
    installThemeToggle();
    applyThemeMode(getThemeMode(), false);
    window.__DOCS_STATE__.manifest = loadManifest();

    const media = window.matchMedia("(prefers-color-scheme: dark)");
    media.addEventListener("change", async () => {
      if (!window.localStorage.getItem("docs-theme-mode")) {
        applyThemeMode(preferredSystemTheme(), false);
        await rerenderInteractiveContent();
      }
    });

    await rerenderInteractiveContent();
  } catch (error) {
    console.error(error);
    document.body.insertAdjacentHTML(
      "afterbegin",
      `<div class="callout"><strong>Documentation error:</strong> ${error.message}</div>`,
    );
  }
});
