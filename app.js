(function () {
  "use strict";
  const data = window.BENCHMARK;
  const $ = (id) => document.getElementById(id);
  if (!data || !Array.isArray(data.metrics) || !Array.isArray(data.pairs)) {
    $("main-context").textContent = "The analysis data could not be loaded. Please refresh, or use the linked results and static figures.";
    return;
  }

  const cohort = $("cohort-select");
  const range = $("range-select");
  const gene = $("gene-search");
  const pairSelect = $("pair-select");
  const svgNS = "http://www.w3.org/2000/svg";
  let visiblePairs = [];
  let currentPairId = null;
  const finite = (value) => value !== null && value !== undefined && Number.isFinite(Number(value));
  const number = (value, digits = 3) => finite(value) ? Number(value).toFixed(digits) : "not available";
  const integer = (value) => finite(value) ? Number(value).toLocaleString("en-US") : "—";
  const percent = (value) => finite(value) ? `${(Number(value) * 100).toFixed(1)}%` : "not reported";
  const signed = (value) => finite(value) ? `${Number(value) >= 0 ? "+" : ""}${Number(value).toFixed(3)}` : "not available";
  const yes = (value) => value === true || value === 1 || value === "TRUE";
  const pairById = new Map(data.pairs.map((pair) => [String(pair.id), pair]));

  function svgElement(tag, attributes, text) {
    const element = document.createElementNS(svgNS, tag);
    Object.entries(attributes || {}).forEach(([key, value]) => element.setAttribute(key, String(value)));
    if (text !== undefined) element.textContent = text;
    return element;
  }

  function chartText(svg, text, x, y, attributes = {}) {
    svg.appendChild(svgElement("text", { x, y, fill: "#666", "font-size": 13, ...attributes }, text));
  }

  function inRange(pair) {
    const d = Number(pair.distance_bp);
    if (range.value === "0–20 kb") return d >= 0 && d < 20000;
    if (range.value === "20–100 kb") return d >= 20000 && d < 100000;
    if (range.value === "100–2,500 kb") return d >= 100000 && d < 2500000;
    return true;
  }

  function filterPairs() {
    return data.pairs.filter((pair) => (cohort.value === "all" || yes(pair.matched)) && inRange(pair));
  }

  function drawPR(metric) {
    const svg = $("pr-chart");
    svg.replaceChildren();
    svg.appendChild(svgElement("title", { id: "pr-title" }, `Precision–recall: ${range.options[range.selectedIndex].text}`));
    svg.appendChild(svgElement("desc", { id: "pr-description" }, `Same ${metric.n} tested pairs. ABC AP ${number(metric.abc_ap)}; distance AP ${number(metric.distance_ap)}. AP is computed separately from drawn line area.`));
    const box = { left: 64, top: 42, width: 506, height: 266 };
    const x = (r) => box.left + Number(r) * box.width;
    const y = (p) => box.top + (1 - Number(p)) * box.height;
    [0, .25, .5, .75, 1].forEach((tick) => {
      svg.appendChild(svgElement("line", { x1: x(0), x2: x(1), y1: y(tick), y2: y(tick), stroke: "#e6e6e6" }));
      chartText(svg, tick === 0 ? "0" : String(tick), x(0) - 10, y(tick) + 5, { "text-anchor": "end" });
      chartText(svg, tick === 0 ? "0" : String(tick), x(tick), y(0) + 24, { "text-anchor": "middle" });
    });
    chartText(svg, "Precision", box.left, 22, { fill: "#333", "font-size": 14 });
    chartText(svg, "Recall", x(.5), 355, { fill: "#333", "text-anchor": "middle", "font-size": 14 });
    if (finite(metric.prevalence)) svg.appendChild(svgElement("line", { x1: x(0), x2: x(1), y1: y(metric.prevalence), y2: y(metric.prevalence), stroke: "#b6b6b6", "stroke-dasharray": "2 5" }));
    const curves = data.curves[`${cohort.value}|${range.value}`] || {};
    [["distance", "#777", "6 4"], ["abc", "#1f7a8c", ""]].forEach(([name, color, dash]) => {
      const points = (curves[name] || []).filter((point) => finite(point[0]) && finite(point[1]));
      if (!points.length) return;
      svg.appendChild(svgElement("polyline", {
        points: points.map((point) => `${x(point[0]).toFixed(2)},${y(point[1]).toFixed(2)}`).join(" "),
        fill: "none", stroke: color, "stroke-width": 2.5, "stroke-dasharray": dash,
        "stroke-linejoin": "round"
      }));
    });
    [["ABC", "#1f7a8c", "", 90], ["Distance", "#777", "6 4", 210], ["Positive fraction", "#b6b6b6", "2 5", 350]].forEach(([label, color, dash, left]) => {
      svg.appendChild(svgElement("line", { x1: left, x2: left + 25, y1: 380, y2: 380, stroke: color, "stroke-width": 2, "stroke-dasharray": dash }));
      chartText(svg, label, left + 33, 384, { "font-size": 12 });
    });
  }

  function updateBenchmark() {
    const metric = data.metrics.find((entry) => entry.cohort === cohort.value && entry.range === range.value);
    visiblePairs = filterPairs();
    if (!metric) {
      $("readout").textContent = "No summary is available for this selection.";
      return;
    }
    drawPR(metric);
    $("metric-abc").textContent = number(metric.abc_ap);
    $("metric-distance").textContent = number(metric.distance_ap);
    $("metric-pairs").textContent = integer(metric.n);
    $("metric-pairs-label").textContent = `${integer(metric.positives)} downregulatory positives`;
    const missingText = cohort.value === "all"
      ? `${integer(metric.missing)} pairs have no matched prediction and are assigned score zero; both methods still rank every pair.`
      : "Both methods are restricted to the same matched subset; this is a sensitivity check, not a replacement for the official cohort.";
    $("readout").textContent = `${integer(metric.n)} pairs · ${integer(metric.positives)} positives (${percent(metric.prevalence)}). ABC minus distance AP: ${signed(metric.delta)}. ${missingText}`;
    const ci = finite(metric.delta_low) && finite(metric.delta_high)
      ? `Paired AP difference: ${signed(metric.delta)}; 95% gene-cluster bootstrap interval ${signed(metric.delta_low)} to ${signed(metric.delta_high)}. ` : "";
    const unresolved = finite(metric.delta_low) && metric.delta_low <= 0 && metric.delta_high >= 0
      ? "The interval includes zero: this stratum does not resolve a precise advantage. " : "";
    const omitted = metric.bootstrap_draws - metric.bootstrap_valid;
    const resampling = omitted > 0 ? `${omitted} of ${metric.bootstrap_draws} resamples had no positives; the interval uses the remaining draws. ` : "";
    $("uncertainty-note").textContent = `${ci}${unresolved}${resampling}These intervals reflect this selected gene set, not independent validation or a guarantee of transfer.`;
    $("chart-caption").textContent = `${range.value === "all" ? "All included distances" : range.value}; ${cohort.value === "all" ? "official zero-fill policy" : "matched subset"}. Ranking thresholds sweep all available scores. AP is a step-weighted summary, not the area under the displayed line segments.`;
    updateGeneOptions();
    updatePairOptions();
  }

  function updateGeneOptions() {
    const genes = [...new Set(visiblePairs.map((pair) => pair.gene))].sort();
    $("gene-options").replaceChildren(...genes.map((name) => {
      const option = document.createElement("option");
      option.value = name;
      return option;
    }));
    if (gene.value && !genes.some((name) => name.toLowerCase() === gene.value.trim().toLowerCase())) gene.value = "";
  }

  function updatePairOptions(requestedId = currentPairId) {
    const query = gene.value.trim().toLowerCase();
    let filtered = visiblePairs.filter((pair) => !query || pair.gene.toLowerCase() === query);
    filtered = filtered.sort((a, b) => a.gene.localeCompare(b.gene) || Number(a.distance_bp) - Number(b.distance_bp) || String(a.id).localeCompare(String(b.id)));
    pairSelect.replaceChildren(...filtered.map((pair) => {
      const option = document.createElement("option");
      option.value = String(pair.id);
      option.textContent = `${pair.gene} · ${pair.chrom}:${integer(pair.start)}–${integer(pair.end)} · ${number(pair.distance_bp / 1000, 1)} kb`;
      return option;
    }));
    $("pair-count").textContent = filtered.length
      ? `${integer(filtered.length)} pairs${query ? ` for ${filtered[0].gene}` : " across all genes"} in the current benchmark selection. Clear the gene field to see all genes.`
      : "No exact gene symbol matches this selection. Choose a suggested symbol, clear the field, or change the cohort filters.";
    $("pair-view").hidden = filtered.length === 0;
    pairSelect.disabled = filtered.length === 0;
    if (!filtered.length) return;
    const choice = filtered.find((pair) => String(pair.id) === String(requestedId)) || filtered[0];
    pairSelect.value = String(choice.id);
    renderPair(choice);
  }

  function fact(label, value) {
    const item = document.createElement("div");
    const key = document.createElement("span");
    const detail = document.createElement("strong");
    key.textContent = label;
    detail.textContent = value;
    item.append(key, detail);
    return item;
  }

  function drawLocus(pair) {
    const svg = $("locus-chart");
    svg.replaceChildren();
    svg.appendChild(svgElement("title", { id: "locus-title" }, `${pair.gene}: element and annotated TSS`));
    svg.appendChild(svgElement("desc", { id: "locus-description" }, `${pair.chrom}, GRCh38. Distance ${integer(pair.distance_bp)} bases. Schematic shows positions, not a measured contact.`));
    const tss = Number(pair.tss);
    const center = (Number(pair.start) + Number(pair.end)) / 2;
    if (!finite(pair.tss)) {
      chartText(svg, "TSS coordinate unavailable; inspect the source output.", 35, 95);
      return;
    }
    const low = Math.min(center, tss);
    const span = Math.max(Math.abs(center - tss), 1000);
    const x = (position) => 55 + 490 * ((position - low + span * .12) / (span * 1.24));
    const ex = x(center), tx = x(tss);
    svg.appendChild(svgElement("line", { x1: 35, x2: 565, y1: 90, y2: 90, stroke: "#999", "stroke-width": 1.5 }));
    svg.appendChild(svgElement("rect", { x: ex - 5, y: 77, width: 10, height: 26, fill: "#1f7a8c" }));
    svg.appendChild(svgElement("line", { x1: tx, x2: tx, y1: 70, y2: 110, stroke: "#222", "stroke-width": 2.5 }));
    chartText(svg, "Tested element", ex, 60, { "text-anchor": "middle", fill: "#1f7a8c" });
    chartText(svg, `${pair.gene} TSS`, tx, 138, { "text-anchor": "middle", fill: "#222" });
    chartText(svg, `${number(pair.distance_bp / 1000, 1)} kb · linear distance, not a measured contact`, 300, 22, { "text-anchor": "middle", "font-size": 12 });
    chartText(svg, `${pair.chrom} · GRCh38 · markers enlarged for visibility`, 300, 178, { "text-anchor": "middle", "font-size": 12 });
  }

  function renderPair(pair) {
    currentPairId = String(pair.id);
    const matched = yes(pair.matched), positive = yes(pair.regulated), above = Number(pair.abc_score) >= .02;
    $("pair-title").textContent = `${pair.gene} · ${pair.chrom}:${integer(pair.start)}–${integer(pair.end)}`;
    const example = (data.examples || []).find((entry) => String(entry.id) === String(pair.id));
    $("example-reason").textContent = example ? example.reason : "";
    $("example-reason").hidden = !example;
    drawLocus(pair);
    $("pair-facts").replaceChildren(
      fact("ABC power-law score", `${number(pair.abc_score, 4)}${matched ? "" : " · zero-filled"}`),
      fact("Source Regulated label", positive ? "TRUE · downregulatory" : "FALSE · not called downregulatory"),
      fact("Reported effect size", finite(pair.effect) ? number(pair.effect, 3) : "not reported"),
      fact("Adjusted p-value", finite(pair.padj) ? Number(pair.padj).toPrecision(3) : "not reported"),
      fact("Power at 25% effect", percent(pair.power)),
      fact("Prediction coverage", matched ? `${integer(pair.overlap_count)} overlapping element(s)` : "No matched prediction")
    );
    let title, interpretation;
    if (!matched) {
      title = positive ? "Observed response, but no matched prediction" : "No matched prediction; the experimental label remains separate";
      interpretation = "The zero comes from the benchmark’s missing-prediction policy, not a model measurement or proof of a negative result. Inspect prediction coverage and the TSS filtering step before interpreting this as a ranking failure.";
    } else if (above && positive) {
      title = "Agreement: above the ABC reference threshold and downregulatory";
      interpretation = "The source reports a downregulatory response and ABC gives this pair a score of at least 0.02. This is within-benchmark agreement, not a new validation experiment or proof that the score measures effect size.";
    } else if (!above && positive) {
      title = "Disagreement: downregulatory response below the ABC threshold";
      interpretation = "The experiment meets the source’s downregulatory label, but the saved score is below 0.02. Compare activity/contact assumptions and interval coverage with the experimental evidence. The supplied threshold is not optimized here.";
    } else if (above) {
      title = "Disagreement: high ABC score without a downregulatory label";
      interpretation = "This is not automatically a false biological connection. Inspect effect direction, adjusted p-value, and power; the assay may not establish the absence of an effect. A significant increase is deliberately not counted as downregulatory.";
    } else {
      title = "Agreement of labels, not proof of no regulation";
      interpretation = "The saved score is below 0.02 and the source does not call a downregulatory response. Assay power and the observed effect remain important; neither result establishes that this element cannot regulate this gene.";
    }
    $("pair-category").textContent = title;
    $("pair-interpretation").textContent = interpretation;
    $("pair-provenance").textContent = `Source: ${pair.reference || "see source table"} · K562 · GRCh38 numeric coordinates. ${matched ? "Scores for multiple overlapping predicted elements are summed." : `Missingness: ${pair.missing_reason || "no eligible overlapping prediction"}.`} Effect size is shown in the source table’s units, without assuming values are comparable across studies.`;
    document.querySelectorAll("#example-buttons button").forEach((button) => button.setAttribute("aria-pressed", String(button.dataset.id === String(pair.id))));
  }

  function selectExample(example) {
    const pair = pairById.get(String(example.id));
    if (!pair) return;
    let reset = false;
    if (!visiblePairs.some((entry) => String(entry.id) === String(pair.id))) {
      cohort.value = "all";
      range.value = "all";
      reset = true;
    }
    gene.value = pair.gene;
    currentPairId = String(pair.id);
    updateBenchmark();
    if (reset) $("pair-count").textContent += " Benchmark filters were reset to all eligible pairs to show this example.";
  }

  function exportCSV() {
    const columns = ["id", "gene", "chrom", "start", "end", "tss", "distance_bp", "abc_score", "matched", "overlap_count", "regulated", "significant", "effect", "padj", "power", "reference", "missing_reason", "category"];
    const escapeCSV = (value) => `"${String(value === undefined || value === null ? "" : value).replace(/"/g, '""')}"`;
    const csv = [columns.join(","), ...visiblePairs.map((pair) => columns.map((key) => escapeCSV(pair[key])).join(","))].join("\r\n");
    const url = URL.createObjectURL(new Blob([csv], { type: "text/csv;charset=utf-8" }));
    const link = document.createElement("a");
    link.href = url;
    link.download = `k562-${cohort.value}-${range.value.replace(/[^a-z0-9]+/gi, "-")}.csv`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  (data.examples || []).forEach((example, index) => {
    if (!pairById.has(String(example.id))) return;
    const button = document.createElement("button");
    button.type = "button";
    button.dataset.id = String(example.id);
    button.textContent = `${index + 1}. ${example.title}`;
    button.setAttribute("aria-pressed", "false");
    button.addEventListener("click", () => selectExample(example));
    $("example-buttons").appendChild(button);
  });
  const overall = data.metrics.find((entry) => entry.cohort === "all" && entry.range === "all");
  if (overall) {
    const direction = Number(overall.delta) > 0 ? "ABC ranks the observed responses better overall" : Number(overall.delta) < 0 ? "Distance ranks the observed responses better overall" : "The overall average precision is equal";
    $("main-finding").textContent = direction;
    $("main-context").textContent = `Across ${integer(overall.n)} eligible pairs, ABC AP is ${number(overall.abc_ap)} versus ${number(overall.distance_ap)} for distance. The useful next question is where that difference persists—and which individual responses remain unexplained.`;
  }
  if (data.official && finite(data.official.abc_auprc) && finite(data.official.distance_auprc)) {
    $("official-result").textContent = `The reproduced official R calculation gives trapezoidal AUPRC ${number(data.official.abc_auprc, 6)} for ABC and ${number(data.official.distance_auprc, 6)} for distance. Those values are separate from the AP values used in the added comparisons above.`;
  }
  cohort.addEventListener("change", updateBenchmark);
  range.addEventListener("change", updateBenchmark);
  $("gene-form").addEventListener("submit", (event) => { event.preventDefault(); updatePairOptions(); });
  gene.addEventListener("change", () => updatePairOptions());
  pairSelect.addEventListener("change", () => {
    const pair = pairById.get(pairSelect.value);
    if (pair) renderPair(pair);
  });
  $("export-cohort").addEventListener("click", exportCSV);
  updateBenchmark();
  if (data.examples && data.examples.length) selectExample(data.examples[0]);
}());
