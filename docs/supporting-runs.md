# Supporting Runs

This page exposes the per-run BattMo supporting data that sits behind the validation figures.

- Select a discharge case to inspect its voltage curve and state evolution.
- Switch between electrolyte, solid potential, and particle-stoichiometry fields.
- Use the heatmap zoom tools to inspect time windows or spatial regions.

<div class="control-row">
  <label for="supporting-case-select">Discharge case</label>
  <select id="supporting-case-select">
    <option>Loading...</option>
  </select>

  <label for="supporting-variable-select">State variable</label>
  <select id="supporting-variable-select">
    <option>Loading...</option>
  </select>
</div>

<div id="supporting-case-meta" class="callout"></div>

<div class="dual-plot-grid">
  <div id="supporting-voltage-plot" class="plot-container"></div>
  <div id="supporting-state-heatmap" class="plot-container"></div>
</div>
