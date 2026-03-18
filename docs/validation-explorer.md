# Validation Explorer

This page exposes the main BattMo-versus-experiment validation dataset in an interactive form.

- Use the rate selector to focus on one discharge case or inspect them all together.
- Hover to read off exact capacity and voltage values.
- Zoom or pan to inspect the curve shape near the end of discharge.

<div class="control-row">
  <label for="validation-case-select">Validation case</label>
  <select id="validation-case-select">
    <option>Loading...</option>
  </select>
</div>

<div id="validation-plot" class="plot-container plot-large"></div>

## Error summary

<div id="validation-summary-table" class="table-container">
  Loading validation metrics...
</div>
