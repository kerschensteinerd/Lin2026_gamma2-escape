# Calcium supplement package

This package presents the completed calcium transformation sweep at publication scale and leaves the current manuscript and artwork unchanged. Figure numbering is provisional.

- `Figure_Sx_calcium_sweep.pdf`, `.svg`, `.png`: eight-panel main figure, 180 × 225 mm.
- `Figure_Sx_calcium_retained_inhibition.pdf`, `.svg`, `.png`: companion sensitivity page, 180 × 181 mm.
- `legend.md`: complete main and companion legends, definitions, and limitations.
- `proposed_manuscript_text.md`: Results and Methods additions and an optional Discussion sentence; not inserted into the manuscript.
- `source_data/`: the exact plotted arrays, tabular map and curve data, and the mapping to the validated numerical archive.
- `build_calcium_supplement.py`: reproducible figure and source-data builder; run with the model repository's Python environment.
- `manifest.json`: input/builder/output hashes, dimensions, fonts, and dependency versions.
- `qa/`: exported-PDF renders, array checks, and visual-review record.

The maps show the entire prescribed parameter domain. The selected transformation preserves preferential fractional enhancement, while the domain-corner example illustrates a contrary regime. Panels A and G show only the selected parameters, with solid control and reduced-inhibition curves. Technical masks are displayed separately from directional reversals. Calcium spread uses the original focal site and nominal 25% retained inhibition; additional electrical, site, and spatial-layout results are separate analyses.

The PDF/SVG exports preserve vector text, axes, curves, contours, and symbols. Exact parameter-cell colors are embedded at 600 dpi to eliminate PDF viewer seam artifacts; the heatmaps are not smoothed. PNG previews use 300 dpi. Both PDF pages were generated with the required artifact-operation marker executed once before the first authoring run.
