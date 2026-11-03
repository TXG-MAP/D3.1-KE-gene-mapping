# Caption-free figures

The same figures as `figures/final/`, drawn from the same scripts in the same run, with the
muted provenance strip under the plot removed — the line carrying the release, the organ
frame, and the test statistics.

These are the ones used in the report, where that text is the figure caption and
would otherwise be printed twice, once as an image and once as a caption. `figures/final/`
keeps the strip and stays the standalone version: it is what
should go into a slide or a Zenodo deposit, where nothing else carries the provenance.

Rebuild with:

```bash
make figures-clean
```

which sets `D31_FIG_CLEAN=1` and re-runs the four figure scripts offline, off the derived
tables. It never writes into `figures/final/`, and it does not restamp
`data/method1/manifest_05_figures.json` — that manifest is the provenance record for the
captioned build.

Eight figures have a strip to remove:

| | figure | script |
| --- | --- | --- |
| FIG01 | method overview | `analysis/figures/dia_01_schematics.py` |
| FIG07 | Method 1 gene Venn, four organs | `analysis/method1/m1_05_figures.py` |
| FIG08 | regex against BERN2 | `analysis/method1/m1_05_figures.py` |
| FIG09 | genes per Key Event | `analysis/method1/m1_05_figures.py` |
| FIG10 | Method 3 coverage | `analysis/method3/m3_15_figures.py` |
| FIG11 | cross-method comparison | `analysis/comparison/cmp_04_figures.py` |
| FIG12 | Method 1 coverage by structure | `analysis/method1/m1_05_figures.py` |
| FIG14 | Method 1 denominator | `analysis/method1/m1_05_figures.py` |

Some of them also lose width, not only height: on FIG09, FIG10 and FIG11 the strip ran wider
than the plot and was what set the figure's bounding box.

Every other figure is identical in both builds and is only in `figures/final/`. FIG04's only
note sits above the tree and states rule 7 of the rubric, so it is part of the diagram rather
than a caption and is not stripped. FIG18's legend and the panel letters on the screenshot
composites stay for the same reason.
