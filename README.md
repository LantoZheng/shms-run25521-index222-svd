# Run 25521 | SHMS angular SVD with 222 indexed clusters

**Start with the [12-page unweighted comparison](comparison/unweighted_highlights.pdf).** The [325-page atlas](comparison/unweighted_complete_atlas.pdf) pairs published and current Run 25521 figures page by page, including old and new optics, the original diagnostic macros, NPS macros, six momentum slices and three foils. Bookmarks identify every section. “Old” and “new” within a panel refer to the starting and emitted matrices for that analysis. Left and right pages compare the two indexing versions. Native plot scales are retained.

## Inputs and outputs

106,058 events across 222 clusters. All available event-count plots include the 147 provisional C199 events under their stored foil. Those events are excluded from both fits. Native gates, a 1,000-event cell cap and equal populations per foil select 73,989 events (24,663 each). Unweighted and weighted runs use identical selected events; the weighted run applies inverse occupancy in foil/column/row/one-percent momentum cells. No inherited event split. Fits use the upstream native C++ ROOT solver, the same starting optics and nominal geometry. Geometric xtar is fixed, and the matrices reconstruct xptar, ytar and yptar. The run is a matrix evaluation, without a full detector replay.

| Version | X′tar RMSE (mrad) | Ytar RMSE (mm) | Y′tar RMSE (mrad) |
|---|---:|---:|---:|
| Published unweighted, 73,188 fitted events | 2.021 | 4.319 | 1.513 |
| Current unweighted, 73,989 fitted events | 1.087 | 3.265 | 1.448 |
| Current weighted, same events, ordinary RMSE | 1.092 | 3.303 | 1.454 |

On the same 64,627 fitted events with unchanged foil and hole targets, historical to current unweighted RMSE is 1.437 to 1.093 mrad, 3.690 to 3.288 mm and 1.472 to 1.458 mrad. The broader overlap is 70,841 events. Evaluating that overlap against the historical targets favors the historical matrix; evaluating against the current targets favors the current matrix. These are fit-sample closure measures. [Common-event results for both target conventions](comparison/common_event_rmse.csv)

## Figures and data

- [Unweighted replay matrix](run25521/unweighted/matrices/nps_shms_newfit_25521_unweighted.dat), [fit report](run25521/unweighted/fit_report.json), [starting versus fitted delta–ytar](run25521/unweighted/plots/comparison/25521_delta_vs_ytar_before_after.pdf), [all unweighted plots](run25521/unweighted/plots)
- [Weighted replay matrix](run25521/weighted/matrices/nps_shms_newfit_25521_weighted.dat), [fit report](run25521/weighted/fit_report.json), [all weighted plots](run25521/weighted/plots)
- [Foil and momentum slice index](run25521/DIAGNOSTICS.md), [NPS legacy diagnostics](run25521/NPS_LEGACY_DIAGNOSTICS.md), [available-event delta](run25521/plots/25521_available_delta_histograms.pdf), [selected-event delta](run25521/plots/25521_delta_histograms.pdf)
- [Cluster-to-hole indexing](indexing/cluster_index.csv), [selection](run25521/selection.json), [plot inventory](plot_inventory.csv), [numerical validation](run25521/validation.json), [HCANA matrix format check](run25521/replay_format_validation.json)

## Reproduction

The historical results and plotting workflow are from [cploen/shms-run25523-svd-results at 914ed10](https://github.com/cploen/shms-run25523-svd-results/tree/914ed103f5850e470d6cc26ed65b8637144ac18d/run25521). The public deliverables include aggregate results and the per-cluster index. Rerunning the fit requires the indexed 106,058-event NPZ and ROOT inputs used to produce them.

1. Check out the pinned upstream commit and apply [the input-selection and count updates](reproduce/index222.patch) using git apply. Use a clean output directory: move the checkout's existing run25521/weighted, run25521/unweighted and run25521/plots outputs aside before running. Keep the top-level weighted/matrices and unweighted/matrices format-check fixtures.
2. Use the pinned [upstream requirements](https://github.com/cploen/shms-run25523-svd-results/blob/914ed103f5850e470d6cc26ed65b8637144ac18d/run25521/reproduce/requirements.txt), a working C++17 compiler and the HCANA reader source at ce5f975288010984a6d66fd3bde4bee95fc36ba2. [Prepare the working inputs](reproduce/prepare_handoff.py) with Python, passing the indexed NPZ path and a new staging directory.
3. From the upstream checkout run bash run25521/reproduce/run.sh INPUT_STAGE WORK_DIR HCANA_CHECKOUT using paths without spaces for the ROOT macro work and checkout directories. The two fits differ only in event weights.
4. Run [compare_atlas.py](reproduce/compare_atlas.py) with the published run25521 directory, the new run25521 directory and an output directory. It asserts matching page counts before pairing pages. [Validation summary](validation_summary.json)
