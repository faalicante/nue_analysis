# Presentation sources

Historical scripts used to prepare the project overview and the H-mu review. Run from the `CNN` project root. Python scripts can be invoked by their path here; they explicitly add the working directory to the import path when needed. JavaScript builders require the OpenAI artifact-tool runtime and the presentation helper paths recorded in their imports. The temporary `node_modules` links and plotting dependency copies were removed; reconnect/install the required runtime before rebuilding a deck.

Final decks, figures and tables remain in `output/`. JSON input snapshots and selected source figures are retained beside these scripts. Rebuilding may also need original CERNBox inputs and regeneration of intermediate previews.

| Code | Purpose |
|---|---|
| [project_overview/build.mjs](project_overview/build.mjs) | Build the original Italian project overview. |
| [project_overview/build_en.mjs](project_overview/build_en.mjs) | Build the English project overview. |
| [project_overview/translate_en.py](project_overview/translate_en.py) | Derive the English builder and notes from the original source. |
| [project_overview/plots_en.py](project_overview/plots_en.py) | Generate the English augmentation and candidate figures. |
| [project_overview/render_en.mjs](project_overview/render_en.mjs) | Render the English deck for visual inspection. |
| [project_overview/render_final.mjs](project_overview/render_final.mjs) | Render the Italian overview for visual inspection. |
| [project_overview/export_keynote.applescript](project_overview/export_keynote.applescript) | Export presentation documents using Keynote. |
| [hmu_review/analyze.py](hmu_review/analyze.py) | Summarize candidate counts, background scans and signal angular performance. |
| [hmu_review/candidate_assets.py](hmu_review/candidate_assets.py) | Generate selected candidate graphics and layer animations. |
| [hmu_review/prepare_deck.py](hmu_review/prepare_deck.py) | Assemble the review builder from overview slides and review-specific results. |
| [hmu_review/build_deck.mjs](hmu_review/build_deck.mjs) | Build the H-mu review presentation. |
| [hmu_review/render_deck.mjs](hmu_review/render_deck.mjs) | Render the retained version-4 review deck. |
| [hmu_review/make_figures.py](hmu_review/make_figures.py) | Generate review performance and schematic figures. |
| [hmu_review/make_bkg_mu_histogram.py](hmu_review/make_bkg_mu_histogram.py) | Plot background-mu distributions by brick. |
| [hmu_review/make_bkg_mu_means.py](hmu_review/make_bkg_mu_means.py) | Plot mean background levels by brick. |
| [hmu_review/make_candidate_rates.py](hmu_review/make_candidate_rates.py) | Compare candidate rates across bricks. |
| [hmu_review/make_cnn_architecture_schematic.py](hmu_review/make_cnn_architecture_schematic.py) | Draw the CNN architecture schematic. |
| [hmu_review/make_roc_pr_test_figure.py](hmu_review/make_roc_pr_test_figure.py) | Prepare the ROC/precision-recall test figure. |
| [hmu_review/make_extra_data_assets.py](hmu_review/make_extra_data_assets.py) | Prepare graphics for additional selected data candidates. |
| [hmu_review/preview_data.py](hmu_review/preview_data.py) | Preview the initial data candidate pool. |
| [hmu_review/preview_extra_data.py](hmu_review/preview_extra_data.py) | Preview the additional selected data candidates. |
