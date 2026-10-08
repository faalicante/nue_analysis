from pathlib import Path
s=Path('visualization/presentations/project_overview/build.mjs').read_text()
r={
'Ricerca di sciami\\ncon una CNN 2+1D':'Shower identification\\nwith a 2+1D CNN',
'Progetto crop_bw\\nSintesi del metodo e dei risultati':'The crop_bw project\\nMethod and results overview',
'15 settembre 2026':'15 September 2026',
'Obiettivo e flusso di analisi':'Objective and analysis workflow',
'Individuare candidati sciame nei volumi e stimarne la direzione':'Find shower candidates in volumes and estimate their direction',
'01   Mappe ROOT':'01   ROOT maps',
'Conteggi grezzi per piano e stima del fondo.':'Raw counts for each layer and a background estimate.',
'02   Dataset e CNN':'02   Dataset and CNN',
'Crop HDF5, classificazione del segnale e regressione delle slope.':'HDF5 crops, signal classification and slope regression.',
'03   Scansione dei volumi':'03   Volume scanning',
'Finestre sovrapposte, raggruppamento dei candidati e proiezioni XZ/YZ.':'Overlapping windows, candidate clustering and XZ/YZ projections.',
'Dataset e controlli di qualità':'Dataset and quality checks',
'volume grezzo per esempio':'raw volume per sample',
'crop all’ingresso della rete':'crop at the network input',
'15.872 esempi nel dataset fitlt5\\n8.656 hard negatives e 7.216 signal':'15,872 samples in the fitlt5 dataset\\n8,656 hard negatives and 7,216 signal samples',
'Separazione per evento di segnale o cella di fondo.\\nEsclusione della cella 58 e controllo del supporto ai bordi.\\nI conteggi originali restano disponibili senza smoothing.':'Splits grouped by signal event or background cell.\\nCell 58 excluded and boundary support checked.\\nOriginal counts are retained without smoothing.',
'Train / validation / test effettivi: 10.186 / 2.377 / 1.996 esempi':'Actual train / validation / test subsets: 10,186 / 2,377 / 1,996 samples',
'CNN 2+1D a due teste':'A two-head 2+1D CNN',
'Convoluzioni separate nello spazio e lungo z':'Separate convolutions in the spatial and z directions',
'Struttura del singolo piano XY':'Structure within each XY layer',
'Evoluzione tra piani consecutivi':'Evolution across consecutive layers',
'Due uscite complementari':'Two complementary outputs',
'Score di segnale e componenti della direzione (sx, sy).':'Signal score and direction components (sx, sy).',
'4 blocchi, canali 16 / 32 / 64 / 128, 196.243 parametri':'4 blocks, channels 16 / 32 / 64 / 128, 196,243 parameters',
'Training e trasformazioni coerenti':'Training and consistent augmentation',
'Jitter ±5 pixel\\nRotazioni e riflessioni\\nSlope trasformate insieme':'Jitter of ±5 pixels\\nRotations and reflections\\nSlopes transformed together',
'Signal >8 mrad per la classificazione. Tutti i signal per la regressione.':'Signal above 8 mrad for classification. All signal samples for regression.',
'Checkpoint separati per classificazione e regressione. Augmentation solo nel training.':'Separate classification and regression checkpoints. Augmentation during training only.',
'Prestazioni sui crop di test':'Performance on test crops',
'Confronto delle tre normalizzazioni fitlt5':'Comparison of the three fitlt5 normalization methods',
'Ingresso CNN':'CNN input', 'Precisione':'Precision',
'Score ≥0,50. Classificazione: 613 signal e 1.117 hard. Regressione: 879 signal.':'Score ≥0.50. Classification: 613 signal and 1,117 hard. Regression: 879 signal samples.',
'Scansione dei volumi completi':'Scanning full volumes',
'Finestre 20 × 20, passo 10 pixel':'20 × 20 windows, 10-pixel stride',
'La scansione conserva tutti i 57 piani e valuta finestre sovrapposte.':'The scan retains all 57 layers and evaluates overlapping windows.',
'Score e direzione':'Score and direction',
'Il checkpoint di classificazione seleziona le finestre. Quello di regressione stima le slope.':'The classification checkpoint selects windows. The regression checkpoint estimates slopes.',
'Candidati da ispezionare':'Candidates for inspection',
'Il raggruppamento delle finestre produce candidati con coordinate, score e proiezioni XZ/YZ.':'Window clustering produces candidates with coordinates, scores and XZ/YZ projections.',
'Script recente sui dati: H − μ, score ≥0,90. Il candidato richiede validazione fisica.':'Recent data scan configuration: H − μ, score ≥0.90. Candidates require physics validation.',
'Efficienza dello scan per angolo':'Scan efficiency by angle',
'H − μ, score ≥0,90, associazione alla verità MC':'H − μ, score ≥0.90, matched to MC truth',
'Efficienza':'Efficiency','θ vero [mrad]':'True θ [mrad]',
'98,7%':'98.7%',
'528 / 535 eventi\\ncon θ >10 mrad':'528 / 535 events\\nwith θ >10 mrad',
'Il limite principale è\\nla regione 5–10 mrad.':'The main limitation is\\nthe 5–10 mrad range.',
'Eventi per bin: 223 / 245 / 253 / 37. Nessun taglio sulla dimensione del candidato.':'Events per bin: 223 / 245 / 253 / 37. No candidate-size cut.',
'Soglia dello score e candidati di fondo':'Score threshold and background candidates',
'H − μ: scansione su 323 mappe di fondo':'H − μ: scan of 323 background maps',
'Score minimo':'Minimum score','Signal associati':'Matched signal','Cluster nel fondo':'Background clusters',
'Da 0,90 a 0,99: 61% di cluster di fondo in meno e 0,75 punti di efficienza in meno.':'From 0.90 to 0.99: 61% fewer background clusters, with a 0.75 percentage-point efficiency loss.',
'Normalizzazione e taglio sulla dimensione':'Normalization and the candidate-size cut',
'Stesso benchmark, score ≥0,90':'Same benchmark, score ≥0.90',
'Il taglio sulla dimensione ha un costo elevato':'The size cut substantially reduces efficiency',
'Con H − μ e almeno 8 voxel: efficienza 84,9% e 124 cluster di fondo.':'H − μ with at least 8 voxels: 84.9% efficiency and 124 background clusters.',
'Ispezione di un candidato nel fondo':'Inspecting a background candidate',
'Cella 271: score 0,99989, θ CNN 8,7 mrad, θ fit 9,5 mrad. Candidato da validare.':'Cell 271: score 0.99989, CNN θ 8.7 mrad, fitted θ 9.5 mrad. Candidate awaiting validation.',
'Stato del progetto e prossimi passi':'Project status and next steps',
'Pipeline operativa':'Operational pipeline',
'Produzione dei dataset, training, scansione e diagnostica dei candidati sono disponibili.':'Dataset production, training, scanning and candidate diagnostics are available.',
'Risultato principale':'Main result',
'Nel benchmark H − μ lo scan recupera il 98.7% del segnale sopra 10 mrad a score 0,90.':'On the H − μ benchmark, the scan finds 98.7% of signal above 10 mrad at score 0.90.',
'Validazione da completare':'Remaining validation',
'Robustezza tra brick e livelli di fondo, piccoli angoli e verifica fisica dei candidati.':'Robustness across bricks and background levels, small angles and candidate physics validation.',
'Prossimo passo proposto: fissare il punto operativo su validazione indipendente e dati reali.':'Proposed next step: select the operating point using independent validation and real data.'
}
for a,b in r.items():
 if a not in s: print('MISSING',a)
 s=s.replace(a,b)
s=s.replace("v.toFixed(n).replace('.',',')","v.toFixed(n)")
s=s.replace("visualization/presentations/project_overview/inputs/qa_02_signal_hdf5_010602.png", "visualization/presentations/project_overview/en_plots/augmentation.png")
old='runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/xz_yz_new_vs_gap5_10/rank_001_event_00271_cell_271_score_0.999890_xz_yz.png'
s=s.replace(old,'visualization/presentations/project_overview/en_plots/candidate.png')
s=s.replace("const B=path.join(W,'visualization/presentations/project_overview')", "const B=path.join(W,'visualization/presentations/project_overview/en')")
s=s.replace("path.join(B,'data.json')","path.join(W,'visualization/presentations/project_overview/data.json')")
s=s.replace('Sintesi_progetto_CNN_SND_v2.pptx','CNN_SND_Project_Overview_EN.pptx')
notes=[
'Project summary based on the local crop_bw repository for SND@LHC nue analysis. Snapshot: 15 September 2026. Main sources: README.md, CNN_DATALOADER.md, cnn21d/model.py, the fitlt5 configurations dated 10 September and their full-run reports. Reported performance is an internal benchmark on the available samples.',
'Sources: README.md, ROOT_TO_HDF5.md, CNN_DATALOADER.md, scanning/scan_cnn21d_volumes.py and scanning/scan_cnn21d_loop2.sh. The pipeline identifies signal-like structures, estimates their direction and clusters adjacent windows into candidates for inspection.',
'Sources: README.md, data_preparation/root_to_hdf5.py, cnn_dataset_signal_p24_l5_fitlt5_hard.selection_summary.json and the fitlt5 full/summary.json reports. The stored dataset contains 15,872 samples: 8,656 hard negatives and 7,216 signal samples, computed by subtraction. The actual training/validation/test subsets contain 10,186/2,377/1,996 samples, after selection and class-fraction sampling. These subset counts need not sum to the stored dataset size. The producer now writes 32x32 crops, superseding an older 20x20 description in ROOT_TO_HDF5.md.',
'Sources: cnn21d/model.py, cnn21d/losses.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml, and the full-run summary. Each block applies Conv3D(1,3,3), GroupNorm, SiLU, Conv3D(3,1,1), GroupNorm, SiLU and MaxPool. Channel counts: 16,32,64,128. Global features concatenate average and maximum pooling. Each head has a hidden layer of 64 units. Parameter count: 196,243.',
'Sources: CNN_DATALOADER.md, training/cnn_dataset.py and training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml. The QA figure illustrates augmentation using HDF5 row 10602 from cnn_dataset.h5; it is not a new fitlt5 performance result. Maximum jitter: 5 pixels. Rotations and reflections transform slopes consistently. Training layer drop: 12% for one layer and 3% for two consecutive layers. Hard negatives are selected using a coherent fit below 5 mrad. Signal above 8 mrad enters binary classification; low-angle signal is excluded from that loss but retained for regression. Classification checkpoints maximize AUPRC. Regression checkpoints minimize angular MAE balanced across angle bins.',
'Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/summary.json. AUPRC, precision and recall use each classification checkpoint at score 0.5. Angular MAE uses each regression checkpoint on 879 signal samples. Classification includes 613 positive signal samples and 1,117 hard negatives. The remaining 266 low-angle signal samples in the test set are excluded from binary classification. MAE means mean absolute error and differs from residual standard deviation. H is the raw count, mu the scalar background estimate, MPV(z) a local background mode for each layer, and sigma(z) its local scale. A small numerical epsilon is omitted from the displayed denominator.',
'Sources: scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh, scanning/produce_data_scan_shards.py, candidates/export_background_scan_candidates.py. Recent script configuration: fitlt5_mu_residual, crop size 20, stride 10, score threshold 0.90, minimum regression score 0.90, and separate classification/regression checkpoints. Window clusters are scan candidates. The CNN score is not a calibrated probability of a neutrino event.',
'Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/signal_scan_summary.json, operating point 0.90. MC truth matching uses recomputed_propagation_limit_40. Angle bins: 5-10,10-20,20-50,50-100 mrad. Event counts: 223,245,253,37. Matched detections: 114,240,252,36. There are no events at or above 100 mrad. The combined efficiency above 10 mrad is 528/535. No size cut is applied.',
'Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/score_size_operating_points/summary.json. Efficiency is truth-matched signal divided by 535 signal events above 10 mrad. Background counts are clusters on 323 maps and 116,603 windows, rather than a per-window probability or a physically normalized rate. Reduction from 184 to 72 clusters: 60.87%. Efficiency loss: 4/535, or 0.75 percentage points.',
'Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/scan_t090/score_size_operating_points/summary.json. The comparison uses the same 535 signal events above 10 mrad and 323 background maps. The size cut requires largest_component_voxels>=8. These internal benchmark results do not establish a final choice for every brick.',
'Source: the rank-1 candidate in runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/gifs_new_vs_gap5_10/manifest.json, with its original XZ/YZ projection manifest. Event271, cell271, candidate0 in background brick b000021. Exact score:0.9998899698. CNN theta:8.7095mrad. Geometric fit theta:9.5240mrad. The plot was regenerated from the same source volume with English labels and the original smoothing and q-hot settings. This background-sample candidate is not a confirmed physics signal.',
'Synthesis of the local project sources. Implemented pipeline, crop AUPRC around0.999 and scan efficiency528/535 above10mrad at score0.90 for H-mu. Proposed next steps: test robustness across bricks and background regimes, select an operating point on independent validation, audit truth matching and unmatched candidates, and perform manual/fit-based physics validation on real data. The source reports do not establish final physics purity or a signal observation.'
]
import json
insert='\nconst englishNotes='+json.dumps(notes)+';\np.slides.items.forEach((s,i)=>s.speakerNotes.textFrame.setText(englishNotes[i]));\n'
s=s.replace("await fs.writeFile(path.join(B,'deck.json')",insert+"await fs.writeFile(path.join(B,'deck.json')")
Path('visualization/presentations/project_overview/en').mkdir(exist_ok=True)
Path('visualization/presentations/project_overview/build_en.mjs').write_text(s)
