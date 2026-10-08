import fs from 'node:fs/promises';
import path from 'node:path';
import { Presentation, PresentationFile } from '@oai/artifact-tool';
import { resolvePresentationFont, applyPresentationChartFont, finalizePresentation } from '/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations/container_tools/artifact_tool_utils.mjs';
const W='/Users/fabioali/SND@LHC/nue_analysis/CNN';
const B=path.join(W,'visualization/presentations/hmu_review/slides'), O=path.join(W,'output/hmu_review_en');
await fs.mkdir(B,{recursive:true});
const A=JSON.parse(await fs.readFile(path.join(O,'analysis_summary.json'),'utf8'));
const E=JSON.parse(await fs.readFile(path.join(O,'candidate_examples_manifest.json'),'utf8'));
const rateLines=(await fs.readFile(path.join(O,'candidate_rates_by_brick.csv'),'utf8')).trim().split(/\r?\n/);
const rateHeader=rateLines[0].split(',');
const rateRows=rateLines.slice(1).map(line=>Object.fromEntries(rateHeader.map((key,index)=>[key,line.split(',')[index]])));
const S='/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/presentations';
const D=JSON.parse(await fs.readFile(path.join(W,'visualization/presentations/project_overview/data.json'),'utf8'));
const family=resolvePresentationFont({fontFamily:'Arial'}); console.log('font',family);
const p=Presentation.create({slideSize:{width:1280,height:720}});
const C={ink:'#102A43',blue:'#126782',muted:'#526777',white:'#FFFFFF',pale:'#EDF3F6',orange:'#C56730'};
function text(s,t,x,y,w,h,size=28,color=C.ink,bold=false){const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});a.text=t;a.text.style={typeface:family,fontSize:size,color,bold,autoFit:'none'};return a;}
function slide(title,notes){const s=p.slides.add();s.background.fill=C.white;text(s,title,64,42,1152,100,43,C.ink,true);text(s,String(p.slides.items.length).padStart(2,'0'),1160,674,55,26,17,C.muted);s.speakerNotes.textFrame.setText(notes);return s;}
function foot(s,t){text(s,t,64,624,1100,50,19,C.muted);}
function para(s,head,body,y,x=64,w=1090){text(s,head,x,y,w,43,28,C.blue,true);text(s,body,x,y+46,w,78,25);}
async function img(s,rel,x,y,w,h){s.images.add({blob:new Uint8Array(await fs.readFile(path.join(W,rel))),contentType:'image/png',alt:rel,fit:'contain',position:{left:x,top:y,width:w,height:h}});}
function table(s,values,widths,y=190,h=330){const t=s.tables.add({rows:values.length,columns:values[0].length,left:64,top:y,width:1152,height:h,columnWidths:widths,values}); t.borders.assign({fill:'#FFFFFF',width:1}); for(let r=0;r<values.length;r++){for(let c=0;c<values[r].length;c++){const a=t.getCell(r,c);a.fill=r===0?C.ink:(r%2?C.pale:C.white);a.text.style={typeface:family,fontSize:r===0?23:25,color:r===0?C.white:C.ink,bold:r===0};}}return t;}
const src=n=>`Fonti: runs/cnn21d_signal_p24_l5/fitlt5_${n}/full/summary.json e scan_t090/score_size_operating_points/summary.json.`;
const fmt=(v,n=1)=>v.toFixed(n);
// 1
{
const s=p.slides.add();s.background.fill=C.ink;text(s,'SND@LHC',72,70,1000,45,27,'#8DCCDF',true);text(s,'Shower identification\nwith a 2+1D CNN',72,209,1120,210,66,C.white,true);text(s,'The crop_bw project\nMethod and results overview',76,466,1080,90,30,'#D4E5EF');text(s,'15 September 2026',76,637,800,32,22,'#B7CEDC');s.speakerNotes.textFrame.setText('Sintesi del progetto presente nella cartella crop_bw, analisi nue SND@LHC. Stato dei file locali al 15 September 2026. Riferimenti principali: README.md, CNN_DATALOADER.md, cnn21d/model.py, configurazioni fitlt5 del 10 settembre e relativi report full. Le prestazioni mostrate sono valutazioni interne sui campioni disponibili.');}
// 2
{
const s=slide('Objective and analysis workflow','Fonti: README.md, ROOT_TO_HDF5.md, CNN_DATALOADER.md, scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh. Il task riconosce strutture di segnale e stima la direzione, poi raggruppa finestre adiacenti in candidati da ispezionare.');
text(s,'Find shower candidates in volumes and estimate their direction',64,151,1100,75,34,C.blue,true);
para(s,'01   ROOT maps','Raw counts for each layer and a background estimate.',260);
para(s,'02   Dataset and CNN','HDF5 crops, signal classification and slope regression.',374);
para(s,'03   Volume scanning','Overlapping windows, candidate clustering and XZ/YZ projections.',488);
}
// 3
{
const s=slide('Dataset and quality checks','Fonti: README.md, data_preparation/root_to_hdf5.py, cnn_dataset_signal_p24_l5_fitlt5_hard.selection_summary.json, configurazioni fitlt5 e summary.json. Dataset fisico: 15872 esempi di cui 8656 hard e 7216 signal (differenza). I subset effettivamente usati dal training sono 10186/2377/1996 e non esauriscono necessariamente il dataset per via dei filtri e delle frazioni di classe. Non sommare i due conteggi come se avessero la stessa selezione. La vecchia descrizione 20x20 in ROOT_TO_HDF5.md è superata dal produttore corrente 32x32.');
text(s,'57 × 32 × 32',64,175,580,85,64,C.blue,true);text(s,'raw volume per sample',68,261,570,45,28);
text(s,'57 × 20 × 20',680,175,530,85,64,C.blue,true);text(s,'crop at the network input',684,261,530,45,28);
text(s,'15,872 samples in the fitlt5 dataset\n8,656 hard negatives and 7,216 signal samples',64,358,1140,95,32,C.ink,true);
text(s,'Splits grouped by signal event or background cell.\nCell 58 excluded and boundary support checked.\nOriginal counts are retained without smoothing.',64,479,1140,119,26);
foot(s,'Actual train / validation / test subsets: 10,186 / 2,377 / 1,996 samples');
}
// 4
{
const s=slide('A two-head 2+1D CNN','Fonti: cnn21d/model.py, cnn21d/losses.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml, full/summary.json. Ogni blocco applica Conv3D(1,3,3), GroupNorm, SiLU, Conv3D(3,1,1), GroupNorm, SiLU e MaxPool. I canali sono 16,32,64,128. Le feature globali concatenano media e massimo. Le teste hanno layer nascosto di dimensione 64. Parametri: 196243.');
text(s,'Separate convolutions in the spatial and z directions',64,155,1120,50,32,C.blue,true);
text(s,'1 × 3 × 3',65,242,500,80,57,C.ink,true);text(s,'Structure within each XY layer',69,326,510,80,28);
text(s,'3 × 1 × 1',686,242,500,80,57,C.ink,true);text(s,'Evolution across consecutive layers',690,326,510,80,28);
para(s,'Two complementary outputs','Signal score and direction components (sx, sy).',446);
foot(s,'4 blocks, channels 16 / 32 / 64 / 128, 196,243 parameters');
}
// 5
{
const s=slide('Training and consistent augmentation','Fonti: CNN_DATALOADER.md, training/cnn_dataset.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml. Figura QA preesistente, illustrativa delle trasformazioni, non un nuovo risultato del run fitlt5. Jitter massimo 5 pixel, rotazioni di 90 gradi e riflessioni coerenti con slope. Layer-drop train: probabilità 12% un piano, 3% due consecutivi. Hard selezionati con fit coerente theta<5 mrad. Nei run fitlt5 la BCE usa signal sopra 8 mrad e hard negativi, esclude signal a piccolo angolo dalla BCE ma li include nella regressione. Checkpoint classificazione sulla AUPRC e regressione sul MAE bilanciato per bin angolari.');
await img(s,'visualization/presentations/project_overview/en_plots/augmentation.png',55,160,758,403);
text(s,'Jitter of ±5 pixels\nRotations and reflections\nSlopes transformed together',858,190,352,155,26,C.ink,true);
text(s,'Signal above 8 mrad for classification. All signal samples for regression.',858,374,352,132,25);
foot(s,'Separate classification and regression checkpoints. Augmentation during training only.');
}
// 6
{
const s=slide('Performance on test crops',src('mu_residual')+' '+src('mpv_z_residual')+' '+src('mpv_z_sigma')+' AUPRC, precision e recall dal checkpoint di classificazione a soglia 0.5. MAE theta dal checkpoint di regressione su 879 signal. La classificazione è valutata su 613 positivi e 1117 hard: 266 signal a basso angolo del test sono esclusi dalla BCE. Il MAE è errore assoluto medio, distinto dalla deviazione standard del residuo.');
text(s,'Comparison of the three fitlt5 normalization methods',64,150,1130,50,29,C.blue);
const names=['mu_residual','mpv_z_residual','mpv_z_sigma'];
table(s,[['CNN input','AUPRC','Precision','Recall','MAE θ'],...names.map((n,i)=>[[ 'H − μ','H − MPV(z)','[H − MPV(z)] / σ(z)'][i],fmt(D[n].test.auprc,4),fmt(D[n].test.precision*100)+'%',fmt(D[n].test.recall*100)+'%',fmt(D[n].reg.theta_mae_mrad,2)+' mrad'])],[395,175,190,180,212],233,295);
foot(s,'Score ≥0.50. Classification: 613 signal and 1,117 hard. Regression: 879 signal samples.');
}

// 7. Exact requested diagram as a source image.
{
 const s=p.slides.add();s.background.fill=C.white;
 await img(s,'output/hmu_review_en/scan_stride_schematic_EN.png',0,0,1280,720);
 s.speakerNotes.textFrame.setText('Source: scanning/scan_cnn21d_volumes.py and current scan settings. Crop20, stride10, all57layers. For a200x200map the number of positions per axis is floor((200-20)/10)+1=19. Adjacent selected windows form8-neighbor connected components.');
}
// 8
{
const s=slide('Scan efficiency by angle',src('mu_residual')+' Source data: output/hmu_review_en/signal_angular_residuals.csv, recomputed from signal_scan_test_t090.h5 with MC matching and propagation limit40. The bins are8-10,10-20,20-50 and50-100mrad. Event counts are78,245,253 and37; detected candidates are72,240,252 and36. The aggregate efficiency uses strict theta>8mrad:600/613. No candidate-size cut.');
text(s,'H − μ, score ≥0.90, matched to MC truth',64,148,1150,46,27,C.blue);
const chart=s.charts.add('bar',{position:{left:60,top:219,width:830,height:359},categories:['8–10','10–20','20–50','50–100'],series:[{name:'Efficiency',values:[92.3,98.0,99.6,97.3],fill:C.blue}],barOptions:{direction:'column',grouping:'clustered',gapWidth:95},hasLegend:false,xAxis:{title:'True θ [mrad]',textStyle:{fontSize:21,fill:C.ink}},yAxis:{min:0,max:110,majorUnit:25,numberFormatCode:'0.0"%"',textStyle:{fontSize:20},majorGridlines:{fill:'#DCE6EC',width:1}},dataLabels:{showValue:true,position:'outEnd',textStyle:{fontSize:23,bold:true,fill:C.ink}}}); applyPresentationChartFont(chart,{fontFamily:family});
text(s,'97.9%',940,256,280,80,58,C.blue,true);text(s,'600 / 613 events\nwith θ >8 mrad',944,350,270,105,28);text(s,'The 8–10 mrad bin\nremains the limiting range.',944,493,270,90,25);
foot(s,'Events per bin: 78 / 245 / 253 / 37. No candidate-size cut.');
}

// 9. Background operating points.
{
const ss=slide('H − μ scan on b21 and b24','Sources: Hmu score_size_operating_points/summary.json for b21 and signal. b24: /Users/fabioali/cernbox/CNN/background_scan_predictions/b000024/fitlt5_mu_residual/summary_t090.json. b24 counts independently recomputed from all13 prediction HDF5 files. Each background sample has323maps and116603windows. Signal uses535events with true theta>10mrad, matched to MC truth. No size cut. b24 means background brick b000024, separate from real-data brick b000224.');
text(ss,'323 maps and 116,603 windows in each background sample',64,148,1150,52,28,C.blue);
table(ss,[['Minimum score','Signal efficiency','b21 candidates','b24 candidates'],...A.b24_scan.map((r,i)=>[r.score_threshold.toFixed(2),(D.mu_residual.op.rows[i].signal_truth_efficiency*100).toFixed(2)+'%',String(D.mu_residual.op.rows[i].background_clusters),String(r.candidate_clusters)])],[250,330,286,286],229,315);
foot(ss,'At score 0.90, b24 has 126 clusters in 106 cells. At score 0.99, clusters decrease by 50.8%.');
}
// 10. Mean background_mu in every brick represented by the converted TXT files.
{
const ss=p.slides.add();ss.background.fill=C.white;
await img(ss,'output/hmu_review_en/mean_bkg_mu_by_brick_EN.png',0,0,1280,720);
ss.speakerNotes.textFrame.setText('Source: output/hmu_review_en/bkg_mu_by_brick_summary.csv. The mean uses every valid scalar background_mu value stored for the data-scan cells in all24bricks with converted candidate TXT files. The x-axis lists every brick together. Blue bars correspond to the Hmu run and orange bars to gap5_10_mu_high. The two files with a conflicted-copy name in b000431 and b000541 are not HDF5 files and are excluded. This plot describes the local background estimate; it is not a candidate rate or a classifier-score distribution.');
}
// 11. Measured angular resolution.
{
const ss=p.slides.add();ss.background.fill=C.white;
await img(ss,'output/hmu_review_en/signal_angular_resolution_EN.png',0,0,1280,720);
ss.speakerNotes.textFrame.setText('Source: output/hmu_review_en/signal_angular_residuals.csv, recomputed from scan_t090/signal_scan_test_t090.h5 and the original signal_scan_volumes HDF5 files. Truth matching uses propagation limit40. The maximum-score representative across truth-matched clusters supplies one direction per detected event. Score>=0.90 and true theta>10mrad, no size cut. 528detections out of535events. Residual sigma5.096584mrad is the population standard deviation, MAE2.816283mrad and bias-0.377704mrad. This full-volume result is distinct from regression-test crops: sigma3.324683mrad, MAE2.234690mrad on879signal crops.');
}
// 12. Candidate rate from the converted TXT files and valid HDF5 cells.
{
const hmuRates=rateRows.filter(row=>row.model==='fitlt5_mu_residual');
const otherRates=rateRows.filter(row=>row.model!=='fitlt5_mu_residual');
const cell=(row)=>row?[row.brick,`${row.candidates} / ${row.cells}`,`${Number(row.candidates_per_100_cells).toFixed(2)}%`]:['','',''];
const rows=Array.from({length:9},(_,i)=>[...cell(hmuRates[i]),...cell(hmuRates[i+9]),...cell(otherRates[i])]);
const ss=slide('Candidate rate by brick','Sources: output/hmu_review_en/candidate_rates_by_brick.csv and bkg_mu_by_brick_summary.csv. The numerator is the number of selected numeric records in the canonical converted TXT file. The denominator is the number of cells across valid source HDF5 files for the same brick. The rate is candidates divided by cells. Hmu rows use fitlt5_mu_residual. Orange-section rows use gap5_10_mu_high. The legacy b121_candidates.txt is excluded to avoid counting the same brick twice.');
text(ss,'201 candidate records across 24 bricks',64,145,1140,42,30,C.blue,true);
text(ss,'H − μ run',64,189,370,26,19,C.blue,true);text(ss,'H − μ run',448,189,370,26,19,C.blue,true);text(ss,'gap5_10_mu_high run',832,189,370,26,19,C.orange,true);
const rateTable=table(ss,[['Brick','Cand. / cells','Rate','Brick','Cand. / cells','Rate','Brick','Cand. / cells','Rate'],...rows],[128,128,128,128,128,128,128,128,128],217,340);
for(let r=0;r<10;r++){for(let c=0;c<9;c++){const isOther=c>=6;rateTable.getCell(r,c).text.style={typeface:family,fontSize:r===0?15:17,color:r===0?C.white:C.ink,bold:r===0};if(r===0&&isOther){rateTable.getCell(r,c).fill=C.orange;}}}
foot(ss,'Rate = selected candidates / valid HDF5 cells.\nTwo non-HDF5 conflicted copies in b000431 and b000541 are excluded from the cell count.');
}
// 13–20. Evidence from the source volumes.
for(const c of E){
const ss=p.slides.add();ss.background.fill=C.white;
await img(ss,path.relative(W,c.projection),30,90,1220,516);
const label=c.category==='signal'?'Signal MC with known direction':c.category==='data'?'Candidate from the selected data TXT files':'Candidate in the background sample';
text(ss,label,64,29,1130,44,31,C.blue,true);
text(ss,'Layer animation: '+c.stem+'_EN.gif',64,640,1130,38,19,C.muted);
ss.speakerNotes.textFrame.setText(`Source volume: ${c.volume_path}, index${c.event_index}, event${c.event_id}, cell${c.cell_id}. Source manifest: ${c.source_manifest||'recomputed signal truth matching'}. XZ/YZ display projections use ROOT TH2 smoothing and positive excess above the Poisson-count threshold alpha=1e-4. These are display transformations. CNN input uses raw H-mu. Projection display crop${c.display_crop_size}pixels at x${c.display_x_start},y${c.display_y_start}. Score${c.presence_score}. English GIF: candidates/${c.stem}_EN.gif,57layers,160msperframe with a fixed color scale within each animation. Examples illustrate morphology and are not a representative sample for measuring rates. Background/data candidates require physics validation.`);
}
// 20. Closing summary.
{
const ss=slide('H − μ results and remaining validation','Sources: Hmu full/summary.json, scan score_size_operating_points/summary.json, output/hmu_review_en/analysis_summary.json and candidate_counts_by_brick.csv. All scan efficiencies use MC truth matching. Candidate counts on background maps and selected TXT lists have different denominators and selection stages.');
para(ss,'Strong crop and scan performance','AUPRC 0.99924 on test crops. The scan finds 528 / 535 signal events above 10 mrad.',160);
para(ss,'b24 background and selected data lists','126 b24 candidates at score 0.90. The H − μ TXT lists contain 177 selected candidates.',298);
para(ss,'Validation priorities','Low-angle efficiency, the angular-error tails and candidate validation across data bricks.',436);
foot(ss,'Figures and all eight layer animations are available in the accompanying English asset folder.');
}

const englishNotes=["Project summary based on the local crop_bw repository for SND@LHC nue analysis. Snapshot: 15 September 2026. Main sources: README.md, CNN_DATALOADER.md, cnn21d/model.py, the fitlt5 configurations dated 10 September and their full-run reports. Reported performance is an internal benchmark on the available samples.", "Sources: README.md, ROOT_TO_HDF5.md, CNN_DATALOADER.md, scanning/scan_cnn21d_volumes.py and scanning/scan_cnn21d_loop2.sh. The pipeline identifies signal-like structures, estimates their direction and clusters adjacent windows into candidates for inspection.", "Sources: README.md, data_preparation/root_to_hdf5.py, cnn_dataset_signal_p24_l5_fitlt5_hard.selection_summary.json and the fitlt5 full/summary.json reports. The stored dataset contains 15,872 samples: 8,656 hard negatives and 7,216 signal samples, computed by subtraction. The actual training/validation/test subsets contain 10,186/2,377/1,996 samples, after selection and class-fraction sampling. These subset counts need not sum to the stored dataset size. The producer now writes 32x32 crops, superseding an older 20x20 description in ROOT_TO_HDF5.md.", "Sources: cnn21d/model.py, cnn21d/losses.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml, and the full-run summary. Each block applies Conv3D(1,3,3), GroupNorm, SiLU, Conv3D(3,1,1), GroupNorm, SiLU and MaxPool. Channel counts: 16,32,64,128. Global features concatenate average and maximum pooling. Each head has a hidden layer of 64 units. Parameter count: 196,243.", "Sources: CNN_DATALOADER.md, training/cnn_dataset.py and training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml. The QA figure illustrates augmentation using HDF5 row 10602 from cnn_dataset.h5; it is not a new fitlt5 performance result. Maximum jitter: 5 pixels. Rotations and reflections transform slopes consistently. Training layer drop: 12% for one layer and 3% for two consecutive layers. Hard negatives are selected using a coherent fit below 5 mrad. Signal above 8 mrad enters binary classification; low-angle signal is excluded from that loss but retained for regression. Classification checkpoints maximize AUPRC. Regression checkpoints minimize angular MAE balanced across angle bins.", "Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/summary.json. AUPRC, precision and recall use each classification checkpoint at score 0.5. Angular MAE uses each regression checkpoint on 879 signal samples. Classification includes 613 positive signal samples and 1,117 hard negatives. The remaining 266 low-angle signal samples in the test set are excluded from binary classification. MAE means mean absolute error and differs from residual standard deviation. H is the raw count, mu the scalar background estimate, MPV(z) a local background mode for each layer, and sigma(z) its local scale. A small numerical epsilon is omitted from the displayed denominator.", "Sources: scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh, scanning/produce_data_scan_shards.py, candidates/export_background_scan_candidates.py. Recent script configuration: fitlt5_mu_residual, crop size 20, stride 10, score threshold 0.90, minimum regression score 0.90, and separate classification/regression checkpoints. Window clusters are scan candidates. The CNN score is not a calibrated probability of a neutrino event.", "Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/signal_scan_summary.json, operating point 0.90. MC truth matching uses recomputed_propagation_limit_40. Angle bins: 5-10,10-20,20-50,50-100 mrad. Event counts: 223,245,253,37. Matched detections: 114,240,252,36. There are no events at or above 100 mrad. The combined efficiency above 10 mrad is 528/535. No size cut is applied.", "Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/score_size_operating_points/summary.json. Efficiency is truth-matched signal divided by 535 signal events above 10 mrad. Background counts are clusters on 323 maps and 116,603 windows, rather than a per-window probability or a physically normalized rate. Reduction from 184 to 72 clusters: 60.87%. Efficiency loss: 4/535, or 0.75 percentage points.", "Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/scan_t090/score_size_operating_points/summary.json. The comparison uses the same 535 signal events above 10 mrad and 323 background maps. The size cut requires largest_component_voxels>=8. These internal benchmark results do not establish a final choice for every brick.", "Source: the rank-1 candidate in runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/gifs_new_vs_gap5_10/manifest.json, with its original XZ/YZ projection manifest. Event271, cell271, candidate0 in background brick b000021. Exact score:0.9998899698. CNN theta:8.7095mrad. Geometric fit theta:9.5240mrad. The plot was regenerated from the same source volume with English labels and the original smoothing and q-hot settings. This background-sample candidate is not a confirmed physics signal.", "Synthesis of the local project sources. Implemented pipeline, crop AUPRC around0.999 and scan efficiency528/535 above10mrad at score0.90 for H-mu. Proposed next steps: test robustness across bricks and background regimes, select an operating point on independent validation, audit truth matching and unmatched candidates, and perform manual/fit-based physics validation on real data. The source reports do not establish final physics purity or a signal observation."];
p.slides.items.slice(0,6).forEach((s,i)=>s.speakerNotes.textFrame.setText(englishNotes[i]));
p.slides.items[7].speakerNotes.textFrame.setText('Source: output/hmu_review_en/signal_angular_residuals.csv, recomputed from scan_t090/signal_scan_test_t090.h5 with MC truth matching and propagation limit40. The displayed bins are8-10,10-20,20-50,50-100mrad. Event counts are78,245,253,37; detected candidates72,240,252,36. The aggregate result uses strict true theta>8mrad:600/613=97.8793%. No candidate-size cut.');
await fs.writeFile(path.join(B,'deck.json'),JSON.stringify(p.toProto()));
await (await PresentationFile.exportPptx(p)).save(path.join(B,'candidate.pptx'));
for(let i=0;i<p.slides.items.length;i++){const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});await fs.writeFile(path.join(B,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await b.arrayBuffer()));}
const result=await finalizePresentation({workspaceDir:W,candidatePath:path.join(B,'candidate.pptx'),finalPath:path.join(O,'CNN_SND_Hmu_Review_EN_v4.pptx'),pythonExecutable:'/Users/fabioali/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',integrityValidatorPath:path.join(S,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(S,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','6','--require-native-table-slide','9','--require-native-table-slide','12'],fontPolicy:{basis:'design',families:[family]},explicitTotalSlideCount:21,requiredNativeTableOwnerSlides:[6,9,12],requiredNativeChartOwnerSlides:[8],materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,receiptPath:path.join(B,'validation_v4.json')});
console.log(result);
