import fs from 'node:fs/promises';
import path from 'node:path';
import { Presentation, PresentationFile } from '@oai/artifact-tool';
import { resolvePresentationFont, applyPresentationChartFont, finalizePresentation } from '/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations/container_tools/artifact_tool_utils.mjs';
const W='/Users/fabioali/SND@LHC/nue_analysis/CNN';
const B=path.join(W,'visualization/presentations/project_overview/en'), O=path.join(W,'output/presentazione_progetto');
await fs.mkdir(B, {recursive:true});
const S='/Users/fabioali/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations';
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
// 7
{
const s=slide('Scanning full volumes','Fonti: scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh, scanning/produce_data_scan_shards.py, candidates/export_background_scan_candidates.py. Configurazione script recente: fitlt5_mu_residual, crop20, stride10, threshold0.90, regressione minima score0.90, due checkpoint separati. I cluster di finestre sono candidati dello scan. Lo score CNN non rappresenta una probabilità fisica calibrata di evento neutrino.');
para(s,'20 × 20 windows, 10-pixel stride','The scan retains all 57 layers and evaluates overlapping windows.',160);
para(s,'Score and direction','The classification checkpoint selects windows. The regression checkpoint estimates slopes.',295);
para(s,'Candidates for inspection','Window clustering produces candidates with coordinates, scores and XZ/YZ projections.',430);
foot(s,'Recent data scan configuration: H − μ, score ≥0.90. Candidates require physics validation.');
}
// 8
{
const s=slide('Scan efficiency by angle',src('mu_residual')+' Fonte aggiuntiva: scan_t090/signal_scan_summary.json, operating_points threshold0.90. Matching con verità MC, recomputed_propagation_limit_40. I bin sono 5-10, 10-20, 20-50, 50-100 mrad con N=223,245,253,37, trovati=114,240,252,36. Nessun evento >=100mrad. Efficiency aggregata theta>10mrad=528/535. Nessun taglio di size.');
text(s,'H − μ, score ≥0.90, matched to MC truth',64,148,1150,46,27,C.blue);
const chart=s.charts.add('bar',{position:{left:60,top:219,width:830,height:359},categories:['5–10','10–20','20–50','50–100'],series:[{name:'Efficiency',values:[51.1,98.0,99.6,97.3],fill:C.blue}],barOptions:{direction:'column',grouping:'clustered',gapWidth:95},hasLegend:false,xAxis:{title:'True θ [mrad]',textStyle:{fontSize:21,fill:C.ink}},yAxis:{min:0,max:110,majorUnit:25,numberFormatCode:'0.0"%"',textStyle:{fontSize:20},majorGridlines:{fill:'#DCE6EC',width:1}},dataLabels:{showValue:true,position:'outEnd',textStyle:{fontSize:23,bold:true,fill:C.ink}}}); applyPresentationChartFont(chart,{fontFamily:family});
text(s,'98.7%',940,256,280,80,58,C.blue,true);text(s,'528 / 535 events\nwith θ >10 mrad',944,350,270,105,28);text(s,'The main limitation is\nthe 5–10 mrad range.',944,493,270,90,25);
foot(s,'Events per bin: 223 / 245 / 253 / 37. No candidate-size cut.');
}
// 9
{
const s=slide('Score threshold and background candidates',src('mu_residual')+' Dati: score_size_operating_points/summary.json. Efficiency truth-matched su 535 signal theta>10 mrad. Il conteggio fondo è il numero di cluster su 323 mappe e 116603 finestre, non una probabilità per finestra o un tasso fisico normalizzato. La riduzione 184 a72 è60.87%.');
text(s,'H − μ: scan of 323 background maps',64,150,1150,48,28,C.blue);
table(s,[['Minimum score','Matched signal','Efficiency','Background clusters'],...D.mu_residual.op.rows.map(r=>[fmt(r.threshold,2),r.signal_truth_matched+' / 535',fmt(r.signal_truth_efficiency*100,2)+'%',String(r.background_clusters)])],[250,315,247,340],221,316);
foot(s,'From 0.90 to 0.99: 61% fewer background clusters, with a 0.75 percentage-point efficiency loss.');
}
// 10
{
const s=slide('Normalization and the candidate-size cut',src('mu_residual')+' '+src('mpv_z_residual')+' '+src('mpv_z_sigma')+' Stesso campione di535 signal theta>10 e323 mappe background. La size è largest_component_voxels>=8. I risultati suggeriscono un compromesso interno al benchmark, non una scelta definitiva per i dati di tutti i brick.');
text(s,'Same benchmark, score ≥0.90',64,149,1150,50,28,C.blue);
table(s,[['CNN input','Efficiency','Background clusters'],...['mu_residual','mpv_z_residual','mpv_z_sigma'].map((n,i)=>[['H − μ','H − MPV(z)','[H − MPV(z)] / σ(z)'][i],fmt(D[n].op.rows[0].signal_truth_efficiency*100,2)+'%',String(D[n].op.rows[0].background_clusters)])],[535,295,322],219,259);
text(s,'The size cut substantially reduces efficiency',64,511,1130,43,31,C.blue,true);
text(s,'H − μ with at least 8 voxels: 84.9% efficiency and 124 background clusters.',64,563,1130,66,27);
}
// 11
{
const s=slide('Inspecting a background candidate','Fonte figura: visualization/presentations/project_overview/en_plots/candidate.png. La figura identifica hard_b000021_hmu_new_vs_gap5_10, rank1, event271 cell271. Score esatto0.999890 dal nome file. Theta CNN8.7mrad e fit9.5mrad dal titolo. Questo esempio illustra l’ispezione di un candidato nel campione background e non costituisce conferma di segnale fisico.');
await img(s,'visualization/presentations/project_overview/en_plots/candidate.png',50,164,1170,445);
foot(s,'Cell 271: score 0.99989, CNN θ 8.7 mrad, fitted θ 9.5 mrad. Candidate awaiting validation.');
}
// 12
{
const s=slide('Project status and next steps','Sintesi basata sui file locali. Risultati: pipeline implementata, test crop AUPRC~0.999, scan H-mu theta>10 528/535 a0.90. Prossimi passi proposti dall’analisi: valutare robustezza per brick e regime di fondo, fissare il punto operativo su validation indipendente, audit truth matching e candidati non associati, validazione manuale/fit su dati reali. Non vi sono in queste fonti una misura finale di purezza fisica sui dati né un risultato di osservazione.');
para(s,'Operational pipeline','Dataset production, training, scanning and candidate diagnostics are available.',162);
para(s,'Main result','On the H − μ benchmark, the scan finds 98.7% of signal above 10 mrad at score 0.90.',296);
para(s,'Remaining validation','Robustness across bricks and background levels, small angles and candidate physics validation.',430);
foot(s,'Proposed next step: select the operating point using independent validation and real data.');
}

const englishNotes=["Project summary based on the local crop_bw repository for SND@LHC nue analysis. Snapshot: 15 September 2026. Main sources: README.md, CNN_DATALOADER.md, cnn21d/model.py, the fitlt5 configurations dated 10 September and their full-run reports. Reported performance is an internal benchmark on the available samples.", "Sources: README.md, ROOT_TO_HDF5.md, CNN_DATALOADER.md, scanning/scan_cnn21d_volumes.py and scanning/scan_cnn21d_loop2.sh. The pipeline identifies signal-like structures, estimates their direction and clusters adjacent windows into candidates for inspection.", "Sources: README.md, data_preparation/root_to_hdf5.py, cnn_dataset_signal_p24_l5_fitlt5_hard.selection_summary.json and the fitlt5 full/summary.json reports. The stored dataset contains 15,872 samples: 8,656 hard negatives and 7,216 signal samples, computed by subtraction. The actual training/validation/test subsets contain 10,186/2,377/1,996 samples, after selection and class-fraction sampling. These subset counts need not sum to the stored dataset size. The producer now writes 32x32 crops, superseding an older 20x20 description in ROOT_TO_HDF5.md.", "Sources: cnn21d/model.py, cnn21d/losses.py, training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml, and the full-run summary. Each block applies Conv3D(1,3,3), GroupNorm, SiLU, Conv3D(3,1,1), GroupNorm, SiLU and MaxPool. Channel counts: 16,32,64,128. Global features concatenate average and maximum pooling. Each head has a hidden layer of 64 units. Parameter count: 196,243.", "Sources: CNN_DATALOADER.md, training/cnn_dataset.py and training/configs/cnn21d_signal_p24_l5_fitlt5_mu_residual.yaml. The QA figure illustrates augmentation using HDF5 row 10602 from cnn_dataset.h5; it is not a new fitlt5 performance result. Maximum jitter: 5 pixels. Rotations and reflections transform slopes consistently. Training layer drop: 12% for one layer and 3% for two consecutive layers. Hard negatives are selected using a coherent fit below 5 mrad. Signal above 8 mrad enters binary classification; low-angle signal is excluded from that loss but retained for regression. Classification checkpoints maximize AUPRC. Regression checkpoints minimize angular MAE balanced across angle bins.", "Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/summary.json. AUPRC, precision and recall use each classification checkpoint at score 0.5. Angular MAE uses each regression checkpoint on 879 signal samples. Classification includes 613 positive signal samples and 1,117 hard negatives. The remaining 266 low-angle signal samples in the test set are excluded from binary classification. MAE means mean absolute error and differs from residual standard deviation. H is the raw count, mu the scalar background estimate, MPV(z) a local background mode for each layer, and sigma(z) its local scale. A small numerical epsilon is omitted from the displayed denominator.", "Sources: scanning/scan_cnn21d_volumes.py, scanning/scan_cnn21d_loop2.sh, scanning/produce_data_scan_shards.py, candidates/export_background_scan_candidates.py. Recent script configuration: fitlt5_mu_residual, crop size 20, stride 10, score threshold 0.90, minimum regression score 0.90, and separate classification/regression checkpoints. Window clusters are scan candidates. The CNN score is not a calibrated probability of a neutrino event.", "Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/signal_scan_summary.json, operating point 0.90. MC truth matching uses recomputed_propagation_limit_40. Angle bins: 5-10,10-20,20-50,50-100 mrad. Event counts: 223,245,253,37. Matched detections: 114,240,252,36. There are no events at or above 100 mrad. The combined efficiency above 10 mrad is 528/535. No size cut is applied.", "Source: runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/score_size_operating_points/summary.json. Efficiency is truth-matched signal divided by 535 signal events above 10 mrad. Background counts are clusters on 323 maps and 116,603 windows, rather than a per-window probability or a physically normalized rate. Reduction from 184 to 72 clusters: 60.87%. Efficiency loss: 4/535, or 0.75 percentage points.", "Sources: runs/cnn21d_signal_p24_l5/fitlt5_{mu_residual,mpv_z_residual,mpv_z_sigma}/full/scan_t090/score_size_operating_points/summary.json. The comparison uses the same 535 signal events above 10 mrad and 323 background maps. The size cut requires largest_component_voxels>=8. These internal benchmark results do not establish a final choice for every brick.", "Source: the rank-1 candidate in runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/gifs_new_vs_gap5_10/manifest.json, with its original XZ/YZ projection manifest. Event271, cell271, candidate0 in background brick b000021. Exact score:0.9998899698. CNN theta:8.7095mrad. Geometric fit theta:9.5240mrad. The plot was regenerated from the same source volume with English labels and the original smoothing and q-hot settings. This background-sample candidate is not a confirmed physics signal.", "Synthesis of the local project sources. Implemented pipeline, crop AUPRC around0.999 and scan efficiency528/535 above10mrad at score0.90 for H-mu. Proposed next steps: test robustness across bricks and background regimes, select an operating point on independent validation, audit truth matching and unmatched candidates, and perform manual/fit-based physics validation on real data. The source reports do not establish final physics purity or a signal observation."];
p.slides.items.forEach((s,i)=>s.speakerNotes.textFrame.setText(englishNotes[i]));
await fs.writeFile(path.join(B,'deck.json'),JSON.stringify(p.toProto()));
await (await PresentationFile.exportPptx(p)).save(path.join(B,'candidate.pptx'));
for(let i=0;i<p.slides.items.length;i++){const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});await fs.writeFile(path.join(B,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await b.arrayBuffer()));}
const result=await finalizePresentation({workspaceDir:W,candidatePath:path.join(B,'candidate.pptx'),finalPath:path.join(O,'CNN_SND_Project_Overview_EN.pptx'),pythonExecutable:'/Users/fabioali/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',integrityValidatorPath:path.join(S,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(S,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','6','--require-native-table-slide','9','--require-native-table-slide','10'],fontPolicy:{basis:'design',families:[family]},explicitTotalSlideCount:12,requiredNativeTableOwnerSlides:[6,9,10],requiredNativeChartOwnerSlides:[8],materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,receiptPath:path.join(B,'validation.json')});
console.log(result);
