from pathlib import Path
s=Path('visualization/presentations/project_overview/build_en.mjs').read_text()
s=s.replace("26.905.11957","26.909.12148").replace("const B=path.join(W,'visualization/presentations/project_overview/en'), O=path.join(W,'output/presentazione_progetto');","const B=path.join(W,'visualization/presentations/hmu_review/slides'), O=path.join(W,'output/hmu_review_en');\nawait fs.mkdir(B,{recursive:true});\nconst A=JSON.parse(await fs.readFile(path.join(O,'analysis_summary.json'),'utf8'));\nconst E=JSON.parse(await fs.readFile(path.join(O,'candidate_examples_manifest.json'),'utf8'));")
preamble=s[:s.index('// 7')]
# Preserve source slides and their original English speaker notes.
notes=s[s.index('const englishNotes='):s.index('p.slides.items.forEach')]
extra=r'''
// 7. Exact requested diagram as a source image.
{
 const s=p.slides.add();s.background.fill=C.white;
 await img(s,'output/hmu_review_en/scan_stride_schematic_EN.png',0,0,1280,720);
 s.speakerNotes.textFrame.setText('Source: scanning/scan_cnn21d_volumes.py and current scan settings. Crop20, stride10, all57layers. For a200x200map the number of positions per axis is floor((200-20)/10)+1=19. Adjacent selected windows form8-neighbor connected components.');
}
'''
old8=s[s.index('// 8'):s.index('// 9')]
extra+=old8+ r'''
// 9. Background operating points.
{
const ss=slide('H − μ scan on b21 and b24','Sources: Hmu score_size_operating_points/summary.json for b21 and signal. b24: /Users/fabioali/cernbox/CNN/background_scan_predictions/b000024/fitlt5_mu_residual/summary_t090.json. b24 counts independently recomputed from all13 prediction HDF5 files. Each background sample has323maps and116603windows. Signal uses535events with true theta>10mrad, matched to MC truth. No size cut. b24 means background brick b000024, separate from real-data brick b000224.');
text(ss,'323 maps and 116,603 windows in each background sample',64,148,1150,52,28,C.blue);
table(ss,[['Minimum score','Signal efficiency','b21 candidates','b24 candidates'],...A.b24_scan.map((r,i)=>[r.score_threshold.toFixed(2),(D.mu_residual.op.rows[i].signal_truth_efficiency*100).toFixed(2)+'%',String(D.mu_residual.op.rows[i].background_clusters),String(r.candidate_clusters)])],[250,330,286,286],229,315);
foot(ss,'At score 0.90, b24 has 126 clusters in 106 cells. At score 0.99, clusters decrease by 50.8%.');
}
// 10. Measured angular resolution.
{
const ss=p.slides.add();ss.background.fill=C.white;
await img(ss,'output/hmu_review_en/signal_angular_resolution_EN.png',0,0,1280,720);
ss.speakerNotes.textFrame.setText('Source: output/hmu_review_en/signal_angular_residuals.csv, recomputed from scan_t090/signal_scan_test_t090.h5 and the original signal_scan_volumes HDF5 files. Truth matching uses propagation limit40. The maximum-score representative across truth-matched clusters supplies one direction per detected event. Score>=0.90 and true theta>10mrad, no size cut. 528detections out of535events. Residual sigma5.096584mrad is the population standard deviation, MAE2.816283mrad and bias-0.377704mrad. This full-volume result is distinct from regression-test crops: sigma3.324683mrad, MAE2.234690mrad on879signal crops.');
}
// 11. Hmu candidates in converted lists.
{
const rows=A.candidate_counts.filter(r=>r.model==='fitlt5_mu_residual');
const ss=slide('H − μ candidates in the converted TXT files','Source: output/hmu_review_en/candidate_counts_by_brick.csv. Each row lists a canonical bNNNNNN_candidates.txt under /Users/fabioali/cernbox/CNN/data_gifs/bNNNNNN/fitlt5_mu_residual. Count one numeric record per candidate, excluding headers and separators. No duplicate records within any canonical TXT. The legacy b121_candidates.txt contains30rows and is excluded to avoid counting the same brick twice. These lists contain selected candidates exported from the converted files and do not enumerate all automatic scan candidates.');
text(ss,'177 selected candidates across 18 bricks',64,149,1140,52,32,C.blue,true);
table(ss,[['Brick','Candidates','Brick','Candidates'],...rows.slice(0,9).map((r,i)=>[r.brick,String(r.candidates),rows[i+9].brick,String(rows[i+9].candidates)])],[340,236,340,236],217,378);
foot(ss,'Canonical TXT files only. The legacy b121_candidates.txt is excluded from the total.');
}
// 12. Other converted candidates.
{
const rows=A.candidate_counts.filter(r=>r.model!=='fitlt5_mu_residual');
const ss=slide('Other converted TXT files','Source: output/hmu_review_en/candidate_counts_by_brick.csv. These six bricks use gap5_10_mu_high, which is a different model from fitlt5_mu_residual. Keep the two model counts distinct. Hmu177 plus other model24 equals201selected candidate records across24bricks. Counts use canonical directory names. Some original manifests use legacy brick identifiers.');
text(ss,'gap5_10_mu_high: 24 candidates across 6 bricks',64,149,1130,52,30,C.blue,true);
table(ss,[['Brick','Candidates'],...rows.map(r=>[r.brick,String(r.candidates)])],[700,452],217,337);
text(ss,'All converted lists: 177 + 24 = 201 candidates',64,590,1120,48,31,C.ink,true);
}
// 13–18. Evidence from the source volumes.
for(const c of E){
const ss=p.slides.add();ss.background.fill=C.white;
await img(ss,path.relative(W,c.projection),30,90,1220,516);
const label=c.category==='signal'?'Signal MC with known direction':c.category==='data'?'Candidate from the selected data TXT files':'Candidate in the background sample';
text(ss,label,64,29,1130,44,31,C.blue,true);
text(ss,'Layer animation: '+c.stem+'_EN.gif',64,640,1130,38,19,C.muted);
ss.speakerNotes.textFrame.setText(`Source volume: ${c.volume_path}, index${c.event_index}, event${c.event_id}, cell${c.cell_id}. Source manifest: ${c.source_manifest||'recomputed signal truth matching'}. XZ/YZ display projections use ROOT TH2 smoothing and positive excess above the Poisson-count threshold alpha=1e-4. These are display transformations. CNN input uses raw H-mu. Projection display crop${c.display_crop_size}pixels at x${c.display_x_start},y${c.display_y_start}. Score${c.presence_score}. English GIF: candidates/${c.stem}_EN.gif,57layers,160msperframe with a fixed color scale within each animation. Examples illustrate morphology and are not a representative sample for measuring rates. Background/data candidates require physics validation.`);
}
// 19. Closing summary.
{
const ss=slide('H − μ results and remaining validation','Sources: Hmu full/summary.json, scan score_size_operating_points/summary.json, output/hmu_review_en/analysis_summary.json and candidate_counts_by_brick.csv. All scan efficiencies use MC truth matching. Candidate counts on background maps and selected TXT lists have different denominators and selection stages.');
para(ss,'Strong crop and scan performance','AUPRC 0.99924 on test crops. The scan finds 528 / 535 signal events above 10 mrad.',160);
para(ss,'b24 background and selected data lists','126 b24 candidates at score 0.90. The H − μ TXT lists contain 177 selected candidates.',298);
para(ss,'Validation priorities','Low-angle efficiency, the angular-error tails and candidate validation across data bricks.',436);
foot(ss,'Figures and all six layer animations are available in the accompanying English asset folder.');
}
'''
finish=s[s.index('await fs.writeFile(path.join(B,\'deck.json\')'):]
finish=finish.replace('CNN_SND_Project_Overview_EN.pptx','CNN_SND_Hmu_Review_EN.pptx').replace("'--require-native-table-slide','10'","'--require-native-table-slide','11','--require-native-table-slide','12'").replace('explicitTotalSlideCount:12','explicitTotalSlideCount:19').replace('requiredNativeTableOwnerSlides:[6,9,10]','requiredNativeTableOwnerSlides:[6,9,11,12]')
# Chart data is intentionally copied from the existing builder's literal measured values.
result=preamble+extra+'\n'+notes+"p.slides.items.slice(0,6).forEach((s,i)=>s.speakerNotes.textFrame.setText(englishNotes[i]));\np.slides.items[7].speakerNotes.textFrame.setText(englishNotes[7]);\n"+finish
Path('visualization/presentations/hmu_review/build_deck.mjs').write_text(result)
