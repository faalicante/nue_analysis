from __future__ import annotations
import ast,csv,glob,json,math,re,sys
from collections import deque
from pathlib import Path
import h5py,numpy as np
sys.path.insert(0,str(Path.cwd()))
from visualization.make_background_candidate_gifs import centered_window_start,estimate_shower_start,fitted_display_crop_size,poisson_hot_threshold,root_th2_smooth
OUT=Path('output/hmu_review_en'); B=Path('visualization/presentations/hmu_review'); BASE=Path('/Users/fabioali/cernbox/CNN')

def load_functions(file,names):
 tree=ast.parse(Path(file).read_text());s='from __future__ import annotations\n'+'\n'.join(ast.unparse(n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names)
 exec(compile(s,file,'exec'),globals())
load_functions('scanning/scan_cnn21d_volumes.py',['connected_components','cluster_contains_truth'])
load_functions('visualization/make_candidate_xzyz_projections.py',['event_lookup_from_volumes','truth_center_bins','nearest_truth_location','window_contains_truth','signal_candidate_row'])
N_LAYERS=57;PLATE_INDEX_OFFSET=3;Z_STEP_UM=1350.

def readtxt(path):
 rows=[]
 for n,line in enumerate(path.read_text().splitlines(),1):
  if re.match(r'^\s*\d+\s*\*',line):
   fields=[x.strip() for x in line.split('*') if x.strip()];assert len(fields)==6,(path,n,fields)
   rows.append(tuple(float(x) for x in fields))
 return rows
counts=[];txtrows={}
for d in sorted(BASE.glob('data_gifs/b*/*')):
 if not d.is_dir():continue
 brick=d.parent.name;model=d.name;p=d/f'{brick}_candidates.txt'
 if not p.exists():continue
 rows=readtxt(p);txtrows[str(d)]=rows
 counts.append({'brick':brick,'model':model,'candidates':len(rows),'unique_rows':len(set(rows)),'unique_cells':len(set(x[0] for x in rows)),'source_txt':str(p)})
with open(OUT/'candidate_counts_by_brick.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(counts[0]));w.writeheader();w.writerows(counts)
legacy=BASE/'data_gifs/b000121/fitlt5_mu_residual/b121_candidates.txt'
print('COUNTS',json.dumps(counts));print('LEGACY',len(readtxt(legacy)))

# Verify b24 operating points directly from the score grids.
b24=[];files=sorted(BASE.glob('background_scan_predictions/b000024/fitlt5_mu_residual/*.predictions.h5'))
for threshold in [.9,.95,.98,.99]:
 nwin=ncl=ncell=nmaps=totwin=0
 for p in files:
  with h5py.File(p) as f:
   assert f.attrs['normalization_mode']=='background_mu_residual'
   scores=f['window_presence_score'][:];nmaps+=len(scores);totwin+=scores.size
   for sc in scores:
    mask=sc>=threshold;nwin+=int(mask.sum());ncl+=len(connected_components(mask));ncell+=int(mask.any())
 b24.append({'score_threshold':threshold,'candidate_clusters':ncl,'selected_windows':nwin,'cells_with_candidates':ncell,'maps':nmaps,'total_windows':totwin})
print('B24',b24)

# Reproduce the scan truth matching and save one angular prediction per detected event.
lookup=event_lookup_from_volumes(str(BASE/'signal_scan_volumes/signal_scan_*.h5'))
scanpath=Path('runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/signal_scan_test_t090.h5')
records=[];candidate_seeds=[]
with h5py.File(scanpath) as f:
 scores=f['window_presence_score'][:];slopes=f['window_slope_xy'][:];ids=f['event_id'][:];gx=f['grid_x_start_bin'][:];gy=f['grid_y_start_bin'][:];crop=int(f.attrs['crop_size'])
 for j,eid in enumerate(ids):
  vp,k=lookup[int(eid)]
  with h5py.File(vp) as src:
   tx,ty,xum,yum,truth=truth_center_bins(src,k,40)
  matched=[]
  for comp in connected_components(scores[j]>=.9):
   if cluster_contains_truth(comp,tx,ty,gx.tolist(),gy.tolist(),crop):
    rc=max(comp,key=lambda rc:float(scores[j][rc]));matched.append((float(scores[j][rc]),rc,comp))
  r={'event_id':int(eid),'theta_true_mrad':truth['theta_true_mrad'],'detected':bool(matched),'theta_pred_mrad':None,'residual_mrad':None}
  if matched:
   score,rc,comp=max(matched,key=lambda x:x[0]);pred=float(1000*np.arctan(np.hypot(*slopes[j][rc])/27));r.update(theta_pred_mrad=pred,residual_mrad=pred-truth['theta_true_mrad'],score=score)
   if 18<truth['theta_true_mrad']<45 and score>.999:
    candidate_seeds.append({'volume_path':str(vp),'source_index':int(k),'row':int(rc[0]),'col':int(rc[1]),'slope':slopes[j][rc].tolist(),'score':score,'truth':truth,'truth_x_um':xum,'truth_y_um':yum})
  records.append(r)
 # Just 12 promising signal candidates for visual selection.
 candidate_seeds.sort(key=lambda x:abs(x['truth']['theta_true_mrad']-28))
 signals=[]
 for seed in candidate_seeds[:14]:
  with h5py.File(seed['volume_path']) as src:
   c=signal_candidate_row(src,seed['source_index'],Path(seed['volume_path']),gx,gy,crop,40,int(src['volumes_raw'].shape[-1]),seed['row'],seed['col'],presence_score=seed['score'],slope=np.array(seed['slope']),candidate_id=0,kind='found',threshold=.9,truth=seed['truth'],nearest_distance_bins=0.,found_window='cluster-max-score',truth_x_um=seed['truth_x_um'],truth_y_um=seed['truth_y_um'])
   c['rank']=len(signals)+1;c['category']='signal';signals.append(c)
with open(OUT/'signal_angular_residuals.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in records))));w.writeheader();w.writerows(records)
metrics={}
for name,lo,hi in [('all_ge5',5,np.inf),('gt10',10,np.inf),('10_to_20',10,20),('20_to_50',20,50),('50_to_100',50,100)]:
 sub=[r for r in records if lo<=r['theta_true_mrad']<hi];det=[r for r in sub if r['detected']];res=np.array([r['residual_mrad'] for r in det]);q=np.quantile(res,[.16,.84])
 metrics[name]={'events':len(sub),'detected':len(det),'efficiency':len(det)/len(sub),'bias_mrad':float(res.mean()),'sigma_mrad':float(res.std()),'mae_mrad':float(np.abs(res).mean()),'central68_halfwidth_mrad':float((q[1]-q[0])/2),'q16_mrad':float(q[0]),'q84_mrad':float(q[1])}
print('ANGULAR',json.dumps(metrics));assert metrics['gt10']['detected']==528

# Rank data examples from the selected TXT lists, using their stored morphology.
data=[]
for d,rows in txtrows.items():
 if '/fitlt5_mu_residual' not in d:continue
 mp=Path(d)/'manifest.json'
 if not mp.exists():continue
 for c in json.load(open(mp)).get('candidates',[]):
  c=dict(c)
  if c.get('start_x_um') is None:continue
  selected=any(int(r[0])==int(c['cell_id']) and abs(r[4]-float(c['start_x_um']))<.01 and abs(r[5]-float(c['start_y_um']))<.01 for r in rows)
  if selected and c.get('start_quality')=='coherent' and 10<c['theta_pred_mrad']<65:
   c['category']='data';c['canonical_brick']=Path(d).parent.name;c['source_manifest']=str(mp);data.append(c)
data.sort(key=lambda c:c.get('start_component_charge',0),reverse=True)
background=[]
for group in ['common_hmu','hmu_exclusive']:
 mp=BASE/f'background_scan_predictions/b000024/comparison_hmu_hmpv_t090/{group}/gifs/manifest.json'
 for c in json.load(open(mp))['candidates']:
  c=dict(c);c['category']='background_b24';c['canonical_brick']='b000024';c['source_manifest']=str(mp);background.append(c)
background.sort(key=lambda c:c.get('start_component_charge',0),reverse=True)
background21=json.load(open('runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/gifs_new_vs_gap5_10/manifest.json'))['candidates']
for c in background21:c['category']='background_b21';c['canonical_brick']='b000021'
background21.sort(key=lambda c:c.get('start_component_charge',0),reverse=True)
for label,arr in [('SIGNAL',signals),('DATA',data[:8]),('B24C',background[:8]),('B21C',background21[:3])]:
 print(label,[(x.get('canonical_brick'),x['event_id'],x['cell_id'],round(x['theta_pred_mrad'],1),x.get('start_component_voxels'),x.get('start_component_layers'),round(x.get('start_component_charge',0))) for x in arr])
(B/'candidates_pool.json').write_text(json.dumps({'signal':signals,'data':data[:10],'background_b24':background[:8],'background_b21':background21[:3]},indent=2))
summary=json.load(open('runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/summary.json'))
report={'model':'fitlt5_mu_residual/full','normalization':'H - mu','classification_test':summary['classification_checkpoint_test'],'regression_test':summary['regression_checkpoint_test'],'angular_scan_metrics':metrics,'b24_scan':b24,'candidate_counts':counts,'txt_total_hmu':sum(r['candidates'] for r in counts if r['model']=='fitlt5_mu_residual'),'txt_total_other_models':sum(r['candidates'] for r in counts if r['model']!='fitlt5_mu_residual'),'legacy_b121_not_counted':str(legacy),'legacy_b121_rows':len(readtxt(legacy))}
(OUT/'analysis_summary.json').write_text(json.dumps(report,indent=2))
