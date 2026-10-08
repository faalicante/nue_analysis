from __future__ import annotations
import ast,io,json,math,sys
from pathlib import Path
import h5py,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from PIL import Image
sys.path.insert(0,str(Path.cwd()))
from visualization.make_background_candidate_gifs import root_th2_smooth,poisson_hot_threshold
O=Path('output/hmu_review_en');(O/'candidates').mkdir(parents=True,exist_ok=True)
selection=json.load(open('visualization/presentations/hmu_review/extra_data_candidates.json'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12})
Z_STEP_UM=1350.
def load_functions(file,names):
 tree=ast.parse(Path(file).read_text());s='from __future__ import annotations\n'+'\n'.join(ast.unparse(n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names);exec(compile(s,file,'exec'),globals())
load_functions('visualization/make_candidate_xzyz_projections.py',['visual_qhot_projections','slope_line'])
manifest=[]
for c in selection:
 if c['category']=='signal':
  stem=f"signal_event{c['event_id']:06d}";title=f"Signal MC · event {c['event_id']}";truth_label=f" · true θ = {c['theta_true_mrad']:.1f} mrad"
 else:
  stem=f"{c['category'].split('_')[0]}_{c['canonical_brick']}_cell{c['cell_id']:03d}";title=f"{'Data' if c['category']=='data' else 'Background'} · {c['canonical_brick']} · cell {c['cell_id']}";truth_label=''
 with h5py.File(c['volume_path']) as f:
  k=c['event_index'];x=c['display_x_start'];y=c['display_y_start'];size=c['display_crop_size'];raw=f['volumes_raw'][k].astype(float);mu=float(f['background_mu'][k]);xe=f['x_edges_um'][k,x:x+size+1]/1000;ye=f['y_edges_um'][k,y:y+size+1]/1000
 projs,threshold=visual_qhot_projections(raw,mu,x,y,size,smooth=True,poisson_tail_probability=1e-4)
 fig,axs=plt.subplots(1,2,figsize=(14,5.9),constrained_layout=True)
 fig.suptitle(f"{title}\nH − μ score = {c['presence_score']:.5f} · CNN θ = {c['theta_pred_mrad']:.1f} mrad{truth_label}",fontsize=16,color='#102A43')
 starts=[c['start_x_um'],c['start_y_um']];cnn=[c['tx_pred_mrad'],c['ty_pred_mrad']];layers=np.arange(1,58)
 for ax,arr,edges,coord,st,sl in zip(axs,projs,[xe,ye],['x','y'],starts,cnn):
  positive=arr[arr>0];vmax=max(float(np.quantile(positive,.995)),1.) if positive.size else 1
  mesh=ax.pcolormesh(np.arange(.5,58.5),edges,arr.T,cmap='viridis',norm=PowerNorm(gamma=.55,vmin=0,vmax=vmax),shading='flat')
  ax.plot(layers,slope_line(st,sl,c['start_layer'])/1000,color='#FF5245',lw=2,label=f'CNN: {sl:+.1f} mrad')
  fit=c.get('selected_component_slope_'+coord)
  if fit is not None:
   fsl=float(fit)*float(np.median(np.diff(edges)))*1e6/1350
   ax.plot(layers,slope_line(st,fsl,c['start_layer'])/1000,'--',color='white',lw=1.8,label=f'Geometric fit: {fsl:+.1f} mrad')
  if c['category']=='signal':
   true=c['t'+coord+'_true_mrad'];vertex=c['truth_vertex_'+coord+'_um'];z=c['truth_vertex_layer']
   ax.plot(layers,slope_line(vertex,true,z)/1000,':',color='#53EF9F',lw=2,label=f'MC: {true:+.1f} mrad')
  ax.scatter([c['start_layer']],[st/1000],c='#FFCE32',marker='x',s=55,lw=2,label='Estimated start')
  ax.set(xlim=(.5,57.5),ylim=(edges[0],edges[-1]),xlabel='XYPseg layer (z)',ylabel=f'{coord} [mm]',title=coord.upper()+'Z projection');ax.legend(loc='best',fontsize=10,framealpha=.88);ax.grid(alpha=.15)
  fig.colorbar(mesh,ax=ax,pad=.015,label='Integrated smoothed q-hot excess')
 output=O/'candidates'/f'{stem}_projection_EN.png';fig.savefig(output,dpi=170);plt.close(fig)
 cc=dict(c);cc.update(stem=stem,title=title,projection=str(output.resolve()),background_mu=mu,qhot_count_threshold=threshold);manifest.append(cc)
 print('Projection',stem,flush=True)
 if '--gifs' in sys.argv:
  smoothed=root_th2_smooth(raw)[:,y:y+size,x:x+size];vmin=float(np.quantile(smoothed,.02));vmax=max(float(np.quantile(smoothed,.998)),vmin+1)
  fig,ax=plt.subplots(figsize=(7.6,7.3),constrained_layout=True)
  im=ax.imshow(smoothed[0],origin='lower',extent=[xe[0],xe[-1],ye[0],ye[-1]],cmap='viridis',vmin=vmin,vmax=vmax,interpolation='nearest',aspect='equal')
  ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]');fig.colorbar(im,ax=ax,shrink=.85,label='Smoothed counts (fixed scale)')
  tt=ax.set_title('');frames=[]
  best=int(np.argmax(np.maximum(smoothed-mu,0).sum(axis=(1,2))))
  for z in range(57):
   im.set_data(smoothed[z]);tt.set_text(f"{title}\nH − μ score {c['presence_score']:.5f} · CNN θ {c['theta_pred_mrad']:.1f} mrad\nLayer {z+1:02d} / 57");b=io.BytesIO();fig.savefig(b,format='png',dpi=110);b.seek(0);frame=Image.open(b).convert('RGB');
   if z in [0,best,56]:frame.save(O/'candidates'/f'{stem}_frame{z+1:02d}_EN.png')
   frames.append(frame.convert('P',palette=Image.Palette.ADAPTIVE,colors=256));b.close()
  plt.close(fig);gp=O/'candidates'/f'{stem}_EN.gif';frames[0].save(gp,save_all=True,append_images=frames[1:],duration=160,loop=0,optimize=False,disposal=2)
  cc['gif']=str(gp.resolve());cc['gif_frames']=57;cc['gif_best_layer']=best+1;cc['gif_duration_ms']=160;cc['gif_vmin']=vmin;cc['gif_vmax']=vmax
  print('GIF',stem,flush=True)
Path('visualization/presentations/hmu_review/new_data_candidates_manifest.json').write_text(json.dumps(manifest,indent=2))
