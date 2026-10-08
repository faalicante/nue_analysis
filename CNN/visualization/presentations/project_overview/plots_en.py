from __future__ import annotations
import ast, json, math, sys
from pathlib import Path
from typing import Any
import h5py, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm,PowerNorm
sys.path.insert(0,str(Path.cwd()))
from visualization.make_background_candidate_gifs import root_th2_smooth,poisson_hot_threshold
out=Path('visualization/presentations/project_overview/en_plots');out.mkdir(exist_ok=True)

def functions_from(file,names,replacements):
 tree=ast.parse(Path(file).read_text())
 code='from __future__ import annotations\n'+'\n\n'.join(ast.unparse(n) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names)
 for a,b in replacements.items(): code=code.replace(a,b)
 exec(compile(code,file,'exec'),globals())
functions_from('visualization/qa_cnn_dataset.py',['_normalization_for','_draw_slope','render_example'],{'Prima: crop centrale, nessuna augmentation':'Before: central crop, no augmentation','Dopo: crop traslato + trasformazione xy':'After: shifted crop + XY transformation'})
CLASS_NAMES={0:'poisson',1:'hard',2:'signal'}
class Tensor:
 def __init__(self,a): self.a=np.asarray(a)
 def __getitem__(self,i): return Tensor(self.a[i])
 def sum(self,dim): return Tensor(self.a.sum(axis=dim))
 def numpy(self): return self.a
with h5py.File('cnn_dataset.h5','r') as f:
 i=10602; raw=f['regions_raw'][i].astype(np.float32); mu=np.float32(np.asarray(f['background_mu'][i]).flat[0]); slope=np.asarray(f['slope_xy'][i])
 before=(raw[:,6:26,6:26]-mu)/np.sqrt(mu+np.float32(1e-6))
 after=(raw[:,7:27,6:26]-mu)/np.sqrt(mu+np.float32(1e-6))
 after=np.rot90(after,k=-2,axes=(-2,-1))[:,:,::-1]
 meta={'sample_type':2,'crop_offset_xy':np.array([0,1]),'rotation_quarter_turns_ccw':2,'reflected_x':True,'reflected_y':False,'positive_visible_fraction':1.013,'positive_may_be_insufficiently_visible':False,'hdf5_index':10602}
 render_example({'volume':Tensor(before[None]),'slope_xy':Tensor(slope)},{'volume':Tensor(after[None]),'slope_xy':Tensor(slope*np.array([1,-1])),'metadata':meta},out/'augmentation.png')
print('QA slope',slope,'mu',mu)
functions_from('visualization/make_candidate_xzyz_projections.py',['visual_qhot_projections','slope_line','render_candidate'],{'inizio stimato':'estimated start','MC vertice originale':'original MC vertex','layer XYPseg (z)':'XYPseg layer (z)','eccesso {display_kind} integrato':'integrated {display_kind} excess','smooth q-hot':'smoothed q-hot'})
Z_STEP_UM=1350.0
manifest=Path('runs/cnn21d_signal_p24_l5/fitlt5_mu_residual/full/scan_t090/background_scan_t090/gifs_new_vs_gap5_10/manifest.json')
c=json.loads(manifest.read_text())['candidates'][0]
r=render_candidate(c,out/'candidate.png','Background brick b000021',smooth_projections=True,projection_poisson_tail_probability=1e-4)
(out/'candidate_metadata.json').write_text(json.dumps(r,indent=2))
print('Candidate',r['theta_fit_mrad'],r['projection_poisson_hot_count_threshold'])
