from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon
O=Path('output/hmu_review_en')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
blue='#126782';ink='#102A43';orange='#D16C34';purple='#7964B0';gray='#526777'
# Exact geometry of the sliding windows.
fig=plt.figure(figsize=(16,9),facecolor='white')
fig.text(.055,.92,'Sliding-window scan with a 10-pixel stride',fontsize=25,weight='bold',color=ink)
fig.text(.055,.864,'20 × 20 XY crop, all 57 layers, overlapping positions',fontsize=17,color=blue)
ax=fig.add_axes([.065,.24,.52,.53]);ax.set_xlim(-3,47);ax.set_ylim(-13,36);ax.set_aspect('equal')
ax.set_xticks(range(0,41,10));ax.set_yticks(range(0,31,10));ax.grid(color='#DCE4EA',lw=.8,zorder=0);ax.set_xlabel('x [pixels]');ax.set_ylabel('y [pixels]')
for x,col,num in [(0,blue,'1'),(10,orange,'2'),(20,purple,'3')]:
 ax.add_patch(Rectangle((x,0),20,20,facecolor=col,edgecolor=col,alpha=.16,zorder=1))
 ax.add_patch(Rectangle((x,0),20,20,fill=False,edgecolor=col,lw=2.7,zorder=2))
 ax.text(x+3,3,num,color=col,weight='bold',fontsize=19)
ax.add_patch(Rectangle((0,10),20,20,fill=False,edgecolor=blue,linestyle='--',lw=2,zorder=3))
ax.text(24,27,'Next row: y + 10',color=blue,fontsize=13)
ax.annotate('',xy=(20,23),xytext=(0,23),arrowprops={'arrowstyle':'<->','lw':1.7,'color':ink});ax.text(10,24,'crop width = 20',ha='center',fontsize=12,color=ink)
ax.annotate('',xy=(10,-6),xytext=(0,-6),arrowprops={'arrowstyle':'<->','lw':2,'color':orange});ax.text(5,-10,'stride = 10',ha='center',fontsize=13,color=orange)
ax.text(30,-10,'10-pixel overlap',ha='center',fontsize=13,color=gray)
ax.set_title('XY view of the first scan positions (schematic)',loc='left',fontsize=15,color=ink,pad=18)
bx=fig.add_axes([.66,.405,.28,.35]);bx.set_xlim(0,10);bx.set_ylim(0,9);bx.axis('off')
for i in range(6):
 dx=i*.32;dy=i*.43
 bx.add_patch(Polygon([(1+dx,1+dy),(6+dx,1+dy),(7.2+dx,2.1+dy),(2.2+dx,2.1+dy)],facecolor='#EDF3F6',edgecolor=blue,lw=1.5))
bx.text(5,7.9,'One CNN input per XY position',ha='center',color=ink,fontsize=16,weight='bold')
bx.text(5,6.8,'57 × 20 × 20 voxels',ha='center',color=blue,fontsize=18,weight='bold')
bx.text(5,.0,'Same XY crop in every layer',ha='center',color=gray,fontsize=13)
fig.text(.66,.31,'For a 200 × 200 XY map',fontsize=17,weight='bold',color=ink)
fig.text(.66,.258,'19 positions per axis',fontsize=17,color=ink)
fig.text(.66,.207,'19 × 19 = 361 CNN windows',fontsize=17,color=blue,weight='bold')
fig.text(.055,.105,'At each position: classify the crop and estimate its slopes. Cluster adjacent selected windows into candidates.',fontsize=14,color=ink)
fig.text(.055,.055,'The stride moves the window in x and y. The scan uses all 57 layers at every position.',fontsize=14,color=gray)
for ext in ['png','svg','pdf']:fig.savefig(O/f'scan_stride_schematic_EN.{ext}',dpi=180)
plt.close(fig)
# Measured resolution on truth-matched scan candidates.
d=json.load(open(O/'analysis_summary.json'));m=d['angular_scan_metrics']['gt10']
rows=list(csv.DictReader(open(O/'signal_angular_residuals.csv')))
res=np.array([float(r['residual_mrad']) for r in rows if r['detected']=='True' and float(r['theta_true_mrad'])>10])
fig,axes=plt.subplots(1,2,figsize=(16,8.4),gridspec_kw={'width_ratios':[1.5,1]},facecolor='white')
fig.subplots_adjust(left=.075,right=.96,bottom=.19,top=.79,wspace=.29)
fig.suptitle('H − μ: angular resolution on detected signal',x=.06,y=.965,ha='left',fontsize=25,weight='bold',color=ink)
fig.text(.06,.892,'Full-volume scan, score ≥ 0.90, true θ > 10 mrad, matched to MC truth',fontsize=17,color=blue)
a=axes[0];bound=max(15,int(np.ceil(np.max(np.abs(res))/5)*5));edges=np.arange(-bound,bound+1.0,1.0)
a.hist(res,bins=edges,color=blue,alpha=.88,edgecolor='white',linewidth=.45)
a.axvline(0,color=ink,ls='--',lw=1.5,label='Zero residual');a.axvline(m['bias_mrad'],color=orange,lw=2,label='Mean residual')
a.axvspan(m['q16_mrad'],m['q84_mrad'],color=orange,alpha=.12,label='Central 68% interval')
a.set_xlabel(r'$\Delta\theta = \theta_{\mathrm{CNN}}-\theta_{\mathrm{true}}$ [mrad]',fontsize=16);a.set_ylabel('Detected events / 1 mrad',fontsize=15);a.grid(axis='y',alpha=.2);a.legend(loc='upper left',fontsize=11)
a.text(.97,.96,f"N = {len(res)}\nσ = {m['sigma_mrad']:.2f} mrad\nBias = {m['bias_mrad']:+.2f} mrad\nMAE = {m['mae_mrad']:.2f} mrad",transform=a.transAxes,ha='right',va='top',fontsize=14,color=ink,bbox={'facecolor':'white','edgecolor':'none','alpha':.92})
b=axes[1];names=['10_to_20','20_to_50','50_to_100'];mm=[d['angular_scan_metrics'][x] for x in names];x=np.arange(3)
b.plot(x,[v['sigma_mrad'] for v in mm],'o-',color=blue,lw=2.2,ms=9,label='Residual standard deviation')
b.plot(x,[v['mae_mrad'] for v in mm],'s--',color=orange,lw=1.8,ms=7,label='Mean absolute error')
for xi,v in zip(x,mm):b.annotate(f"N = {v['detected']}",(xi,v['sigma_mrad']),xytext=(0,15),textcoords='offset points',ha='center',fontsize=12,color=ink)
b.set_xticks(x,['10–20','20–50','50–100']);b.set_xlabel('True θ [mrad]',fontsize=16);b.set_ylabel('Angular error [mrad]',fontsize=15);b.set_ylim(0,max(v['sigma_mrad'] for v in mm)*1.35);b.grid(alpha=.2);b.legend(loc='upper left',fontsize=11)
fig.text(.06,.105,'Resolution σ is the population standard deviation of the residuals, not a Gaussian-fit width.',fontsize=13,color=gray)
fig.text(.06,.062,f"One maximum-score representative per matched event. Detection efficiency: {m['detected']}/{m['events']} = {m['efficiency']*100:.2f}%. No size cut.",fontsize=13,color=gray)
for ext in ['png','pdf','svg']:fig.savefig(O/f'signal_angular_resolution_EN.{ext}',dpi=180)
plt.close(fig)
print('Created scan schematic and angular-resolution figures')
