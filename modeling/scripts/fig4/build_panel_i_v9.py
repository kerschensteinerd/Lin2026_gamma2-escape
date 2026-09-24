"""Draw only Figure 4I using new 25%-retained outputs; F-H never regenerated."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'fig3'))
from build_figure_v3 import GRAY,ORANGE,SUMMARY_MARKER_SIZE_PT,SUMMARY_MARKER_EDGE_PT,SUMMARY_ERROR_WIDTH_PT
PROJECT=Path(__file__).resolve().parents[3]
DATA=PROJECT/'vg3_cid3_rf_simulation-main/outputs/fig4_v9_retained_inhibition_20260916/nominal25'
OUT=PROJECT/'figures/main/assets/Fig4_v9'
W,H=135*72/25.4-16,384.

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    model=json.loads((DATA/'figure_data.json').read_text())
    experiment_path=PROJECT/'analysis/fig3_model/experimental_fig3.json'
    data=json.loads(experiment_path.read_text())
    plt.rcParams.update({'font.family':'Arial','font.size':7.3,'axes.labelsize':7.3,'xtick.labelsize':6.4,'ytick.labelsize':6.4,'axes.linewidth':.5,'axes.spines.top':False,'axes.spines.right':False,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
    fig=plt.figure(figsize=(W/72,H/72))
    def text(x,top,label,**kwargs):return fig.text(x/W,1-top/H,label,va='top',**kwargs)
    text(2*W/3,256,'I',fontsize=10,fontweight='bold')
    text(14+2*W/3,256.4,'OFF comparison',fontsize=8)
    text(33+2*W/3,266,'Lines: fixed-parameter model',fontsize=5.5,color='#555555')
    text(33+2*W/3,272,'Points: imaging',fontsize=5.5,color='#555555')
    ax=fig.add_axes([(33+2*W/3)/W,1-(280+80)/H,82/W,80/H])
    for key,genotype,color in [('ctrl','Ctrl',GRAY),('reduced','cKO',ORANGE)]:
        curve=model['contrast'][key]
        ax.plot(curve['contrast_pct'],curve['auc_shape'],color=color,lw=1.05)
        group=data['groups'][genotype+'_OFF']['contrasts']
        means=[group[str(c)]['normalized_mean'] for c in (.15,.3,.7)]
        sems=[group[str(c)]['normalized_sem'] for c in (.15,.3,.7)]
        ax.errorbar([15,30,70],means,yerr=sems,fmt='o',color=color,mfc=color,ms=SUMMARY_MARKER_SIZE_PT,mew=SUMMARY_MARKER_EDGE_PT,elinewidth=SUMMARY_ERROR_WIDTH_PT,capsize=0,clip_on=False,zorder=5)
    ax.set(xlim=(12,73),ylim=(-.03,1.2),xticks=[15,30,70],yticks=[0,.5,1],xlabel='Contrast (%)',ylabel='Normalized AUC')
    ax.xaxis.labelpad=ax.yaxis.labelpad=1.5
    for ext in ('pdf','svg','png'):fig.savefig(OUT/f'Fig4_panel_I.{ext}',dpi=300)
    plt.close(fig)
    (OUT/'panel_I_manifest.json').write_text(json.dumps({'model_source_sha256':hashlib.sha256((DATA/'figure_data.json').read_bytes()).hexdigest(),'imaging_source_sha256':hashlib.sha256(experiment_path.read_bytes()).hexdigest(),'imaging_points_unchanged':True,'model_retained_effective_inhibition':.25,'calcium_parameters_fixed':True,'new_fit_performed':False,'source_canvas_points':[W,H],'panel_box_points':[2*W/3,256,W/3,128],'label':'Lines: fixed-parameter model'},indent=2)+'\n')
    print(OUT/'Fig4_panel_I.pdf')
if __name__=='__main__':main()
