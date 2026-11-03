"""Plot proved parameter ranges and the analytic local profile (not simulation)."""
from pathlib import Path
import math
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/entropy-matplotlib-cache')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
out=root/"figures"; out.mkdir(exist_ok=True)
plum="#613D68"; ochre="#B57732"; jade="#3D7B72"; ink="#29252B"; grey="#777078"
paper="#FCFAF7"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                     "axes.spines.right":False,"text.color":ink,"axes.labelcolor":ink,
                     "xtick.color":ink,"ytick.color":ink})
fig=plt.figure(figsize=(8.4,6.5),facecolor=paper)
gs=fig.add_gridspec(2,1,height_ratios=[1.25,3.0],hspace=.65)
ax=fig.add_subplot(gs[0],facecolor=paper)
qc=4/7
ax.plot([0,qc],[0,0],color=jade,lw=14,solid_capstyle="butt")
ax.plot([qc,1],[0,0],color=plum,lw=14,solid_capstyle="butt")
ax.scatter([qc],[0],s=115,facecolor=paper,edgecolor=ink,lw=1.6,zorder=4)
ax.text(.23,.22,"Local concavity",ha="center",color=jade,weight="bold",fontsize=12)
ax.text(.79,.22,"Local minima possible",ha="center",color=plum,weight="bold",fontsize=12)
ax.annotate("The boundary value remains open",xy=(qc,-.02),xytext=(qc,-.43),
            ha="center",fontsize=9,color=grey,arrowprops={"arrowstyle":"-","color":grey})
ax.set_xlim(0,1);ax.set_ylim(-.57,.48)
ax.set_yticks([]);ax.set_xticks([0,qc,1],["0","4/7","1"])
ax.spines["left"].set_visible(False);ax.spines["bottom"].set_visible(False)
ax.tick_params(axis="x",length=0,pad=5)
ax.set_xlabel("Retained signal variance fraction, q",labelpad=7)
ax.set_title("An exact local threshold:  sup C = 4/7",loc="left",weight="bold",fontsize=15,pad=12)

ax2=fig.add_subplot(gs[1],facecolor=paper)
q=3/5;gamma=21*math.sqrt(5)*q/(49*q-4)
s=np.linspace(-1.05,1.05,501)
gain=1-np.exp(-4*gamma*s*s)
inflection=1/math.sqrt(8*gamma)
ax2.axvspan(-inflection,inflection,color=ochre,alpha=.12,zorder=0)
ax2.plot(s,gain,color=plum,lw=2.8)
ax2.scatter([0],[0],color=ochre,s=45,zorder=5)
ax2.axhline(1,color=grey,linestyle=(0,(3,4)),lw=.7,alpha=.7)
ax2.set_xlim(-1.05,1.05);ax2.set_ylim(-.04,1.07)
ax2.set_xticks([-1,-.5,0,.5,1]);ax2.set_yticks([0,.25,.5,.75,1])
ax2.set_xlabel(r"Scaled distance from equal weighting,  $s=\sqrt{n}(t-\frac{1}{2})$",labelpad=9)
ax2.set_ylabel("Entropy gain / midpoint deficit",labelpad=9)
ax2.set_title("The limiting entropy dip above the threshold",loc="left",weight="bold",fontsize=12,pad=10)
ax2.text(.5,.94,r"$q=3/5,\quad 1-e^{-4\gamma(q)s^2}$",transform=ax2.transAxes,
         ha="center",va="top",fontsize=11)
ax2.annotate("Strict local minimum",xy=(0,.012),xytext=(.42,.16),
             arrowprops={"arrowstyle":"-","lw":.8,"color":ochre},
             color=ochre,fontsize=10)
ax2.grid(axis="y",color=grey,alpha=.13,lw=.6)
fig.subplots_adjust(left=.12,right=.96,top=.92,bottom=.19)
fig.text(.12,.025,"Symmetric, variance-one inputs near the Gaussian. The lower curve is an analytic asymptotic profile;\n"
                     "it does not show a measured finite-density entropy effect.",fontsize=8.5,color=grey)
fig.savefig(out/"entropy_transition.png",dpi=190,facecolor=paper)
fig.savefig(out/"entropy_transition.pdf",facecolor=paper)
plt.close(fig)
print("Saved analytic overview figure.")

