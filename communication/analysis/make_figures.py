"""Generate vector and raster figures for the communication working manuscript.
Run: python communication/analysis/make_figures.py
"""
from pathlib import Path
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE = Path(__file__).resolve().parent
FIG = HERE.parent/'figures'
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})

def save(fig,name):
    fig.savefig(FIG/f'{name}.svg',bbox_inches='tight')
    fig.savefig(FIG/f'{name}.png',bbox_inches='tight')
    plt.close(fig)

def geometry():
    fig,ax=plt.subplots(figsize=(8,5.1))
    R=6371.; L=5000.; y=math.sqrt(R*R-(L/2)**2)
    ax.add_patch(Circle((0,0),R,facecolor='#eaf2f7',edgecolor='#31566b',lw=1.6))
    ax.plot([-L/2,L/2],[y,y],color='#ba4d3c',lw=2.5)
    ax.scatter([-L/2,L/2],[y,y],s=75,c=['#234663','#ad5036'],zorder=5)
    ax.annotate('Source',(-L/2,y),(-4300,7100),arrowprops={'arrowstyle':'->'},ha='center')
    ax.annotate('Receiver',(L/2,y),(4300,7100),arrowprops={'arrowstyle':'->'},ha='center')
    ax.text(0,y-1350,'Through-Earth chord  L = 5,000 km',ha='center',color='#9e4132')
    ax.annotate('',xy=(L/2,y+400),xytext=(-L/2,y+400),arrowprops={'arrowstyle':'<->','color':'#ba4d3c'})
    ax.text(0,-4500,'Earth cross-section (to scale)',ha='center',color='#31566b')
    ax.set(xlim=(-6900,6900),ylim=(-6900,8000),aspect='equal')
    ax.axis('off');fig.tight_layout();save(fig,'01_earth_chord')

def benchmark():
    n=np.arange(1,13);mu=.81*n;ber=.5*np.exp(-mu)
    fig,ax=plt.subplots(figsize=(7.5,4.6))
    input_mean=2*1402/3454
    input_se=np.sqrt(1402)/(3454/2)
    ax.fill_between(n,.5*np.exp(-(input_mean+1.96*input_se)*n),
                    .5*np.exp(-(input_mean-1.96*input_se)*n),
                    color='#7aa6bb',alpha=.22,label='Approx. 95% input-mean interval')
    ax.semilogy(n,ber,'o-',lw=2,color='#215d77',label='Ideal Poisson OOK model')
    ax.scatter([1,5,9],[.22318,.00893,.000355],marker='s',s=50,color='#c96140',zorder=5,label='Fixed-seed simulation')
    ax.axhline(.01,ls='--',color='#666',lw=1);ax.annotate('1% raw BER',xy=(11,.01),xytext=(7,.016),arrowprops={'arrowstyle':'->'},color='#555')
    ax.set(xlabel='Independent pulses pooled per raw bit',ylabel='Uncoded bit-error probability',xticks=n,ylim=(1e-5,.4))
    ax.grid(alpha=.25,which='both');ax.legend(frameon=False,loc='upper right');fig.tight_layout();save(fig,'02_pooled_pulses_ber')

def sensitivity():
    l=np.linspace(1000,12000,400)
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    for mass,color in [(1,'#ba4d3c'),(10,'#255e78'),(40,'#467b59')]:
        mw=.9760362708238478*(l/1000)**2*10/mass
        ax.plot(l,mw,lw=2.2,label=f'{mass} kt receiver',color=color)
    ax.set(yscale='log',xlabel='Through-Earth chord length (km)',ylabel='Peak neutrino-carried power during on slots (MW)',xlim=(1000,12000),ylim=(.05,3000))
    ax.grid(alpha=.25,which='both');ax.legend(frameon=False,loc='upper left')
    ax.text(.02,.02,'3 GeV, 1 mrad, 50% selected efficiency; 1 slot/s (1 s slots), 1% BER\nAverage power at equiprobable OOK is half the peak; not electrical power',transform=ax.transAxes,fontsize=8,va='bottom',bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
    fig.tight_layout();save(fig,'03_energy_vs_chord')

def design_grid():
    masses=np.array([1,10,40]);angles=np.array([.1,1,10])
    values=np.array([[24.400906770596194 * 10/m * (angle/1.)**2 for angle in angles] for m in masses])
    fig,ax=plt.subplots(figsize=(7.5,4.2))
    im=ax.imshow(np.log10(values),cmap='viridis',aspect='auto',vmin=-2,vmax=5)
    ax.set(xticks=np.arange(3),xticklabels=['0.1','1','10'],yticks=np.arange(3),yticklabels=['1','10','40'],xlabel='Assumed beam half-angle (mrad)',ylabel='Receiver fiducial mass (kt)')
    for i in range(3):
        for j in range(3):
            ax.text(j,i,f'{values[i,j]:.3g}',ha='center',va='center',color='white' if np.log10(values[i,j])<2 else 'black',weight='bold')
    c=fig.colorbar(im,ax=ax,pad=.02);c.set_label('log10(peak on-slot neutrino-carried power, MW)')
    ax.set_title('5,000 km · 3 GeV · 1 slot/s · 1% ideal BER')
    fig.tight_layout();save(fig,'04_mass_divergence_grid')

def operating_points():
    from tradeoffs import detector_kt
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    power=np.geomspace(.1,10000,300)
    for length,color in [(1000,'#467b59'),(5000,'#255e78'),(12000,'#ba4d3c')]:
        ax.loglog(power,[detector_kt(length,p) for p in power],color=color,lw=2,label=f'{length:,} km chord')
    for p in (1,10,100):
        ax.axvline(p,color='#888',alpha=.25,ls=':')
    ax.set(xlabel='Peak neutrino-carried power during an on slot (MW)',ylabel='Required fiducial detector mass (kt)',xlim=(.1,10000),ylim=(.001,100000))
    ax.grid(alpha=.25,which='both');ax.legend(frameon=False)
    ax.text(.02,.025,'3 GeV · 1 mrad · 50% efficiency · 1 s slot · 1% ideal raw BER\nPower is peak in the on slot; no conversion or background included',transform=ax.transAxes,fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    fig.tight_layout();save(fig,'05_detector_power_tradeoff')

def finance_latency():
    from tradeoffs import required_mw,latency_window_ms
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    power=np.geomspace(1,100000,300)
    delay=np.array([required_mw(5000,10)/p*1000 for p in power])
    ax.loglog(power,delay,color='#255e78',lw=2.5,label='Ideal 1%-BER count-collection time')
    limit=latency_window_ms(5000)
    ax.axhline(limit,color='#ba4d3c',ls='--',lw=2,label=f'Ideal surface-fiber advantage: {limit:.2f} ms')
    ax.fill_between(power,limit,100000,color='#ba4d3c',alpha=.08)
    ax.set(xlabel='Peak neutrino-carried power during an on slot (MW)',ylabel='Ideal on-slot count-collection time (ms)',xlim=(1,100000),ylim=(.1,100000))
    ax.grid(alpha=.25,which='both');ax.legend(frameon=False,loc='upper right')
    ax.text(.02,.025,'5,000 km chord · 10 kt · 1 mrad · 3 GeV\nFiber: shortest surface arc, group index 1.4677; no routing or processing',transform=ax.transAxes,fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    fig.tight_layout();save(fig,'06_ideal_finance_latency')

def payload_budget():
    from throughput import raw_slots_per_s,payload_bits_per_s
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    length=np.linspace(1000,12000,300)
    for power,color in [(1,'#467b59'),(10,'#255e78'),(100,'#ba4d3c')]:
        rates=[payload_bits_per_s(l,10,power) for l in length]
        ax.plot(length,rates,color=color,lw=2,label=f'{power} MW peak on-slot power')
    ax.set(yscale='log',xlabel='Through-Earth chord length (km)',ylabel='Illustrative payload rate (bit/s)',xlim=(1000,12000))
    ax.grid(alpha=.25,which='both');ax.legend(frameon=False)
    ax.text(.02,.025,'10 kt · 1 mrad · 3 GeV · 1% ideal uncoded BER\nPeak power during on slot; assumed FEC 1/2 and frame payload 0.8\nNominal bookkeeping only; no code gain, retries or acquisition modelled',transform=ax.transAxes,fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    fig.tight_layout();save(fig,'07_illustrative_payload_rate')

if __name__=='__main__':
    geometry();benchmark();sensitivity();design_grid();operating_points();finance_latency();payload_budget()
