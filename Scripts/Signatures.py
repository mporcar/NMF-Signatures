from os.path import dirname, abspath
from os.path import join
from warnings import warn
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

from multiprocessing import Pool
from functools import partial

import pandas as pd
import nimfa
from collections import defaultdict, Counter
import urllib

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch

from matplotlib.pyplot import savefig, imshow, set_cmap



def reorder(C):
    """
    Reorder consensus matrix.
    :param C: Consensus matrix.
    :type C: `numpy.ndarray`
    """
    Y = 1 - C
    Z = linkage(squareform(Y), method='average')
    ivl = leaves_list(Z)
    ivl = ivl[::-1]
    return {'M':C[:, ivl][ivl, :],'leaves':ivl, 'Z':Z}

def coph_corr(D):
    """
    Given a consensus matrix D returns the cophenetic coefficient
    """
    # upper diagonal elements of consensus
    avec = np.array([D[i, j] for i in range(D.shape[0] - 1)
                    for j in range(i + 1, D.shape[1])])
    # consensus entries are similarities, conversion to distances
    Y = 1 - avec
    Z = linkage(Y, method='average')
        # cophenetic correlation coefficient of a hierarchical clustering
        # defined by the linkage matrix Z and matrix Y from which Z was
        # generated
    return cophenet(Z, Y)[0]

def dispersion(D):
    dispersion = np.sum(4 * np.multiply(D - 0.5, D - 0.5))/D.size
    return dispersion

def consensus(D):
    rank=D['rank']
    nruns=D['nruns']
    n=D[0]['coef'].shape[1]
    consensus = np.zeros((n, n))
    for i in range(nruns):
        consensus+=D[i]['summary']['connectivity']
    consensus/=nruns
    return consensus
    
    
def silhouette(D):
    rank=D['rank']
    nruns=D['nruns']
    nmfres = pd.DataFrame([])
    for i in range(nruns):
        A=pd.DataFrame(D[i]['basis'])
        nmfres=pd.concat([nmfres,A],axis=1)
    model = KMeans(n_clusters=rank).fit(nmfres.T)
    labels = model.labels_
    General=metrics.silhouette_score(nmfres.T, labels, metric='cosine')
    samples= metrics.silhouette_samples(nmfres.T, labels, metric='cosine')
    return {'AVG': General, 'samples': samples, 'labels':labels}

def get_processes(D):
    rank=D['rank']
    nruns=D['nruns']
    nmfres = pd.DataFrame([])
    for i in range(nruns):
        A=pd.DataFrame(D[i]['basis'])
        nmfres=pd.concat([nmfres,A],axis=1)
    model = KMeans(n_clusters=rank).fit(nmfres.T)
    processes = model.cluster_centers_.T
    return processes



def run_one_object_method_C(df, rank, num, num_core,type_nmf):
    """
    Aquesta funció fará un diccionari complet amb les dades de una k en concret,
    calcularà totes les metriques, aixi quan tinguem els resultats ja no fara falta fer més calculs.
    """
    results = {}
    l = list()
    step=num//num_core
    for i in range(num_core):
        l.append(range(i*step,(i+1)*step))
    f = partial(one_run, df=df, rank=rank, type_nmf=type_nmf)
    with Pool(num_core) as pool:
        for D in pool.map(f, l):
            results.update(D)
    # results['samples']=df.index.values
    # results['types']=df.columns.values
    results['nruns']=num
    results['rank']=rank
    summary=eval_nmf(results)
    results['summary']=summary
    results['type_nmf']=type_nmf
    return {'results': results, 'summary':summary }

def get_metrics(D):
    nruns = D['nruns']
    M={}
    euclidean = 0
    evar=0
    kl = 0
    rss = 0
    sparseness_0=0
    sparseness_1=0

    for i in range(nruns):
        euclidean += D[i]['summary']['euclidean']
        evar += D[i]['summary']['evar']
        kl += D[i]['summary']['kl']
        rss += D[i]['summary']['rss']
        sparseness_0 += D[i]['summary']['sparseness'][0]
        sparseness_1 += D[i]['summary']['sparseness'][1]
    return {'euclidean':euclidean/nruns, 'evar':evar/nruns, 'kl' : kl/nruns, 'rss': rss/nruns, 'sparseness': [sparseness_0/nruns,sparseness_1/nruns]}
       

def eval_nmf(D):
    summary={}
    summary['consensus']=consensus(D)
    summary['dispersion']=dispersion(summary['consensus'])
    summary['coph_corr']=coph_corr(summary['consensus'])
    summary['silhouette']=silhouette(D)
    summary['processes']=get_processes(D)
    summary.update(get_metrics(D))
    return summary

def one_run(iterador, df, rank,type_nmf):
    V=np.array(df)
    V=np.transpose(V)
    D={}
    for i in iterador:
        d={}
        if(type_nmf=="nmf"):
            nmf = nimfa.Nmf(V, rank=rank, seed="random_vcol", max_iter=1000000, update='divergence', objective='conn', conn_change=40)
        if(type_nmf=="snmf"):
            nmf = nimfa.Snmf(V, rank=rank, seed="random_vcol", max_iter=1000000, conn_change=40, version = 'l')
        if(type_nmf=="nsnmf"):
            nmf = nimfa.Nsnmf(V, rank=rank, seed="random_vcol", max_iter=1000000, objective='conn', conn_change=40)
        fit = nmf()
        S=fit.summary()
        SS={}
        SS['connectivity']=S['connectivity']
        SS['euclidean']=S['euclidean']
        SS['evar']=S['evar']
        SS['kl']=S['kl']
        SS['rss']=S['rss']
        SS['sparseness']=S['sparseness']
        d['summary']=SS
        d['n_iter']=fit.n_iter
        d['distance']=fit.distance
        H=pd.DataFrame(fit.basis())
        E=extract_norm(H)
        d['basis']=E['P']
        d['coef']=pd.DataFrame(E['R']*fit.coef())
        D[i]=d
    return D

            
def clean_axis(ax):
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
        
def extract_norm(P):
    diag=P.sum(0)
    n=P.shape[1]
    R=np.zeros(n*n).reshape(n,n)
    for i in range(n):
        P.ix[:,i]=P.ix[:,i]/diag[i]
        R[i,i]=diag[i]
    D={'P':P,'R':R}
    return D

def heatplot_consensus(C):  
    fig = plt.figure(figsize=(13.9, 10))
    heatmapGS = gridspec.GridSpec(1, 2, wspace=.01, hspace=0., width_ratios=[0.25,1])

    C = 1 - C
    Y = sch.linkage(C, method='average')

    denAX = fig.add_subplot(heatmapGS[0,0])
    denD = sch.dendrogram(Y, orientation='right', link_color_func=lambda k: 'black')
    clean_axis(denAX)

    heatmapAX = fig.add_subplot(heatmapGS[0,1])
    D = C[denD['leaves'], :][:, denD['leaves']]
    axi = heatmapAX.imshow(D, interpolation='nearest', aspect='equal', origin='lower', cmap='RdBu') 
    clean_axis(heatmapAX)

    cb = fig.colorbar(axi, fraction=0.046, pad=0.04, aspect=10) 
    cb.set_label('Distance', fontsize=20)