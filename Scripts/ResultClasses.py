import pandas as pd
import nimfa
from collections import defaultdict, Counter
import urllib
import numpy as np
from Signatures import *
from transform import *
from distances import *
from normalize import *
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch
from sklearn.cluster import KMeans
import seaborn as sns
import plotly.plotly as py
import dill as pickle
import scipy.spatial.distance as spdist
from plotly.tools import FigureFactory as FF
import plotly.graph_objs as go
import math

def minor_tick_labels():
    major_labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    flanks = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 
              'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    minor_labels = []
    for subs in major_labels:
        for flank in flanks:
            minor_labels.append(flank[0]+subs[0]+flank[1])
    return minor_labels

def plot_signature(signature, plot_title, ax, ymax=0.3):
    """
    Args:
        signature: signature-like object in lexicographic order
        plot_title: string
        ymax: float
    Returns:
        produces the signature bar plot
    """
    vector = np.array([signature[k] for k in sorted(signature)])
    barlist = ax.bar(range(96), vector)
    color_list = ['#72bcd4', 'k', 'r', '#7e7e7e', 'g', '#e6add8']
    for category in range(6):
        for i in range(16):
            barlist[category * 16 + i].set_color(color_list[category])
    ax.set_title(plot_title, weight='bold')
    ax.set_xlim([0, 96])
    ax.set_ylim([0, ymax])
    labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    major_ticks = np.arange(8, 8 + 16 * 5 + 1, 16)
    minor_ticks = np.arange(0.5, 96.5, 1)
    ax.tick_params(length=0, which='major', pad=15, labelsize=12)
    ax.tick_params(length=0, which='minor', pad=4, labelsize=4)
    ax.set_xticks(major_ticks, minor=False)
    ax.set_xticklabels(labels, minor=False)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(minor_tick_labels(), minor=True, rotation=90)

class summaryNMF(object):
    def __init__(self, data):
        self.data_dict=data
        self.cosmic_signatures = None
        
    def load_cosmic(self):
        file_name = 'data/cosmic_signatures.pickle.gz'
        with gzip.open(file_name) as f:
            self.cosmic_signatures = pickle.load(f)
            
    def load_brca(self):
        file_name = 'data/BRCA_results.pickle.gz'
        with gzip.open(file_name) as f:
            self.brca_signatures = pickle.load(f)
            
    def load_tnt_counts(self):
        file_name = 'data/tnt_counts_dict.pickle.gz'
        with gzip.open(file_name) as f:
            self.tnt_counts = pickle.load(f)
          
        
    
    def min_k(self):
        return self.data_dict['mink']
    
    def max_k(self):
        return self.data_dict['maxk']
        
    def error(self):
        min_k=self.min_k()
        max_k=self.max_k()
        error = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            error[i]=self.data_dict[i+min_k]['euclidean']
        return error
    
    def coph(self):
        min_k=self.min_k()
        max_k=self.max_k()
        coph = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            coph[i]=self.data_dict[i+min_k]['coph_corr']
        return coph
    
    def evar(self):
        min_k=self.min_k()
        max_k=self.max_k()
        evar = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            evar[i]=self.data_dict[i+min_k]['evar']
        return evar
    
    def kl(self):
        min_k=self.min_k()
        max_k=self.max_k()
        kl = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            kl[i]=self.data_dict[i+min_k]['kl']
        return kl
    
    def rss(self):
        min_k=self.min_k()
        max_k=self.max_k()
        rss = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            rss[i]=self.data_dict[i+min_k]['rss']
        return rss
    
    def dispersion(self):
        min_k=self.min_k()
        max_k=self.max_k()
        dispersion = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            dispersion[i]=self.data_dict[i+min_k]['disersion']
        return dispersion
    
    def silhouette(self):
        min_k=self.min_k()
        max_k=self.max_k()
        silhouette = np.zeros(max_k-min_k+1)
        for i in range(max_k-min_k+1):
            silhouette[i]=self.data_dict[i+min_k]['silhouette']['AVG']
        return silhouette
    
    def sparseness(self):
        min_k=self.min_k()
        max_k=self.max_k()
        sparseness = np.zeros((max_k-min_k+1,2))
        for i in range(max_k-min_k+1):
            sparseness[i]=self.data_dict[i+min_k]['sparseness']
        return sparseness
    
    def plot_metrics(self):
        min_k=self.min_k()
        max_k=self.max_k()
        rank_cands=range(min_k,max_k+1)
        dpi=250
        spar_h,spar_w=zip(*self.sparseness())
        
        plt.plot(rank_cands, self.coph(), 'o-', label='Cophenetic correlation', linewidth=2)
        #plt.plot(rank_cands, self.dispersion(),'o-', label='Dispersion', linewidth=2)
        #plt.plot(rank_cands, spar_h, 'o-', label='Sparsity (Basis)', linewidth=2)
        #plt.plot(rank_cands, spar_w, 'o-', label='Sparsity (Mixture)', linewidth=2)
        plt.plot(rank_cands, self.evar(), 'o-', label='Explained variance', linewidth=2)
        plt.plot(rank_cands, self.silhouette(), 'o-', label='AVG Silhouette', linewidth=2)
        plt.legend(loc=3,mode = 'expand', bbox_to_anchor=(0., 1.02, 1., .102), ncol=3, numpoints=1)
        plt.xlabel('Rank candidates')
        plt.ylabel('Value')
        plt.savefig('Images/Summary.png',dpi=dpi);
       
    
    def distance_cosmic(self, dist, rank):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.cosmic_signatures is None:
            self.load_cosmic()
        procs = pd.DataFrame(self.get_processes(rank)).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(self.cosmic_signatures.keys()),
                         dtype=np.float64)
        P=self.get_processes(rank)
        types= self.get_substitutions()
        for ind in df.index:
            vec1 = np.array([self.cosmic_signatures[ind][k] for k in sorted(self.cosmic_signatures[ind])])
            for col in df.columns:
                proc=dict(zip(types,P[:,col]))
                vec2 = np.array([proc[k] for k in sorted(proc)])
                d = dist(*(vec1, vec2))
                df.loc[ind, col] = d
        return df
    
    def distance_brca(self, dist, rank1,rank2):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.brca_signatures is None:
            self.load_brca()
        
        if self.tnt_counts is None:
            self.load_tnt_counts()
        P=self.get_processes(rank1)
        Pbraca=self.brca_signatures[rank2]['processes']
        procs = pd.DataFrame(P).to_dict()
        procsBrca=pd.DataFrame(Pbraca).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(procsBrca.keys()),
                         dtype=np.float64)
        types= self.get_substitutions()
        
        for ind in df.index:
            proc1=dict(zip(types,Pbraca[:,ind]))
            proc1=normalize_profile(proc1,self.tnt_counts['exon'])
            vec1 = np.array([proc1[k] for k in sorted(proc1)])
            for col in df.columns:
                proc2=dict(zip(types,P[:,col]))
                proc2=normalize_profile(proc2,self.tnt_counts['full'])
                vec2 = np.array([proc2[k] for k in sorted(proc2)])
                d = dist(*(vec1, vec2))
                df.loc[ind, col] = d
        return df
    
    
    def heatmap_brca(self,dist,rank1,rank2):
        df = self.distance_brca(dist,rank1,rank2)
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_title(str(dist))
        sns.heatmap(df, annot=True, fmt='g', cmap='viridis')
        plt.show()
    
    
    def heatmap_cosmic(self,dist,rank):
        df = self.distance_cosmic(dist,rank)
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_title(str(dist))
        sns.heatmap(df, annot=True, fmt='g', cmap='viridis')
        plt.show()
    
    def plot_cophenetic(self):
        min_k=self.min_k()
        max_k=self.max_k()
        rank_cands=range(min_k,max_k+1)
        plt.plot(rank_cands, self.coph(), 'o-', label='Cophenetic correlation', linewidth=2)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, numpoints=1);
        
    
    def plot_error(self):
        min_k=self.min_k()
        max_k=self.max_k()
        rank_cands=range(min_k,max_k+1)
        plt.plot(rank_cands, self.rss(), 'o-', label='Reconstruction error', linewidth=2)
        plt.plot(rank_cands, self.kl(), 'o-', label='Kullback-Leiber divergence', linewidth=2)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, numpoints=1);
        
    
    def get_samples(self):
        return self.data_dict['samples']
    
    def get_substitutions(self):
        return self.data_dict['types']
        
    def get_processes(self, rank):
        processes=self.data_dict[rank]['processes']
        processes[processes < 0] = 0
        return processes
        
    def plot_processes(self,rank,yMax=0.2):
        processes = self.get_processes(rank)
        processes = pd.DataFrame(processes)
        # Image size
        width, height = 600, 1000 
        # Pixel border around image
        border = 1
        dpi = 72.0
        figsize= (width+2*border)/float(dpi), (height+2*border)/float(dpi) 
        fig, axes = plt.subplots(ncols=1, nrows=len(processes.columns), figsize=figsize)
        i = 0
        for col in processes.columns:
            sig=dict(zip(self.get_substitutions(), np.array(processes[[col]])))
            plot_signature(sig, col, ax=axes[i], ymax=yMax)
            i += 1
        plt.show()
        
        
    def get_consensus(self,rank):
        return self.data_dict[rank]['consensus']
    
    
class ResultsNMF(object):
    
    def __init__(self, file_name):
        self.file_name = file_name
        with gzip.open(file_name, 'rb') as rs:
            self.data_dict = pickle.load(rs)
        self.cosmic_signatures = None
        
    def get_nruns(self):
        return self.data_dict['nruns']
    
    def get_rank(self):
        return self.data_dict['rank']
    
    def get_samples(self):
        return self.data_dict['samples']
    
    def get_contexts(self):
        return self.data_dict['typs']
    
    def get_mean_error(self):
        n = self.get_nruns()
        error = 0
        for i in range(n):
            error+=self.data_dict[i]['summary']['euclidean']
        return error/n
    
    def get_consensus(self):
        consensus = np.zeros((self.get_samples().shape[0], self.get_samples().shape[0]))
        n = self.get_nruns()
        for i in range(n):
            consensus+=self.data_dict[i]['summary']['connectivity']
        return consensus/n
    
    def plot_consensus(self):
        """
        Plot reordered consensus matrix.
        :param C: Reordered consensus matrix.
        :type C: numpy.ndarray`
        :param rank: Factorization rank.
        :type rank: `int`
        """
        C=self.get_consensus()
        C=reorder(C)
        imshow(C)
        set_cmap("RdBu_r")
        #plt.show()
        savefig("consensus_%d.png")
    
    def get_all_processes(self):
        d={}
        n = self.get_nruns()
        for i in range(n):
            d[i]=self.data_dict[i]['basis']
        return d
    
    def get_all_weights(self):
        d={}
        n = self.get_nruns()
        for i in range(n):
            d[i]=self.data_dict[i]['coef']
        return d
    
    def get_cluster(self):
        d=self.get_all_processes()
        nmfres = pd.DataFrame([])
        for i in d.keys():
            A=pd.DataFrame(d[i])
            nmfres=pd.concat([nmfres,A],axis=1)
        random.seed(1)
        model = KMeans(n_clusters=self.get_rank(),random_state=1234).fit(nmfres.T)
        return model
    
    def get_silhouette(self):
        cluster = self.get_cluster()
        labels = cluster.labels_
        d=self.get_all_processes()
        nmfres = pd.DataFrame([])
        for i in d.keys():
            A=pd.DataFrame(d[i])
            nmfres=pd.concat([nmfres,A],axis=1)
        return metrics.silhouette_score(nmfres.T, labels, metric='euclidean')
    
    def get_processes(self):
        processes = self.get_cluster().cluster_centers_.T
        processes[processes < 0] = 0
        return processes
    
    def load_cosmic(self):
        file_name = 'data/cosmic_signatures.pickle.gz'
        with gzip.open(file_name) as f:
            self.cosmic_signatures = pickle.load(f)
            
    def cosine_distance(self):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.cosmic_signatures is None:
            self.load_cosmic()
        procs = pd.DataFrame(self.get_processes()).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(self.cosmic_signatures.keys()),
                         dtype=np.float64)
        for ind in df.index:
            vec1 = np.array([self.cosmic_signatures[ind][k] for k in sorted(self.cosmic_signatures[ind])])
            for col in df.columns:
                vec2 = np.array([procs[col][k] for k in sorted(procs[col])])
                d = spdist.cosine(vec1, vec2)
                df.loc[ind, col] = d
        return df
    
    
    def distance_cosmic(self, dist):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.cosmic_signatures is None:
            self.load_cosmic()
        procs = pd.DataFrame(self.get_processes()).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(self.cosmic_signatures.keys()),
                         dtype=np.float64)
        for ind in df.index:
            vec1 = np.array([self.cosmic_signatures[ind][k] for k in sorted(self.cosmic_signatures[ind])])
            for col in df.columns:
                vec2 = np.array([procs[col][k] for k in sorted(procs[col])])
                d = dist(*vec1, *vec2)
                df.loc[ind, col] = d
        return df
    
    def kl_distance(self):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.cosmic_signatures is None:
            self.load_cosmic()
        procs = pd.DataFrame(self.get_processes()).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(self.cosmic_signatures.keys()),
                         dtype=np.float64)
        for ind in df.index:
            vec1 = np.array([self.cosmic_signatures[ind][k] for k in sorted(self.cosmic_signatures[ind])])
            for col in df.columns:
                vec2 = np.array([procs[col][k] for k in sorted(procs[col])])
                d = kullback_leibler(vec1, vec2)
                df.loc[ind, col] = d
        return df
    
    def battacharyya_distance(self):
        """heatmap with 1 - cos(u, v) for each process u and cosmic signature v"""
        
        if self.cosmic_signatures is None:
            self.load_cosmic()
        procs = pd.DataFrame(self.get_processes()).to_dict()
        df = pd.DataFrame(columns=sorted(procs.keys()), index=sorted(self.cosmic_signatures.keys()),
                         dtype=np.float64)
        for ind in df.index:
            vec1 = np.array([self.cosmic_signatures[ind][k] for k in sorted(self.cosmic_signatures[ind])])
            for col in df.columns:
                vec2 = np.array([procs[col][k] for k in sorted(procs[col])])
                d = battacharyya(vec1, vec2)
                df.loc[ind, col] = d
        return df

    def heatmap_cosine_distance(self):
        df = self.cosine_distance()
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_title('Cosine Distances')
        sns.heatmap(df, annot=True, fmt='g', cmap='viridis')
        plt.show()
        
    def heatmap_kl_distance(self):
        df = self.kl_distance()
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_title('Kullback Leibler Distances')
        sns.heatmap(df, annot=True, fmt='g', cmap='viridis')
        plt.show()
        
    def heatmap_battacharyya_distance(self):
        df = self.battacharyya_distance()
        fig, ax = plt.subplots(figsize=(8,8))
        ax.set_title('Battacharyya Distances')
        sns.heatmap(df, annot=True, fmt='g', cmap='viridis')
        plt.show()
        
    def plot_processes(self,yMax=0.2):#Falta arreglar el dataset amb noms y keys
        processes = self.get_processes()
        processes = pd.DataFrame(processes)
        # Image size
        width, height = 600, 1000 
        # Pixel border around image
        border = 1
        dpi = 72.0
        figsize= (width+2*border)/float(dpi), (height+2*border)/float(dpi) 
        fig, axes = plt.subplots(ncols=1, nrows=len(processes.columns), figsize=figsize)
        i = 0
        for col in processes.columns:
            sig = processes[[col]].to_dict()[col]
            plot_signature(sig, col, ax=axes[i], ymax=yMax)
            i += 1
        plt.show()
        
    