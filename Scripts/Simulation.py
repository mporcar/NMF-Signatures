import dill as pickle
import numpy as np
import matplotlib.pyplot as plt
import nimfa
from scipy.spatial import distance
import scipy
from Signatures import *
from simFunctions import *
import urllib
import gzip


n_samples=200
n_processes=3
l=50
L=90
mink=2
maxk=10
runsPerSim=50
numberSim=50
numcores=10
NoiseNum=3
spar=0.3
sep=0.3
ndim=np.array([30,40])
nprocesses=np.array([3,5,7])
type_nmf="snmf"
resultsdir="Results"
for q in ndim:
    for j in nprocesses:
        inputs={'n_samples':n_samples,'n_processes':j,'spar':spar,'dim':q,'sep': sep, 'l':l, 'L':L,'runsPerSim':runsPerSim,'numberSim':numberSim,'mink':mink, 'maxk':maxk, 'Noise': NoiseNum,'type_nmf':type_nmf}
        D_total={}
        for i in range(numberSim):
            print("============ Iteration number = %d =================" % i)
            E=generate_exposures(l,L,j,n_samples)
            P=generate_process(spar,j,q,sep)
            X=poisonMatrix(E*P)
            Noise=np.random.randint(NoiseNum, size=(n_samples, q))
            X=X+Noise
            D={}
            D['TrueProcesses']=P
            D['TrueExposures']=E
            for k in range(mink,maxk):
                print("====== rank number = %d =========" % k)
                d=run_one_object_method_C(X,k,runsPerSim,numcores,type_nmf)
                D[k]=d['summary']
            D_total[i]=D

        D_total['inputs']=inputs
        summary_res=resultsdir+'/'+type_nmf+str(q)+'simdata'+str(j)+'.pickle.gz'
        with gzip.open(summary_res, 'wb') as rs:
            pickle.dump(D_total, rs,pickle.HIGHEST_PROTOCOL) 

