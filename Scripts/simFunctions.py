import numpy as np
import matplotlib.pyplot as plt
import nimfa
from scipy.spatial import distance
import scipy
from Signatures import *


def poisonMatrix(X):
    M=X
    n=M.shape[0]
    m=M.shape[1]
    for i in range(n):
        for j in range(m):
            M[i,j]=np.random.poisson(M[i,j], 1)
    return M

def arrange_sign(FP, P,dist):
    n=FP.shape[1]
    D=dist_sign(FP, P,dist)
    pos = np.ones((n,2))
    for i in range(n):
        index=np.unravel_index(D.argmin(), D.shape)
        pos[i]=index
        D[index[0],:]=1
        D[:,index[1]]=1
    return pos

def dist_sign(FP, P,dist):
    n=FP.shape[1]
    D=np.ones((n,n))
    for i in range(n):
        for j in range(n):
            D[i][j]=dist(*(FP[:,i],P[j]))
    return D

def calculate_dist(FP, P,dist):
    n=FP.shape[0]
    D=np.ones((n,n))
    pos = np.ones((n,2))
    totaldist=0
    for i in range(n):
        for j in range(n):
            D[i][j]=dist(*(FP[i],P[j]))
    for i in range(n):
        index=np.unravel_index(D.argmin(), D.shape)
        pos[i]=index
        totaldist+=D[index[0],index[1]]
        D[index[0],:]=1
        D[:,index[1]]=1
    return totaldist

def plot_sign(FP,P,i,dist):
    width, height = 2000, 1000
    border = 1
    dpi = 200.0
    figsize= (width+2*border)/float(dpi), (height+2*border)/float(dpi) 
    fig, ax = plt.subplots(figsize=figsize)
    bar_width = 0.35
    opacity = 0.4
    dim=FP.shape[0]
    labels = np.arange(1,dim+1)
    index = np.arange(dim)
    major_ticks = [k+0.15 for k in index]
    ax.tick_params(length=0, which='major', pad=10, labelsize=12)
    ax.set_xticks(major_ticks, minor=False)
    ax.set_xticklabels(labels, minor=False)
    
    indx=arrange_sign(FP,P,dist)
    rects1 = plt.bar(index, P[int(indx[i][1])].T, bar_width,
                     alpha=opacity,
                     color='b',
                     label='True values')

    rects2 = plt.bar(index + bar_width,np.matrix(FP[:,int(indx[i][0])]).T, bar_width,
                     alpha=opacity,
                     color='r',
                     label='Fitted values')

    plt.xlabel('Type Mutation')
    plt.ylabel('Freq')
    plt.title('Comparison between real and fitted values')
    plt.tight_layout()
    plt.legend()
    plt.savefig('Images/TrueVsFitted{0}.png'.format(i),dpi=dpi)
    plt.show()
    
    
    
def produce_exposure(l,L,rank):
    n=np.random.randint(l,L)
    exp=np.random.random_sample(rank)
    exp=exp/sum(exp)*n
    return(exp)

def produce_process(spar,dim):
    sumP=0
    while(sumP==0):
        process=np.random.binomial(1, 1-spar,dim)*np.random.random_sample(dim)
        sumP=sum(process)
    return(process/sumP)

def generate_exposures(l,L,rank, num_samples):
    E=[produce_exposure(l,L,rank)]
    for i in range(num_samples-1):
        E.append(produce_exposure(l,L,rank))
    return np.matrix(E)
    
def generate_process(spar,rank,dim,sep):
    k=1
    P=[produce_process(spar,dim)]
    while(k<rank):
        p=produce_process(spar,dim)
        if(is_ok(sep,P,p)):
            k+=1
            P.append(p)
    return np.matrix(P)
            
def is_ok(sep,P,p):
    for i in range(len(P)):
        if(distance.euclidean(P[i],p)<sep):
            return False
    return True