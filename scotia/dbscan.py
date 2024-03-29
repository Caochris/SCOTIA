import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from forest_fire_clustering.forest_fire_clustering import FFC
from sklearn.metrics.cluster import adjusted_rand_score
from scipy.spatial import distance_matrix
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

def dbscan_ff_cell(X, X_index_arr, min_cluster_size = 10, eps_range = list(range(10,150,5)), ftem_l = list(range(1,60,1)), ff_sed = 123, unclustered_ratio = 0.2, maxcluster_ratio = 0.8):
    """DBSCAN cell clustering

    Dynamic determination of the eps parameter of DBSCAN by 
    finding the most consistent clustering results between DBSCAN and Foreset Fire clustering (FFC)

    ---------------
    Required inputs: X, cell coordinates array (2D)
                     X_index_arr, cell index/label array (1D)
    ---------------
    Parameters:
    -min_cluster_size: minimum size of clusters identified by DBSCAN
    -eps_range: range of eps to search (DBSCAN)
    -ftem_l: range of fire temperature to search (FFC)
    -unclustered_ratio: the maximum ratio of unclustered cell to total cells 
                        (exclude the cases that most cells are not clustered)
    -maxcluster_ratio: the maximum ratio of cells belonging to the largest clusters to total cells 
                        (exclude the cases that most cells are clustered into the same cluster)
    ----------------
    Returns: list of cell clusters, final eps used for DBSCAN clustering
    ----------------
    Example:
    >>> pos_arr = np.array(pd.read_csv('./input_files/position.csv',index_col=0))
    >>> idx_l, eps = dbscan_ff_cell(pos_arr,np.array(range(pos_arr.shape[0])),min_cluster_size=5,eps_range = list(range(10,50,1)))
    >>> idx_l
    [array([ 0,  1,  2,  3,  5,  6,  7,  8, 14]), 
    array([ 4,  9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26])]
    >>> eps
    35

    """
    fi_eps = None
    cell_index_l = []

    #run dbscan
    db_label_l = []
    for eps in eps_range:
        dbscan_tmp = DBSCAN(eps=eps).fit(X)
        db_label_l.append(dbscan_tmp.labels_)
        
    #run ffc
    ff_label_l = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cluster_obj = FFC(n_jobs = 5)
        np.random.seed(ff_sed)
        for ftem in ftem_l:
            cluster_obj.preprocess(X)
            cluster_obj.fit(fire_temp=ftem)
            ff_label_l.append(cluster_obj.cluster_labels)
    
    
    #compare dbscan with ffc
    score_t1 = {}
    score_t2 = {}
    score_t3 = {}
    for eps,label_db in zip(eps_range,db_label_l):
        score_t1[str(eps)]=0
        score_t2[str(eps)]=0
        score_t3[str(eps)]=0
        for ftem,label_ff in zip(ftem_l,ff_label_l):
            label_db_new = label_db[label_db!=-1]
            label_ff_new = label_ff[label_db!=-1]
            #check whether more than one cluster
            if len(set(label_db_new))>1 and len(set(label_ff_new))>1:
                #check the ratio of unclustered cells
                if label_db[label_db==-1].shape[0] < np.max([X.shape[0]*unclustered_ratio,20]):
                    #check the ratio of the cells belonging to the largest-size cluster
                    if max(Counter(label_db_new).values()) < X.shape[0]*maxcluster_ratio:
                        score_t1[str(eps)] += adjusted_rand_score(label_db_new,label_ff_new)
                    elif max(Counter(label_db_new).values()) < X.shape[0]:
                        score_t2[str(eps)] += adjusted_rand_score(label_db_new,label_ff_new)
            elif len(set(label_db_new))>=1 and len(set(label_ff_new))>=1:
                if label_db[label_db==-1].shape[0] < X.shape[0]:
                    if  max(Counter(label_db_new).values()) < X.shape[0]:
                        score_t3[str(eps)] += adjusted_rand_score(label_db_new,label_ff_new)

    #select the eps that results in the highest clustering consistency
    if max(score_t1.values()) >0:
        fi_eps = int([k for k,v in score_t1.items() if v == max(score_t1.values())][-1])
    elif max(score_t2.values()) >0:
        fi_eps = int([k for k,v in score_t2.items() if v == max(score_t2.values())][-1])
    elif max(score_t3.values()) >0:
        fi_eps = int([k for k,v in score_t3.items() if v == max(score_t3.values())][-1])
      
    if fi_eps:
        print('eps: '+str(fi_eps))
        labels = db_label_l[eps_range.index(fi_eps)]
        labels_unique = np.unique(labels)
        for l in labels_unique:
            my_members = labels == l
            if l >-1 and X[my_members, :].shape[0]>=min_cluster_size:
                cell_index_l.append(X_index_arr[my_members])
        
    return cell_index_l, fi_eps
   
