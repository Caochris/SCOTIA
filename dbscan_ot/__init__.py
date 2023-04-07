import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from forest_fire_clustering.forest_fire_clustering import FFC
from sklearn.metrics.cluster import adjusted_rand_score
from matplotlib.collections import LineCollection
from scipy.spatial import distance_matrix
from collections import Counter
from random import shuffle
import sys
import warnings
import pickle
import ot
import bezier

def bdscan_ff_cell(X, X_index_arr = None, min_cluster_size = 10, eps_l = list(range(15,150,5)), ftem_l = list(range(1,60,1)), ff_sed = 123, unclustered_ratio = 0.2, maxcluster_ratio = 0.8):
    """DBSCAN cell clustering
    Dynamic determination of the eps parameters of DBSCAN by finding the most consistent clustering results between DBSCAN and FFC
    Required inputs: X, cell coordinates array
    Key parameters:
    -X_index_arr: cell index array
    -min_cluster_size: minimum size of clusters identified by DBSCAN
    -eps_l: range of eps to search (DBSCAN)
    -ftem_l: range of fire temperature to search (FFC)
    -unclustered_ratio: the maximum ratio of unclustered cell to total cells (exclude the cases that most cells are not clustered)
    -maxcluster_ratio: the maximum ratio of cells belonging to the largest clusters to total cells (exclude the cases that most cells are clustered into the same cluster)
    """
    cluster_pos_l = []
    cell_index_l = []
    fi_eps = None
    
    #run dbscan
    db_label_l = []
    for eps in eps_l:
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
    for eps,label_db in zip(eps_l,db_label_l):
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
        labels = db_label_l[eps_l.index(fi_eps)]
        labels_unique = np.unique(labels)
        for l in labels_unique:
            my_members = labels == l
            if l >-1 and X[my_members, :].shape[0]>=min_cluster_size:
                cluster_pos_l.append(X[my_members, :])
                if X_index_arr.any():
                    cell_index_l.append(X_index_arr[my_members])
        
    if X_index_arr.any(): 
        return cluster_pos_l, cell_index_l, fi_eps
    else:
        return cluster_pos_l, fi_eps

def sel_pot_inter_cluster_pairs(S_all_arr,cluster_cell_df,effect_range=50):
    """select potentially communicating cell cluster pairs (spatially proximal)
    Required inputs:
    -S_all_arr: cell by cell spatial distance array without any filtering
    -cluster_cell_df: cell clustering results
    Key parameters:
    -effect_range: used for checking whether two cell clusters are spatially proximal to each other, also a normalization factor
    """
    S_all_arr_new = np.zeros_like(S_all_arr)
    S_all_arr_new[:] = np.inf
    anno_t = list(set(cluster_cell_df['cell_type']))
    #print(anno_t)
    for l_cell_type in anno_t:
        for r_cell_type in anno_t:
            l_cell_clusters_df = cluster_cell_df[cluster_cell_df['cell_type']==l_cell_type]
            r_cell_clusters_df = cluster_cell_df[cluster_cell_df['cell_type']==r_cell_type]
            for j in range(l_cell_clusters_df.shape[0]):
                for k in range(r_cell_clusters_df.shape[0]):
                    if j!=k or l_cell_type != r_cell_type:
                        l_cell_ids_index = np.array(l_cell_clusters_df.iloc[j]['cell_idx'])
                        r_cell_ids_index = np.array(r_cell_clusters_df.iloc[k]['cell_idx'])

                        dis_mtx = S_all_arr[l_cell_ids_index,:][:,r_cell_ids_index]
                        #make sure cell clusters are spatially proximal to each other
                        #(at least one cell-cell pair that is within the effect_range distance)
                        mask_close = np.where(dis_mtx<=effect_range)
                        if dis_mtx[mask_close].shape[0] >= 1:
                            #print(l_cell_type,r_cell_type)
                            S_all_arr_new[l_cell_ids_index[:,None], r_cell_ids_index] = S_all_arr[l_cell_ids_index[:,None], r_cell_ids_index]
    return S_all_arr_new

def source_target_ot(dis_arr, exp_df, meta_df, known_lr_pairs, reg = 1, reg_m = 2, dist_cutoff = 50, min_likeli_cutoff = 0.01):
    """unbalanced ot between source and target cells
    Required inputs:
    -dis_arr: cell by cell spatial distance array (get from function sel_pot_inter_cluster_pairs, Inf for excluded cluster pairs)
    -exp_df: ligand and receptor expression dataframe
    -known_lr_pairs: ligand-receptor pairs
    -meta_df: meata data, including annotation information
    Key parameters:
    -reg: entropy regularization parameter for penalizing one-to-one transport cases
    -reg_m: parameter for relaxing the unbalanced distribution between source and target cells
    -dist_cutoff: distance cutoff for the post-processing of ot results
    -min_likeli_cutoff: likelihood cutoff for the post-processing of ot results
    -effect_range: the normalization factor
    """
    ga_df_final = pd.DataFrame([])
    mask = np.where((dis_arr==0)|(np.isnan(dis_arr)==True)|(dis_arr>dist_cutoff))
    dis_arr[mask] = np.inf
    dis_arr = dis_arr/dist_cutoff

    if not (dis_arr==np.inf).all().all():
        l_gene_list = list(known_lr_pairs['l_gene'])
        r_gene_list = list(known_lr_pairs['r_gene'])

        for l_gene, r_gene in zip(l_gene_list,r_gene_list):
            #print(l_gene,r_gene)
            ##filter cells with no expression
            cell_l_exp = np.array(exp_df[l_gene])
            cell_r_exp = np.array(exp_df[r_gene])

            cell_l_exp_left = cell_l_exp[np.where(cell_l_exp > 0)[0]]
            cell_r_exp_left = cell_r_exp[np.where(cell_r_exp > 0)[0]]
            
            cell_id_all = np.array(range(dis_arr.shape[0]))
            source_cell_id_left = cell_id_all[np.where(cell_l_exp > 0)[0]]
            target_cell_id_left = cell_id_all[np.where(cell_r_exp > 0)[0]]

            S_all_arr_left = dis_arr[np.where(cell_l_exp > 0)[0],:][:,np.where(cell_r_exp > 0)[0]]
            exp2 = np.outer(cell_l_exp_left,cell_r_exp_left)
            S_all_arr_left = S_all_arr_left/exp2

            ###remove all Inf rows and columns
            S_all_arr_left2 = S_all_arr_left[~(S_all_arr_left==np.inf).all(1)]
            S_all_arr_left2 = S_all_arr_left2[:,~(S_all_arr_left2==np.inf).all(0)]

            source_cell_id_left2 = source_cell_id_left[~(S_all_arr_left==np.inf).all(1)]
            target_cell_id_left2 = target_cell_id_left[~(S_all_arr_left==np.inf).all(0)]

            cell_l_exp_left2 = cell_l_exp_left[~(S_all_arr_left==np.inf).all(1)]
            cell_r_exp_left2 = cell_r_exp_left[~(S_all_arr_left==np.inf).all(0)]
            cell_l_exp_left2 = np.append(cell_l_exp_left2, 1)
            cell_r_exp_left2 = np.append(cell_r_exp_left2, 1)

            S_all_arr_left2 = np.insert(S_all_arr_left2, S_all_arr_left2.shape[0], np.ones(S_all_arr_left2.shape[1])*np.Inf, axis = 0)
            S_all_arr_left2 = np.insert(S_all_arr_left2, S_all_arr_left2.shape[1], np.ones(S_all_arr_left2.shape[0])*np.Inf, axis = 1)

            S_all_arr_left2[-1,-1] = 0.1

            if S_all_arr_left2.shape[0]>1:
                #print(l_gene,r_gene)
                ga = ot.sinkhorn_unbalanced(cell_l_exp_left2, cell_r_exp_left2, S_all_arr_left2, reg, reg_m)
                ga = ga/ga[-1,-1]
                ga = ga[:-1,:-1]
                ###post_process
                ga_df = pd.DataFrame(ga)
                ga_df.index = list(source_cell_id_left2)
                ga_df.columns = list(target_cell_id_left2)
                #print('bbbbbbb')
                ga_df = ga_df.stack().reset_index()
                ga_df.columns = ['l_id','r_id','likelihood_norm']
                ga_df = ga_df[ga_df['likelihood_norm']>min_likeli_cutoff]

                if ga_df.shape[0]>0:
                    #print(l_gene,r_gene)
                    ga_df['lr_pairs'] = [l_gene+"_"+r_gene for x in range(ga_df.shape[0])]
                    ga_df['l_anno'] = list(meta_df.loc[list(ga_df['l_id']),'annotation'])
                    ga_df['r_anno'] = list(meta_df.loc[list(ga_df['r_id']),'annotation'])
                    ga_df_final = pd.concat([ga_df_final,ga_df])
                else:
                    ga_df = pd.DataFrame([])
    return ga_df_final

def post_ot(ot_data_df, label,it_n_label = None):
    """post-processing of ot results by calculating the 
       averaged likelihoods of each LR pair for each cell type pari
    Required inputs:
    -ot_data: cell-cell interaction likelihood results from ot analysis
    -label: sample_id
    Optimal inputs:
    -it_n_label: labels for permutation tests
    """
    df_all = {}
    c_t_l = list(set(ot_data_df['cell_pairs']))

    for c_t_p in c_t_l:
        df_tmp_sel = ot_data_df[ot_data_df['cell_pairs']==c_t_p ]
        lr_pair_l = list(set(df_tmp_sel['ligand_recptor']))
        for lr_pair in lr_pair_l:
            df_tmp_sel2 = df_tmp_sel[df_tmp_sel['ligand_recptor']==lr_pair]
            if df_tmp_sel2.shape[0]>=8:
                groupby_source = df_tmp_sel2.groupby('source_cell_idx')['likelihood'].max()
                if it_n_label == None:
                    df_all[lr_pair+"|"+label+"|"+c_t_p] = groupby_source.mean()
                else:
                    df_all[lr_pair+"|"+label+"|"+c_t_p+"|"+str(it_n_label)] = groupby_source.mean()

    df_sum = pd.DataFrame({'label':df_all.keys(),'ave_likelihood':df_all.values()})
    return df_sum

def permutation_test(X_all,it_n=50,random_range=20):
    """permutation test, shuffle expression and randomize cooridnates
    Required inputs:
    -X_all: cell coordinates array
    Key parameters:
    -it_n: the total number of permutations
    -random_range: the range to select a random number for randomization of cell coordinates
    """

    shuffle_exp_idx = []
    new_pos = np.array([])
    for it_i in range(it_n):
        cell_idx = list(range(X_all.shape[0]))
        shuffle(cell_idx)
        shuffle_exp_idx.append(cell_idx)
        random_x = np.random.choice(range(-1*random_range,random_range+1), size=X_all.shape[0], replace=True)
        random_y = np.random.choice(range(-1*random_range,random_range+1), size=X_all.shape[0], replace=True)
        new_pos_tmp = X_all + np.array([list(random_x),list(random_y)]).T
        
        if new_pos.shape[0]==0:
            new_pos = new_pos_tmp
        else:
            new_pos = np.hstack((new_pos,new_pos_tmp))
    shuffle_idx_df = pd.DataFrame(shuffle_exp_idx).T
    tmp_pos_df = pd.DataFrame(np.round(new_pos,2))
    return tmp_pos_df,shuffle_idx_df

def curved_edges(Xs, Xt):
    """calculate curve for LR visualization
    Xs, Xt are the coordinates array of the source and target cells.
    """
    dist_ratio = 0.2
    bezier_precision=20
    rnd = np.where(np.random.randint(2, size=Xs.shape[0])==0, -1, 1)

    coords_node1 = Xs
    coords_node2 = Xt


    # Swap node1/node2 allocations to make sure the directionality works correctly
    should_swap = coords_node1[:,0] > coords_node2[:,0]
    coords_node1[should_swap], coords_node2[should_swap] = coords_node2[should_swap], coords_node1[should_swap]

    # Distance for control points
    dist = dist_ratio * np.sqrt(np.sum((coords_node1-coords_node2)**2, axis=1))

    # Gradients of line connecting node & perpendicular
    m1 = (coords_node2[:,1]-coords_node1[:,1])/(coords_node2[:,0]-coords_node1[:,0])
    m2 = -1/m1

    # Temporary points along the line which connects two nodes
    t1 = dist/np.sqrt(1+m1**2)
    v1 = np.array([np.ones(Xs.shape[0]),m1])
    coords_node1_displace = coords_node1 + (v1*t1).T
    coords_node2_displace = coords_node2 - (v1*t1).T

    # Control points, same distance but along perpendicular line
    # rnd gives the 'polarity' to determine which side of the line the curve should arc
    t2 = dist/np.sqrt(1+m2**2)
    v2 = np.array([np.ones(Xs.shape[0]),m2])
    coords_node1_ctrl = coords_node1_displace + (rnd*v2*t2).T
    coords_node2_ctrl = coords_node2_displace + (rnd*v2*t2).T

    # Combine all these four (x,y) columns into a 'node matrix'
    node_matrix = np.array([coords_node1, coords_node1_ctrl, coords_node2_ctrl, coords_node2])
    # Create the Bezier curves and store them in a list
    curveplots = []
    for i in range(Xs.shape[0]):
        nodes = node_matrix[:,i,:].T
        #print(nodes)
        curveplots.append(bezier.Curve.from_nodes(nodes).evaluate_multi(np.linspace(0,1,bezier_precision)).T)

    # Return an array of these curves
    curves = np.array(curveplots)
    return curves
