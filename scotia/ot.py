import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
import ot
from .dbscan import dbscan_ff_cell

def sel_pot_inter_cluster_pairs(S_all_arr,cluster_cell_df,effect_range=50):
    """select potentially communicating cell pairs (spatially proximal)
       filtering out cell pair that exceed the effect range

    ---------------
    Required inputs:
    -S_all_arr: cell by cell spatial distance array without any filtering
    -cluster_cell_df: dbscan clustering results

    ---------------
    Key parameters:
    -effect_range: used for checking whether two cell clusters are spatially proximal 
                   to each other, also a normalization factor
    ---------------
    Returns: modified cell by cell spatial distance array
            (filtered cell pairs were marked with Inf)
    ---------------
    Example:
    >>> pos_arr = np.array(pd.read_csv('./input_files/position.csv',index_col=0))
    >>> S_all_arr = distance_matrix(pos_arr,pos_arr)
    >>> S_all_arr.shape
    (30, 30)
    >>> idx_l, eps = dbscan_ff_cell(pos_arr,np.array(range(pos_arr.shape[0])),min_cluster_size=5,eps_l = list(range(10,50,1)))
    >>> cluster_df = pd.DataFrame({'cell_type':['ct1' for x in idx_l],'cell_idx':idx_l})
    >>> S_all_arr_new = sel_pot_inter_cluster_pairs(S_all_arr,cluster_df)
    >>> S_all_arr_new
    >>> S_all_arr_new.shape
    (30, 30)
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

    ---------------
    Required inputs:
    -dis_arr: cell by cell spatial distance array (get from function sel_pot_inter_cluster_pairs, Inf for excluded cell pairs)
    -exp_df: ligand and receptor expression dataframe
    -known_lr_pairs: ligand-receptor pairs
    -meta_df: meata data, including annotation information

    ---------------
    Key parameters:
    -reg: entropy regularization parameter for penalizing one-to-one transport cases
    -reg_m: parameter for relaxing the unbalanced distribution between source and target cells
    -dist_cutoff: distance cutoff for the post-processing of ot results
    -min_likeli_cutoff: likelihood cutoff for the post-processing of ot results
    -effect_range: the normalization factor

    ---------------
    Returns: cell-cell interaction likelihood dataframe:
             Columns: source_cell_idx, target_cell_idx, interaction likelihood, LR pair, source_celltype, target_celltype

    ---------------
    Example:
    >>> known_lr_pairs = pd.read_csv('./input_files/known_lr_pairs.csv',index_col=0)
    >>> S_all_arr_new = np.load('./input_files/S_all_arr_new.npy',allow_pickle=True)
    >>> exp_df = pd.read_csv('./input_files/normed_exp_mtx.csv',index_col=0)
    >>> meta_df = pd.read_csv('./input_files/metadata.csv',index_col=0)
    >>> ot_results = source_target_ot(S_all_arr_new, exp_df, meta_df, known_lr_pairs)
    >>> ot_results.shape
    (2709, 6)
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
       averaged likelihoods of each LR pair for each cell type pair
    ---------------
    Required inputs:
    -ot_data: cell-cell interaction likelihood results from ot analysis
    -label: sample_id
    ---------------
    Optimal inputs:
    -it_n_label: labels for permutation tests
    ---------------
    Returns: Summary dataframe for each LR pair for each cell type pair
             Columns: 'LR|sampleid|sourcecelltype_targetcelltype', averaged_interaction_likelihood
    ---------------
    Example:
    # doctest: +NORMALIZE_WHITESPACE
    >>> ot_results =  pd.read_csv('./input_files/ot_results.csv',index_col=0)
    >>> post_ot(ot_results,label='test').head(1) 
                                        label  ave_likelihood
    0  Dll4_Notch2|test|Erythroid_Erythroidpro        0.098149
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
