import pandas as pd
import numpy as np
import anndata as ad
from .dbscan import dbscan_ff_cell
from scipy.spatial import distance_matrix
from .ot import sel_pot_inter_cluster_pairs, source_target_ot, post_ot
import os

def lr_score(adata, lr_list, sample_col, fov_col, celltype_col, output_path):
    """---------------
    Required inputs: 
        adata: AnnData object with spatial data in adata.obsm['spatial']
        lr_list: DataFrame/array of ligand-receptor pairs with 2 columns
        sample_col: Sample column name in adata.obs
        fov_col: FOV column name in adata.obs
        celltype_col: Cell type column name in adata.obs
        output_path: Output directory path
     -----------------
     """
    
    # Validate inputs
    required_keys = {
        'obsm': ['spatial'],
        'obs': [sample_col, fov_col, celltype_col]
    }
    
    for attr, keys in required_keys.items():
        for key in keys:
            if not hasattr(adata, attr) or (attr == 'obsm' and key not in getattr(adata, attr)) or \
               (attr == 'obs' and key not in adata.obs.columns):
                raise ValueError(f"Missing required {attr} key: {key}")
    
    if not os.path.exists(output_path):
        raise ValueError(f"Output path does not exist: {output_path}")
    
    # Convert lr_list to DataFrame if needed
    if isinstance(lr_list, np.ndarray) and lr_list.shape[1] == 2:
        lr_list = pd.DataFrame(lr_list, columns=['l_gene', 'r_gene'])
    elif not isinstance(lr_list, pd.DataFrame) or lr_list.shape[1] != 2:
        raise ValueError("lr_list must have exactly 2 columns")

    # Add spatial coordinates to obs
    adata.obs['x_pos'] = adata.obsm['spatial'][:, 0]
    adata.obs['y_pos'] = adata.obsm['spatial'][:, 1]
    adata.obs['annotation'] = adata.obs[celltype_col]

    # Create output directories
    output_dirs = [f"{output_path}/{subdir}" for subdir in ['clustering', 'ot', 'ot/summary']]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)
    
    # Process expression data
    exp_df_all = pd.DataFrame(adata.X.todense(), columns=adata.var.index, index=adata.obs.index)
    df_quantile = exp_df_all[exp_df_all > 0].quantile(q=0.99, axis=0)

    for sample_id in set(adata.obs[sample_col]):
        adata_sample = adata.obs[adata.obs[sample_col]==sample_id]
        for fov in set(adata_sample[fov_col]):
            
            adata_fov = adata_sample[adata_sample[fov_col]==fov]
            adata_fov['index'] = range(adata_fov.shape[0])
            cell_type_l = list(set(adata_fov[celltype_col]))

            #get the size of fov
            fov_size_x = adata_fov['x_pos'].max()-adata_fov['x_pos'].min()
            fov_size_y = adata_fov['y_pos'].max()-adata_fov['y_pos'].min()
            fov_size = np.max([fov_size_x,fov_size_y])
            if fov_size/4 >50:
                search_range = list(range(10, int(fov_size/4), 5))
            else:
                search_range = list(range(10, int(fov_size/4), 1))
            
            #Run DBSCAN clustering
            cluster_results = run_dbscan(adata_fov, celltype_col, mcell=5, sr=search_range)
            np.save(output_path+'/clustering/'+sample_id+'_fov_'+str(fov)+'_dbscan.cell.clusters',cluster_results)
            
            #Run Optimal Transport
            ot_result, summary_ot_results = run_ot_analysis(adata_fov, exp_df_all, df_quantile, 
                                            lr_list, cluster_results, sample_id, fov, output_path)
            ot_result.to_csv(output_path+'/ot/'+sample_id+'_fov_'+str(fov)+".ot.csv",header = True, index = False, sep = "\t")
            summary_ot_results.to_csv(output_path+'/ot/summary/'+sample_id+'_fov_'+str(fov)+".ot.csv",header = True, index = False, sep = "\t")
            
            
            print(f"Processed {sample_id} FOV {fov}")
    print(f"Analysis complete. Results saved in {output_path}/ot")

def run_dbscan(fov_data, celltype_col, mcell, sr):
    """Run DBSCAN clustering on FOV data"""
    celltype = []
    cell_idx = []
    for cell_type in set(fov_data[celltype_col]):
        print(cell_type)
        cells = fov_data[fov_data[celltype_col] == cell_type]
        coords = cells[['x_pos', 'y_pos']].values
        if coords.shape[0] >= mcell:
            idx_list, _ = dbscan_ff_cell(coords, X_index_arr=np.arange(len(cells)),
                            min_cluster_size=mcell, eps_range = sr)
            if idx_list:
                celltype += [cell_type for x in idx_list]
                cell_idx += idx_list
    return pd.DataFrame({'cell_type': celltype, 'cell_idx': cell_idx})
    
    
def run_ot_analysis(fov_data, exp_df, quantiles, lr_list, cluster_df, sample_id, fov, output_path):
    """Run Optimal Transport analysis on clustered data"""
    ga_df_final = pd.DataFrame({})
    final_summary = pd.DataFrame({})
    
    #coordinates
    cell_id_all = np.array(range(fov_data.shape[0]))
    coord = np.array(fov_data[['x_pos','y_pos']])
    S_all_arr = distance_matrix(coord,coord)
    
    #expression
    exp_df_fov = exp_df.loc[list(fov_data.index)]
    exp_df_fov = exp_df_fov/quantiles
    exp_df_fov[exp_df_fov>1]=1
    exp_df_fov.index = cell_id_all
    
    fov_data.index = range(fov_data.shape[0])

    #select potentially communicating cell cluster pairs (spatially adjacent)
    S_all_arr_new = sel_pot_inter_cluster_pairs(S_all_arr,cluster_df)

    #optimal transport between source and target cells
    ga_df_final = source_target_ot(S_all_arr_new, exp_df_fov, fov_data, lr_list)
    
    if ga_df_final.shape[0]>0:
        ga_df_final.columns = ['source_cell_idx','receptor_cell_idx','likelihood','ligand_recptor','source_cell_type','target_cell_type']
        ga_df_final['cell_pairs'] = ga_df_final['source_cell_type']+"_"+ga_df_final['target_cell_type']
        final_summary = post_ot(ga_df_final,label=sample_id)
        
    return ga_df_final.iloc[:,:-1], final_summary
