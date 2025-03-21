import pandas as pd
import numpy as np
import anndata as ad
from .dbscan import dbscan_ff_cell
from scipy.spatial import distance_matrix
from .ot import sel_pot_inter_cluster_pairs, source_target_ot, post_ot
import os

def lr_score(adata, lr_list, sample_col, fov_col, celltype_col, output_path, sel_sample_id=[], sel_fov_id=[]):
    """---------------
    Required inputs: 
        adata: AnnData object with spatial data in adata.obsm['spatial']
        lr_list: DataFrame/array of ligand-receptor pairs with 2 columns
        sample_col: str, Sample column name in adata.obs
        fov_col: str, FOV column name in adata.obs
        celltype_col: str, Cell type column name in adata.obs
        output_path: str, Output directory path
        sel_sample_id: list, Select certain samples, default is [] which using all samples
        sel_fov_id: list, Select certain fovs, applied when sel_sample_id is not [], default is [] which using all fovs
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
    elif isinstance(lr_list, pd.DataFrame) and lr_list.shape[1] == 2:
        lr_list.columns = ['l_gene','r_gene']
    else:
        raise ValueError("lr_list must have exactly 2 columns, e.g., l_gene, r_gene")
    
    if len(sel_sample_id) >0:
        invalid_samples = [s_id for s_id in sel_sample_id if s_id not in set(adata.obs[sample_col])]
        if len(invalid_samples) >0:
            raise ValueError(f"Selected sample id(s) {invalid_samples} do not exist in the {sample_col} column")
    
    if len(sel_fov_id) >0:
        invalid_fovs = [f_id for f_id in sel_fov_id if f_id not in set(adata.obs[fov_col])]
        if len(invalid_fovs) >0:
            raise ValueError(f"Selected fov id(s) {invalid_fovs} do not exist in the {fov_col} column")
    
    for sample_id in sel_sample_id:
        sample_fovs = set(adata.obs[adata.obs[sample_col]==sample_id][fov_col])
        missing_fovs = [f_id for f_id in sel_fov_id if f_id not in sample_fovs]
            
        if missing_fovs:
            raise ValueError(f"FOV ID(s) {missing_fovs} not found in sample {sample_id}. Each FOV must exist in all selected samples.")

    # Add spatial coordinates to obs
    adata.obs['x_pos'] = adata.obsm['spatial'][:, 0].copy()
    adata.obs['y_pos'] = adata.obsm['spatial'][:, 1].copy()
    adata.obs['annotation'] = adata.obs[celltype_col]

    # Create output directories
    output_dirs = [f"{output_path}/{subdir}" for subdir in ['clustering', 'ot', 'ot/summary']]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)
    
    # Process expression data
    exp_df_all = pd.DataFrame(adata.X.todense(), columns=adata.var.index, index=adata.obs.index)
    df_quantile = exp_df_all[exp_df_all > 0].quantile(q=0.99, axis=0)
    
    if len(sel_sample_id)>0:
        adata_sel = adata[adata.obs[sample_col].isin(sel_sample_id)]
    if len(sel_fov_id)>0:
        adata_sel = adata_sel[adata_sel.obs[fov_col].isin(sel_fov_id)]

    for sample_id in set(adata_sel.obs[sample_col]):
        adata_sample = adata_sel[adata_sel.obs[sample_col]==sample_id]
        for fov in set(adata_sample.obs[fov_col]):
            
            adata_fov = adata_sample[adata_sample.obs[fov_col]==fov]
            adata_fov_obs = adata_fov.obs.copy()
            adata_fov_obs['index'] = range(adata_fov_obs.shape[0])
            cell_type_l = list(set(adata_fov_obs[celltype_col]))

            #get the size of fov
            fov_size_x = adata_fov_obs['x_pos'].max()-adata_fov_obs['x_pos'].min()
            fov_size_y = adata_fov_obs['y_pos'].max()-adata_fov_obs['y_pos'].min()
            fov_size = np.min([fov_size_x,fov_size_y])
            if fov_size/4 >50:
                search_range = list(range(15, int(fov_size/4), 5))
            else:
                search_range = list(range(15, int(fov_size/4), 1))
            
            #Run DBSCAN clustering
            cluster_results = run_dbscan(adata_fov_obs, celltype_col, mcell=10, sr=search_range)
            np.save(output_path+'/clustering/'+sample_id+'_fov_'+str(fov)+'_dbscan.cell.clusters',cluster_results)
            

            #Run Optimal Transport
            exp_df_tmp = pd.DataFrame(adata_fov.X.todense(), columns=adata_fov.var.index, index=adata_fov.obs.index)
            ot_result, summary_ot_results = run_ot_analysis(adata_fov_obs, exp_df_tmp, df_quantile, 
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
            idx_list, _ = dbscan_ff_cell(coords, X_index_arr=np.array(cells['index']),
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
