# Scotia
Scotia is a Python package for inferring ligand–receptor interactions from spatial imaging data.

# Installation
This package requires Python >=3.6.
Before installing scotia, it is highly recommended to create a new environment using [conda](https://docs.anaconda.com/free/anaconda/install/index.html).
```
conda create -n scotia_env python=3.9
conda activate scotia_env
```
In this new conda environment, you can git clone this repository to your workspace, and then install it with
```
python setup.py install
```
After installation, you can test it by
```
python
import scotia
```
# Usage/Example
Check out this [notebook](https://github.com/Caochris/SCOTIA/blob/master/notebook/scotia_example.ipynb) for a quick start tutorial. The example data used in the tutorial were included in the [example](https://github.com/Caochris/SCOTIA/tree/master/example) folder.

# Main functions
## DBSCAN cell clustering
```
scotia.dbscan_ff_cell(X_pos)
```
Required inputs:

X_pos: cell coordinates array.

Key parameters:

-min_cluster_size: the minimum number of cells in the clusters identified by DBSCAN, default: 10.

-eps_l: the maximum distance between neighboring cells within each cluster, default: range(15,150,5).

## Select adjacent cluster pairs
```
scotia.sel_pot_inter_cluster_pairs(dis_arr,cluster_cell_df)
```
Required inputs: 

dis_arr: cell by cell spatial distance array.

cluster_cell_df: cell clustering results.

## OT transport
```
scotia.source_target_ot(filtered_dis_arr, exp_df, meta_df, known_lr_pairs)
```
Required inputs:

dis_arr: cell by cell spatial distance array (get from function `sel_pot_inter_cluster_pairs`, Inf for excluded cluster pairs).

exp_df:  gene expression dataframe.

meta_df: metadata, including annotation information.

known_lr_pairs: ligand-receptor pairs.

## Summarize OT results
```
scotia.post_ot(ot_data_df, label)
```
Required inputs:

ot_data_df: cell-cell interaction scores from OT analysis.

label: sample_id for output.

## Permutation test
```
scotia.permutation_test(X_pos)
```
Required inputs:

X_pos: cell coordinates array.

Key parameters:

-it_n: the total number of permutations.
