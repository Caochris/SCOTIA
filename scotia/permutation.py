import numpy as np
import pandas as pd
from random import shuffle
import pickle

def permutation_test(X_all,it_n=50,random_range=20):
    """permutation test, shuffle expression and randomize coordinates

    ---------------
    Required inputs:
    -X_all: cell coordinates array

    ---------------
    Key parameters:
    -it_n: the total number of permutations
    -random_range: the range to select a random number for randomization of cell coordinates

    ---------------
    Returns: randomized cell coordinates, idx dataframe for shuffled expression matrix
    
    ---------------
    Example:
    >>> import numpy as np
    >>> import pandas as pd
    >>> pos_arr = np.array(pd.read_csv('./input_files/position.csv',index_col=0))
    >>> it_n = 50 #permuation times
    >>> random_pos, shuffled_exp = permutation_test(pos_arr,int(it_n))
    >>> print(random_pos.shape)
    (30, 100)
    >>> print(shuffled_exp.shape)
    (30, 50)
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
