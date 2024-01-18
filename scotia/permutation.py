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
    >>> pos_arr = np.array(pd.read_csv('/your_work_space/input_files/position.csv',index_col=0))
    >>> it_n = 50 #permuation times
    >>> random_pos, shuffled_exp = scotia.permutation_test(pos_arr,int(it_n))
    >>> random_pos.loc[:5,:3]
           0       1      2      3
    0  -4.93   21.25   9.07  40.25
    1  25.36   49.69  -8.64  34.69
    2   1.14   64.98  26.14  79.98
    3  23.46   41.57  30.46  12.57
    4  12.66  112.57   0.66  90.57
    5  21.38   36.74  22.38  45.74
    >>> shuffled_exp.loc[:5,:3]
        0   1   2   3
    0   4   1  29   5
    1  15  16  10  11
    2  29   4  15  18
    3  13   7  23  14
    4   0   9   1  29
    5  27   8  16  19
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
