import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from scipy.spatial import distance_matrix
from collections import Counter
import pickle
import bezier

def curved_edges(Xs, Xt):
    """calculate curve for LR visualization
    ---------------
    Required inputs:
    -Xs, the coordinates array of source cells.
    -Xt, the coordinates array of target cells.

    ---------------
    Returns: visualization curves indicating the interaction strength

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
