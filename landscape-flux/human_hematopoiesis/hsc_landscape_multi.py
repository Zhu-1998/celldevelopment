# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 02:03:39 2022

@author: USER
"""



import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import os
import dynamo as dyn
dyn.dynamo_logger.main_silence()

import warnings
warnings.filterwarnings('ignore')


#帮助您调试与版本相关的错误（如果有的话）
dyn.get_all_dependencies_version()

#模拟带有白色背景的 ggplot2 绘图样式
dyn.configuration.set_figure_params('dynamo', background='white')


adata_labeling = anndata.read("/home/wj/datadisk/zlg/singlecell/HSC/hematopoiesis_v1.h5ad")

VecFld = adata_labeling.uns['VecFld_umap']

def vector_field_function(x, VecFld=VecFld, dim=None, kernel="full", X_ctrl_ind=None):
    """Learn an analytical function of vector field from sparse single cell samples on the entire space robustly.

		Reference: Regularized vector field learning with sparse approximation for mismatch removal, Ma, Jiayi, etc. al, Pattern Recognition
		"""

    x = np.array(x).reshape((1, -1))
    if np.size(x) == 1:
        x = x[None, :]
    K = dyn.vf.utils.con_K(x, VecFld["X_ctrl"], VecFld["beta"])
    
    if X_ctrl_ind is not None:
        C = np.zeros_like(VecFld["C"])
        C[X_ctrl_ind, :] = VecFld["C"][X_ctrl_ind, :]
    else:
        C = VecFld["C"]

    K = K.dot(C)
    return K




import time

import math
# from scipy.integrate import odeint
from scipy.interpolate import griddata

time_start = time.time() #开始计时

VecFnc = vector_field_function
x_lim=[0, 20]
y_lim=[0, 20]
xyGridSpacing=1
dt=2e-1
tol=1e-2   #Tolerance to test for convergence.    if abs(Pot - Pot_old) > tol: print(1, "Warning: not converged!\n")
numTimeSteps=2000000
starttime = 100000
#(7,2) stabletime=2.5e5
Tra_grid = 200
# Calculate total no. of paths for defined grid spacing
numPaths = int(np.diff(x_lim) / xyGridSpacing + 1) * int(np.diff(y_lim) / xyGridSpacing + 1)

# Initialize "path" variable matrices
x_path = np.zeros((numPaths, numTimeSteps))  # x-coord. along path
y_path = np.zeros((numPaths, numTimeSteps))  # y-coord. along path

F_x = np.zeros((numPaths, numTimeSteps))
F_y = np.zeros((numPaths, numTimeSteps))

num_tra = np.zeros((Tra_grid, Tra_grid))
total_Fx = np.zeros((Tra_grid, Tra_grid))
total_Fy = np.zeros((Tra_grid, Tra_grid))

path_tag = np.ones((numPaths, 1), dtype="int")  # tag for given path (to denote basin of attraction)
# ** initialized to 1 for all paths **

# Initialize "Path counter" to 1
path_counter = 0

# Assign array to keep track of no. of paths per attractor
numPaths_att = None


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loop over x-y grid

# sigma = 1.0  # Standard deviation.
# mu = 0.0  # Mean.
# tau = 0.05  # Time constant.
D=0.2  #noise size
for i in np.arange(x_lim[0], x_lim[1] + xyGridSpacing, xyGridSpacing):
    for j in np.arange(y_lim[0], y_lim[1] + xyGridSpacing, xyGridSpacing):

        # *** Init conds for given (x,y) ***
        # Initialize coords.
        x0 = i
        y0 = j
        
        
        # Initialize "path" variables
        x_p = x0
        y_p = y0

        # Initialize global arrays (time t = 0 counts as "time step #1")
        x_path[path_counter, 0] = x_p
        y_path[path_counter, 0] = y_p
        
        dxdt, dydt = VecFnc([x_p, y_p])
        
        F_x[path_counter,  0] = dxdt
        F_y[path_counter,  0] = dydt        

        # Evaluate potential (Integrate) over trajectory from init cond to  stable steady state
        for n_steps in np.arange(1, numTimeSteps):
            print(path_counter, n_steps)

            # update dxdt, dydt

            dxdt, dydt = VecFnc([x_p, y_p])
            

            # update x, y
            # dx = dxdt * dt + dt * np.sqrt(2*D) * ((-(x_p - mu) / tau) + sigma * np.sqrt(2. / tau) * np.sqrt(dt) * np.random.randn())
            # dy = dydt * dt + dt * np.sqrt(2*D) * ((-(y_p - mu) / tau) + sigma * np.sqrt(2. / tau) * np.sqrt(dt) * np.random.randn())
            
            dx = dxdt * dt + np.sqrt(2*D) * np.sqrt(dt) * np.random.randn()
            dy = dydt * dt + np.sqrt(2*D) * np.sqrt(dt) * np.random.randn()
            
            # dx = dxdt * dt
            # dy = dydt * dt

            x_p = x_p + dx
            y_p = y_p + dy
            
            # if x_p < x_lim[0]:
            #     x_p = x_lim[0] 
            # if y_p < x_lim[0]:
            #     y_p = y_lim[0]
            
            if x_p < x_lim[0]:
                x_p = 2 * x_lim[0] - x_p
            if y_p < x_lim[0]:
                y_p = 2 * y_lim[0] - y_p
                
            if x_p > x_lim[1]:
                x_p = 2 * x_lim[1] - x_p
            if y_p > y_lim[1]:
                y_p = 2 * y_lim[1] - y_p
                

                
            dxdt, dydt = VecFnc([x_p, y_p])    
            F_x[path_counter,  n_steps] = dxdt
            F_y[path_counter,  n_steps] = dydt

            x_path[path_counter, n_steps] = x_p
            y_path[path_counter, n_steps] = y_p
            
            if n_steps > starttime:
                A = int((x_p - x_lim[0]) * Tra_grid / (x_lim[1] - x_lim[0]))   
                B = int((y_p - y_lim[0]) * Tra_grid / (y_lim[1] - y_lim[0]))
                if A < Tra_grid and B<Tra_grid:
        # 			A = math.floor((x_p - x_lim[0]) * Tra_grid / (x_lim[1] - x_lim[0]))
        #           B = math.floor((y_p - y_lim[0]) * Tra_grid / (y_lim[1] - y_lim[0]))
                    num_tra[A, B] = num_tra[A, B] + 1;
                    total_Fx[A, B] = total_Fx[A, B] + dxdt
                    total_Fy[A, B] = total_Fy[A, B] + dydt

        path_counter = path_counter + 1

# p_tra = num_tra / (numPaths * (numTimeSteps - starttime))
p_tra = num_tra / (sum(sum(num_tra)))
pot_U = -np.log(p_tra)
mean_Fx = total_Fx / num_tra
mean_Fy = total_Fy / num_tra

xlin = np.linspace(x_lim[0], x_lim[1], Tra_grid)
ylin = np.linspace(y_lim[0], y_lim[1], Tra_grid)
Xgrid, Ygrid = np.meshgrid(xlin, ylin)

np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/num_tra.csv", num_tra, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/total_Fx.csv", total_Fx, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/total_Fy.csv", total_Fy, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/p_tra.csv", p_tra, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/pot_U.csv", pot_U, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/mean_Fx.csv", mean_Fx, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/mean_Fy.csv", mean_Fy, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/Xgrid.csv", Xgrid, delimiter=",") 
np.savetxt("/home/wj/datadisk/zlg/singlecell/HSC/0.2_1_2000000_0.2_200/minmax_fantan/Ygrid.csv", Ygrid, delimiter=",") 

time_end = time.time()    #结束计时
time_c= time_end - time_start   #运行所花时间
print('time cost', time_c, 's')
