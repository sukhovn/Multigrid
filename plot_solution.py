from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import os
import math

#This function checks if number is float

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#This function reads a line of floats separated by spaces from a file

def read_coordinates(s):
	file = open(s,mode='r')
	data = file.read().split('\n')
	x = np.array([float(i) for i in data[0].split(' ') if is_float(i) == True])
	y = np.array([float(i) for i in data[1].split(' ') if is_float(i) == True])
	return x, y

def read_file(s):
	file = open(s,mode='r')
	strng = file.read().split(' ')
	return np.array([float(i) for i in strng if is_float(i) == True])

data_folder = "result"

xfunc, yfunc = read_coordinates("./" + data_folder + "/x.dat")
nx = len(xfunc)
ny = len(yfunc)

yfunc, xfunc = np.meshgrid(yfunc, xfunc)
# yfunc, xfunc = xfunc*np.sin(yfunc), xfunc*np.cos(yfunc)

zfunc = read_file("./" + data_folder + "/u.dat").reshape((nx, ny))
# zfunc = zfunc.reshape((nx-2, ny-2))
# zfunc = np.array([np.concatenate(([0.0], zfunc[i], [0.0])) for i in range(zfunc.shape[0])])
# zfunc = np.concatenate((np.zeros((1, zfunc.shape[1])), zfunc, np.zeros((1, zfunc.shape[1]))))

ax = plt.axes(projection='3d')
ax.plot_surface(xfunc, yfunc, zfunc, rstride=1, cstride=1,
                cmap='jet', edgecolor='none')
ax.view_init(45, -130)
plt.show(block=True)