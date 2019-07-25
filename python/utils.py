import numpy as np
import matplotlib.pyplot as plt
import pylab as pil

def DrawSquareLattice_Problem(L,h,J, outfile=None):
    
    # coordinates of the spins
    fig = plt.figure(figsize = (5,5))
    spins_x = np.array([[i for j in range(L)] for i in range(L)]).flatten()
    spins_y = np.array([[i for i in range(L)] for j in range(L)]).flatten()
    
    
    # now draw the lines
    for i in range(L*L):
        for j in range(L*L):
            if J[i,j] != 0 :
                if J[i,j] > 0:
                    color = 'red'
                    lw = 1
                    ls = '-'
                    z = 2
                else:
                    if (J[i,j] <= -1):
                        lw = 1.5
                        ls = '-'
                    else: 
                        lw = 1
                        ls = '--'
                        
                    color = 'black'
                    z = 1
                
                points_x = [spins_x[i],spins_x[j]]
                points_y = [spins_y[i],spins_y[j]]
                
                plt.plot(points_x, points_y, color=color, linewidth = lw, linestyle = ls, zorder=z)
                
    colors = [ 'grey' for i in range(L*L)]
    for i in range(L*L):
        if h[i] > 0:
            colors[i] = 'red'
        elif h[i] < 0:
            colors[i] = 'magenta'
    
    plt.scatter(spins_x, spins_y, marker = 'o', color = colors, zorder = 3, s = 50)
                
    if outfile is not None:
        pil.savefig(outfile, bbox_inches='tight')
    plt.show()
    
    
    
def DrawSquareLattice_Solution(L,spins, outfile=None):
    fig = plt.figure(figsize = (5,5))
    spins_x = np.array([[i for j in range(L)] for i in range(L)]).flatten()
    spins_y = np.array([[i for i in range(L)] for j in range(L)]).flatten()
    
    for i in range(spins.shape[0]):
        if spins[i]>0:
            color = 'C0'
            marker = '^'
        else:
            color = 'red'
            marker = 'v'
            
        x = spins_x[i]
        y = spins_y[i]
        plt.scatter(x, y, marker = marker, color = color, s = 50 )
        
    if outfile is not None:
        pil.savefig(outfile, bbox_inches='tight')
        
    plt.show()
    
    
def DrawChimeraGraph_Problem(L,h,J,outfile=None):
    # number of spins
    n_spins = 8*(L**2)
    fig = plt.figure(figsize = (5,5))
    centers_x = np.array([[1.5* i for i in range(L)] for j in range(L)]).flatten()
    centers_y = np.array([[1.5* j for i in range(L)] for j in range(L)]).flatten()
    spins_x = np.array([[-0.5+centers_x[i],
                        -1./6.+centers_x[i],
                        1./6.+centers_x[i],
                        0.5+centers_x[i],
                        centers_x[i],
                        centers_x[i],
                        centers_x[i],
                        centers_x[i]] for i in range(L**2)]).flatten()
    spins_y = np.array([[centers_y[i],
                        centers_y[i],
                        centers_y[i],
                        centers_y[i],
                        -0.5+centers_y[i],
                        -1./6.+centers_y[i],
                        1./6.+centers_y[i],
                        0.5+centers_y[i]] for i in range(L**2)]).flatten()
        
    # now draw the lines
    for i in range(8*L*L):
        for j in range(8*L*L):
            if J[i,j] != 0 :
                if J[i,j] > 0:
                    color = 'red'
                    lw = 1
                    ls = '-'
                    z = 2
                else:
                    if (J[i,j] <= -1):
                        lw = 1.5
                        ls = '-'
                    else:
                        lw = 1
                        ls = '--'
                                                                                    
                    color = 'black'
                    z = 1
                                                                                            
                points_x = [spins_x[i],spins_x[j]]
                points_y = [spins_y[i],spins_y[j]]
                                                                                            
                plt.plot(points_x, points_y, color=color, linewidth = lw, linestyle = ls, zorder=z)

    colors = [ 'grey' for i in range(L*L)]
    for i in range(8*L*L):
        if h[i] > 0:
            colors[i] = 'red'
        elif h[i] < 0:
            colors[i] = 'magenta'

    plt.scatter(spins_x, spins_y, marker = 'o', color = colors, zorder = 3, s = 50)

    if outfile is not None:
        pil.savefig(outfile, bbox_inches='tight')
    plt.show()
    
    
