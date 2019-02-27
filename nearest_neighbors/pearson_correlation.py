
# PREAMBLE:

import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt

def pearson_correlation(trajectory_data):
        """ Calculates the Pearson Correlation Matrix for node pairs
        Usage: node_correlation, node_variance, node_average = pearson_correlation(trajectory_data)
        Arguments:
        trajectory_data: multidimensional numpy array; first index (rows) correspond to timestep, second index correspond to positions of each node; 

        Returns:
        node_correlation: a nNodes x nNodes square matrix (numpy array) filled with the Pearson Correlation Coefficients for all node pairs
        node_variance: one dimensional numpy array containing the variances of the data
        node_average: one dimensional numpy array containing the averages of the data

        """

        # ----------------------------------------
        # CALCULATING THE AVERAGE OF TRAJECTORY DATA
        # ----------------------------------------
        nSteps = len(trajectory_data)
        nSteps_range = range(nSteps)
        nNodes = len(trajectory_data[0])
        nNodes_range = range(nNodes)
        for ts in nSteps_range:
                # removing center of geometry translational motion
                center_of_geometry = np.mean(trajectory_data[ts])
                trajectory_data[ts] -= center_of_geometry
                # no rotations to worry about...

        node_average = np.sum(trajectory_data,axis=0)/nSteps 
       
        # ----------------------------------------
        # PREPARE NUMPY ARRAYS
        # ----------------------------------------
        node_variance = np.zeros(nNodes,dtype=np.float64)
        node_covariance = np.zeros((nNodes,nNodes),dtype=np.float64)
        node_correlation = np.zeros((nNodes,nNodes),dtype=np.float64)

        # ----------------------------------------
        # CALCULATING PEARSON CORRELATION COEFFICIENT MATRIX
        # ----------------------------------------
        for ts in nSteps_range:
                for i in nNodes_range:
                        node_variance[i] += trajectory_data[ts,i]**2
                        for j in nNodes_range[i:]:
                                node_covariance[i,j] += trajectory_data[ts,i]*trajectory_data[ts,j]

        node_variance /= nSteps
        node_variance -= node_average**2
        node_covariance /= nSteps

        for i in nNodes_range:
                for j in nNodes_range[i:]:
                        node_covariance[i,j] -= node_average[i]*node_average[j]
                        node_correlation[i,j] = node_covariance[i,j]/np.sqrt(node_variance[i]*node_variance[j])
                        node_correlation[j,i] = node_correlation[i,j]

        # ----------------------------------------
        # OUTPUT OF AVERAGE, VARIANCE, COVARIANCE, AND CORRELATION MATRICES
        # ----------------------------------------
        np.savetxt('node_positional_average.dat',node_average)
        np.savetxt('node_positional_variance.dat',node_variance)
        np.savetxt('node_positional_covariance.dat',node_covariance)
        np.savetxt('node_positional_correlation.dat',node_correlation)

        # ----------------------------------------
        # PLOTTING VARIANCE
        # ----------------------------------------
        plt.plot(nNodes_range,node_variance,'k')
        plt.xlabel('Node Index',size=14)
        plt.ylabel(r'Node Variance ($\AA^{2}$)',size=14)
        plt.tight_layout()
        plt.savefig('node_positional_variance.png',dpi=600,transparent=True)
        plt.close()

        # ----------------------------------------
        # PLOTTING COVARIANCE
        # ----------------------------------------
        fig, ax = plt.subplots()
        temp = plt.pcolormesh(range(nNodes+1),range(nNodes+1),node_covariance,cmap='bwr')
        cb1 = plt.colorbar()
        cb1.set_label(r'Node-Node Covariance ($\AA^{2}$)')
        
        xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((-0.5,nNodes+0.5))
        plt.ylim((-0.5,nNodes+0.5))
        plt.xlabel('Node Index',size=14)
        plt.ylabel('Node Index',size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig('node_positional_covariance.png',dpi=600,transparent=True)
        plt.close()
        
        # ----------------------------------------
        # PLOTTING CORRELATION
        # ----------------------------------------
        fig, ax = plt.subplots()
        temp = plt.pcolormesh(range(nNodes+1),range(nNodes+1),node_correlation,cmap='bwr',vmin=-1.0,vmax=1.0)
        cb1 = plt.colorbar()
        cb1.set_label('Node-Node Correlation')
        
        xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((-0.5,nNodes+0.5))
        plt.ylim((-0.5,nNodes+0.5))
        plt.xlabel('Node Index',size=14)
        plt.ylabel('Node Index',size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig('node_positional_correlation.png',dpi=600,transparent=True)
        plt.close()
        
        return node_correlation, node_variance, node_average
