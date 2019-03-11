
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.cm as cm
import numpy as np

# define function for making figures
def define_figure(xlabel="X",ylabel="Y"):
    # setup plot parameters
    fig = plt.figure(figsize=(10,8), dpi= 80, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
    ax.set_xlabel('Nodes',size=20)
    ax.set_ylabel('Nodes',size=20)
    plt.tick_params(axis='both',labelsize=20)
    return ax

def plot_square_matrix(square_matrix,figure_name,axes_label='',cbar_label='',plotting_cmap='bwr',v_range=None,minor_ticks=10,major_ticks=100):
        """
        """
        nNodes = len(square_matrix)
        node_range = range(nNodes+1)
        fig, ax = plt.subplots()
        ax.tick_params(which='major',length=6,width=2)
        ax.tick_params(which='minor',length=3,width=1)
        ax.xaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.xaxis.set_major_locator(MultipleLocator(major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.yaxis.set_major_locator(MultipleLocator(major_ticks))

        if v_range != None:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap,vmin=v_range[0],vmax=v_range[1])
        else:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap)
        cb1 = plt.colorbar()
        cb1.set_label(r'%s'%(cbar_label))

        xlabels = [str(int(x)+1) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)+1) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((0,nNodes))
        plt.ylim((0,nNodes))
        plt.xlabel(axes_label,size=14)
        plt.ylabel(axes_label,size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.plot(np.arange(0,nNodes,0.01),np.arange(0,nNodes,0.01),'--',lw=2,color='k')
        plt.savefig(figure_name,dpi=600,transparent=True)
        plt.close()

def plot_square_matrix_inline(square_matrix,figure_name,axes_label='',cbar_label='',plotting_cmap='bwr',v_range=None,minor_ticks=10,major_ticks=100):
        """
        """
        nNodes = len(square_matrix)
        node_range = range(nNodes+1)
        fig, ax = plt.subplots()
        ax.tick_params(which='major',length=6,width=2)
        ax.tick_params(which='minor',length=3,width=1)
        ax.xaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.xaxis.set_major_locator(MultipleLocator(major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.yaxis.set_major_locator(MultipleLocator(major_ticks))

        if v_range != None:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap,vmin=v_range[0],vmax=v_range[1])
        else:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap)
        cb1 = plt.colorbar()
        cb1.set_label(r'%s'%(cbar_label))

        xlabels = [str(int(x)+1) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)+1) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((0,nNodes))
        plt.ylim((0,nNodes))
        plt.xlabel(axes_label,size=14)
        plt.ylabel(axes_label,size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.plot(np.arange(0,nNodes,0.01),np.arange(0,nNodes,0.01),'--',lw=2,color='k')


