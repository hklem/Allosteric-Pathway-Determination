import numpy as np
import matplotlib.pyplot as plt

def autocorrelation_particle1(N,dt):
    # stystem description
    positions = np.loadtxt("positions.dat")
    nSteps = len(positions)
    nSteps_range = range(nSteps)
    nNodes = len(positions[0])

    # center of geometry motion
    for ts in nSteps_range:
        center_of_geometry = np.mean(positions[ts])
        positions[ts] -= center_of_geometry


    x_avg = np.sum(positions, axis=0)/nSteps

    x1_average = np.full(nSteps,x_avg[0])

    x_t = np.zeros(nSteps)
    for i in range(nSteps):
        x_t[i] = positions[i,0]

    auto_correlationfs1 = np.zeros(nSteps)

    for framestep in range(nSteps):
        for frame in range(nSteps-framestep):
            auto_correlationfs1[framestep] += (x_t[frame]-x1_average[frame])*(x_t[frame+framestep]-x1_average[frame+framestep])
        auto_correlationfs1[framestep] /= nSteps-framestep

    # plot
    framestep = np.arange(nSteps)
    plt.plot(framestep*400,auto_correlationfs1)
    plt.xlabel("time(fs)")
    plt.xlim(0,7500)
    plt.ylabel("auto-correlation")
    plt.grid()
    plt.tight_layout()
    plt.savefig("auto_correlation.pdf",dpi=600)





