"""

    1D electrostatic particle-in-cell solver for electron two-stream instability.


    for the course Introduction to Plasma Dynamics by Giovanni Lapenta, KULeuven.

"""

import os
import time

#start_time = time.clock()

import numpy as np   #array syntax

import pylab as plt  #plot

import matplotlib

import matplotlib.patches as mpatches   #plot

import scipy

import scipy.fftpack

from scipy import sparse   #special functions, optimization, linear algebra

from scipy.sparse import linalg


# Output folder

path = './Results'

if not os.path.exists(path):

    os.makedirs(path)



# Set plotting parameters

params = {'axes.labelsize': 'large',

              'xtick.labelsize': 'medium',

              'ytick.labelsize': 'medium',

              'font.size': 15,

              'font.family': 'sans-serif',

              'text.usetex': False,

              'mathtext.fontset': 'stixsans',}

plt.rcParams.update(params)

## Switch on interactive plotting mode

plt.ion()

# Simulation parameters

L = 3*2*np.pi/(np.sqrt(3.0/2)/2.) # Domain size

DT = 0.02 # Time step

NT = 5000  # Number of time steps

TOut = round(NT/100) # Output period

verbose = True

NG = 384 # Number of grid cells

N = NG * 20 # Number of particles

WP = 1. # Plasma frequency

QM = -1. # Charge/mass ratio

V0 = .9 # Stream velocity

VT = 0.00000001 # Thermal speed



# perturbation

XP1 = 0.0001

mode = 3



Q = WP**2 / (QM*N/L)  # rho0*L/N: charge carried by a single particle?

rho_back = -Q*N/L  # Background charge density?

dx = L / NG # Grid step



# Auxilliary vectors

p = np.concatenate([np.arange(N), np.arange(N)])  # Some indices up to N

Poisson = sparse.spdiags(([1, -2, 1] * np.ones((1, NG-1), dtype=int).T).T, [-1, 0, 1], NG-1, NG-1)

Poisson = Poisson.tocsc()


# Cell center coordinates

xg = np.linspace(0, L-dx, NG) + dx/2


# electrons

xp = np.linspace(0, L-L/N, N).T   # Particle positions

vp = VT * np.random.randn(N) # particle thermal spread

pm = np.arange(N)

pm = 1 - 2 * np.mod(pm+1, 2)

vp += pm * V0 # Drift plus thermal spread



# Add electron perturbation to excite the desired mode

xp += XP1 * np.cos(2 * np.pi * mode / L * xp)

xp[np.where(xp < 0)] += L

xp[np.where(xp >= L)] -= L



histEnergy, histPotE, histKinE, t = [], [], [], []



if verbose:

    plt.figure(1, figsize=(16,9))



# Main cycle

for it in range(NT+1):

    # update particle position xp

    xp += vp * DT

    # Periodic boundary condition

    xp[np.where(xp < 0)] += L

    xp[np.where(xp >= L)] -= L



    # Project particles->grid

    g1 = np.floor(xp/dx - 0.5)

    g = np.concatenate((g1, g1+1))

    fraz1 = 1 - np.abs(xp/dx - g1 - 0.5)

    fraz = np.concatenate((fraz1, 1-fraz1))

    g[np.where(g < 0)] += NG

    g[np.where(g > NG-1)] -= NG



    mat = sparse.csc_matrix((fraz, (p, g)), shape=(N, NG))

    rho = Q / dx * mat.toarray().sum(axis=0) + rho_back



    # Compute electric field potential

    Phi = linalg.spsolve(Poisson, -dx**2 * rho[0:NG-1])

    Phi = np.concatenate((Phi,[0]))

    # Electric field on the grid

    Eg = (np.roll(Phi, 1) - np.roll(Phi, -1)) / (2*dx)

    # Electric field fft

    #ft = abs(scipy.fft(Eg))
    #k = scipy.fftpack.fftfreq(Eg.size,xg[1]-xg[0])

    # interpolation grid->particle and velocity update

    vp += mat * QM * Eg * DT

    bins,edges=np.histogram(vp,bins=40,range=(-3.2,3.2))
    left,right = edges[:-1],edges[1:]
    vc = np.array([left,right]).T.flatten()
    fv = np.array([bins,bins]).T.flatten()

    Etot = 0.5 * (Eg**2).sum() * dx

    histEnergy.append(Etot+0.5 * Q/QM * (vp**2).sum())

    histPotE.append(0.5 * (Eg**2).sum() * dx)

    histKinE.append(0.5 * Q/QM * (vp**2).sum())


    t.append(it*DT)



    if (np.mod(it, TOut) == 0) and verbose:

        # Phase space

        plt.clf()

        plt.subplot(2, 2, 1)

        plt.scatter(xp[0:-1:2], vp[0:-1:2], s=0.5, marker='.', color='blue')

        plt.scatter(xp[1:-1:2], vp[1:-1:2], s=0.5, marker='.', color='red')

        plt.xlim(0, L)

        plt.ylim(-4, 4)

        plt.xlabel('x')

        plt.ylabel('v')

        plt.legend((mpatches.Patch(color='w'), ), (r'$\omega_{pe}t=$' + str(DT*it), ), loc=1, frameon=False)



        # Electric field

        plt.subplot(2, 2, 2)

        plt.xlim(0, L)

        plt.ylim(-1.5, 1.5)

        plt.xlabel('x')

        plt.plot(xg, Eg, label='E', linewidth=2)

        plt.legend(loc=1)




        # Energies

        plt.subplot(2, 2, 3)

        plt.xlim(0, NT*DT)

        plt.ylim(1e-12, 100)

        plt.xlabel('time')

        plt.yscale('log')

        plt.plot(t, histPotE, label='Potential', linewidth=2)

        plt.plot(t, histKinE, label='Kinetic', linewidth=2)

        plt.plot(t, histEnergy, label='Total Energy', linestyle='--', linewidth=2)

        plt.legend(loc=4)



        # Electron distribution function

        plt.subplot(2, 2, 4)

        plt.xlim(-4, 4)

        #plt.ylim(0, N/2)

        plt.xlabel('v')

        plt.plot(vc,fv, label='f(v)', linewidth=2)

        plt.legend(loc=1)

        plt.pause(0.000000000000001)

        print(it)

        plt.savefig(os.path.join(path, 'twostream%3.3i' % (it/TOut,) + '.png'))


np.savetxt('field_ene.txt',(t,histPotE))

#print ('Time elapsed: ', time.clock() - start_time)

# Comment this line if you want the figure to automatically close at the end of the simulation
raw_input('Press enter...')
