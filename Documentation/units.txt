Units in iPic3D

---------------------------------------------------------------------------------------------------
1) The input file.

Everything is normalized to ion plasma parameters.
So that in the input file Lx, Ly, Lz are in units of ion skin depth:

[L] == d_i.

The time unit is reversed ion plasma frequency:

[t] == 1 / f_{p,i}.

Hence, the light speed is unity:

c == 1 (given in the inputfile).

Charge is normalized to e:

[e] == 1.

Mass:

m_i == 1

Density:

[rho] = 1/4/pi (normalization is done in EMfields Init* method)

Energy:

[Energy] = m_i * c^2

---------------------------------------------------------------------------------------------------
2) Pressures.

If you need to compute, e.g., plasma beta you have to compute the gas and the magnetic pressures:

p_{gas} = u_{th,i}^2 * \rho_i / (4*pi)   +   (m_e/m_i) * u_{th,e}^2 * \rho_e / (4*pi)
(where the first term in rhs is the ion pressure, and the second term is the electron pressure).
Note, that you should take the value of a single (x, y, or z) component of the thermal velocity!

p_{mag} = B^2 / (8*pi)


For those who forget, I remind that

\{beta} = p_{gas} / p_{mag}.

---------------------------------------------------------------------------------------------------
3) Magnetic permittivity (permeability) of vacuum

\mu_0 = 1. 

---------------------------------------------------------------------------------------------------
4) Alfven speed equals to the B0 in the inputfile

v_A = B_0
---------------------------------------------------------------------------------------------------

5) Currents in iPIC3D (exactly Faraday's law):

j = curl(B)*c/FourPI
---------------------------------------------------------------------------------------------------

6) Velocity distribution

kT/m = vth
where k is Boltzmann constant, T is temperature, m is particle mass for specie; 
vth is initial thermal speed for specie (one component!).

So the distribution of the absolute value of velocity (u) should follow the Maxwellian:
F = 4.*pi * (0.5/pi/uth**2)**(1.5) * u**2 * exp(-0.5*u**2/uth**2) * Ns * du,
where uth is thermal speed for the specie, Ns=npcelx*npcely*npcelz, du is velocity bin size.
---------------------------------------------------------------------------------------------------