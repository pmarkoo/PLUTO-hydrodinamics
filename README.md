# PLUTO-hydrodinamics
Hydrodynamic Simulations of Supernova Remnants
We modelled the dynamical evolution of SNRs by numerically solving the time-dependent Euler partial differential equations (PDEs) of fluid dynamics, also known as hyperbolic conservation laws, that we write as:

\begin{equation}\label{eq:CL}
  \pd{\vec{U}}{t} + \nabla\cdot\vec{F} = 0.
\end{equation}

Here $\vec{U}$ and $\vec{F}$ represent a state and flux vectors, respectively, 
which can be written in the form:
