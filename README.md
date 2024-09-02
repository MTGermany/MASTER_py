# MASTER_py - Python code for generic longitudinal macroscopic models
----------------------------------------------------------------------



## Description
--------------

MASTER_py is a very basic Python code for a generic longitudinal macroscopic model with relaxation and constant anticipation distance

Features:

- Basis for future 2d macroscopic models

- As is, the model poroduces a png image Density_Distribution.png as main output

- Upwind finite differences. The necessary downwind info is obtained by the anticipation term (d_t+V d_x)V=(Ve(rho_a)-V)/tau with rho_a(x,t)=rho(x+gamma*Vf*timegap,t)

- Important: Constant anticipation distance since, otherwise, instabilities near rho_max due to low or vanishing anticipation

- By formulating the relaxation as the local solution of the ODE rather than first-order Euler, LWR-type models can be simulated for relaxation times of 0.001 s, or so (Euler would lead to a local relaxational numerical instability for tau<0.5*dt_update)

- Together with catching density values >rho_max or <SMALL_VAL, the above features ensure that the simulation is stable for any initial conditions and a wide range of model and simulation parameters




## References 
-------------

[1] M. Treiber and A. Kesting. [Traffic Flow Dynamics, Data, Models and Simulation](http://www.traffic-flow-dynamics.org). Springer 2013. [Link](http://www.springer.com/physics/complexity/book/978-3-642-32459-8)

[2] D. Helbing, A. Hennecke, V. Shvetsov, and M. Treiber (2001)
MASTER: Macroscopic traffic simulation based on a gas-kinetic,
non-local traffic model. Transportation Research B 35, 183-211.

[3] M. Treiber, A. Hennecke, and D. Helbing (1999)
Derivation, Properties, and Simulation of a Gas-Kinetic-Based, Non-Local Traffic Model,
Phys. Rev. E 59, 239-253. 
