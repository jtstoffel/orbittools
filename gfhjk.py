import math_tools as mt
from casadi import *
import constants


mu = constants.planets['earth']['gravitational_constant']
a0 = 24443 # km 
e0 = 0.725 
i0 = 7 # deg
w0 = 0 # deg
raas0 = 0 # deg
v0 = 0 # deg
x0_classical = [a0,e0,i0,w0,raas0,v0]
x0_equinoctial = mt.classical_to_equinoctial(x0_classical)

af = 42164
ef = 0
i_f = 0

control_int = 100
Nstates = 6
Ncontrols = 3
opti = casadi.Opti()
X = opti.variable(Nstates,control_int+1)
U = opti.variable(Ncontrols,control_int)
p = X[0,:]
f = X[1,:]
g = X[2,:]
h = X[3,:]
k = X[4,:]
L = X[5,:]

def system_dynamics(x,u):
    p = x[0]
    f = x[1]
    g = x[2]
    h = x[3]
    k = x[4]
    L = x[5]
    ur = u[0]
    ut = u[1]
    un = u[2]
    w = 1 + (f*cos(L)) + (g*sin(L))
    r = p / w
    alpha = h**2 - k**2
    beta = 1 + h**2 + k**2
    pdot = (2*p/w)*sqrt(p/mu)*ut 
    fdot = ((sqrt(p/mu)*sin(L))*ur) + ((sqrt(p/mu)*(1/w)*(f+(w+1)*cos(L)))*ut) + ((-sqrt(p/mu)*(g/w)*(h*sin(L)-k*cos(L)))*un)
    gdot = ((-sqrt(p/mu)*cos(L))*ur) + ((sqrt(p/mu)*(g+(w+1)*sin(L)))*ut) + ((sqrt(p/mu)*(f/w)*(h*sin(L)-k*cos(L)))*un)
    hdot = un*sqrt(p/mu)*beta*cos(L)/(2*w)
    kdot = un*sqrt(p/mu)*beta*sin(L)/(2*w)
    Ldot = (un*sqrt(p/mu)*(h*sin(L)-k*cos(L))) + sqrt(mu*p)*((w/p)**2)
    xdot = vertcat(pdot,fdot,gdot,hdot,kdot,Ldot)
    return xdot

tf = opti.variable()
opti.minimize(tf)
dt = tf / control_int
print(range(control_int))
for K in range(control_int):
    k1 = system_dynamics(X[:,K] , U[:,K])
    k2 = system_dynamics(X[:,K]+dt/2*k1 , U[:,K])
    k3 = system_dynamics(X[:,K]+dt/2*k2, U[:,K])
    k4 = system_dynamics(X[:,K]+dt*k3,   U[:,K])
    x_next = X[:,K] + dt/6*(k1+2*k2+2*k3+k4)
    opti.subject_to(X[:,K+1]==x_next) 

# 
#opti.subject_to(X[:,0]==x0_equinoctial)

# Final orbit state constraints
#opti.subject_to(p[control_int]==af*(1-ef**2))
#opti.subject_to( sqrt(f[control_int]**2 + g[control_int]**2) == ef)
#opti.subject_to( sqrt(h[control_int]**2 + k[control_int]**2) == tan(i_f/2))
#opti.subject_to(dot(U,U) == 1)
#opti.subject_to(tf>=0)

#opti.solver('ipopt')
#sol = opti.solve()


