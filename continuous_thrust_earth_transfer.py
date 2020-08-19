from casadi import *
import constants

mu = constants.planets['earth']['gravitational_constant']
T = 20
m = 10
t_max = 10
n_int = 5
p = MX.sym('p')
f = MX.sym('f')
g = MX.sym('g')
h = MX.sym('h')
k = MX.sym('k')
L = MX.sym('L')

x = vertcat(p,f,g,h,k,L)

ur = MX.sym('ur')
ut = MX.sym('ut')
un = MX.sym('un')
u = vertcat(ur,ut,un)

w = 1 + (f*cos(L)) + (g*sin(L))
r = p / w
alpha = h**2 - k**2
beta = 1 + h**2 + k**2

P = T * u / m # no non-spherical Earth gravity perturbation or atmospheric drag added yet

pdot = (2*p/w)*sqrt(p/mu)*ut 
fdot = ((sqrt(p/mu)*sin(L))*ur) + ((sqrt(p/mu)*(1/w)*(f+(w+1)*cos(L)))*ut) + ((-sqrt(p/mu)*(g/w)*(h*sin(L)-k*cos(L)))*un)
gdot = ((-sqrt(p/mu)*cos(L))*ur) + ((sqrt(p/mu)*(g+(w+1)*sin(L)))*ut) + ((sqrt(p/mu)*(f/w)*(h*sin(L)-k*cos(L)))*un)
hdot = un*sqrt(p/mu)*beta*cos(L)/(2*w)
kdot = un*sqrt(p/mu)*beta*sin(L)/(2*w)
Ldot = (un*sqrt(p/mu)*(h*sin(L)-k*cos(L))) + sqrt(mu*p)*((w/p)**2)
xdot = vertcat(pdot,fdot,gdot,hdot,kdot,Ldot)

# Objective term
Ll = dot(u,u)
# Formulate discrete time dynamics

dae = {'x':x, 'p':u, 'ode':xdot, 'quad':Ll}
opts = {'tf':t_max/n_int}
F = integrator('F', 'cvodes', dae, opts)

