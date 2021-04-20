from casadi import *
import math_tools as mt
import constants as con

#################################################################
p_i = 21837080.052835
f_i = 0
g_i = 0
h_i = -0.25396764647494
k_i = 0
L_i = pi
w_i = 1

e_f = 0.73550320568829 # final eccentricity 
Isp = 450 # specific impulse

tau = -50 # -99 < tau < 0
mu = 1.407645794e16
J2 = 1082.639e-6
J4 = -1.608e-6
p_f = 40007346.015232 # final p
tanif2 = 0.6176125876099
g0 = 32.174 # ft/s^2
T = 4.446618e-3 # max thrust
Re = 20925662.73

#################################################################



# States
Nstates = 7
p = MX.sym('p')
f = MX.sym('f')
g = MX.sym('g')
h = MX.sym('h')
k = MX.sym('k')
L = MX.sym('L')
w = MX.sym('w')
y = vertcat(p,f,g,h,k,L,w)

# Controls
Ncontrols = 3
ur = MX.sym('ur')
ut = MX.sym('ut')
un = MX.sym('un')
u = vertcat(ur,ut,un)

# Equations of Motion
q = 1 + (f*cos(L)) + (g*sin(L))
chi = sqrt(h**2 + k**2)
ss = 1 + chi**2
alpha = h**2 - k**2

delta_T = g0 * T * (1 + 0.01*tau)  * u / w
delta = delta_T

pdot = (2 * p * sqrt(p/mu) / q) * delta_T[1]
fdot = ((sqrt(p/mu)*sin(L))*delta[0]) + ((sqrt(p/mu)*(((q + 1) * cos(L)) + f) / q)*delta[1]) + ((-sqrt(p/mu) * (g/q) * ((h*sin(L)) - (k*cos(L))))*delta[2])
gdot = ((-sqrt(p/mu)*cos(L))*delta[0]) + ((sqrt(p/mu)*(((q + 1) * sin(L)) + g) / q)*delta[1]) + ((sqrt(p/mu) * (f/q) * ((h*sin(L)) - (k*cos(L))))*delta[2]) 
hdot = (sqrt(p/mu) * ss * cos(L) / (2*q))*delta[2]
kdot = (sqrt(p/mu) * ss * sin(L) / (2*q))*delta[2]
Ldot = (sqrt(p/mu) * (1/q) * ((h*sin(L)) - (k*cos(L))))*delta[2] + (sqrt(mu*p) * ((q/p)*(q/p)))
wdot = -T * (1 + 0.01*tau) / Isp
ydot = vertcat(pdot,fdot,gdot,hdot,kdot,Ldot,wdot)


r = p/q
r_vec = MX.sym('r_vec',3)
r_vec[0] = (r/ss) * ((cos(L)) + (alpha*cos(L)) + (2*h*k*sin(L)))
r_vec[1] = (r/ss) * ((sin(L)) - (alpha*sin(L)) + (2*h*k*cos(L)))
r_vec[2] = (2*r/ss) * ((h*sin(L)) - (k*cos(L)))

v_vec = MX.sym('v_vec',3)
v_vec[0] = (-1/ss) * sqrt(mu/p) * (sin(L) + (alpha*sin(L)) - (2*h*k*cos(L)) + g - (2*f*h*k) + (alpha*g))
v_vec[1] = (-1/ss) * sqrt(mu/p) * (-cos(L) + (alpha*cos(L)) + (2*h*k*sin(L)) - f + (2*g*h*k) + (alpha*f))
v_vec[2] = (2/ss) * sqrt(mu/p) * ((h*cos(L)) + (k*sin(L)) + (f*h) + (g*k))


# Continous time dynamics
F = Function('F', [y, u], [ydot], ['y', 'u'], ['ydot'])

# Discritization
Ngridpoints = 1000

# Problem definition
opti = Opti()
tf = opti.variable()
dt = tf/Ngridpoints
Y = opti.variable(Nstates,Ngridpoints)
Ymid = opti.variable(Nstates,Ngridpoints-1)
U = opti.variable(Ncontrols,Ngridpoints)
Umid = opti.variable(Ncontrols,Ngridpoints-1)

opti.set_initial(Y[0,:], p_i)
opti.set_initial(Y[1,:], f_i)
opti.set_initial(Y[2,:], g_i)
opti.set_initial(Y[3,:], h_i)
opti.set_initial(Y[4,:], L_i)
opti.set_initial(Y[5,:], w_i)

# State Constraints
opti.subject_to(vec(Y[0,:]) >= 0.1 * p_i)
opti.subject_to(vec(Ymid[0,:]) >= 0.1 * p_i)
opti.subject_to(vec(Y[0,:]) <= 5 * p_f)
opti.subject_to(vec(Ymid[0,:]) <= 5 * p_f)
opti.subject_to(vec(Y[1:5,:]) >= -1)
opti.subject_to(vec(Ymid[1:5,:]) >= -1)
opti.subject_to(vec(Y[1:5,:]) <= 1)
opti.subject_to(vec(Ymid[1:5,:]) <= 1)
opti.subject_to(vec(Y[5,:]) >= pi)
opti.subject_to(vec(Ymid[5,:]) >= pi)
opti.subject_to(vec(Y[5,:]) <= 34*pi)
opti.subject_to(vec(Ymid[5,:]) <= 34*pi)
opti.subject_to(vec(Y[6,:]) >= 0.001)
opti.subject_to(vec(Ymid[6,:]) >= 0.001)
opti.subject_to(vec(Y[6,:]) <= 1.01)
opti.subject_to(vec(Ymid[6,:]) <= 1.01)

# Initial Conditions
opti.subject_to(Y[0,0] == p_i)
opti.subject_to(Y[1,0] == f_i)
opti.subject_to(Y[2,0] == g_i)
opti.subject_to(Y[3,0] == h_i)
opti.subject_to(Y[4,0] == k_i)
opti.subject_to(Y[5,0] == L_i)
opti.subject_to(Y[6,0] == w_i)

# Final Conditions
opti.subject_to(Y[0,-1] == p_f)
opti.subject_to(sqrt(Y[1,-1]**2 + Y[2,-1]) == e_f) # final eccentricity constraint
opti.subject_to(sqrt(Y[3,-1]**2 + Y[4,-1])**2 == tanif2) # final inclination constraint
opti.subject_to(Y[1,-1]*Y[3,-1] + Y[2,-1]*Y[4,-1] == 0)
opti.subject_to(Y[2,-1]*Y[3,-1] - Y[4,-1]*Y[1,-1] <= 0)

# Control Constraints
opti.subject_to(vec(U) >= -1.1)
opti.subject_to(vec(Umid) >= -1.1)
opti.subject_to(vec(U) <= 1.1)
opti.subject_to(vec(Umid) <= 1.1)


# Hermite Simpson Constraints

for k in range(1,Ngridpoints-1):
    opti.subject_to(vec(Ymid[:,k-1] - (0.5 * (Y[:,k] - Y[:,k-1])) - ((dt/8)*(F(Y[:,k-1],U[:,k-1]) - F(Y[:,k],U[:,k])))) == 0)          # Hermite
    opti.subject_to(vec(Y[:,k] - Y[:,k-1] - ((dt/6)*(F(Y[:,k],U[:,k]) - (4 * F(Ymid[:,k-1],Umid[:,k-1])) - F(Y[:,k-1],U[:,k-1])))) == 0 ) # Simpson



# Objective
opti.minimize(-Y[6,-1])

# Solver
opti.solver('ipopt')

# Solution
sol = opti.solve()
