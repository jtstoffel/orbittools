import numpy as np

def classical_to_equinoctial(x):
    semimajor_axis = a = x[0]
    eccentricity = e = x[1]
    inclination = i = x[2]
    argument_of_perigee = w = x[3]
    right_ascension_of_ascending_node = raas = x[4]
    true_anomaly = v = x[5]
    p = a * (1 - e**2)
    f = e * np.cos(w + raas)
    g = e * np.sin(w + raas)
    h = np.tan(i / 2.0) * np.sin(raas)
    k = np.tan(i / 2.0) * np.cos(raas)
    L = raas + w + v
    return np.array([p, f, g, h, k, L])

def equinoctial_to_ECI(x,mu):
    p = x[0]
    f = x[1]
    g = x[2]
    h = x[3]
    k = x[4]
    L = x[5]
    alpha = h**2 - k**2
    beta = 1 + h**2 + k**2
    w = 1 + (f * np.cos(L)) + (g * np.sin(L))
    r = p/w
    x = (r / beta) * (np.cos(L) + (alpha * np.cos(L)) + (2 * h * k * np.sin(L)))
    y = (r / beta) * (np.sin(L) - (alpha * np.sin(L)) + (2 * h * k * np.cos(L)))
    z = (2 * r / beta) * ((h * np.sin(L)) - (k * np.cos(L)))
    vx = (-1 / beta) * np.sqrt(mu / p) * (np.sin(L) + (alpha * np.sin(L)
                                                       ) - (2 * h * k * np.cos(L)) + g - (2 * f * h * k) + (alpha * g))
    vy = (-1 / beta) * np.sqrt(mu / p) * (-np.cos(L) + (alpha * np.cos(L)
                                                        ) + (2 * h * k * np.cos(L)) - f + (2 * g * h * k) + (alpha * f))
    vz = (2 / beta) * np.sqrt(mu / p) * \
        ((h * np.cos(L)) + (k * np.sin(L)) + (f * h) + (g * k))
    return np.array([x, y, z , vx, vy, vz])

def ROR(phi):
    return np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])

def ROI(h,hx,hy):
    hz = np.sqrt(h**2 - hx**2 - hy**2)
    r11 = hz/np.sqrt(h**2-hy**2)
    r12 = -hx*hy/(h*np.sqrt(h**2-hy**2))
    r13 = hx/h
    r21 = 0
    r22 = (h**2-hy**2)/(h*np.sqrt(h**2-hy**2))
    r23 = hy/h
    r31 = -hx/np.sqrt(h**2-hy**2)
    r32 = -hy*hz/(h*np.sqrt(h**2-hy**2))
    r33 = hz/h
    return np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])

def xyz_position(h,hx,hy,ex,ey,phi,mu):
    r = h**2/(mu*(1+(ex*np.cos(phi))+(ey*np.sin(phi))))
    ror = ROR(phi)
    roi = ROI(h,hx,hy)
    return roi @ ror @ np.array([[r],[0],[0]])

def classical_to_momentum_eccentricity(x,mu):
    a = x[0] # km 
    e = x[1] 
    h = np.sqrt(a*mu*(1-e**2))
    xyz = equinoctial_to_ECI(classical_to_equinoctial(x))
    R= xyz[:2]
    V = xyz[-3:]
    H = np.cross(R,V)
    hx = H[0]
    hy = H[1]
    E = np.transpose(ROI(h,hx,hy)) @ ((np.cross(V,H)/mu)-(R/np.norm(R)))
    ex = E[0]
    ey = E[1]
    return np.array([h,hx,hy,ex,ey])
