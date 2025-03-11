# MIT License
#
# Copyright (c) 2025 Jacob Nuttall
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# 
#   CoolMPM : A Total Lagrangian MPM Code
# 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # 

import numpy as np
from copy import deepcopy
from tlmpm.fem_shape import XI, N, B, X, dXI_dx

################ INITIALIZATION #############################
def initialize_grid(params,n_cells=3,n_ppc=3,n_en=3):
    L = params['L']
    
    n_nodes = (n_en-1)*n_cells + 1
    xnodes, dx = np.linspace(0,L,n_nodes, retstep=True)

    gridn = np.array([np.arange(n*(n_en-1),(n_en-1)*n+n_en) for n in range(n_cells)])
    cellbounds = gridn[:,[0,-1]]
    xnodes, xnodes[gridn], gridn
    
    xbounds = xnodes[cellbounds]
    dx = xbounds[:, 1] - xbounds[:, 0]
    h = L/n_cells

    grid = {
        'n_cells': n_cells,
        'n_el': n_cells,
        'n_en': n_en, 
        'n_ppc': n_ppc, 
        'n_nodes': n_nodes,
        'cellids': np.arange(n_cells),
        'nodeids': np.arange(n_nodes),
        'cells': np.copy(gridn),
        'cellbounds': np.copy(cellbounds),
        'h': h,
        'dx': dx,
        'X': np.copy(xnodes), 
        'mass': xnodes*0, 
        'int_force': xnodes*0, 
        'ext_force': xnodes*0,
        'drag_force': xnodes*0,
        'force': xnodes*0, 
        'v': xnodes*0, 
        'v_tilde': xnodes*0,
        'v_new': xnodes*0,
        'momentum': xnodes*0,
        'momentum_tilde': xnodes*0,
        'momentum_new': xnodes*0,
        'd': xnodes*0, 
        'el_length': np.full(n_cells, dx),
    }
    return grid 


def initialize_particles(params, grid):
    rho = params['rho']
    A = params['A']
    n_ppc = grid['n_ppc']
    
    particles = {
        'ids': [],
        'n_mp': 0,
        'X': [],
        'dX': [],
        'x': [],
        'dx': [],
        'u': [],
        'v': [],
        'momentum': [],
        'l_grad': [],
        'defrate': [],
        'defgrad': [],
        'stress': [],
        'FPK_stress': [],
        'SPK_stress': [],
        'volume': [],
        'ref_volume': [],
        'mass': [],
        'rho': [],
        'e11': [],
        'e22': [],
        'e33': [],
        'le11': [],
        'E11': [],
        't': [],
    }

    id = -1
    
    for cid in grid['cellids']:

        for _ in range(n_ppc):
            id += 1
            particles['ids'].append(id)
            particles['n_mp'] += 1
        
        bounds = grid['cellbounds'][cid]
        x_n = grid['X'][bounds]
        x_p, dx = np.linspace(*x_n, 2*n_ppc+1, retstep=True)
        x_p = x_p[1:-1:2]
        dx = 2*dx 

        zeros = np.zeros_like(x_p)
        ones = np.ones_like(x_p)
        dx_ = ones*dx 
        dX_ = ones*dx 
        mass_ = ones*rho*A*dx 
        rho_ = ones*rho 
        vol_ = ones*A*dx 
        Vol_ = ones*A*dx

        particles['x'].extend(x_p.copy())
        particles['X'].extend(x_p.copy()) 
        particles['u'].extend(np.zeros_like(x_p))
        particles['dX'].extend(dX_)
        particles['dx'].extend(dx_)
        particles['v'].extend(zeros.copy())
        particles['l_grad'].extend(zeros.copy())
        particles['defgrad'].extend(ones.copy())
        particles['defrate'].extend(zeros.copy())
        particles['mass'].extend(mass_.copy())
        particles['rho'].extend(rho_.copy())
        particles['ref_volume'].extend(Vol_.copy())
        particles['volume'].extend(vol_.copy())
        particles['FPK_stress'].extend(zeros.copy())
        particles['SPK_stress'].extend(zeros.copy())
        particles['stress'].extend(zeros.copy())
        particles['e11'].extend(zeros.copy())
        particles['e22'].extend(zeros.copy())
        particles['e33'].extend(zeros.copy())
        particles['le11'].extend(zeros.copy())
        particles['E11'].extend(zeros.copy())
        particles['t'].extend(zeros.copy())

    for key, item in particles.items():
        particles[key] = np.array(item)
    
    return particles

def initialize_history(particles):
    particles_history = {}
    for key in particles.keys():
        particles_history[key] = [deepcopy(particles[key])]
    return particles_history

def update_history(particles, particles_history):
    for key in particles_history.keys():
        particles_history[key].append(deepcopy(particles[key]))

def calc_cid(pid, particles, grid):
    h = grid['h']
    Xp = particles['X'][pid]
    e = int(np.floor(Xp/h))
    return e

def particle_in_cell(pid, cid, particles, grid):
    cellbound = grid['cellbounds'][cid]
    x1, x2 = grid['X'][cellbound]
    xp = particles['X'][pid]
    return x1 <= xp <= x2 

def right_particle_pid(particles):
    return np.argmax(particles['x'])

############ PARTICLES TO GRID ##########################################
def reset_grid(grid):
    grid['mass'] *= 0
    grid['int_force'] *= 0
    grid['ext_force'] *= 0 
    grid['drag_force'] *= 0
    grid['momentum_new'] *= 0 
    grid['momentum_tilde'] *= 0
    grid['momentum'] *= 0
    grid['v'] *= 0 
    grid['v_new'] *= 0

def p2g_momentum(pid, cid, particles, grid):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    shape = N(n_en, xi)

    m_p = particles['mass'][pid]
    v_p = particles['v'][pid]     

    grid['momentum'][cell] += m_p*v_p*shape 

def p2g_ext_force(pid, cid, particles, grid, params, t):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    shape = N(n_en, xi)

    X_p = particles['X'][pid]
    m_p = particles['mass'][pid]

    # Calc concentrated force / traction
    conc_force = 0
    if pid == right_particle_pid(particles):

        # Calculate the rightmost edge of the rightmost particle volume
        traction_bound = X_p + 0.5*particles['dX'][pid]
        xi_bound = XI(traction_bound, *x_cell)
        shape_bound = N(n_en, xi_bound)
        conc_force = params['tau'](t)*params['A']*shape_bound

    # Calc body force
    body_force = shape*m_p*params['b'](t)

    grid['ext_force'][cell] += body_force + conc_force

def p2g_int_force(pid, cid, particles, grid):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    je = dXI_dx(*x_cell)
    shapegrad = B(n_en, xi)*je

    Vol_p = particles['ref_volume'][pid]
    FPK_p = particles['FPK_stress'][pid]
    grid['int_force'][cell] += -Vol_p*FPK_p*shapegrad

def p2g_mass(pid, cid, particles, grid):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    shape = N(n_en, xi)

    m_p = particles['mass'][pid]
    grid['mass'][cell] += shape*m_p 
    

    ########## GRID TO PARTICLES #################
    # Routines described in "grid-to-particles" and shown in TLMPM algorithm

def g2p_velocity(pid, cid, particles, grid, alpha):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    shape = N(n_en, xi)
    
    v_p = particles['v'][pid]
    v_g = np.dot(grid['v'][cell], shape) 
    v_tilde = np.dot(grid['v_tilde'][cell], shape)
    v_p = alpha*(v_p + (v_tilde - v_g)) + (1 - alpha)*v_tilde
    particles['v'][pid] = v_p

def g2p_defrate(pid, cid, particles, grid):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    je = dXI_dx(*x_cell)
    shapegrad = B(n_en, xi)*je

    v_g = grid['v'][cell]
    defgrad_rate_p = np.dot(shapegrad, v_g)
    particles['defrate'][pid] = defgrad_rate_p

def g2p_position(pid, cid, particles, grid, dt):
    n_en = grid['n_en']
    cell = grid['cells'][cid]
    cellbound = grid['cellbounds'][cid]
    x_cell = grid['X'][cellbound]
    X_p = particles['X'][pid]
    xi = XI(X_p, *x_cell)
    shape = N(n_en, xi)

    v_n = grid['v'][cell]    
    x_p_inc = dt*np.dot(v_n, shape)
    particles['x'][pid] += x_p_inc


############## EXTRAPOLATE GRID QUANTITIES ####################
# Weight material quantities by their mass to extrapolate to grid nodes
def extrapolate_quantity_p2g(Qp, particles, grid):
    n_en = grid['n_en']
    n_nodes = grid['n_nodes']
    m_n = np.zeros(n_nodes)
    Q_n = np.zeros(n_nodes)
    for pid in particles['ids']:

        cid = calc_cid(pid, particles, grid)
        cellbound = grid['cellbounds'][cid]
        cell = grid['cells'][cid]
        X_cell = grid['X'][cellbound]
        X_p = particles['X'][pid]
        xi = XI(X_p, *X_cell)
        shape = N(n_en, xi)
        Q_p = Qp[pid]
        m_p = particles['mass'][pid]
        Q_n[cell] += shape*Q_p*m_p
        m_n[cell] += shape*m_p
        
    return Q_n / m_n 


############### MAIN TLMMPM BODY ########################
def TotalLagrangianMPM(
        particles, grid, params, calc_SPK_stress, 
        alpha, dtfactor, T,
        plot=False, plot_func=None, nplots=10):
    """ The Total Lagrangian Algorithm for Axially-Loaded bar with zero lateral displacement in 1D

    Args:
        particles (dict): dictionary containing the current particle data
        grid (dict): dictionary containing the grid data
        params (dict): material property definitions and external force 
        calc_SPK_stress (callable): a function with signature (defgrad, params) returning the consitutive form of SPK stress
            For Neo-Hookean, calc_SPK_stress = lambda defgrad, params : mu + (dlambda*np.log(J) - mu)/defgrad^T.defgrad
        alpha (float): 0 < alpha < 1 : Newmark time integration parameter. Set to 0.5 for explicit central difference. 
        dtfactor (float): Multiplier for CFL condition. Recommended is dtfactor < 0.2.
        T (float): Total time to run simulation.
        plot (bool, optional): If we should plot the configuration. Defaults to False.
        plot_func (callable, optional): Function hook with signature (t,k) that manages the plotting. Defaults to None.
        nplots (int, optional): The mnumber of times to call plot_func. Defaults to 10.

    Returns:
        particles, grid, particle_history: particle_history contains all computed material point quantities over time.
    """
    
    # Initialization
    particle_history = initialize_history(particles)
    t = 0

    # Calculate dt size
    h = grid['h']
    wavespeed = params['wavespeed']
    dt = dtfactor*h/wavespeed
    k = 0

    Nt = int(np.ceil(T/dt))
    if plot: kplot = Nt//nplots

    while k < Nt:
        
        reset_grid(grid)

        if plot:
            if k % kplot == 0:
                plot_func(t,k)

        # ~~~ Particles to grid ~~~
        for pid in particles['ids']:
            cid = calc_cid(pid, particles, grid)
            p2g_mass(pid, cid, particles, grid)
            p2g_momentum(pid, cid, particles, grid)
            p2g_ext_force(pid, cid, particles, grid, params, t)
            p2g_int_force(pid, cid, particles, grid)

        #   Force-Balance Equation
        F_net = grid['force'] = grid['ext_force'] + grid['int_force']

        # ~~~ Update velocities (Double mapping) ~~~
        #   First momentum update
        m_n = grid['mass']
        phi_n = grid['momentum']
        phi_tilde_n = phi_n + dt * F_net 

        #   Set Dirichlet BC
        phi_n[0] = 0
        phi_tilde_n[0] = 0
        grid['momentum'] = phi_n 
        grid['momentum_tilde'] = phi_tilde_n

        #   Momentum update: Particlesphi
        grid['v'] = phi_n/m_n
        grid['v_tilde'] = phi_tilde_n/m_n
        
        for pid in particles['ids']:
            cid = calc_cid(pid, particles, grid)
            g2p_velocity(pid, cid, particles, grid, alpha)

        #   Momentum update: grid double mapping
        grid['momentum'][:] = 0
        for pid in particles['ids']:
            cid = calc_cid(pid, particles, grid)
            p2g_momentum(pid, cid, particles, grid)
            
        phi_n = grid['momentum']
        phi_n[0] = 0
        grid['momentum'] = phi_n 
        grid['v'] = phi_n/m_n
        
        # ~~~ Particles-to-Grid ~~~ 
        for pid in particles['ids']:
            cid = calc_cid(pid, particles, grid)
            g2p_position(pid, cid, particles, grid, dt)
            g2p_defrate(pid, cid, particles, grid)

            # Update def gradient F
            defrate = particles['defrate'][pid]
            defgrad = particles['defgrad'][pid] = particles['defgrad'][pid] + dt*defrate 

            J = defgrad # in 1D only

            right_CG = defgrad*defgrad
            b = defgrad*defgrad 
            b_inv = 1/b

            SPK_stress = calc_SPK_stress(defgrad, params)
            FPK_stress = defgrad*SPK_stress 
            stress = FPK_stress # in 1D only

            # Strain measures 
            E11 = 0.5*(right_CG - 1.0)
            e11 = 0.5*(1.0 - b_inv)
            le11 = np.log(defgrad)

            # Volume update
            A = params['A']
            ref_volume = particles['ref_volume'][pid]
            volume = J*ref_volume
            dx = volume/A

            # Update particle quantities
            particles['defgrad'][pid] = defgrad 
            particles['l_grad'][pid] = defrate/defgrad 
            particles['SPK_stress'][pid] = SPK_stress
            particles['FPK_stress'][pid] = FPK_stress
            particles['stress'][pid] = stress 
            particles['e11'][pid] = e11 
            particles['le11'][pid] = le11
            particles['E11'][pid] = E11 
            particles['t'][pid] = t
            particles['volume'][pid] = volume
            particles['dx'][pid] = dx 

        update_history(particles, particle_history)

        k+=1
        t+=dt

    else:
        if plot: plot_func(t, k)

    for key, item in particle_history.items():
        particle_history[key] = np.array(item)
        
    return grid, particles, particle_history

def select_data(particle_history, pid=0):
    """ Convience function to switch between particle data from a history.

    Args:
        particle_history (dict): dictionary containing time series of material properties
        pid (int, optional): particle id. Defaults to 0.
    """
    global E11, le11, e11, sig11, S11, ts, xs, Xs, dXs, dxs

    for key, item in particle_history.items():
        particle_history[key] = np.array(item)

    E11 = particle_history['E11'][:,pid]
    le11 = particle_history['le11'][:,pid]
    e11 = particle_history['e11'][:,pid]
    sig11 = particle_history['stress'][:,pid]
    S11 = particle_history['SPK_stress'][:,pid]
    ts = particle_history['t'][:,pid]
    xs = particle_history['x'][:,pid]
    Xs = particle_history['X'][:,pid]
    dXs = particle_history['dX'][:,pid]
    dxs = particle_history['dx'][:,pid]