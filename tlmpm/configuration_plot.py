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

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np  

def draw_rectangle(ax,x0,x1,y0,y1,**kwargs):
    ax.plot([x0,x0,x1,x1,x0],[y0,y1,y1,y0,y0], **kwargs)


def plot_configuration(ax,particles, grid, 
                fig, hdisplay, 
                params, 
                use_colormap=False, 
                cmin=None, cmax=None, 
                Qp=None, cbar=None, 
                cmap=None, label=None,
                show_labels=False,
                draw_rect=True,
                draw_nodes=True):

    # Collect needed quantities for plotting
    xnodes = grid['X']
    xparts = particles['x']
    Xparts = particles['X']
    dxcells = grid['dx']
    cellbounds = grid['cellbounds']
    xbounds = np.unique(grid['X'][cellbounds].flatten())
    xbounds = np.sort(xbounds)
    dx = grid['dx']
    nodeids = grid['nodeids']
    cellids = grid['cellids']
    partids = particles['ids']
    nnodes = len(nodeids)
    ncells = len(cellids)
    A = params['A']
    cross_len = np.sqrt(A)
    dxparts = particles['dx']
    dXparts = particles['dX']

    ones_ = np.ones(2)
    range_ = np.array([-1.0,1.0])

    ax.set_xticks(xbounds)
    ax.set_yticks([])

    # Plot cell quantities and text
    for cid in cellids:
        cellbound = grid['cellbounds'][cid]
        x0, x1 = grid['X'][cellbound]
        dx_cell = grid['dx'][cid]
        if show_labels: ax.text(x0 + 0.5*dx_cell, -0.85*cross_len, 'c'+str(cid+1))    

    # Cell node lines: vertical
    for ci in xbounds:
        if draw_nodes: ax.plot(ci*ones_, cross_len*range_,color='k')
    # Cell node lines: horizontal
    for x0, x1 in zip(xnodes[:-1],xnodes[1:]):
        if draw_nodes: ax.plot([x0,x1],-ones_*cross_len,color='k')
        if draw_nodes: ax.plot([x0,x1],ones_*cross_len,color='k')
    
    # Plot node quantities and text
    for xi in xnodes:
        if draw_nodes: ax.plot(xi*ones_,cross_len*range_,color='k',linestyle=':')
    # Node labels
    for nodeid, xnode in zip(nodeids, xnodes[:]):
        if show_labels: plt.text(xnode, 0.35*cross_len, 'n'+str(nodeid+1))

    current_color='blue'
    if use_colormap:
        norm = mcolors.Normalize(cmin, cmax)
        colormap = cm.get_cmap(cmap)
        rgba_values = [colormap(norm(v)) for v in Qp]
        
    # Plot particle quantities and text
    for pid, xpart, Xpart, dxpart, dXpart in zip(partids, xparts, Xparts, dxparts, dXparts):
        # plt.text(xpart, -0.93*cross_len, np.round(xpart,4))

        # Method to color based on stress
        if use_colormap:
            current_color = rgba_values[pid]
            sm = cm.ScalarMappable(cmap=colormap, norm=norm)
            sm.set_array([])
            
        # Shade in reference volume
        ax.fill_between([Xpart - 0.5*dXpart, Xpart + 0.5*dXpart], -0.5*cross_len, 0.5*cross_len, color='gray', alpha=0.075)

        # Draw rectangle around reference volume
        if draw_rect:
            draw_rectangle(ax, Xpart - 0.5*dXpart, Xpart + 0.5*dXpart, *(np.array([-0.5,0.5])*cross_len),
                        color='grey', alpha=0.4, linewidth=1, linestyle=':')
        if show_labels: plt.text(Xpart, 0.05*cross_len, 'p'+str(pid+1), color='grey')
        # plt.plot(Xpart*ones_, 0.5*(range_ - 1)*cross_len, linestyle=':', color='grey', alpha=0.75)

        # Shade in current volume
        ax.fill_between([xpart - 0.5*dxpart, xpart + 0.5*dxpart], -0.5*cross_len, 0.5*cross_len, color=current_color, alpha=0.75)

        # Draw rectangle around current volume
        if draw_rect:
            draw_rectangle(ax, xpart - 0.5*dxpart, xpart + 0.5*dxpart, *(np.array([-0.5,0.5])*cross_len),
                        color='k', linewidth=0.7, linestyle='--')
        if show_labels: plt.text(xpart, 0.05*cross_len, 'p'+str(pid+1))
        if draw_nodes: ax.plot(xpart*ones_, 0.5*(range_ - 1)*cross_len, linestyle=':', color='b')

    # Draw circles for nodes, reference particle positions, and current position
    if draw_nodes: ax.scatter(xnodes,xnodes*0+0.3*cross_len,color='r',label='Nodes')
    if draw_nodes: ax.scatter(Xparts, Xparts*0, marker='o', facecolors='none', edgecolors='b')
    if draw_nodes: ax.scatter(xparts, xparts*0, color='b', facecolors='b')
