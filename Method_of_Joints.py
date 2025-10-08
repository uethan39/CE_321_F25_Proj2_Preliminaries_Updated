#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 12:37:32 2021

@author: kendrick shepherd
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 12:37:32 2021

@author: kendrick shepherd
Updated to fix single- and two-bar force computations
"""

import sys
import numpy as np
import Geometry_Operations as geom

# Determine unknown bars at a node
def UnknownBars(node):
    unknown_bars = []
    for bar in node.bars:
        if not bar.is_computed:
            unknown_bars.append(bar)
    return unknown_bars

# Determine if a node is solvable
def NodeIsViable(node):
    unknown_bars = UnknownBars(node)
    if len(unknown_bars) == 0:
        return False
    elif len(unknown_bars) == 1:
        return True
    elif len(unknown_bars) == 2:
        # Node is solvable if it has at least one known reaction or one bar already computed
        has_known_force = False
        if (0 in node.ConstraintType() and not np.isnan(node.xforce_reaction)) or \
           (1 in node.ConstraintType() and not np.isnan(node.yforce_reaction)):
            has_known_force = True
        connected_computed = any(bar.is_computed for bar in node.bars if bar not in unknown_bars)
        return has_known_force or connected_computed
    else:
        return False

# Compute unknown force for a single bar at a node
def SumOfForcesInLocalX(node, bar):
    vec = geom.BarNodeToVector(node, bar)
    length = np.linalg.norm(vec)
    if length == 0:
        sys.exit("Zero-length bar detected")
    
    cos_theta = vec[0] / length
    sin_theta = vec[1] / length

    sum_fx = node.xforce_external + (node.xforce_reaction if not np.isnan(node.xforce_reaction) else 0)
    sum_fy = node.yforce_external + (node.yforce_reaction if not np.isnan(node.yforce_reaction) else 0)

    for other_bar in node.bars:
        if other_bar.is_computed and other_bar != bar:
            vec_o = geom.BarNodeToVector(node, other_bar)
            length_o = np.linalg.norm(vec_o)
            cos_o = vec_o[0] / length_o
            sin_o = vec_o[1] / length_o
            sum_fx += other_bar.axial_load * cos_o
            sum_fy += other_bar.axial_load * sin_o

    bar.axial_load = -(sum_fx * cos_theta + sum_fy * sin_theta)
    bar.is_computed = True

# Compute unknown forces for two unknown bars at a node using 2x2 system
def SumOfForcesInLocalY(node, unknown_bars):
    if len(unknown_bars) != 2:
        sys.exit("SumOfForcesInLocalY requires exactly two unknown bars")

    bar1, bar2 = unknown_bars

    vec1 = geom.BarNodeToVector(node, bar1)
    vec2 = geom.BarNodeToVector(node, bar2)

    # Cosine and sine along global x and y
    cos1 = vec1[0] / np.linalg.norm(vec1)
    sin1 = vec1[1] / np.linalg.norm(vec1)
    cos2 = vec2[0] / np.linalg.norm(vec2)
    sin2 = vec2[1] / np.linalg.norm(vec2)

    sum_fx = node.xforce_external + (node.xforce_reaction if not np.isnan(node.xforce_reaction) else 0)
    sum_fy = node.yforce_external + (node.yforce_reaction if not np.isnan(node.yforce_reaction) else 0)

    for other_bar in node.bars:
        if other_bar.is_computed:
            vec_o = geom.BarNodeToVector(node, other_bar)
            length_o = np.linalg.norm(vec_o)
            cos_o = vec_o[0] / length_o
            sin_o = vec_o[1] / length_o
            sum_fx += other_bar.axial_load * cos_o
            sum_fy += other_bar.axial_load * sin_o

    A = np.array([[cos1, cos2],
                  [sin1, sin2]])
    b = np.array([-sum_fx, -sum_fy])

    try:
        F = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        sys.exit("Cannot solve forces at node: bars may be collinear or degenerate")

    bar1.axial_load = F[0]
    bar2.axial_load = F[1]
    bar1.is_computed = True
    bar2.is_computed = True

# Main iterative solver
def IterateUsingMethodOfJoints(nodes, bars):
    max_iterations = 5000
    iteration = 0

    while any(not bar.is_computed for bar in bars):
        iteration += 1
        if iteration > max_iterations:
            sys.exit("Exceeded maximum iterations, possible truss inconsistency")

        progress_made = False

        for node in nodes:
            if NodeIsViable(node):
                unknown_bars = UnknownBars(node)
                if len(unknown_bars) == 1:
                    SumOfForcesInLocalX(node, unknown_bars[0])
                elif len(unknown_bars) == 2:
                    SumOfForcesInLocalY(node, unknown_bars)
                progress_made = True

        if not progress_made:
            sys.exit("No solvable nodes found this iteration. Check truss geometry or constraints.")