#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 12:37:32 2021

@author: kendrick shepherd
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
        # More than 2 unknowns → cannot solve directly
        return False

# Compute unknown force for a single bar at a node
def SumOfForcesInLocalX(node, bar):
    sum_x = node.xforce_external + (node.xforce_reaction if not np.isnan(node.xforce_reaction) else 0)
    sum_y = node.yforce_external + (node.yforce_reaction if not np.isnan(node.yforce_reaction) else 0)

    for other_bar in node.bars:
        if other_bar.is_computed:
            sum_x += other_bar.axial_load * geom.CosineBars(bar, other_bar)
            sum_y += other_bar.axial_load * geom.SineBars(bar, other_bar)

    bar.axial_load = -sum_x
    bar.is_computed = True

# Compute unknown forces for two unknown bars at a node
def SumOfForcesInLocalY(node, unknown_bars):
    if len(unknown_bars) != 2:
        sys.exit("SumOfForcesInLocalY requires exactly two unknown bars")

    bar1, bar2 = unknown_bars

    # Compute sine and cosine between bars
    sin_theta = geom.SineBars(bar1, bar2)
    cos_theta = geom.CosineBars(bar1, bar2)

    if sin_theta == 0:
        sys.exit("Cannot compute forces: bars are collinear")

    # Sum of known forces at the node
    sum_x = node.xforce_external + (node.xforce_reaction if not np.isnan(node.xforce_reaction) else 0)
    sum_y = node.yforce_external + (node.yforce_reaction if not np.isnan(node.yforce_reaction) else 0)

    for bar in node.bars:
        if bar.is_computed:
            sum_x += bar.axial_load * geom.CosineBars(bar1, bar)
            sum_y += bar.axial_load * geom.SineBars(bar1, bar)

    # Solve for unknown bar forces
    bar2.axial_load = (sum_y - sum_x * cos_theta) / sin_theta
    bar1.axial_load = -(sum_x + bar2.axial_load * cos_theta)
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