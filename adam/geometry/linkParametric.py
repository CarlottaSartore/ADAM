import logging
from dataclasses import dataclass, field
from os import error

import casadi as cs
import numpy as np
from prettytable import PrettyTable
from urdf_parser_py.urdf import URDF
from enums import *
import math
from adam.geometry import utils

class linkParametric:

    def __init__(
        self,
        link_name: str,
        length_multiplier,
        density, 
        robot
    ) -> None:
        self.link_name = link_name
        self.density_multiplier = density
        self.length_multiplier = length_multiplier
        link_list = [link_name for corresponding_link in robot.links if corresponding_link.name == link_name]
        if len(link_list) != 0:
            self.link =link_list[0]
        else:
            print(f"Error modifying link {link_name} it does not exist in the urdf")

        self.compute_volume()
        self.compute_mass()

    """Function that starting from a multiplier (casadi variable) and link visual characteristics computes the link volume"""
    def compute_volume(self):
        self.volume = cs.SX.zeros(1)
        """Modifies a link's volume by a given multiplier, in a manner that is logical with the link's geometry"""
        if (self.link.geometry_type == Geometry.BOX):
            self.visual_data = cs.SX.zeros(3)
            if(self.link.dimension == Side.WIDTH):
                self.visual_data_new[0]= self.link.visual_data.size[0] *self.length_multiplier
            elif (self.link.dimension == Side.HEIGHT):
                self.visual_data_new[1] = self.link.visual_data.size[1]*self.length_multiplier
            elif (self.link.dimension == Side.DEPTH):
                self.visual_data_new[2] = self.link.visual_data.size[2]*self.length_multiplier
            self.volume = self.visual_data_new[0] * self.visual_data_new[1] * self.visual_data_new[2]
        elif (self.link.geometry_type == Geometry.CYLINDER):
            self.visual_data_new = cs.SX.zeros(1)
            self.visual_data_new = self.link.visual_data.length *self.length_multiplier
            self.volume = math.pi * self.link.visual_data.radius ** 2 * self.visual_data_new
        elif (self.link.geometry_type == Geometry.SPHERE):
            self.visual_data_new = cs.SX.zeros(1)
            self.visual_data_new = self.link.visual_data.radius *self.length_multiplier** (1./3)
            self.volume = 4 * math.pi * self.visual_data_new ** 3 / 3 

    """Function that computes the mass starting from the density, the length multiplier and the link"""
    def compute_mass(self):
        """Changes the mass of a link by preserving a given density."""
        self.mass = cs.SX.zeros(1)
        self.mass = self.volume * self.density    

    def modify_origin(self):
        """Modifies the position of the origin by a given amount"""
        visual_obj = self.get_visual()
        geometry_type, visual_data = self.get_geometry(visual_obj)
        xyz_rpy = matrix_to_xyz_rpy(visual_obj.origin)
        if (geometry_type == Geometry.BOX):
            if (self.dimension is not None):
                if (self.dimension == Side.WIDTH):
                    index_to_change = 0
                if (self.dimension == Side.HEIGHT):
                    index_to_change = 1
                if (self.dimension == Side.DEPTH):
                    index_to_change = 2
                if (self.calculate_origin_from_dimensions):
                    xyz_rpy[index_to_change] = (visual_data.size[index_to_change] if not self.flip_direction else -visual_data.size[index_to_change]) / 2
                else:
                    xyz_rpy[index_to_change] = 0
                xyz_rpy[index_to_change] += self.origin_modifier 
            else:
                print(f"Error modifying link {self.element.name}'s origin: Box geometry with no dimension")
        elif (geometry_type == Geometry.CYLINDER):
            xyz_rpy[2] = -visual_data.length / 2 + self.origin_modifier
            visual_obj.origin = xyz_rpy_to_matrix(xyz_rpy)
        elif (geometry_type == Geometry.SPHERE):
            return

    "TODO"        
    def compute_inertia_parametric(self):
        ixx = cs.SX.zeros(1)
        iyy = cs.SX.zeros(1)
        izz = cs.SX.zeros(1)
        I = [ixx, iyy, izz]

        """Calculates inertia (ixx, iyy and izz) with the formula that corresponds to the geometry
        Formulas retrieved from https://en.wikipedia.org/wiki/List_of_moments_of_inertia"""
        geometry_type, visual_data = self.get_geometry(self.get_visual())
        mass = self.get_mass()
        if (geometry_type == Geometry.BOX):
            return mass / 12 * np.array([visual_data.size[1] ** 2 + visual_data.size[2] ** 2, 
                                visual_data.size[0] ** 2 + visual_data.size[2] ** 2,
                                visual_data.size[0] ** 2 + visual_data.size[1] ** 2])
        elif (geometry_type == Geometry.CYLINDER):
            i_xy_incomplete = (3 ** visual_data.radius ** 2 + visual_data.length ** 2) / 12
            return mass * np.array([i_xy_incomplete, i_xy_incomplete, visual_data.radius ** 2 / 2])
        elif (geometry_type == Geometry.SPHERE):
            inertia = 2 * mass * visual_data.radius ** 2 / 5
            return np.array([inertia, inertia, inertia])   
 
