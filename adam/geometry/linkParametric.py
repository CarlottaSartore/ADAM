import logging
from dataclasses import dataclass, field
from os import error

import casadi as cs
import numpy as np
from prettytable import PrettyTable
from urdf_parser_py.urdf import URDF
from urdfpy import xyz_rpy_to_matrix, matrix_to_xyz_rpy
from enums import *
import math
from adam.geometry import utils

"definition of I parametric. Note that the off diagonal element are considered to be zero"

class I_parametric:
    def __init__(self) -> None:
        self.ixx =cs.SX.zeros(1)
        self.ixy = 0 
        self.ixz = 0
        self.iyy= cs.SX.zeros(1)
        self.iyz = 0
        self.izz = cs.SX.zeros(1)

class linkParametric:

    def __init__(
        self,
        link_name: str,
        length_multiplier,
        density, 
        robot
    ) -> None:
        self.link_name = link_name
        self.density = density
        self.length_multiplier = length_multiplier
        link_list = [link_name for corresponding_link in robot.links if corresponding_link.name == link_name]
        if len(link_list) != 0:
            self.link =link_list[0]
        else:
            print(f"Error modifying link {link_name} it does not exist in the urdf")

        (self.volume, self.visual_data_new) = self.compute_volume()
        self.mass = self.compute_mass()
        self.I = self.compute_inertia_parametric()

    """Function that starting from a multiplier (casadi variable) and link visual characteristics computes the link volume"""
    def compute_volume(self):
        volume = cs.SX.zeros(1)
        """Modifies a link's volume by a given multiplier, in a manner that is logical with the link's geometry"""
        if (self.link.geometry_type == Geometry.BOX):
            visual_data_new = cs.SX.zeros(3)
            if(self.link.dimension == Side.WIDTH):
                visual_data_new[0]= self.link.visual_data.size[0] *self.length_multiplier
            elif (self.link.dimension == Side.HEIGHT):
                visual_data_new[1] = self.link.visual_data.size[1]*self.length_multiplier
            elif (self.link.dimension == Side.DEPTH):
                visual_data_new[2] = self.link.visual_data.size[2]*self.length_multiplier
            volume = self.visual_data_new[0] * self.visual_data_new[1] * self.visual_data_new[2]
        elif (self.link.geometry_type == Geometry.CYLINDER):
            visual_data_new = cs.SX.zeros(2)
            visual_data_new[0] = self.link.visual_data.length *self.length_multiplier 
            visual_data_new[1] = self.link.visual_data.radius
            volume = math.pi * self.link.visual_data.radius ** 2 * self.visual_data_new
        elif (self.link.geometry_type == Geometry.SPHERE):
            visual_data_new = cs.SX.zeros(1)
            visual_data_new = self.link.visual_data.radius *self.length_multiplier** (1./3)
            volume = 4 * math.pi * self.visual_data_new ** 3 / 3 

        return volume, visual_data_new

    """Function that computes the mass starting from the density, the length multiplier and the link"""
    def compute_mass(self):
        """Changes the mass of a link by preserving a given density."""
        mass = cs.SX.zeros(1)
        mass = self.volume * self.density    
        return mass 
    
    def modify_origin(self):
        origin = cs.SX.casadi(6)
        visual = self.get_visual()
        """Modifies the position of the origin by a given amount"""
        xyz_rpy = matrix_to_xyz_rpy(visual.origin)
        
        if (self.link.geometry_type == Geometry.BOX):
            "For now hardcoded could be then changed Correspictive line in ergocub gazebo simulator --> check for modifyOrigin"
            "if (self.link.dimension == Side.DEPTH):"
            index_to_change = 2
            "if (self.calculate_origin_from_dimensions):"
            xyz_rpy[index_to_change] = (self.visual_data_new.size[index_to_change]) / 2
            "else:"
            "xyz_rpy[index_to_change] = 0"
            "xyz_rpy[index_to_change] += self.origin_modifier" 
        elif (self.link.geometry_type == Geometry.CYLINDER):
            xyz_rpy[2] = -self.visual_data_new[0] / 2
        elif (self.link.geometry_type == Geometry.SPHERE):
            "in case of a sphere the origin of the frame does not change"
            xyz_rpy = xyz_rpy 
        
        origin = xyz_rpy_to_matrix(xyz_rpy)
        
        return origin

    def compute_inertia_parametric(self):
        I = I_parametric
        """Calculates inertia (ixx, iyy and izz) with the formula that corresponds to the geometry
        Formulas retrieved from https://en.wikipedia.org/wiki/List_of_moments_of_inertia"""
        if (self.link.geometry_type == Geometry.BOX):
            I.ixx = self.mass / 12 * np.array([self.visual_data_new.size[1] ** 2 + self.visual_data_new.size[2] ** 2, 
                                self.visual_data_new.size[0] ** 2 + self.visual_data_new.size[2] ** 2,
                                self.visual_data_new.size[0] ** 2 + self.visual_data_new.size[1] ** 2])
        elif (self.link.geometry_type == Geometry.CYLINDER):
            i_xy_incomplete = (3 ** self.visual_data_new[0] ** 2 + self.visual_data_new[1] ** 2) / 12
            I.ixx =  self.mass * np.array([i_xy_incomplete, i_xy_incomplete, self.visual_data_new[1] ** 2 / 2])
        elif (self.link.geometry_type == Geometry.SPHERE):
            I.ixx= 2 * self.mass * self.visual_data_new[0] ** 2 / 5 
        I.iyy = I.ixx
        I.izz = I.ixx
        return I
 
    def get_visual(self):
        """Returns the visual object of a link"""
        return self.link.element.visuals[0]

    def get_principal_length(self):
        if (self.link.geometry_type == Geometry.CYLINDER):
            return self.visual_data_new[0]
        elif (self.link.geometry_type ==  Geometry.BOX):
            return self.visual_data_new[2]
        else:
            "TODO understand why in case of cylinder it is zero, most likely because things does not change"
            return 0
        
class jointParametric: 
    def __init__(self, joint_name, parent_link_name,parent_link, joint) -> None:
        self.jointName  = joint_name
        self.parent_link_name = parent_link
        self.joint = joint
        self.parent_link = parent_link
        self.origin = self.modify()
    def modify(self):
        length = self.parent_link.get_principal_length()
        xyz_rpy = matrix_to_xyz_rpy(self.joint.origin)
        xyz_rpy[2] = length 
        return xyz_rpy
        