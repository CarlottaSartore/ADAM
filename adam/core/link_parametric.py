from urdfpy import matrix_to_xyz_rpy
from enum import Enum
import math

class I_parametric():
    def __init__(self) -> None:
        self.ixx = 0.0
        self.ixy = 0.0
        self.ixz = 0.0
        self.iyy = 0.0
        self.iyz = 0.0
        self.izz = 0.0

class LinkCharacteristics: 
    def __init__(self, offset =0.0, dimension = None, flip_direction = False, calculate_origin_from_dimension = True) -> None:
        
        self.offset = offset
        self.flip_direction = flip_direction
        self.dimension = dimension
        self.calculate_origin_from_dimension = calculate_origin_from_dimension

class JointCharacteristics: 
    def __init__(self, offset = 0.0, take_half_length = False, flip_direction = False, modify_origin = True) -> None:
        self.offset = offset
        self.take_half_length = take_half_length
        self.flip_direction = flip_direction
        self.modify_origin = modify_origin

class Geometry(Enum):
    """The different types of geometries that constitute the URDF"""
    BOX = 1
    CYLINDER = 2
    SPHERE = 3

class Side(Enum):
    """The possible sides of a box geometry"""
    WIDTH = 1
    HEIGHT = 2
    DEPTH = 3

class linkParametric():
    def __init__(self, link_name: str, length_multiplier, density, robot, link, link_characteristics = LinkCharacteristics()) -> None:
        self.name = link_name
        self.density = density
        self.length_multiplier = length_multiplier
        self.link = link
        self.link_characteristic = link_characteristics
        self.geometry_type, self.visual_data = self.get_geometry(self.get_visual())
        (self.volume, self.visual_data_new) = self.compute_volume()
        self.mass = self.compute_mass()
        self.I = self.compute_inertia_parametric()
        self.origin = self.modify_origin()
        self.inertial = self.I
        parent_joint_list = [corresponding_joint for corresponding_joint in robot.joints if corresponding_joint.child == link.name]
        self.parent_joint = (parent_joint_list[0] if parent_joint_list else None)
        child_joint_list = [corresponding_joint for corresponding_joint in robot.joints if corresponding_joint.parent == link.name]
        self.child_joint = (child_joint_list[0] if child_joint_list else None)
        link_offset, joint_offset = self.compute_offsets()
        link_characteristics.offset = link_offset

    def compute_offsets(self): 
        # Taking the principal direction i.e. the length 
        visual = self.get_visual()
        xyz_rpy = matrix_to_xyz_rpy(visual.origin) 
        v_l= 0.0
        if self.geometry_type == Geometry.CYLINDER:
            if(xyz_rpy[3] < 0.0 or xyz_rpy[4] > 0.0):
                v_l = 2*self.visual_data.radius # returning the diameter, since the orientation of the shape is such that the radius is the principal lenght 
            else: 
                v_l=self.visual_data.length # returning the lenght, since the orientation of the shape is such that the radius is the principal lenght 
        j_0 = matrix_to_xyz_rpy(self.child_joint.origin)[2]
        v_o = xyz_rpy[2]
        link_offset = (v_o - v_l/2)* math.copysign(1,v_o)
        joint_offset_temp = v_o + v_l*math.copysign(1,j_0)/2 - j_0
        joint_offset = v_o-joint_offset_temp
        return link_offset, joint_offset
        
        

    def get_visual(self):
        """Returns the visual object of a link"""
        return self.link.visuals[0]

    @staticmethod
    def get_geometry(visual_obj):
        """Returns the geometry type and the corresponding geometry object for a given visual"""
        if visual_obj.geometry.box is not None:
            return [Geometry.BOX, visual_obj.geometry.box]
        if visual_obj.geometry.cylinder is not None:
            return [Geometry.CYLINDER, visual_obj.geometry.cylinder]
        if visual_obj.geometry.sphere is not None:
            return [Geometry.SPHERE, visual_obj.geometry.sphere]

    """Function that starting from a multiplier and link visual characteristics computes the link volume"""
    def compute_volume(self):
        volume = 0.0
        """Modifies a link's volume by a given multiplier, in a manner that is logical with the link's geometry"""
        if self.geometry_type == Geometry.BOX:
            visual_data_new =[0.0, 0.0, 0.0]
            visual_data_new[0] = self.visual_data.size[0] * self.length_multiplier[0]
            visual_data_new[1] = self.visual_data.size[1] * self.length_multiplier[1]
            visual_data_new[2] = self.visual_data.size[2] * self.length_multiplier[2]
            volume = visual_data_new[0] * visual_data_new[1] * visual_data_new[2]
        elif self.geometry_type == Geometry.CYLINDER:
            visual_data_new = [0.0, 0.0]
            visual_data_new[0] = self.visual_data.length * self.length_multiplier[0]
            visual_data_new[1] = self.visual_data.radius * self.length_multiplier[1]
            volume = math.pi * visual_data_new[1] ** 2 * visual_data_new[0]
        elif self.geometry_type == Geometry.SPHERE:
            visual_data_new = 0.0
            visual_data_new = self.visual_data.radius * self.length_multiplier[0]
            volume = 4 * math.pi * visual_data_new ** 3 / 3
        return volume, visual_data_new

    """Function that computes the mass starting from the density, the length multiplier and the link"""
    def compute_mass(self):
        """Changes the mass of a link by preserving a given density."""
        mass = 0.0
        mass = self.volume * self.density
        return mass

    def modify_origin(self):
        origin = [0.0,0.0,0.0,0.0,0.0,0.0]
        visual = self.get_visual()
        """Modifies the position of the origin by a given amount"""
        xyz_rpy = matrix_to_xyz_rpy(visual.origin) 
        length = self.get_principal_length()/2
        if self.geometry_type == Geometry.BOX:
            index_to_change = 2 
            if (self.link_characteristic.dimension == Side.DEPTH):
                index_to_change = 2
            elif(self.link_characteristic.dimension == Side.WIDTH):
                index_to_change = 0
            elif(self.link_characteristic.dimension == Side.HEIGHT):
                index_to_change = 1

            if(self.link_characteristic.calculate_origin_from_dimension):
                temp = self.visual_data_new[index_to_change]
                xyz_rpy[index_to_change] = temp 
                origin[2] = temp 
            origin[2] += self.link_characteristic.offset
        elif self.geometry_type == Geometry.CYLINDER:
            origin[2] =(length+ self.link_characteristic.offset)*math.copysign(1,xyz_rpy[2])
            origin[0] = xyz_rpy[0]
            origin[1] = xyz_rpy[1]
            origin[3] = xyz_rpy[3]
            origin[4] = xyz_rpy[4]
            origin[5] = xyz_rpy[5]
        elif self.geometry_type == Geometry.SPHERE:
            "in case of a sphere the origin of the framjoint_name_list[0]:link_parametric.JointCharacteristics(0.0295),e does not change"
            origin = xyz_rpy
        return origin

    def compute_inertia_parametric(self):
        I = I_parametric
        """Calculates inertia (ixx, iyy and izz) with the formula that corresponds to the geometry
        Formulas retrieved from https://en.wikipedia.org/wiki/List_of_moments_of_inertia"""
        if self.geometry_type == Geometry.BOX:
            I.ixx = self.mass * (self.visual_data_new[1] ** 2+ self.visual_data_new[2] ** 2)/12
            I.iyy = self.mass* (self.visual_data_new[0]**2 + self.visual_data_new[2]**2)/12
            I.izz = self.mass * (self.visual_data_new[0]**2 + self.visual_data_new[1]**2)/12
        elif self.geometry_type == Geometry.CYLINDER:
            i_xy_incomplete = (
                3 *(self.visual_data_new[1] ** 2) + self.visual_data_new[0] ** 2
            ) / 12
            I.ixx = self.mass * i_xy_incomplete
            I.iyy = self.mass * i_xy_incomplete
            I.izz = self.mass * self.visual_data_new[1] ** 2 / 2
            return I
        elif self.geometry_type == Geometry.SPHERE:
            I.ixx = 2 * self.mass * self.visual_data_new** 2 / 5
            I.iyy = I.ixx
            I.izz = I.ixx
        return I

    def get_principal_length(self):   
        visual = self.get_visual()
    #    """Modifies the position of the origin by a given amount"""
        xyz_rpy = matrix_to_xyz_rpy(visual.origin) 
        if self.geometry_type == Geometry.CYLINDER:
            if(xyz_rpy[3] < 0.0 or xyz_rpy[4] > 0.0):
                return 2*self.visual_data_new[1] # returning the diameter, since the orientation of the shape is such that the radius is the principal lenght 
            else: 
                return self.visual_data_new[0] # returning the lenght, since the orientation of the shape is such that the radius is the principal lenght 

        elif self.geometry_type == Geometry.BOX:
            index = 2
            if (self.link_characteristic.dimension == Side.DEPTH):
                index = 2
            elif(self.link_characteristic.dimension == Side.WIDTH):
                index = 0
            elif(self.link_characteristic.dimension == Side.HEIGHT):
                index = 1
            return self.visual_data_new[index]
        
        else:
        
            return 0

class jointParametric:
    def __init__(self, joint_name, parent_link, joint, joint_characteristic) -> None:
        self.jointName = joint_name
        self.parent_link_name = parent_link
        self.joint = joint
        self.parent_link = parent_link
        self.jointCharacteristic  = joint_characteristic
        self.link_offset, joint_offset = self.parent_link.compute_offsets()
        joint_characteristic.offset = joint_offset
        self.origin = self.modify()
        
    def modify(self):
        length = self.parent_link.get_principal_length()
        # Ack for avoiding depending on casadi 
        xyz_rpy = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        xyz_rpy[0] = self.joint.origin[0,3]
        xyz_rpy[1] = self.joint.origin[1,3]
        xyz_rpy[2] = self.joint.origin[2,3]
        xyz_rpy_temp=  matrix_to_xyz_rpy(self.joint.origin)
        xyz_rpy[3] = xyz_rpy_temp[3]
        xyz_rpy[4] = xyz_rpy_temp[4]
        xyz_rpy[5] = xyz_rpy_temp[5]
        xyz_rpy[2] = length/2*math.copysign(1,xyz_rpy[2]) + self.jointCharacteristic.offset
        return xyz_rpy
