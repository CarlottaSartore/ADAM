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
    def __init__(self, link_name: str, length_multiplier, density, link) -> None:
        self.name = link_name
        self.density = density
        self.length_multiplier = length_multiplier
        self.link = link
        self.geometry_type, self.visual_data = self.get_geometry(self.get_visual())
        link_offset = self.compute_offset()
        self.offset = link_offset
        
        (self.volume, self.visual_data_new) = self.compute_volume()
        self.mass = self.compute_mass()
        self.I = self.compute_inertia_parametric()
        self.origin = self.modify_origin()
        self.inertial = self.I
       
    def get_principal_lenght(self): 
        visual = self.get_visual()
        xyz_rpy = [*visual.origin.xyz, *visual.origin.rpy] 
        if self.geometry_type == Geometry.CYLINDER:
            if(xyz_rpy[3] < 0.0 or xyz_rpy[4] > 0.0):
                v_l = 2*self.visual_data.radius # returning the diameter, since the orientation of the shape is such that the radius is the principal lenght 
            else: 
                v_l=self.visual_data.length # returning the lenght, since the orientation of the shape is such that the radius is the principal lenght 
        elif(self.geometry_type == Geometry.SPHERE): 
            v_l = self.visual_data.radius
        elif(self.geometry_type == Geometry.BOX): 
            v_l= self.visual_data.size[2]
        else:
            raise Exception(f"THE GEOMETRY IS NOT SPECIFIED")
        return v_l 

    def get_principal_lenght_parametric(self): 
        visual = self.get_visual()
        xyz_rpy = [*visual.origin.xyz, *visual.origin.rpy] 
        if self.geometry_type == Geometry.CYLINDER:
            if(xyz_rpy[3] < 0.0 or xyz_rpy[4]):
                v_l = 2*self.visual_data_new[1] # returning the diameter, since the orientation of the shape is such that the radius is the principal lenght 
            else: 
                v_l=self.visual_data_new[0] # returning the lenght, since the orientation of the shape is such that the radius is the principal lenght 
        elif(self.geometry_type == Geometry.SPHERE): 
            v_l = self.visual_data_new
        elif(self.geometry_type == Geometry.BOX): 
            v_l= self.visual_data_new[2]
        else:
            raise Exception(f"THE GEOMETRY IS NOT SPECIFIED")
        return v_l 
   
    def compute_offset(self): 
        visual = self.get_visual()
        xyz_rpy = [*visual.origin.xyz, *visual.origin.rpy] 
        v_l=  self.get_principal_lenght()
        v_o = xyz_rpy[2]
        if(v_o<0):
            link_offset = v_l/2 + v_o
        else:
            link_offset = (v_o - v_l/2)
        return link_offset

    def compute_joint_offset(self,joint_i, parent_offset): 
         # Taking the principal direction i.e. the length 
        visual = self.get_visual()
        xyz_rpy = [*visual.origin.xyz, *visual.origin.rpy] 
        v_l= self.get_principal_lenght()
        j_0 = joint_i.origin.xyz[2]
        v_o =  xyz_rpy[2]
        if(j_0<0):
            joint_offset_temp = -(v_l + j_0-parent_offset)
            joint_offset = joint_offset_temp
        else:
            joint_offset_temp = v_l + parent_offset - j_0
            joint_offset = joint_offset_temp
        return joint_offset
    
    def get_visual(self):
        """Returns the visual object of a link"""
        return self.link.visuals[0]

    @staticmethod
    def get_geometry(visual_obj):
        if hasattr(visual_obj.geometry, "size"):
            return [Geometry.BOX, visual_obj.geometry]
        elif hasattr(visual_obj.geometry, "length"):
            return [Geometry.CYLINDER, visual_obj.geometry]
        elif hasattr(visual_obj.geometry,"radius"):
            return [Geometry.SPHERE, visual_obj.geometry]

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
        xyz_rpy = [*visual.origin.xyz, *visual.origin.rpy] 
        v_o = visual.origin.xyz[2] 
        length = self.get_principal_lenght_parametric()
        if(v_o<0):
            origin[2] = self.offset-length/2
            origin[0] = xyz_rpy[0]
            origin[1] = xyz_rpy[1]
            origin[3] = xyz_rpy[3]
            origin[4] = xyz_rpy[4]
            origin[5] = xyz_rpy[5]
        else:
            origin[2] = length/2 +self.offset
            origin[0] = xyz_rpy[0]
            origin[1] = xyz_rpy[1]
            origin[3] = xyz_rpy[3]
            origin[4] = xyz_rpy[4]
            origin[5] = xyz_rpy[5]
        if self.geometry_type == Geometry.SPHERE:
            origin = xyz_rpy
        return origin

    def compute_inertia_parametric(self):
        I = I_parametric()
        
        visual = self.get_visual() 
        xyz_rpy =  [*visual.origin.xyz, *visual.origin.rpy]
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

            if(xyz_rpy[3]>0 and xyz_rpy[4] == 0.0 and xyz_rpy[5] == 0.0):
                itemp = I.izz
                I.iyy = itemp
                I.izz = I.ixx
            elif(xyz_rpy[4]>0.0):
                itemp = I.izz
                I.ixx = itemp
                I.izz = I.iyy
            return I
        elif self.geometry_type == Geometry.SPHERE:
            I.ixx = 2 * self.mass * self.visual_data_new** 2 / 5
            I.iyy = I.ixx
            I.izz = I.ixx
        return I

class jointParametric:
    def __init__(self, joint_name:str, parent_link:linkParametric, joint) -> None:
        self.jointName = joint_name
        self.parent_link_name = parent_link
        self.joint = joint
        self.parent_link = parent_link
        joint_offset = self.parent_link.compute_joint_offset(joint, self.parent_link.offset)
        self.offset = joint_offset
        self.origin = self.modify(self.parent_link.offset)
        
    def modify(self, parent_joint_offset):
        length = self.parent_link.get_principal_lenght_parametric()
        # Ack for avoiding depending on casadi 
        vo = self.parent_link.origin[2]
        xyz_rpy = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        xyz_rpy[0] = self.joint.origin.xyz[0]
        xyz_rpy[1] = self.joint.origin.xyz[1]
        xyz_rpy[2] = self.joint.origin.xyz[2]
        xyz_rpy_temp=  [*self.joint.origin.xyz, *self.joint.origin.rpy]
        xyz_rpy[3] = xyz_rpy_temp[3]
        xyz_rpy[4] = xyz_rpy_temp[4]
        xyz_rpy[5] = xyz_rpy_temp[5]
        if(xyz_rpy[2]<0): 
            xyz_rpy[2] = -length +parent_joint_offset - self.offset   
        else:
            xyz_rpy[2] = vo+ length/2 - self.offset
        return xyz_rpy