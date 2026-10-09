import pygmsh
import numpy as np
import meshio
import pandas as pd

#mi_r = 0.1; # Mesh resolution

# Add melt inclusion parameters here


# create class objects for different crystal systems...


# Sphere - melt inclusion

# Olivine

# Feldspars

# Define functions for points, faces, surfaces, volumes etc. 

# Add scaling function in here as well? - yes


class sphere(object):
    
    
    def __init__(self, geom, mi_cp, mi_x, mi_y, mi_z, mi_r):
        self.geom = geom
        self.mi_cp = mi_cp # Melt inclusion central point
        self.mi_x = mi_x # Melt inclusion x axis
        self.mi_y = mi_y # Melt inclusion y axis
        self.mi_z = mi_z # Melt inclusion z axis
        self.mi_r = mi_r # Mesh resolution

    def surface(self):
        mi_p0 = self.geom.add_point(self.mi_cp, self.mi_r)
        mi_p1 = self.geom.add_point([self.mi_cp[0] - self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p2 = self.geom.add_point([self.mi_cp[0] + self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p3 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] - self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p4 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] + self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p5 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] - self.mi_z], self.mi_r)
        mi_p6 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] + self.mi_z], self.mi_r)
        
        #points = [mi_p0, mi_p1, mi_p2, mi_p3, mi_p4, mi_p5, mi_p6]
        
        mi_c1 = self.geom.add_ellipse_arc(mi_p2, mi_p0, mi_p0, mi_p4)
        mi_c2 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p4)
        mi_c3 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p3) 
        mi_c4 = self.geom.add_ellipse_arc(mi_p3, mi_p0, mi_p0, mi_p2) 
        mi_c5 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p2) 
        mi_c6 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p1) 
        mi_c7 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p6) 
        mi_c8 = self.geom.add_ellipse_arc(mi_p6, mi_p0, mi_p0, mi_p3) 
        mi_c9 = self.geom.add_ellipse_arc(mi_p4, mi_p0, mi_p0, mi_p5) 
        mi_c10 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p3) 
        mi_c11 = self.geom.add_ellipse_arc(mi_p3, mi_p0, mi_p0, mi_p6)
        mi_c12 = self.geom.add_ellipse_arc(mi_p6, mi_p0, mi_p0, mi_p4)
        mi_c13 = self.geom.add_ellipse_arc(mi_p2, mi_p0, mi_p0, mi_p6)
        
        #arcs = [mi_c1, mi_c2, mi_c3, mi_c4, mi_c5, mi_c6, mi_c7, mi_c8, mi_c9, mi_c10, mi_c11, mi_c12, mi_c13]
        
        mi_ll1 = self.geom.add_line_loop([mi_c9, mi_c5, mi_c1])
        mi_ll2 = self.geom.add_line_loop([mi_c5, -mi_c4, -mi_c10])
        mi_ll3 = self.geom.add_line_loop([mi_c6, mi_c3, -mi_c10])
        mi_ll4 = self.geom.add_line_loop([mi_c2, mi_c9, mi_c6])
        mi_ll5 = self.geom.add_line_loop([mi_c12, -mi_c2, mi_c7])
        mi_ll6 = self.geom.add_line_loop([mi_c7, mi_c8, -mi_c3])
        mi_ll7 = self.geom.add_line_loop([mi_c12, -mi_c1, mi_c13])
        mi_ll8 = self.geom.add_line_loop([mi_c13, mi_c8, mi_c4])
        
        mi_s1 = self.geom.add_surface(mi_ll1)
        mi_s2 = self.geom.add_surface(mi_ll2)
        mi_s3 = self.geom.add_surface(mi_ll3)
        mi_s4 = self.geom.add_surface(mi_ll4)
        mi_s5 = self.geom.add_surface(mi_ll5)
        mi_s6 = self.geom.add_surface(mi_ll6)
        mi_s7 = self.geom.add_surface(mi_ll7)
        mi_s8 = self.geom.add_surface(mi_ll8)
        
        surfaces = [mi_s1, mi_s2, mi_s3, mi_s4, mi_s5, mi_s6, mi_s7, mi_s8]

        return surfaces
        
    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())

 
class olivine(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution
        
        
        ol_x = 1.0 # Distance along x (a)
        ol_y = 0.46 # Distance along y (b)
        ol_z = 0.92 # Half-length of the crystal along z, which is parallel to the c axis
        
        self.sf = self.z/ol_z
        
        
        	# Create olivine points
    def surfaces(self):
        ol_p1 = self.geom.add_point([0.7061357248994403*self.sf, -0.0000000165845153*self.sf, 0.6865771927157608*self.sf], self.m_r)
        ol_p2 = self.geom.add_point([0.4187895557184330*self.sf, 0.4599999885112464*self.sf, 0.4469191173484135*self.sf], self.m_r)
        ol_p3 = self.geom.add_point([0.6266660391977698*self.sf, 0.1704586163692148*self.sf, 0.7865511922113492*self.sf], self.m_r)
        ol_p4 = self.geom.add_point([0.5645677927000951*self.sf, 0.3036563522524822*self.sf, 0.6303102254605011*self.sf], self.m_r)
        ol_p5 = self.geom.add_point([0.5205871277453623*self.sf, 0.0566915522734537*self.sf, 0.9199999857404487*self.sf], self.m_r)
        ol_p6 = self.geom.add_point([-0.7061357248994403*self.sf, -0.0000000161873903*self.sf, 0.6865771927157608*self.sf], self.m_r)
        ol_p7 = self.geom.add_point([-0.4187895559771333*self.sf, 0.4599999887467703*self.sf, 0.4469191173484135*self.sf], self.m_r)
        ol_p8 = self.geom.add_point([-0.6266660392936343*self.sf, 0.1704586167216467*self.sf, 0.7865511922113492*self.sf], self.m_r)
        ol_p9 = self.geom.add_point([-0.5645677928708690*self.sf, 0.3036563525699906*self.sf, 0.6303102254605011*self.sf], self.m_r)
        ol_p10 = self.geom.add_point([-0.5205871277772453*self.sf, 0.0566915525662277*self.sf, 0.9199999857404486*self.sf], self.m_r)
        ol_p11 = self.geom.add_point([0.4187895559771333*self.sf, -0.4599999887467703*self.sf, 0.4469191173484134*self.sf], self.m_r)
        ol_p12 = self.geom.add_point([0.6266660392936343*self.sf, -0.1704586167216468*self.sf, 0.7865511922113491*self.sf], self.m_r)
        ol_p13 = self.geom.add_point([0.5645677928708690*self.sf, -0.3036563525699906*self.sf, 0.6303102254605010*self.sf], self.m_r)
        ol_p14 = self.geom.add_point([0.5205871277772451*self.sf, -0.0566915525662277*self.sf, 0.9199999857404487*self.sf], self.m_r)
        ol_p15 = self.geom.add_point([-0.4187895557184330*self.sf, -0.4599999885112465*self.sf, 0.4469191173484135*self.sf], self.m_r)
        ol_p16 = self.geom.add_point([-0.6266660391977699*self.sf, -0.1704586163692148*self.sf, 0.7865511922113491*self.sf], self.m_r)
        ol_p17 = self.geom.add_point([-0.5645677927000951*self.sf, -0.3036563522524822*self.sf, 0.6303102254605010*self.sf], self.m_r)
        ol_p18 = self.geom.add_point([-0.5205871277453623*self.sf, -0.0566915522734537*self.sf, 0.9199999857404486*self.sf], self.m_r)
        ol_p19 = self.geom.add_point([0.7061357248994403*self.sf, -0.0000000165845154*self.sf, -0.6865771927157607*self.sf], self.m_r)
        ol_p20 = self.geom.add_point([0.4187895557184330*self.sf, 0.4599999885112465*self.sf, -0.4469191173484135*self.sf], self.m_r)
        ol_p21 = self.geom.add_point([0.6266660391977699*self.sf, 0.1704586163692148*self.sf, -0.7865511922113491*self.sf], self.m_r)
        ol_p22 = self.geom.add_point([0.5645677927000951*self.sf, 0.3036563522524822*self.sf, -0.6303102254605010*self.sf], self.m_r)
        ol_p23 = self.geom.add_point([0.5205871277453623*self.sf, 0.0566915522734537*self.sf, -0.9199999857404486*self.sf], self.m_r)
        ol_p24 = self.geom.add_point([-0.7061357248994404*self.sf, -0.0000000161873904*self.sf, -0.6865771927157607*self.sf], self.m_r)
        ol_p25 = self.geom.add_point([-0.4187895559771333*self.sf, 0.4599999887467703*self.sf, -0.4469191173484134*self.sf], self.m_r)
        ol_p26 = self.geom.add_point([-0.6266660392936343*self.sf, 0.1704586167216468*self.sf, -0.7865511922113491*self.sf], self.m_r)
        ol_p27 = self.geom.add_point([-0.5645677928708690*self.sf, 0.3036563525699906*self.sf, -0.6303102254605010*self.sf], self.m_r)
        ol_p28 = self.geom.add_point([-0.5205871277772451*self.sf, 0.0566915525662277*self.sf, -0.9199999857404487*self.sf], self.m_r)
        ol_p29 = self.geom.add_point([0.4187895559771333*self.sf, -0.4599999887467703*self.sf, -0.4469191173484135*self.sf], self.m_r)
        ol_p30 = self.geom.add_point([0.6266660392936343*self.sf, -0.1704586167216467*self.sf, -0.7865511922113492*self.sf], self.m_r)
        ol_p31 = self.geom.add_point([0.5645677928708690*self.sf, -0.3036563525699906*self.sf, -0.6303102254605011*self.sf], self.m_r)
        ol_p32 = self.geom.add_point([0.5205871277772453*self.sf, -0.0566915525662277*self.sf, -0.9199999857404486*self.sf], self.m_r)
        ol_p33 = self.geom.add_point([-0.4187895557184330*self.sf, -0.4599999885112464*self.sf, -0.4469191173484135*self.sf], self.m_r)
        ol_p34 = self.geom.add_point([-0.6266660391977698*self.sf, -0.1704586163692148*self.sf, -0.7865511922113492*self.sf], self.m_r)
        ol_p35 = self.geom.add_point([-0.5645677927000951*self.sf, -0.3036563522524822*self.sf, -0.6303102254605011*self.sf], self.m_r)
        ol_p36 = self.geom.add_point([-0.5205871277453623*self.sf, -0.0566915522734537*self.sf, -0.9199999857404487*self.sf], self.m_r)
        
        #return [ol_p1, ol_p2, ol_p3, ol_p4, ol_p5, ol_p6, ol_p7, ol_p8, ol_p9, ol_p10, ol_p11, ol_p12, ol_p13, ol_p14, ol_p15, ol_p16, ol_p17, ol_p18, ol_p19, ol_p20, ol_p21, ol_p22, ol_p23, ol_p24, ol_p25, ol_p26, ol_p27, ol_p28, ol_p29, ol_p30, ol_p31, ol_p32, ol_p33, ol_p34, ol_p35, ol_p36]
        
        ol_l1 = self.geom.add_line(ol_p1, ol_p3)
        ol_l2 = self.geom.add_line(ol_p1, ol_p12)
        ol_l3 = self.geom.add_line(ol_p1, ol_p19)
        ol_l4 = self.geom.add_line(ol_p2, ol_p4)
        ol_l5 = self.geom.add_line(ol_p2, ol_p7)
        ol_l6 = self.geom.add_line(ol_p2, ol_p20)
        ol_l7 = self.geom.add_line(ol_p3, ol_p4)
        ol_l8 = self.geom.add_line(ol_p3, ol_p5)
        ol_l9 = self.geom.add_line(ol_p4, ol_p22)
        ol_l10 = self.geom.add_line(ol_p5, ol_p10)
        ol_l11 = self.geom.add_line(ol_p5, ol_p14)
        ol_l12 = self.geom.add_line(ol_p6, ol_p8)
        ol_l13 = self.geom.add_line(ol_p6, ol_p16)
        ol_l14 = self.geom.add_line(ol_p6, ol_p24)
        ol_l15 = self.geom.add_line(ol_p7, ol_p9)
        ol_l16 = self.geom.add_line(ol_p7, ol_p25)
        ol_l17 = self.geom.add_line(ol_p8, ol_p9)
        ol_l18 = self.geom.add_line(ol_p8, ol_p10)
        ol_l19 = self.geom.add_line(ol_p9, ol_p27)
        ol_l20 = self.geom.add_line(ol_p10, ol_p18)
        ol_l21 = self.geom.add_line(ol_p11, ol_p13)
        ol_l22 = self.geom.add_line(ol_p11, ol_p15)
        ol_l23 = self.geom.add_line(ol_p11, ol_p29)
        ol_l24 = self.geom.add_line(ol_p12, ol_p13)
        ol_l25 = self.geom.add_line(ol_p12, ol_p14)
        ol_l26 = self.geom.add_line(ol_p13, ol_p31)
        ol_l27 = self.geom.add_line(ol_p14, ol_p18)
        ol_l28 = self.geom.add_line(ol_p15, ol_p17)
        ol_l29 = self.geom.add_line(ol_p15, ol_p33)
        ol_l30 = self.geom.add_line(ol_p16, ol_p17)
        ol_l31 = self.geom.add_line(ol_p16, ol_p18)
        ol_l32 = self.geom.add_line(ol_p17, ol_p35)
        ol_l33 = self.geom.add_line(ol_p19, ol_p21)
        ol_l34 = self.geom.add_line(ol_p19, ol_p30)
        ol_l35 = self.geom.add_line(ol_p20, ol_p22)
        ol_l36 = self.geom.add_line(ol_p20, ol_p25)
        ol_l37 = self.geom.add_line(ol_p21, ol_p22)
        ol_l38 = self.geom.add_line(ol_p21, ol_p23)
        ol_l39 = self.geom.add_line(ol_p23, ol_p28)
        ol_l40 = self.geom.add_line(ol_p23, ol_p32)
        ol_l41 = self.geom.add_line(ol_p24, ol_p26)
        ol_l42 = self.geom.add_line(ol_p24, ol_p34)
        ol_l43 = self.geom.add_line(ol_p25, ol_p27)
        ol_l44 = self.geom.add_line(ol_p26, ol_p27)
        ol_l45 = self.geom.add_line(ol_p26, ol_p28)
        ol_l46 = self.geom.add_line(ol_p28, ol_p36)
        ol_l47 = self.geom.add_line(ol_p29, ol_p31)
        ol_l48 = self.geom.add_line(ol_p29, ol_p33)
        ol_l49 = self.geom.add_line(ol_p30, ol_p31)
        ol_l50 = self.geom.add_line(ol_p30, ol_p32)
        ol_l51 = self.geom.add_line(ol_p32, ol_p36)
        ol_l52 = self.geom.add_line(ol_p33, ol_p35)
        ol_l53 = self.geom.add_line(ol_p34, ol_p35)
        ol_l54 = self.geom.add_line(ol_p34, ol_p36)
        
        ol_ll1 = self.geom.add_line_loop([ol_l45, ol_l46, -ol_l54, -ol_l42, ol_l41])
        ol_ll2 = self.geom.add_line_loop([ol_l19, -ol_l44, -ol_l41, -ol_l14, ol_l12, ol_l17])
        ol_ll3 = self.geom.add_line_loop([ol_l42, ol_l53, -ol_l32, -ol_l30, -ol_l13, ol_l14])
        ol_ll4 = self.geom.add_line_loop([ol_l52, -ol_l32, -ol_l28, ol_l29])
        ol_ll5 = self.geom.add_line_loop([ol_l43, -ol_l19, -ol_l15, ol_l16])
        ol_ll6 = self.geom.add_line_loop([ol_l13, ol_l31, -ol_l20, -ol_l18, -ol_l12])
        ol_ll7 = self.geom.add_line_loop([ol_l39, ol_l46, -ol_l51, -ol_l40])
        ol_ll8 = self.geom.add_line_loop([ol_l39, -ol_l45, ol_l44, -ol_l43, -ol_l36, ol_l35, -ol_l37, ol_l38])
        ol_ll9 = self.geom.add_line_loop([ol_l51, -ol_l54, ol_l53, -ol_l52, -ol_l48, ol_l47, -ol_l49, ol_l50])
        ol_ll10 = self.geom.add_line_loop([ol_l40, -ol_l50, -ol_l34, ol_l33, ol_l38])
        ol_ll11 = self.geom.add_line_loop([ol_l49, -ol_l26, -ol_l24, -ol_l2, ol_l3, ol_l34])
        ol_ll12 = self.geom.add_line_loop([ol_l33, ol_l37, -ol_l9, -ol_l7, -ol_l1, ol_l3])
        ol_ll13 = self.geom.add_line_loop([ol_l26, -ol_l47, -ol_l23, ol_l21])
        ol_ll14 = self.geom.add_line_loop([ol_l35, -ol_l9, -ol_l4, ol_l6])
        ol_ll15 = self.geom.add_line_loop([ol_l1, ol_l8, ol_l11, -ol_l25, -ol_l2])
        ol_ll16 = self.geom.add_line_loop([ol_l11, ol_l27, -ol_l20, -ol_l10])
        ol_ll17 = self.geom.add_line_loop([ol_l30, -ol_l28, -ol_l22, ol_l21, -ol_l24, ol_l25, ol_l27, -ol_l31])
        ol_ll18 = self.geom.add_line_loop([ol_l8, ol_l10, -ol_l18, ol_l17, -ol_l15, -ol_l5, ol_l4, -ol_l7])
        ol_ll19 = self.geom.add_line_loop([ol_l29, -ol_l48, -ol_l23, ol_l22])
        ol_ll20 = self.geom.add_line_loop([ol_l16, -ol_l36, -ol_l6, ol_l5])
        
        ol_s1 = self.geom.add_plane_surface(ol_ll1, holes=None)
        ol_s2 = self.geom.add_plane_surface(ol_ll2, holes=None)
        ol_s3 = self.geom.add_plane_surface(ol_ll3, holes=None)
        ol_s4 = self.geom.add_plane_surface(ol_ll4, holes=None)
        ol_s5 = self.geom.add_plane_surface(ol_ll5, holes=None)
        ol_s6 = self.geom.add_plane_surface(ol_ll6, holes=None)
        ol_s7 = self.geom.add_plane_surface(ol_ll7, holes=None)
        ol_s8 = self.geom.add_plane_surface(ol_ll8, holes=None)
        ol_s9 = self.geom.add_plane_surface(ol_ll9, holes=None)
        ol_s10 = self.geom.add_plane_surface(ol_ll10, holes=None)
        ol_s11 = self.geom.add_plane_surface(ol_ll11, holes=None)
        ol_s12 = self.geom.add_plane_surface(ol_ll12, holes=None)
        ol_s13 = self.geom.add_plane_surface(ol_ll13, holes=None)
        ol_s14 = self.geom.add_plane_surface(ol_ll14, holes=None)
        ol_s15 = self.geom.add_plane_surface(ol_ll15, holes=None)
        ol_s16 = self.geom.add_plane_surface(ol_ll16, holes=None)
        ol_s17 = self.geom.add_plane_surface(ol_ll17, holes=None)
        ol_s18 = self.geom.add_plane_surface(ol_ll18, holes=None)
        ol_s19 = self.geom.add_plane_surface(ol_ll19, holes=None)
        ol_s20 = self.geom.add_plane_surface(ol_ll20, holes=None)
        
        
        surfaces = [ol_s1, ol_s2, ol_s3, ol_s4, ol_s5, ol_s6, ol_s7, ol_s8, ol_s9, ol_s10, ol_s11, ol_s12, ol_s13, ol_s14, ol_s15, ol_s16, ol_s17, ol_s18, ol_s19, ol_s20]
        
        return surfaces
        
    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())


# plagioclase geometries

# Scaling convention for all crystals: the mesh z axis is parallel to the crystal c axis, and
# z is the half-length of the crystal along it (the crystal spans -z..+z along c).
#
# feldspar_656: points checked - identical to the corners in feldspar_656_orientated.txt (c along z).

class feldspar_656(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution
        
        
        #fsp_x = 0.581002 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.450000 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.102766 # Half-length of the crystal along z, which is parallel to the c axis
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5810018256315345*self.sf, -0.0000000217394182*self.sf, 0.8514804479211366*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.5810018256315342*self.sf, -0.0000000217394182*self.sf, -0.5616907996479438*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.5810018067177953*self.sf, 0.0000000217394174*self.sf, -0.8514804571459848*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.5810018067177950*self.sf, 0.0000000217394174*self.sf, 0.5616908257741478*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.3146822408283196*self.sf, 0.4499999866084922*self.sf, 0.9813731882015383*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([0.3146822408283198*self.sf, 0.4499999866084922*self.sf, -0.9295672896694439*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.3146822221676567*self.sf, 0.4500000101575011*self.sf, -0.9813732183109790*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.3146822221676567*self.sf, 0.4500000101575012*self.sf, 0.9295672949110551*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([0.0657897987457684*self.sf, 0.4499999959213320*self.sf, 1.1027661528813291*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.0657897989988447*self.sf, 0.4500000008446622*self.sf, -1.1027661737659222*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-0.3112542369170341*self.sf, 0.4500000100292358*self.sf, 0.9343024898370291*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.3112542366639577*self.sf, 0.4499999867367582*self.sf, -0.9343025107216223*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.3146822410813961*self.sf, -0.4499999867854670*self.sf, 0.9813731882015383*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([0.3146822410813961*self.sf, -0.4499999867854670*self.sf, -0.9295672896694439*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.3146822219145803*self.sf, -0.4500000099805264*self.sf, -0.9813732183109790*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.3146822219145803*self.sf, -0.4500000099805265*self.sf, 0.9295672949110549*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.0657897989988447*self.sf, -0.4499999959583315*self.sf, 1.1027661528813291*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-0.0657897987457684*self.sf, -0.4500000008076626*self.sf, -1.1027661737659222*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.3112542366639578*self.sf, -0.4500000098541890*self.sf, 0.9343024898370291*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.3112542369170341*self.sf, -0.4499999869118051*self.sf, -0.9343025107216223*self.sf], self.m_r)
        
        fsp_l21 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l22 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l23 = self.geom.add_line(fsp_p1, fsp_p13)
        fsp_l24 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l25 = self.geom.add_line(fsp_p2, fsp_p14)
        fsp_l26 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l27 = self.geom.add_line(fsp_p3, fsp_p7)
        fsp_l28 = self.geom.add_line(fsp_p3, fsp_p15)
        fsp_l29 = self.geom.add_line(fsp_p4, fsp_p8)
        fsp_l30 = self.geom.add_line(fsp_p4, fsp_p16)
        fsp_l31 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l32 = self.geom.add_line(fsp_p5, fsp_p9)
        fsp_l33 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l34 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l35 = self.geom.add_line(fsp_p7, fsp_p10)
        fsp_l36 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l37 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l38 = self.geom.add_line(fsp_p9, fsp_p17)
        fsp_l39 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l40 = self.geom.add_line(fsp_p10, fsp_p18)
        fsp_l41 = self.geom.add_line(fsp_p11, fsp_p19)
        fsp_l42 = self.geom.add_line(fsp_p12, fsp_p20)
        fsp_l43 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l44 = self.geom.add_line(fsp_p13, fsp_p17)
        fsp_l45 = self.geom.add_line(fsp_p14, fsp_p20)
        fsp_l46 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l47 = self.geom.add_line(fsp_p15, fsp_p18)
        fsp_l48 = self.geom.add_line(fsp_p16, fsp_p19)
        fsp_l49 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l50 = self.geom.add_line(fsp_p18, fsp_p20)
        
        fsp_ll1 = self.geom.add_line_loop([fsp_l38, fsp_l49, -fsp_l41, -fsp_l37])
        fsp_ll2 = self.geom.add_line_loop([fsp_l44, -fsp_l38, -fsp_l32, -fsp_l22, fsp_l23])
        fsp_ll3 = self.geom.add_line_loop([fsp_l40, fsp_l50, -fsp_l42, -fsp_l39])
        fsp_ll4 = self.geom.add_line_loop([fsp_l21, fsp_l25, -fsp_l43, -fsp_l23])
        fsp_ll5 = self.geom.add_line_loop([fsp_l24, -fsp_l31, -fsp_l22, fsp_l21])
        fsp_ll6 = self.geom.add_line_loop([fsp_l27, fsp_l34, -fsp_l29, -fsp_l26])
        fsp_ll7 = self.geom.add_line_loop([fsp_l28, fsp_l46, -fsp_l30, -fsp_l26])
        fsp_ll8 = self.geom.add_line_loop([fsp_l40, -fsp_l47, -fsp_l28, fsp_l27, fsp_l35])
        fsp_ll9 = self.geom.add_line_loop([fsp_l33, -fsp_l39, -fsp_l35, fsp_l34, fsp_l36, -fsp_l37, -fsp_l32, fsp_l31])
        fsp_ll10 = self.geom.add_line_loop([fsp_l41, -fsp_l48, -fsp_l30, fsp_l29, fsp_l36])
        fsp_ll11 = self.geom.add_line_loop([fsp_l49, -fsp_l48, -fsp_l46, fsp_l47, fsp_l50, -fsp_l45, -fsp_l43, fsp_l44])
        fsp_ll12 = self.geom.add_line_loop([fsp_l25, fsp_l45, -fsp_l42, -fsp_l33, -fsp_l24])

        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)

        
        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12]
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())


# feldspar_657: points checked - identical to the corners in feldspar_657_orientated.txt (c along z).

class feldspar_657(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution
        
        
        #fsp_x = 0.581002 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.450000 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.092127 # Half-length of the crystal along z, which is parallel to the c axis
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5810018256315344*self.sf, -0.0000000217394182*self.sf, 0.8292284088079225*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.5810018256315342*self.sf, -0.0000000217394182*self.sf, -0.5957970189823630*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.5810018067177952*self.sf, 0.0000000217394175*self.sf, -0.8292284180327706*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.5810018067177950*self.sf, 0.0000000217394174*self.sf, 0.5957970451085674*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.3146822408283197*self.sf, 0.4499999866084922*self.sf, 0.9591211490883241*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([0.3146822408283198*self.sf, 0.4499999866084922*self.sf, -0.8375135102955085*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.3821811468061060*self.sf, 0.3359471343711228*self.sf, -0.8704349263823588*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.3146822221676567*self.sf, 0.4500000101575011*self.sf, -0.9591211791977649*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.3146822221676567*self.sf, 0.4500000101575012*self.sf, 0.8375135038230509*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.3821811343486344*self.sf, 0.3359471523814530*self.sf, 0.8704349282599010*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.1439093175986418*self.sf, 0.4499999929983240*self.sf, 1.0424126755984344*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-0.1439093242307899*self.sf, 0.4500000037676702*self.sf, -1.0424126933717457*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.0419789504971035*self.sf, 0.2777683334580627*self.sf, 1.0921274453617813*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.0419789506533182*self.sf, 0.2777683258208139*self.sf, -1.0921274582530671*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.3477496793953546*self.sf, 0.2777683480405877*self.sf, 0.9179962925655165*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([0.3477496792391401*self.sf, 0.2777683112382886*self.sf, -0.9179963054568023*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.3146822410813961*self.sf, -0.4499999867854670*self.sf, 0.9591211490883240*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([0.3146822410813961*self.sf, -0.4499999867854670*self.sf, -0.8375135102955086*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([0.3821811469950399*self.sf, -0.3359471345860585*self.sf, -0.8704349263823586*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-0.3146822219145803*self.sf, -0.4500000099805264*self.sf, -0.9591211791977649*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.3146822219145804*self.sf, -0.4500000099805265*self.sf, 0.8375135038230509*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-0.3821811341597006*self.sf, -0.3359471521665173*self.sf, 0.8704349282599010*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([0.1439093178517182*self.sf, -0.4499999930792575*self.sf, 1.0424126755984344*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-0.1439093239777136*self.sf, -0.4500000036867369*self.sf, -1.0424126933717457*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([0.0419789506533182*self.sf, -0.2777683334816715*self.sf, 1.0921274453617813*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([-0.0419789504971035*self.sf, -0.2777683257972053*self.sf, -1.0921274582530671*self.sf], self.m_r)
        fsp_p27 = self.geom.add_point([-0.3477496792391400*self.sf, -0.2777683478450161*self.sf, 0.9179962925655166*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([0.3477496793953548*self.sf, -0.2777683114338602*self.sf, -0.9179963054568023*self.sf], self.m_r)

        
        fsp_l29 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l30 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l31 = self.geom.add_line(fsp_p1, fsp_p17)
        fsp_l32 = self.geom.add_line(fsp_p2, fsp_p7)
        fsp_l33 = self.geom.add_line(fsp_p2, fsp_p19)
        fsp_l34 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l35 = self.geom.add_line(fsp_p3, fsp_p8)
        fsp_l36 = self.geom.add_line(fsp_p3, fsp_p20)
        fsp_l37 = self.geom.add_line(fsp_p4, fsp_p10)
        fsp_l38 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l39 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l40 = self.geom.add_line(fsp_p5, fsp_p11)
        fsp_l41 = self.geom.add_line(fsp_p6, fsp_p7)
        fsp_l42 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l43 = self.geom.add_line(fsp_p7, fsp_p16)
        fsp_l44 = self.geom.add_line(fsp_p8, fsp_p9)
        fsp_l45 = self.geom.add_line(fsp_p8, fsp_p12)
        fsp_l46 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l47 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l48 = self.geom.add_line(fsp_p10, fsp_p15)
        fsp_l49 = self.geom.add_line(fsp_p11, fsp_p13)
        fsp_l50 = self.geom.add_line(fsp_p12, fsp_p14)
        fsp_l51 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p13, fsp_p25)
        fsp_l53 = self.geom.add_line(fsp_p14, fsp_p16)
        fsp_l54 = self.geom.add_line(fsp_p14, fsp_p26)
        fsp_l55 = self.geom.add_line(fsp_p15, fsp_p27)
        fsp_l56 = self.geom.add_line(fsp_p16, fsp_p28)
        fsp_l57 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l58 = self.geom.add_line(fsp_p17, fsp_p23)
        fsp_l59 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l60 = self.geom.add_line(fsp_p18, fsp_p24)
        fsp_l61 = self.geom.add_line(fsp_p19, fsp_p28)
        fsp_l62 = self.geom.add_line(fsp_p20, fsp_p21)
        fsp_l63 = self.geom.add_line(fsp_p20, fsp_p24)
        fsp_l64 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l65 = self.geom.add_line(fsp_p21, fsp_p23)
        fsp_l66 = self.geom.add_line(fsp_p22, fsp_p27)
        fsp_l67 = self.geom.add_line(fsp_p23, fsp_p25)
        fsp_l68 = self.geom.add_line(fsp_p24, fsp_p26)
        fsp_l69 = self.geom.add_line(fsp_p25, fsp_p27)
        fsp_l70 = self.geom.add_line(fsp_p26, fsp_p28)

        
        fsp_ll1 = self.geom.add_line_loop([fsp_l55, -fsp_l66, -fsp_l38, fsp_l37, fsp_l48])
        fsp_ll2 = self.geom.add_line_loop([fsp_l51, fsp_l55, -fsp_l69, -fsp_l52])
        fsp_ll3 = self.geom.add_line_loop([fsp_l69, -fsp_l66, -fsp_l64, fsp_l65, fsp_l67])
        fsp_ll4 = self.geom.add_line_loop([fsp_l52, -fsp_l67, -fsp_l58, -fsp_l31, fsp_l30, fsp_l40, fsp_l49])
        fsp_ll5 = self.geom.add_line_loop([fsp_l51, -fsp_l48, -fsp_l46, fsp_l47, fsp_l49])
        fsp_ll6 = self.geom.add_line_loop([fsp_l37, -fsp_l46, -fsp_l44, -fsp_l35, fsp_l34])
        fsp_ll7 = self.geom.add_line_loop([fsp_l47, -fsp_l40, fsp_l39, fsp_l42, -fsp_l45, fsp_l44])
        fsp_ll8 = self.geom.add_line_loop([fsp_l30, fsp_l39, fsp_l41, -fsp_l32, -fsp_l29])
        fsp_ll9 = self.geom.add_line_loop([fsp_l31, fsp_l57, fsp_l59, -fsp_l33, -fsp_l29])
        fsp_ll10 = self.geom.add_line_loop([fsp_l58, -fsp_l65, -fsp_l62, fsp_l63, -fsp_l60, -fsp_l57])
        fsp_ll11 = self.geom.add_line_loop([fsp_l59, fsp_l61, -fsp_l70, -fsp_l68, -fsp_l60])
        fsp_ll12 = self.geom.add_line_loop([fsp_l33, fsp_l61, -fsp_l56, -fsp_l43, -fsp_l32])
        fsp_ll13 = self.geom.add_line_loop([fsp_l56, -fsp_l70, -fsp_l54, fsp_l53])
        fsp_ll14 = self.geom.add_line_loop([fsp_l42, fsp_l50, fsp_l53, -fsp_l43, -fsp_l41])
        fsp_ll15 = self.geom.add_line_loop([fsp_l54, -fsp_l68, -fsp_l63, -fsp_l36, fsp_l35, fsp_l45, fsp_l50])
        fsp_ll16 = self.geom.add_line_loop([fsp_l38, -fsp_l64, -fsp_l62, -fsp_l36, fsp_l34])

        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)
        fsp_s13 = self.geom.add_plane_surface(fsp_ll13, holes=None)
        fsp_s14 = self.geom.add_plane_surface(fsp_ll14, holes=None)
        fsp_s15 = self.geom.add_plane_surface(fsp_ll15, holes=None)
        fsp_s16 = self.geom.add_plane_surface(fsp_ll16, holes=None)

        
        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16] 
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())

# anorthite_743: points checked - identical to the corners in An_743_orientated.txt (c along z).

class anorthite_743(object):
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution
        
        
        #fsp_x = 0.560000 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.511164 # Half-width of the crystal along y (for reference only)
        fsp_z = 0.704836 # Half-length of the crystal along z, which is parallel to the c axis
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5599999975490981*self.sf, 0.0627419733626839*self.sf, 0.0826665360731313*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.5599999975490981*self.sf, 0.0627419733626842*self.sf, -0.1202585312403163*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([0.5599999975763512*self.sf, -0.0341760025805563*self.sf, 0.0880002573664282*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.5599999975763511*self.sf, -0.0341760025805564*self.sf, -0.1149248099470194*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([-0.5599999975490981*self.sf, -0.0627419733626839*self.sf, -0.0826665360731313*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.5599999975490981*self.sf,-0.0627419733626842*self.sf, 0.1202585312403163*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.5599999975763512*self.sf, 0.0341760025805563*self.sf, -0.0880002573664282*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.5599999975763511*self.sf, 0.0341760025805564*self.sf, 0.1149248099470194*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([0.5446587691904488*self.sf, 0.0903567025366361*self.sf, 0.1477651615710680*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([0.3749858924993268*self.sf, 0.3957736110431138*self.sf, 0.2139169452348454*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.3749858924993267*self.sf, 0.3957736110431135*self.sf, -0.1972294602197983*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.4529022750631763*self.sf, 0.2555214909120938*self.sf, -0.4912044431366064*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.4433904909632132*self.sf, 0.2726430235391399*self.sf, -0.4874960043082677*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.5446587691904488*self.sf, -0.0903567025366361*self.sf, -0.1477651615710680*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.3749858924993268*self.sf, -0.3957736110431138*self.sf, -0.2139169452348454*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.3749858924993267*self.sf, -0.3957736110431135*self.sf, 0.1972294602197983*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([-0.4529022750631763*self.sf, -0.2555214909120938*self.sf, 0.4912044431366064*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-0.4433904909632132*self.sf, -0.2726430235391399*self.sf, 0.4874960043082677*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.5313458377570518*self.sf, 0.0828832209160659*self.sf, -0.2151097357660006*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-0.3109662933435937*self.sf, 0.4574911191799947*self.sf, -0.3434780968876873*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.3109662933435938*self.sf, 0.4574911191799951*self.sf, 0.1718006994967579*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-0.4854326880684028*self.sf, 0.1609277987525359*self.sf, 0.3588353476639533*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([-0.3331305972616552*self.sf, 0.4198155566438884*self.sf, 0.2701212619519472*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([0.5313458377570518*self.sf, -0.0828832209160659*self.sf, 0.2151097357660006*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([0.3109662933435937*self.sf, -0.4574911191799947*self.sf, 0.3434780968876873*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([0.3109662933435938*self.sf, -0.4574911191799951*self.sf, -0.1718006994967579*self.sf], self.m_r)
        fsp_p27 = self.geom.add_point([0.4854326880684028*self.sf, -0.1609277987525359*self.sf, -0.3588353476639533*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([0.3331305972616552*self.sf, -0.4198155566438884*self.sf, -0.2701212619519472*self.sf], self.m_r)
        fsp_p29 = self.geom.add_point([0.2274382018244973*self.sf, 0.4892321344583517*self.sf, 0.1849240451820411*self.sf], self.m_r)
        fsp_p30 = self.geom.add_point([0.2274382018244972*self.sf, 0.4892321344583518*self.sf, -0.3197557633349201*self.sf], self.m_r)
        fsp_p31 = self.geom.add_point([-0.2103065149353340*self.sf, 0.5111638347351615*self.sf, -0.2860129884019268*self.sf], self.m_r)
        fsp_p32 = self.geom.add_point([-0.2103065149353340*self.sf, 0.5111638347351616*self.sf, 0.2489277529940093*self.sf], self.m_r)
        fsp_p33 = self.geom.add_point([0.1809658262454565*self.sf, 0.4915604737239746*self.sf, -0.4762429502633959*self.sf], self.m_r)
        fsp_p34 = self.geom.add_point([-0.1716869844313617*self.sf, 0.5092289353853449*self.sf, 0.3789719268449100*self.sf], self.m_r)
        fsp_p35 = self.geom.add_point([-0.2274382018244973*self.sf, -0.4892321344583517*self.sf, -0.1849240451820411*self.sf], self.m_r)
        fsp_p36 = self.geom.add_point([-0.2274382018244972*self.sf, -0.4892321344583518*self.sf, 0.3197557633349201*self.sf], self.m_r)
        fsp_p37 = self.geom.add_point([0.2103065149353340*self.sf, -0.5111638347351615*self.sf, 0.2860129884019268*self.sf], self.m_r)
        fsp_p38 = self.geom.add_point([0.2103065149353340*self.sf, -0.5111638347351616*self.sf, -0.2489277529940093*self.sf], self.m_r)
        fsp_p39 = self.geom.add_point([-0.1809658262454565*self.sf, -0.4915604737239746*self.sf, 0.4762429502633959*self.sf], self.m_r)
        fsp_p40 = self.geom.add_point([0.1716869844313617*self.sf, -0.5092289353853449*self.sf, -0.3789719268449100*self.sf], self.m_r)
        fsp_p41 = self.geom.add_point([0.5145789627524026*self.sf, 0.0392261264210874*self.sf, 0.2811987897357633*self.sf], self.m_r)
        fsp_p42 = self.geom.add_point([0.5145789627782524*self.sf, -0.0527022796766373*self.sf, 0.2862579182647106*self.sf], self.m_r)
        fsp_p43 = self.geom.add_point([0.3048903786275969*self.sf, 0.4166726597823603*self.sf, 0.3629518319679624*self.sf], self.m_r)
        fsp_p44 = self.geom.add_point([0.3168995630060078*self.sf, -0.3887237769342264*self.sf, 0.4014037256662981*self.sf], self.m_r)
        fsp_p45 = self.geom.add_point([-0.1066626426732671*self.sf, 0.4372921148336893*self.sf, 0.5630419157949236*self.sf], self.m_r)
        fsp_p46 = self.geom.add_point([-0.2538873800542000*self.sf, -0.3601264514204584*self.sf, 0.6789106357680711*self.sf], self.m_r)
        fsp_p47 = self.geom.add_point([-0.3203844642672524*self.sf, 0.0740012176973576*self.sf, 0.6875322356518778*self.sf], self.m_r)
        fsp_p48 = self.geom.add_point([-0.3203844641788359*self.sf, -0.2404294541507459*self.sf, 0.7048364101703807*self.sf], self.m_r)
        fsp_p49 = self.geom.add_point([-0.5145789627524026*self.sf, -0.0392261264210874*self.sf, -0.2811987897357633*self.sf], self.m_r)
        fsp_p50 = self.geom.add_point([-0.5145789627782524*self.sf, 0.0527022796766373*self.sf, -0.2862579182647106*self.sf], self.m_r)
        fsp_p51 = self.geom.add_point([-0.3048903786275969*self.sf, -0.4166726597823603*self.sf, -0.3629518319679624*self.sf], self.m_r)
        fsp_p52 = self.geom.add_point([-0.3168995630060078*self.sf, 0.3887237769342264*self.sf, -0.4014037256662981*self.sf], self.m_r)
        fsp_p53 = self.geom.add_point([0.1066626426732671*self.sf, -0.4372921148336893*self.sf, -0.5630419157949236*self.sf], self.m_r)
        fsp_p54 = self.geom.add_point([0.2538873800542000*self.sf, 0.3601264514204584*self.sf, -0.6789106357680711*self.sf], self.m_r)
        fsp_p55 = self.geom.add_point([0.3203844642672524*self.sf, -0.0740012176973576*self.sf, -0.6875322356518778*self.sf], self.m_r)
        fsp_p56 = self.geom.add_point([0.3203844641788359*self.sf, 0.2404294541507459*self.sf, -0.7048364101703807*self.sf], self.m_r)
        fsp_p57 = self.geom.add_point([0.3149623313684918*self.sf, 0.4337932734841546*self.sf, 0.3182727812669887*self.sf], self.m_r)
        fsp_p58 = self.geom.add_point([0.2510532692061323*self.sf, 0.4742740601045964*self.sf, 0.2652015262273856*self.sf], self.m_r)
        fsp_p59 = self.geom.add_point([-0.3149623313684918*self.sf, -0.4337932734841546*self.sf, -0.3182727812669887*self.sf], self.m_r)
        fsp_p60 = self.geom.add_point([-0.2510532692061323*self.sf, -0.4742740601045964*self.sf, -0.2652015262273856*self.sf], self.m_r)
        fsp_p61 = self.geom.add_point([-0.3033127516314823*self.sf, 0.4615720577018110*self.sf, -0.3572757818103139*self.sf], self.m_r)
        fsp_p62 = self.geom.add_point([0.3033127516314823*self.sf, -0.4615720577018110*self.sf, 0.3572757818103139*self.sf], self.m_r)
        fsp_p63 = self.geom.add_point([0.2946575583867223*self.sf, -0.4459924175442804*self.sf, 0.3940030353574674*self.sf], self.m_r)
        fsp_p64 = self.geom.add_point([-0.2946575583867223*self.sf, 0.4459924175442804*self.sf, -0.3940030353574674*self.sf], self.m_r)
        fsp_p65 = self.geom.add_point([-0.1380376161518490*self.sf, 0.4937681267406804*self.sf, 0.4543709013458120*self.sf], self.m_r)
        fsp_p66 = self.geom.add_point([0.1380376161518490*self.sf, -0.4937681267406804*self.sf, -0.4543709013458120*self.sf], self.m_r)
        fsp_p67 = self.geom.add_point([-0.2206230984092732*self.sf, -0.4201760444605870*self.sf, 0.6445237325129128*self.sf], self.m_r)
        fsp_p68 = self.geom.add_point([0.2206230984092732*self.sf, 0.4201760444605870*self.sf, -0.6445237325129128*self.sf], self.m_r)
        fsp_p69 = self.geom.add_point([-0.1702597150921688*self.sf, 0.5066598023707853*self.sf, 0.3949154095337355*self.sf], self.m_r)
        fsp_p70 = self.geom.add_point([0.1702597150921688*self.sf, -0.5066598023707853*self.sf, -0.3949154095337355*self.sf], self.m_r)
        fsp_p71 = self.geom.add_point([-0.2379438194026194*self.sf, -0.4027754741170599*self.sf, 0.6581026809596990*self.sf], self.m_r)
        fsp_p72 = self.geom.add_point([0.2379438194026194*self.sf, 0.4027754741170599*self.sf, -0.6581026809596990*self.sf], self.m_r)
        fsp_p73 = self.geom.add_point([-0.3815174507653690*self.sf, 0.0414045495222755*self.sf, 0.6406912114151304*self.sf], self.m_r)
        fsp_p74 = self.geom.add_point([-0.4201061958642431*self.sf, 0.0433379064747741*self.sf, 0.5851015525764195*self.sf], self.m_r)
        fsp_p75 = self.geom.add_point([0.3815174507653690*self.sf, -0.0414045495222755*self.sf, -0.6406912114151304*self.sf], self.m_r)
        fsp_p76 = self.geom.add_point([0.4201061958642431*self.sf, -0.0433379064747741*self.sf, -0.5851015525764195*self.sf], self.m_r)
        fsp_p77 = self.geom.add_point([-0.3815174506970069*self.sf, -0.2017070680136205*self.sf, 0.6540704591602384*self.sf], self.m_r)
        fsp_p78 = self.geom.add_point([-0.4201061957958807*self.sf, -0.1997737110611223*self.sf, 0.5984808003215281*self.sf], self.m_r)
        fsp_p79 = self.geom.add_point([0.3815174506970069*self.sf, 0.2017070680136205*self.sf, -0.6540704591602384*self.sf], self.m_r)
        fsp_p80 = self.geom.add_point([0.4201061957958807*self.sf, 0.1997737110611223*self.sf, -0.5984808003215281*self.sf], self.m_r)


        
        fsp_l81 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l82 = self.geom.add_line(fsp_p1, fsp_p3)
        fsp_l83 = self.geom.add_line(fsp_p1, fsp_p9)
        fsp_l84 = self.geom.add_line(fsp_p2, fsp_p4)
        fsp_l85 = self.geom.add_line(fsp_p2, fsp_p12)
        fsp_l86 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l87 = self.geom.add_line(fsp_p3, fsp_p24)
        fsp_l88 = self.geom.add_line(fsp_p4, fsp_p27)
        fsp_l89 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l90 = self.geom.add_line(fsp_p5, fsp_p7)
        fsp_l91 = self.geom.add_line(fsp_p5, fsp_p14)
        fsp_l92 = self.geom.add_line(fsp_p6, fsp_p8)
        fsp_l93 = self.geom.add_line(fsp_p6, fsp_p17)
        fsp_l94 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l95 = self.geom.add_line(fsp_p7, fsp_p19)
        fsp_l96 = self.geom.add_line(fsp_p8, fsp_p22)
        fsp_l97 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l98 = self.geom.add_line(fsp_p9, fsp_p41)
        fsp_l99 = self.geom.add_line(fsp_p10, fsp_p11)
        fsp_l100 = self.geom.add_line(fsp_p10, fsp_p57)
        fsp_l101 = self.geom.add_line(fsp_p11, fsp_p13)
        fsp_l102 = self.geom.add_line(fsp_p11, fsp_p30)
        fsp_l103 = self.geom.add_line(fsp_p12, fsp_p13)
        fsp_l104 = self.geom.add_line(fsp_p12, fsp_p80)
        fsp_l105 = self.geom.add_line(fsp_p13, fsp_p72)
        fsp_l106 = self.geom.add_line(fsp_p14, fsp_p15)
        fsp_l107 = self.geom.add_line(fsp_p14, fsp_p49)
        fsp_l108 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l109 = self.geom.add_line(fsp_p15, fsp_p59)
        fsp_l110 = self.geom.add_line(fsp_p16, fsp_p18)
        fsp_l111 = self.geom.add_line(fsp_p16, fsp_p36)
        fsp_l112 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l113 = self.geom.add_line(fsp_p17, fsp_p78)
        fsp_l114 = self.geom.add_line(fsp_p18, fsp_p71)
        fsp_l115 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l116 = self.geom.add_line(fsp_p19, fsp_p50)
        fsp_l117 = self.geom.add_line(fsp_p20, fsp_p21)
        fsp_l118 = self.geom.add_line(fsp_p20, fsp_p61)
        fsp_l119 = self.geom.add_line(fsp_p21, fsp_p23)
        fsp_l120 = self.geom.add_line(fsp_p21, fsp_p32)
        fsp_l121 = self.geom.add_line(fsp_p22, fsp_p23)
        fsp_l122 = self.geom.add_line(fsp_p22, fsp_p74)
        fsp_l123 = self.geom.add_line(fsp_p23, fsp_p69)
        fsp_l124 = self.geom.add_line(fsp_p24, fsp_p25)
        fsp_l125 = self.geom.add_line(fsp_p24, fsp_p42)
        fsp_l126 = self.geom.add_line(fsp_p25, fsp_p26)
        fsp_l127 = self.geom.add_line(fsp_p25, fsp_p62)
        fsp_l128 = self.geom.add_line(fsp_p26, fsp_p28)
        fsp_l129 = self.geom.add_line(fsp_p26, fsp_p38)
        fsp_l130 = self.geom.add_line(fsp_p27, fsp_p28)
        fsp_l131 = self.geom.add_line(fsp_p27, fsp_p76)
        fsp_l132 = self.geom.add_line(fsp_p28, fsp_p70)
        fsp_l133 = self.geom.add_line(fsp_p29, fsp_p30)
        fsp_l134 = self.geom.add_line(fsp_p29, fsp_p34)
        fsp_l135 = self.geom.add_line(fsp_p29, fsp_p58)
        fsp_l136 = self.geom.add_line(fsp_p30, fsp_p33)
        fsp_l137 = self.geom.add_line(fsp_p31, fsp_p32)
        fsp_l138 = self.geom.add_line(fsp_p31, fsp_p33)
        fsp_l139 = self.geom.add_line(fsp_p31, fsp_p61)
        fsp_l140 = self.geom.add_line(fsp_p32, fsp_p34)
        fsp_l141 = self.geom.add_line(fsp_p33, fsp_p68)
        fsp_l142 = self.geom.add_line(fsp_p34, fsp_p69)
        fsp_l143 = self.geom.add_line(fsp_p35, fsp_p36)
        fsp_l144 = self.geom.add_line(fsp_p35, fsp_p40)
        fsp_l145 = self.geom.add_line(fsp_p35, fsp_p60)
        fsp_l146 = self.geom.add_line(fsp_p36, fsp_p39)
        fsp_l147 = self.geom.add_line(fsp_p37, fsp_p38)
        fsp_l148 = self.geom.add_line(fsp_p37, fsp_p39)
        fsp_l149 = self.geom.add_line(fsp_p37, fsp_p62)
        fsp_l150 = self.geom.add_line(fsp_p38, fsp_p40)
        fsp_l151 = self.geom.add_line(fsp_p39, fsp_p67)
        fsp_l152 = self.geom.add_line(fsp_p40, fsp_p70)
        fsp_l153 = self.geom.add_line(fsp_p41, fsp_p42)
        fsp_l154 = self.geom.add_line(fsp_p41, fsp_p43)
        fsp_l155 = self.geom.add_line(fsp_p42, fsp_p44)
        fsp_l156 = self.geom.add_line(fsp_p43, fsp_p45)
        fsp_l157 = self.geom.add_line(fsp_p43, fsp_p57)
        fsp_l158 = self.geom.add_line(fsp_p44, fsp_p46)
        fsp_l159 = self.geom.add_line(fsp_p44, fsp_p63)
        fsp_l160 = self.geom.add_line(fsp_p45, fsp_p47)
        fsp_l161 = self.geom.add_line(fsp_p45, fsp_p65)
        fsp_l162 = self.geom.add_line(fsp_p46, fsp_p48)
        fsp_l163 = self.geom.add_line(fsp_p46, fsp_p71)
        fsp_l164 = self.geom.add_line(fsp_p47, fsp_p48)
        fsp_l165 = self.geom.add_line(fsp_p47, fsp_p73)
        fsp_l166 = self.geom.add_line(fsp_p48, fsp_p77)
        fsp_l167 = self.geom.add_line(fsp_p49, fsp_p50)
        fsp_l168 = self.geom.add_line(fsp_p49, fsp_p51)
        fsp_l169 = self.geom.add_line(fsp_p50, fsp_p52)
        fsp_l170 = self.geom.add_line(fsp_p51, fsp_p53)
        fsp_l171 = self.geom.add_line(fsp_p51, fsp_p59)
        fsp_l172 = self.geom.add_line(fsp_p52, fsp_p54)
        fsp_l173 = self.geom.add_line(fsp_p52, fsp_p64)
        fsp_l174 = self.geom.add_line(fsp_p53, fsp_p55)
        fsp_l175 = self.geom.add_line(fsp_p53, fsp_p66)
        fsp_l176 = self.geom.add_line(fsp_p54, fsp_p56)
        fsp_l177 = self.geom.add_line(fsp_p54, fsp_p72)
        fsp_l178 = self.geom.add_line(fsp_p55, fsp_p56)
        fsp_l179 = self.geom.add_line(fsp_p55, fsp_p75)
        fsp_l180 = self.geom.add_line(fsp_p56, fsp_p79)
        fsp_l181 = self.geom.add_line(fsp_p57, fsp_p58)
        fsp_l182 = self.geom.add_line(fsp_p58, fsp_p65)
        fsp_l183 = self.geom.add_line(fsp_p59, fsp_p60)
        fsp_l184 = self.geom.add_line(fsp_p60, fsp_p66)
        fsp_l185 = self.geom.add_line(fsp_p61, fsp_p64)
        fsp_l186 = self.geom.add_line(fsp_p62, fsp_p63)
        fsp_l187 = self.geom.add_line(fsp_p63, fsp_p67)
        fsp_l188 = self.geom.add_line(fsp_p64, fsp_p68)
        fsp_l189 = self.geom.add_line(fsp_p65, fsp_p69)
        fsp_l190 = self.geom.add_line(fsp_p66, fsp_p70)
        fsp_l191 = self.geom.add_line(fsp_p67, fsp_p71)
        fsp_l192 = self.geom.add_line(fsp_p68, fsp_p72)
        fsp_l193 = self.geom.add_line(fsp_p73, fsp_p74)
        fsp_l194 = self.geom.add_line(fsp_p73, fsp_p77)
        fsp_l195 = self.geom.add_line(fsp_p74, fsp_p78)
        fsp_l196 = self.geom.add_line(fsp_p75, fsp_p76)
        fsp_l197 = self.geom.add_line(fsp_p75, fsp_p79)
        fsp_l198 = self.geom.add_line(fsp_p76, fsp_p80)
        fsp_l199 = self.geom.add_line(fsp_p77, fsp_p78)
        fsp_l200 = self.geom.add_line(fsp_p79, fsp_p80)
        
        fsp_ll1 = self.geom.add_line_loop([fsp_l168, fsp_l171, -fsp_l109, -fsp_l106, fsp_l107])
        fsp_ll2 = self.geom.add_line_loop([fsp_l107, fsp_l167, -fsp_l116, -fsp_l95, -fsp_l90, fsp_l91])
        fsp_ll3 = self.geom.add_line_loop([fsp_l89, fsp_l92, -fsp_l94, -fsp_l90])
        fsp_ll4 = self.geom.add_line_loop([fsp_l93, fsp_l113, -fsp_l195, -fsp_l122, -fsp_l96, -fsp_l92])
        fsp_ll5 = self.geom.add_line_loop([fsp_l93, fsp_l112, -fsp_l110, -fsp_l108, -fsp_l106, -fsp_l91, fsp_l89])
        fsp_ll6 = self.geom.add_line_loop([fsp_l170, fsp_l174, fsp_l178, -fsp_l176, -fsp_l172, -fsp_l169, -fsp_l167, fsp_l168])
        fsp_ll7 = self.geom.add_line_loop([fsp_l170, fsp_l175, -fsp_l184, -fsp_l183, -fsp_l171])
        fsp_ll8 = self.geom.add_line_loop([fsp_l145, fsp_l184, fsp_l190, -fsp_l152, -fsp_l144])
        fsp_ll9 = self.geom.add_line_loop([fsp_l143, -fsp_l111, -fsp_l108, fsp_l109, fsp_l183, -fsp_l145])
        fsp_ll10 = self.geom.add_line_loop([fsp_l110, fsp_l114, -fsp_l191, -fsp_l151, -fsp_l146, -fsp_l111])
        fsp_ll11 = self.geom.add_line_loop([fsp_l95, fsp_l115, fsp_l117, fsp_l119, -fsp_l121, -fsp_l96, -fsp_l94])
        fsp_ll12 = self.geom.add_line_loop([fsp_l115, fsp_l118, fsp_l185, -fsp_l173, -fsp_l169, -fsp_l116])
        fsp_ll13 = self.geom.add_line_loop([fsp_l172, fsp_l177, -fsp_l192, -fsp_l188, -fsp_l173])
        fsp_ll14 = self.geom.add_line_loop([fsp_l185, fsp_l188, -fsp_l141, -fsp_l138, fsp_l139])
        fsp_ll15 = self.geom.add_line_loop([fsp_l141, fsp_l192, -fsp_l105, -fsp_l101, fsp_l102, fsp_l136])
        fsp_ll16 = self.geom.add_line_loop([fsp_l138, -fsp_l136, -fsp_l133, fsp_l134, -fsp_l140, -fsp_l137])
        fsp_ll17 = self.geom.add_line_loop([fsp_l139, -fsp_l118, fsp_l117, fsp_l120, -fsp_l137])
        fsp_ll18 = self.geom.add_line_loop([fsp_l119, fsp_l123, -fsp_l142, -fsp_l140, -fsp_l120])
        fsp_ll19 = self.geom.add_line_loop([fsp_l121, fsp_l123, -fsp_l189, -fsp_l161, fsp_l160, fsp_l165, fsp_l193, -fsp_l122])
        fsp_ll20 = self.geom.add_line_loop([fsp_l189, -fsp_l142, -fsp_l134, fsp_l135, fsp_l182])
        fsp_ll21 = self.geom.add_line_loop([fsp_l161, -fsp_l182, -fsp_l181, -fsp_l157, fsp_l156])
        fsp_ll22 = self.geom.add_line_loop([fsp_l99, fsp_l102, -fsp_l133, fsp_l135, -fsp_l181, -fsp_l100])
        fsp_ll23 = self.geom.add_line_loop([fsp_l101, -fsp_l103, -fsp_l85, -fsp_l81, fsp_l83, fsp_l97, fsp_l99])
        fsp_ll24 = self.geom.add_line_loop([fsp_l97, fsp_l100, -fsp_l157, -fsp_l154, -fsp_l98])
        fsp_ll25 = self.geom.add_line_loop([fsp_l82, fsp_l87, fsp_l125, -fsp_l153, -fsp_l98, -fsp_l83])
        fsp_ll26 = self.geom.add_line_loop([fsp_l154, fsp_l156, fsp_l160, fsp_l164, -fsp_l162, -fsp_l158, -fsp_l155, -fsp_l153])
        fsp_ll27 = self.geom.add_line_loop([fsp_l124, fsp_l127, fsp_l186, -fsp_l159, -fsp_l155, -fsp_l125])
        fsp_ll28 = self.geom.add_line_loop([fsp_l81, fsp_l84, -fsp_l86, -fsp_l82])
        fsp_ll29 = self.geom.add_line_loop([fsp_l85, fsp_l104, -fsp_l198, -fsp_l131, -fsp_l88, -fsp_l84])
        fsp_ll30 = self.geom.add_line_loop([fsp_l105, -fsp_l177, fsp_l176, fsp_l180, fsp_l200, -fsp_l104, fsp_l103])
        fsp_ll31 = self.geom.add_line_loop([fsp_l158, fsp_l163, -fsp_l191, -fsp_l187, -fsp_l159])
        fsp_ll32 = self.geom.add_line_loop([fsp_l151, -fsp_l187, -fsp_l186, -fsp_l149, fsp_l148])
        fsp_ll33 = self.geom.add_line_loop([fsp_l149, -fsp_l127, fsp_l126, fsp_l129, -fsp_l147])
        fsp_ll34 = self.geom.add_line_loop([fsp_l144, -fsp_l150, -fsp_l147, fsp_l148, -fsp_l146, -fsp_l143])
        fsp_ll35 = self.geom.add_line_loop([fsp_l124, fsp_l126, fsp_l128, -fsp_l130, -fsp_l88, -fsp_l86, fsp_l87])
        fsp_ll36 = self.geom.add_line_loop([fsp_l130, fsp_l132, -fsp_l190, -fsp_l175, fsp_l174, fsp_l179, fsp_l196, -fsp_l131])
        fsp_ll37 = self.geom.add_line_loop([fsp_l128, fsp_l132, -fsp_l152, -fsp_l150, -fsp_l129])
        fsp_ll38 = self.geom.add_line_loop([fsp_l198, -fsp_l200, -fsp_l197, fsp_l196])
        fsp_ll39 = self.geom.add_line_loop([fsp_l162, fsp_l166, fsp_l199, -fsp_l113, fsp_l112, fsp_l114, -fsp_l163])
        fsp_ll40 = self.geom.add_line_loop([fsp_l178, fsp_l180, -fsp_l197, -fsp_l179])
        fsp_ll41 = self.geom.add_line_loop([fsp_l164, fsp_l166, -fsp_l194, -fsp_l165])
        fsp_ll42 = self.geom.add_line_loop([fsp_l195, -fsp_l199, -fsp_l194, fsp_l193])
        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)
        fsp_s13 = self.geom.add_plane_surface(fsp_ll13, holes=None)
        fsp_s14 = self.geom.add_plane_surface(fsp_ll14, holes=None)
        fsp_s15 = self.geom.add_plane_surface(fsp_ll15, holes=None)
        fsp_s16 = self.geom.add_plane_surface(fsp_ll16, holes=None)
        fsp_s17 = self.geom.add_plane_surface(fsp_ll17, holes=None)
        fsp_s18 = self.geom.add_plane_surface(fsp_ll18, holes=None)
        fsp_s19 = self.geom.add_plane_surface(fsp_ll19, holes=None)
        fsp_s20 = self.geom.add_plane_surface(fsp_ll20, holes=None)
        fsp_s21 = self.geom.add_plane_surface(fsp_ll21, holes=None)
        fsp_s22 = self.geom.add_plane_surface(fsp_ll22, holes=None)
        fsp_s23 = self.geom.add_plane_surface(fsp_ll23, holes=None)
        fsp_s24 = self.geom.add_plane_surface(fsp_ll24, holes=None)
        fsp_s25 = self.geom.add_plane_surface(fsp_ll25, holes=None)
        fsp_s26 = self.geom.add_plane_surface(fsp_ll26, holes=None)
        fsp_s27 = self.geom.add_plane_surface(fsp_ll27, holes=None)
        fsp_s28 = self.geom.add_plane_surface(fsp_ll28, holes=None)
        fsp_s29 = self.geom.add_plane_surface(fsp_ll29, holes=None)
        fsp_s30 = self.geom.add_plane_surface(fsp_ll30, holes=None)
        fsp_s31 = self.geom.add_plane_surface(fsp_ll31, holes=None)
        fsp_s32 = self.geom.add_plane_surface(fsp_ll32, holes=None)
        fsp_s33 = self.geom.add_plane_surface(fsp_ll33, holes=None)
        fsp_s34 = self.geom.add_plane_surface(fsp_ll34, holes=None)
        fsp_s35 = self.geom.add_plane_surface(fsp_ll35, holes=None)
        fsp_s36 = self.geom.add_plane_surface(fsp_ll36, holes=None)
        fsp_s37 = self.geom.add_plane_surface(fsp_ll37, holes=None)
        fsp_s38 = self.geom.add_plane_surface(fsp_ll38, holes=None)
        fsp_s39 = self.geom.add_plane_surface(fsp_ll39, holes=None)
        fsp_s40 = self.geom.add_plane_surface(fsp_ll40, holes=None)
        fsp_s41 = self.geom.add_plane_surface(fsp_ll41, holes=None)
        fsp_s42 = self.geom.add_plane_surface(fsp_ll42, holes=None)

        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16, fsp_s17, fsp_s18, fsp_s19, fsp_s20, fsp_s21, fsp_s22, fsp_s23, fsp_s24, fsp_s25, fsp_s26, fsp_s27, fsp_s28, fsp_s29, fsp_s30, fsp_s31, fsp_s32, fsp_s33, fsp_s34, fsp_s35, fsp_s36, fsp_s37, fsp_s38, fsp_s39, fsp_s40, fsp_s41, fsp_s42]
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())
  
#mesh = pygmsh.generate_mesh(geom)

#meshio.write("spherical.vtk", mesh)


  
#mesh = pygmsh.generate_mesh(geom)

#meshio.write("spherical.vtk", mesh)





