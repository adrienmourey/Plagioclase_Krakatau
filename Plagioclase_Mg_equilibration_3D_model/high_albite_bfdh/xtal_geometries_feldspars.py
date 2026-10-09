"""
Crystal geometries for pygmsh, in the same style as xtal_geometries2.py.

Classes
  microcline_732                 SHAPE 732.SHP, corners/edges taken directly from Microcline_732_orientated.txt
  albite_733                     feldspar habit on the albite cell (733.SHP)
  albite_733_bfdh                BFDH morphology on the albite C-1 cell (733_BFDH.SHP)
  anorthite_habit                feldspar habit on the anorthite cell, c = 14.17 A (Anorthite_habit.SHP)
  anorthite_bfdh_I1bar           BFDH morphology, I-1 lattice (Anorthite_BFDH_I1bar.SHP)
  anorthite_bfdh_P1bar           BFDH morphology, P-1 lattice (Anorthite_BFDH_P1bar.SHP)
  anorthite_bfdh_C1bar_subcell   BFDH morphology on the 7 A C-1 subcell (Anorthite_BFDH_C1bar_subcell.SHP)
  plagioclase_An86_mints         An86 plagioclase with the habit of the MINTS drawings (An86_MINTS_habit.SHP)
  bytownite_dhz_fig210e          bytownite with the habit of Deer, Howie and Zussman Fig. 210e (Bytownite_DHZ_fig210e.SHP)
  high_albite_bfdh               BFDH morphology of high albite, C-1 (HighAlbite_BFDH_C1bar.SHP)

Coordinates use SHAPE's Cartesian frame (c along z, b in the yz plane), with the crystal
centre at the origin. The mesh z axis is parallel to the crystal c axis, and sf = z/fsp_z, where
fsp_z is the half-length of the crystal along z (its largest z coordinate). A class called with z
therefore spans -z to +z along c (total length 2z along the c axis).

Every face's line loop runs anticlockwise when viewed from outside the crystal.
The generated classes were checked for closure (each edge used by exactly two faces,
Euler characteristic V - E + F = 2) and positive enclosed volume.
"""
import pygmsh
import numpy as np
import meshio



class microcline_732(object):
    # Microcline, SHAPE file 732.SHP (title 731Mic), point group 2/m
    # a, b, c = 8.58, 12.97, 7.22 A; alpha, beta, gamma = 90.0000, 116.0000, 90.0000 deg
    # Forms (distance): (1 1 0) 1.2, (1 3 0) 1.08, (0 1 0) 0.8, (0 0 1) 2.1, (-1 0 1) 2.03, (-2 0 1) 1.77, (-1 1 1) 2
    # 20 faces, 36 corners, 54 edges; corners and edges as listed by SHAPE V7.3.0

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.396091 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.800000 # Half-width of the crystal along y (for reference only)
        fsp_z = 2.278310 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([1.3960909763952816*self.sf, -0.0000000522377111*self.sf, 1.6555449870890684*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([1.3960909763952816*self.sf, -0.0000000522377110*self.sf, -1.0900624864940640*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-1.3960909308461023*self.sf, 0.0000000522377094*self.sf, -1.6555450093048893*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-1.3960909308461023*self.sf, 0.0000000522377094*self.sf, 1.0900625495686693*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.9898810391466366*self.sf, 0.6831924018556138*self.sf, 1.8536668110514836*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([0.9898810391466364*self.sf, 0.6831924018556138*self.sf, -1.6342043591278137*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.9996865511230051*self.sf, 0.6667008018811011*self.sf, -1.6389868268606362*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.9898810061912144*self.sf, 0.6831924553977624*self.sf, -1.8536668590194179*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.9898810061912151*self.sf, 0.6831924553977624*self.sf, 1.6342043656683609*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.9996865223878080*self.sf, 0.6667008490435680*self.sf, 1.6389868362294271*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.7815280707160475*self.sf, 0.7999999678827943*self.sf, 1.9552873490130871*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.7815280707160474*self.sf, 0.7999999678827943*self.sf, -1.6626301700741539*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([-0.7815280011977103*self.sf, 0.8000000263678594*self.sf, -1.9552874202670516*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.7815280011977099*self.sf, 0.8000000263678592*self.sf, 1.6626301761909159*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([0.3453245692128166*self.sf, 0.7999999842042753*self.sf, 2.1680380288172572*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.3453245823338381*self.sf, 0.8000000100463812*self.sf, -2.1680380597653093*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.1192341983122505*self.sf, 0.4197453068766941*self.sf, 2.2783096882886991*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-0.1192341985483117*self.sf, 0.4197452944883356*self.sf, -2.2783097077691652*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.8528526691875550*self.sf, 0.4197453432493893*self.sf, 1.8423162414349046*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.8528526689514930*self.sf, 0.4197452581156381*self.sf, -1.8423162609153700*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([0.9898810395308585*self.sf, -0.6831924024123150*self.sf, 1.8536668110514836*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([0.9898810395308584*self.sf, -0.6831924024123149*self.sf, -1.6342043591278137*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([0.9996865514979522*self.sf, -0.6667008024433169*self.sf, -1.6389868268606365*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-0.9898810058069926*self.sf, -0.6831924548410614*self.sf, -1.8536668590194179*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([-0.9898810058069931*self.sf, -0.6831924548410614*self.sf, 1.6342043656683609*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([-0.9996865220128609*self.sf, -0.6667008484813524*self.sf, 1.6389868362294271*self.sf], self.m_r)
        fsp_p27 = self.geom.add_point([0.7815280711659607*self.sf, -0.7999999683223193*self.sf, 1.9552873490130873*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([0.7815280711659609*self.sf, -0.7999999683223192*self.sf, -1.6626301700741537*self.sf], self.m_r)
        fsp_p29 = self.geom.add_point([-0.7815280007477967*self.sf, -0.8000000259283342*self.sf, -1.9552874202670514*self.sf], self.m_r)
        fsp_p30 = self.geom.add_point([-0.7815280007477962*self.sf, -0.8000000259283343*self.sf, 1.6626301761909159*self.sf], self.m_r)
        fsp_p31 = self.geom.add_point([0.3453245696627299*self.sf, -0.7999999843984831*self.sf, 2.1680380288172572*self.sf], self.m_r)
        fsp_p32 = self.geom.add_point([-0.3453245818839245*self.sf, -0.8000000098521736*self.sf, -2.1680380597653093*self.sf], self.m_r)
        fsp_p33 = self.geom.add_point([0.1192341985483119*self.sf, -0.4197453069437505*self.sf, 2.2783096882886991*self.sf], self.m_r)
        fsp_p34 = self.geom.add_point([-0.1192341983122504*self.sf, -0.4197452944212792*self.sf, -2.2783097077691652*self.sf], self.m_r)
        fsp_p35 = self.geom.add_point([-0.8528526689514935*self.sf, -0.4197453427697517*self.sf, 1.8423162414349048*self.sf], self.m_r)
        fsp_p36 = self.geom.add_point([0.8528526691875544*self.sf, -0.4197452585952757*self.sf, -1.8423162609153700*self.sf], self.m_r)

        fsp_l37 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l38 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l39 = self.geom.add_line(fsp_p1, fsp_p21)
        fsp_l40 = self.geom.add_line(fsp_p2, fsp_p7)
        fsp_l41 = self.geom.add_line(fsp_p2, fsp_p23)
        fsp_l42 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l43 = self.geom.add_line(fsp_p3, fsp_p8)
        fsp_l44 = self.geom.add_line(fsp_p3, fsp_p24)
        fsp_l45 = self.geom.add_line(fsp_p4, fsp_p10)
        fsp_l46 = self.geom.add_line(fsp_p4, fsp_p26)
        fsp_l47 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l48 = self.geom.add_line(fsp_p5, fsp_p11)
        fsp_l49 = self.geom.add_line(fsp_p6, fsp_p7)
        fsp_l50 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l51 = self.geom.add_line(fsp_p7, fsp_p20)
        fsp_l52 = self.geom.add_line(fsp_p8, fsp_p9)
        fsp_l53 = self.geom.add_line(fsp_p8, fsp_p13)
        fsp_l54 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l55 = self.geom.add_line(fsp_p9, fsp_p14)
        fsp_l56 = self.geom.add_line(fsp_p10, fsp_p19)
        fsp_l57 = self.geom.add_line(fsp_p11, fsp_p12)
        fsp_l58 = self.geom.add_line(fsp_p11, fsp_p15)
        fsp_l59 = self.geom.add_line(fsp_p12, fsp_p16)
        fsp_l60 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l61 = self.geom.add_line(fsp_p13, fsp_p16)
        fsp_l62 = self.geom.add_line(fsp_p14, fsp_p15)
        fsp_l63 = self.geom.add_line(fsp_p15, fsp_p17)
        fsp_l64 = self.geom.add_line(fsp_p16, fsp_p18)
        fsp_l65 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l66 = self.geom.add_line(fsp_p17, fsp_p33)
        fsp_l67 = self.geom.add_line(fsp_p18, fsp_p20)
        fsp_l68 = self.geom.add_line(fsp_p18, fsp_p34)
        fsp_l69 = self.geom.add_line(fsp_p19, fsp_p35)
        fsp_l70 = self.geom.add_line(fsp_p20, fsp_p36)
        fsp_l71 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l72 = self.geom.add_line(fsp_p21, fsp_p27)
        fsp_l73 = self.geom.add_line(fsp_p22, fsp_p23)
        fsp_l74 = self.geom.add_line(fsp_p22, fsp_p28)
        fsp_l75 = self.geom.add_line(fsp_p23, fsp_p36)
        fsp_l76 = self.geom.add_line(fsp_p24, fsp_p25)
        fsp_l77 = self.geom.add_line(fsp_p24, fsp_p29)
        fsp_l78 = self.geom.add_line(fsp_p25, fsp_p26)
        fsp_l79 = self.geom.add_line(fsp_p25, fsp_p30)
        fsp_l80 = self.geom.add_line(fsp_p26, fsp_p35)
        fsp_l81 = self.geom.add_line(fsp_p27, fsp_p28)
        fsp_l82 = self.geom.add_line(fsp_p27, fsp_p31)
        fsp_l83 = self.geom.add_line(fsp_p28, fsp_p32)
        fsp_l84 = self.geom.add_line(fsp_p29, fsp_p30)
        fsp_l85 = self.geom.add_line(fsp_p29, fsp_p32)
        fsp_l86 = self.geom.add_line(fsp_p30, fsp_p31)
        fsp_l87 = self.geom.add_line(fsp_p31, fsp_p33)
        fsp_l88 = self.geom.add_line(fsp_p32, fsp_p34)
        fsp_l89 = self.geom.add_line(fsp_p33, fsp_p35)
        fsp_l90 = self.geom.add_line(fsp_p34, fsp_p36)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l38, fsp_l37, fsp_l40, -fsp_l49, -fsp_l47])  # face 1 (1 1 0)
        fsp_ll2 = self.geom.add_line_loop([-fsp_l43, fsp_l42, fsp_l45, -fsp_l54, -fsp_l52])  # face 2 (-1 1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l73, -fsp_l41, -fsp_l37, fsp_l39, fsp_l71])  # face 3 (1 -1 0)
        fsp_ll4 = self.geom.add_line_loop([fsp_l78, -fsp_l46, -fsp_l42, fsp_l44, fsp_l76])  # face 4 (-1 -1 0)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l48, fsp_l47, fsp_l50, -fsp_l57])  # face 5 (1 3 0)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l53, fsp_l52, fsp_l55, -fsp_l60])  # face 6 (-1 3 0)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l74, -fsp_l71, fsp_l72, fsp_l81])  # face 7 (1 -3 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l79, -fsp_l76, fsp_l77, fsp_l84])  # face 8 (-1 -3 0)
        fsp_ll9 = self.geom.add_line_loop([fsp_l62, -fsp_l58, fsp_l57, fsp_l59, -fsp_l61, fsp_l60])  # face 9 (0 1 0)
        fsp_ll10 = self.geom.add_line_loop([fsp_l85, -fsp_l83, -fsp_l81, fsp_l82, -fsp_l86, -fsp_l84])  # face 10 (0 -1 0)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l87, -fsp_l82, -fsp_l72, -fsp_l39, fsp_l38, fsp_l48, fsp_l58, fsp_l63, fsp_l66])  # face 11 (0 0 1)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l88, -fsp_l85, -fsp_l77, -fsp_l44, fsp_l43, fsp_l53, fsp_l61, fsp_l64, fsp_l68])  # face 12 (0 0 -1)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l89, -fsp_l66, fsp_l65, fsp_l69])  # face 13 (-1 0 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l68, fsp_l67, fsp_l70, -fsp_l90])  # face 14 (1 0 -1)
        fsp_ll15 = self.geom.add_line_loop([-fsp_l56, -fsp_l45, fsp_l46, fsp_l80, -fsp_l69])  # face 15 (-2 0 1)
        fsp_ll16 = self.geom.add_line_loop([-fsp_l51, -fsp_l40, fsp_l41, fsp_l75, -fsp_l70])  # face 16 (2 0 -1)
        fsp_ll17 = self.geom.add_line_loop([-fsp_l62, -fsp_l55, fsp_l54, fsp_l56, -fsp_l65, -fsp_l63])  # face 17 (-1 1 1)
        fsp_ll18 = self.geom.add_line_loop([-fsp_l59, -fsp_l50, fsp_l49, fsp_l51, -fsp_l67, -fsp_l64])  # face 18 (1 1 -1)
        fsp_ll19 = self.geom.add_line_loop([fsp_l89, -fsp_l80, -fsp_l78, fsp_l79, fsp_l86, fsp_l87])  # face 19 (-1 -1 1)
        fsp_ll20 = self.geom.add_line_loop([fsp_l90, -fsp_l75, -fsp_l73, fsp_l74, fsp_l83, fsp_l88])  # face 20 (1 -1 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16, fsp_s17, fsp_s18, fsp_s19, fsp_s20]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class albite_733(object):
    # Albite, feldspar habit (forms and distances adapted from microcline 732), point group -1
    # Source: 733.SHP (title 733Alb)
    # a, b, c = 8.14, 12.8, 7.16 A; alpha, beta, gamma = 94.3333, 116.5667, 87.6500 deg
    # Forms (distance): (0 1 0) 0.8, (0 0 1) 2.1, (1 1 0) 1.2, (1 -1 0) 1.2, (1 3 0) 1.08, (1 -3 0) 1.08, (-1 0 1) 2.03, (-2 0 1) 1.77, (-1 1 1) 2, (-1 -1 1) 2
    # 20 faces, 36 corners, 54 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.381494 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.806204 # Half-width of the crystal along y (for reference only)
        fsp_z = 2.327536 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.7640240877332332*self.sf, 0.8062036725010842*self.sf, 1.9100583919024197*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.2410164671792897*self.sf, 0.8019748902831294*self.sf, 2.1715891843955730*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.7791025679768174*self.sf, 0.7937267098235474*self.sf, -2.0237630862040050*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.3868948100231912*self.sf, 0.7968979087265797*self.sf, -2.2198871993688059*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.7640240877332332*self.sf, 0.8062036725010843*self.sf, -1.6635044513203874*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.7791025679768174*self.sf, 0.7937267098235473*self.sf, 1.6784383349108030*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.7791025679768173*self.sf, -0.7937267098235473*self.sf, 2.0237630862040055*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([0.3868948100231912*self.sf, -0.7968979087265797*self.sf, 2.2198871993688059*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.7640240877332330*self.sf, -0.8062036725010842*self.sf, -1.9100583919024197*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.2410164671792897*self.sf, -0.8019748902831294*self.sf, -2.1715891843955730*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-0.7640240877332334*self.sf, -0.8062036725010843*self.sf, 1.6635044513203876*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.7791025679768174*self.sf, -0.7937267098235473*self.sf, -1.6784383349108030*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([1.3814943291570807*self.sf, 0.0027421577268755*self.sf, 1.6625524017689604*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([0.9985115723151570*self.sf, 0.6710695065314568*self.sf, 1.8031860735300986*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([1.0053942936184357*self.sf, -0.6596563587342419*self.sf, 1.9005849503678665*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([0.0937375724137114*self.sf, 0.5425829408539519*self.sf, 2.2648015411178273*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.0937375724137114*self.sf, -0.2853213536617616*self.sf, 2.3275363767582475*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-1.3814943291570809*self.sf, -0.0027421577268756*self.sf, -1.6625524017689599*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.9985115723151570*self.sf, -0.6710695065314568*self.sf, -1.8031860735300984*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-1.0053942936184357*self.sf, 0.6596563587342420*self.sf, -1.9005849503678665*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.0937375724137114*self.sf, -0.5425829408539519*self.sf, -2.2648015411178273*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-0.0937375724137114*self.sf, 0.2853213536617616*self.sf, -2.3275363767582475*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([1.3814943291570809*self.sf, 0.0027421577268756*self.sf, -1.1189373677286802*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([0.9985115723151571*self.sf, 0.6710695065314569*self.sf, -1.6166348316894905*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([1.0569663912159728*self.sf, 0.5690624228105596*self.sf, -1.6380998080292664*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([-1.3814943291570809*self.sf, -0.0027421577268756*self.sf, 1.1189373677286802*self.sf], self.m_r)
        fsp_p27 = self.geom.add_point([-0.9985115723151569*self.sf, -0.6710695065314570*self.sf, 1.6166348316894905*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([-1.0569663912159728*self.sf, -0.5690624228105596*self.sf, 1.6380998080292664*self.sf], self.m_r)
        fsp_p29 = self.geom.add_point([1.0053942936184355*self.sf, -0.6596563587342420*self.sf, -1.6206758365452769*self.sf], self.m_r)
        fsp_p30 = self.geom.add_point([-1.0053942936184355*self.sf, 0.6596563587342420*self.sf, 1.6206758365452769*self.sf], self.m_r)
        fsp_p31 = self.geom.add_point([0.9805546635943507*self.sf, -0.6743730166624762*self.sf, -1.6560131770100703*self.sf], self.m_r)
        fsp_p32 = self.geom.add_point([-0.9805546635943505*self.sf, 0.6743730166624761*self.sf, 1.6560131770100706*self.sf], self.m_r)
        fsp_p33 = self.geom.add_point([-0.9004264673679501*self.sf, 0.5345446197130359*self.sf, 1.7841979803083656*self.sf], self.m_r)
        fsp_p34 = self.geom.add_point([-0.9004264673679501*self.sf, -0.2933596748026777*self.sf, 1.8469328159487859*self.sf], self.m_r)
        fsp_p35 = self.geom.add_point([0.9004264673679502*self.sf, -0.5345446197130359*self.sf, -1.7841979803083658*self.sf], self.m_r)
        fsp_p36 = self.geom.add_point([0.9004264673679501*self.sf, 0.2933596748026777*self.sf, -1.8469328159487859*self.sf], self.m_r)

        fsp_l37 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l38 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l39 = self.geom.add_line(fsp_p1, fsp_p14)
        fsp_l40 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l41 = self.geom.add_line(fsp_p2, fsp_p16)
        fsp_l42 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l43 = self.geom.add_line(fsp_p3, fsp_p6)
        fsp_l44 = self.geom.add_line(fsp_p3, fsp_p20)
        fsp_l45 = self.geom.add_line(fsp_p4, fsp_p5)
        fsp_l46 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l47 = self.geom.add_line(fsp_p5, fsp_p24)
        fsp_l48 = self.geom.add_line(fsp_p6, fsp_p32)
        fsp_l49 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l50 = self.geom.add_line(fsp_p7, fsp_p12)
        fsp_l51 = self.geom.add_line(fsp_p7, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l53 = self.geom.add_line(fsp_p8, fsp_p17)
        fsp_l54 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l55 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l56 = self.geom.add_line(fsp_p9, fsp_p19)
        fsp_l57 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l58 = self.geom.add_line(fsp_p10, fsp_p21)
        fsp_l59 = self.geom.add_line(fsp_p11, fsp_p27)
        fsp_l60 = self.geom.add_line(fsp_p12, fsp_p31)
        fsp_l61 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l62 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l63 = self.geom.add_line(fsp_p13, fsp_p23)
        fsp_l64 = self.geom.add_line(fsp_p14, fsp_p24)
        fsp_l65 = self.geom.add_line(fsp_p15, fsp_p29)
        fsp_l66 = self.geom.add_line(fsp_p16, fsp_p17)
        fsp_l67 = self.geom.add_line(fsp_p16, fsp_p33)
        fsp_l68 = self.geom.add_line(fsp_p17, fsp_p34)
        fsp_l69 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l70 = self.geom.add_line(fsp_p18, fsp_p20)
        fsp_l71 = self.geom.add_line(fsp_p18, fsp_p26)
        fsp_l72 = self.geom.add_line(fsp_p19, fsp_p27)
        fsp_l73 = self.geom.add_line(fsp_p20, fsp_p30)
        fsp_l74 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l75 = self.geom.add_line(fsp_p21, fsp_p35)
        fsp_l76 = self.geom.add_line(fsp_p22, fsp_p36)
        fsp_l77 = self.geom.add_line(fsp_p23, fsp_p25)
        fsp_l78 = self.geom.add_line(fsp_p23, fsp_p29)
        fsp_l79 = self.geom.add_line(fsp_p24, fsp_p25)
        fsp_l80 = self.geom.add_line(fsp_p25, fsp_p36)
        fsp_l81 = self.geom.add_line(fsp_p26, fsp_p28)
        fsp_l82 = self.geom.add_line(fsp_p26, fsp_p30)
        fsp_l83 = self.geom.add_line(fsp_p27, fsp_p28)
        fsp_l84 = self.geom.add_line(fsp_p28, fsp_p34)
        fsp_l85 = self.geom.add_line(fsp_p29, fsp_p31)
        fsp_l86 = self.geom.add_line(fsp_p30, fsp_p32)
        fsp_l87 = self.geom.add_line(fsp_p31, fsp_p35)
        fsp_l88 = self.geom.add_line(fsp_p32, fsp_p33)
        fsp_l89 = self.geom.add_line(fsp_p33, fsp_p34)
        fsp_l90 = self.geom.add_line(fsp_p35, fsp_p36)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l40, -fsp_l37, fsp_l38, -fsp_l45, -fsp_l42, fsp_l43])  # face 1 (0 1 0)
        fsp_ll2 = self.geom.add_line_loop([fsp_l57, -fsp_l50, fsp_l49, fsp_l52, -fsp_l55, fsp_l54])  # face 2 (0 -1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l51, -fsp_l62, fsp_l61, -fsp_l39, fsp_l37, fsp_l41, fsp_l66, -fsp_l53, -fsp_l49])  # face 3 (0 0 1)
        fsp_ll4 = self.geom.add_line_loop([fsp_l56, -fsp_l69, fsp_l70, -fsp_l44, fsp_l42, fsp_l46, -fsp_l74, -fsp_l58, -fsp_l54])  # face 4 (0 0 -1)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l61, fsp_l63, fsp_l77, -fsp_l79, -fsp_l64])  # face 5 (1 1 0)
        fsp_ll6 = self.geom.add_line_loop([fsp_l83, -fsp_l81, -fsp_l71, fsp_l69, fsp_l72])  # face 6 (-1 -1 0)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l78, -fsp_l63, fsp_l62, fsp_l65])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l70, fsp_l71, fsp_l82, -fsp_l73])  # face 8 (-1 1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l38, fsp_l39, fsp_l64, -fsp_l47])  # face 9 (1 3 0)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l72, -fsp_l56, fsp_l55, fsp_l59])  # face 10 (-1 -3 0)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l65, -fsp_l51, fsp_l50, fsp_l60, -fsp_l85])  # face 11 (1 -3 0)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l48, -fsp_l43, fsp_l44, fsp_l73, fsp_l86])  # face 12 (-1 3 0)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l68, -fsp_l66, fsp_l67, fsp_l89])  # face 13 (-1 0 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l75, fsp_l74, fsp_l76, -fsp_l90])  # face 14 (1 0 -1)
        fsp_ll15 = self.geom.add_line_loop([-fsp_l88, -fsp_l86, -fsp_l82, fsp_l81, fsp_l84, -fsp_l89])  # face 15 (-2 0 1)
        fsp_ll16 = self.geom.add_line_loop([-fsp_l80, -fsp_l77, fsp_l78, fsp_l85, fsp_l87, fsp_l90])  # face 16 (2 0 -1)
        fsp_ll17 = self.geom.add_line_loop([-fsp_l67, -fsp_l41, fsp_l40, fsp_l48, fsp_l88])  # face 17 (-1 1 1)
        fsp_ll18 = self.geom.add_line_loop([-fsp_l60, -fsp_l57, fsp_l58, fsp_l75, -fsp_l87])  # face 18 (1 -1 -1)
        fsp_ll19 = self.geom.add_line_loop([-fsp_l83, -fsp_l59, -fsp_l52, fsp_l53, fsp_l68, -fsp_l84])  # face 19 (-1 -1 1)
        fsp_ll20 = self.geom.add_line_loop([-fsp_l76, -fsp_l46, fsp_l45, fsp_l47, fsp_l79, fsp_l80])  # face 20 (1 1 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16, fsp_s17, fsp_s18, fsp_s19, fsp_s20]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class albite_733_bfdh(object):
    # Albite, BFDH morphology on the C-1 cell (distance = 1/d_eff, (001) = 1), point group -1
    # Source: 733_BFDH.SHP (title 733AlbBFDH)
    # a, b, c = 8.14, 12.8, 7.16 A; alpha, beta, gamma = 94.3333, 116.5667, 87.6500 deg
    # Forms (distance): (0 0 1) 1, (0 1 0) 1.00147, (1 1 0) 1.00708, (1 -1 0) 1.01412, (-1 -1 1) 1.08048, (-1 1 1) 1.14451, (0 -2 1) 1.36919
    # 14 faces, 24 corners, 36 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.163463 # Half-width of the crystal along x (for reference only)
        #fsp_y = 1.006225 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.245525 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5841006922698811*self.sf, 1.0062254859629440*self.sf, 0.7523788516643328*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.2182307393235995*self.sf, 1.0032672415888666*self.sf, 0.9353327035453581*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([1.1634625756785151*self.sf, -0.0047949441206472*self.sf, 0.5396335322299260*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.6733540218397439*self.sf, -0.8679886380610127*self.sf, 0.8498217192458764*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([-0.2718846618208865*self.sf, 0.1400614879709722*self.sf, 1.2455252241999051*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([0.3074703742820325*self.sf, -0.8709469931629781*self.sf, 1.0327824191391279*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.5841006922698811*self.sf, -1.0062254859629440*self.sf, -0.7523788516643327*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.2182307393235995*self.sf, -1.0032672415888666*self.sf, -0.9353327035453581*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-1.1634625756785151*self.sf, 0.0047949441206472*self.sf, -0.5396335322299260*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.6733540218397439*self.sf, 0.8679886380610127*self.sf, -0.8498217192458764*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.2718846618208866*self.sf, -0.1400614879709723*self.sf, -1.2455252241999053*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-0.3074703742820325*self.sf, 0.8709469931629781*self.sf, -1.0327824191391279*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.5841006922698810*self.sf, 1.0062254859629438*self.sf, -0.5396350273654484*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.6003025891210341*self.sf, 0.9966489839419437*self.sf, 0.5396333857767255*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.6003025891210341*self.sf, 0.9966489839419437*self.sf, -0.7523670249259181*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.2344189415633227*self.sf, 0.9996073390439090*self.sf, -0.9353277248191698*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([-0.5841006922698810*self.sf, -1.0062254859629440*self.sf, 0.5396350273654487*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([0.6003025891210341*self.sf, -0.9966489839419437*self.sf, -0.5396333857767255*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([0.6003025891210341*self.sf, -0.9966489839419437*self.sf, 0.7523670249259181*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.2344189415633227*self.sf, -0.9996073390439090*self.sf, 0.9353277248191698*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([1.1634625756785151*self.sf, -0.0047949441206472*self.sf, -0.7523803467998552*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-1.1634625756785151*self.sf, 0.0047949441206472*self.sf, 0.7523803467998552*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([1.0904179902655202*self.sf, -0.1334432303240492*self.sf, -0.8498259064312725*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-1.0904179902655202*self.sf, 0.1334432303240492*self.sf, 0.8498259064312725*self.sf], self.m_r)

        fsp_l25 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l26 = self.geom.add_line(fsp_p1, fsp_p3)
        fsp_l27 = self.geom.add_line(fsp_p1, fsp_p13)
        fsp_l28 = self.geom.add_line(fsp_p2, fsp_p5)
        fsp_l29 = self.geom.add_line(fsp_p2, fsp_p14)
        fsp_l30 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l31 = self.geom.add_line(fsp_p3, fsp_p21)
        fsp_l32 = self.geom.add_line(fsp_p4, fsp_p6)
        fsp_l33 = self.geom.add_line(fsp_p4, fsp_p19)
        fsp_l34 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l35 = self.geom.add_line(fsp_p5, fsp_p24)
        fsp_l36 = self.geom.add_line(fsp_p6, fsp_p20)
        fsp_l37 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l38 = self.geom.add_line(fsp_p7, fsp_p9)
        fsp_l39 = self.geom.add_line(fsp_p7, fsp_p17)
        fsp_l40 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l41 = self.geom.add_line(fsp_p8, fsp_p18)
        fsp_l42 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l43 = self.geom.add_line(fsp_p9, fsp_p22)
        fsp_l44 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l45 = self.geom.add_line(fsp_p10, fsp_p15)
        fsp_l46 = self.geom.add_line(fsp_p11, fsp_p12)
        fsp_l47 = self.geom.add_line(fsp_p11, fsp_p23)
        fsp_l48 = self.geom.add_line(fsp_p12, fsp_p16)
        fsp_l49 = self.geom.add_line(fsp_p13, fsp_p16)
        fsp_l50 = self.geom.add_line(fsp_p13, fsp_p21)
        fsp_l51 = self.geom.add_line(fsp_p14, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p14, fsp_p24)
        fsp_l53 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l54 = self.geom.add_line(fsp_p17, fsp_p20)
        fsp_l55 = self.geom.add_line(fsp_p17, fsp_p22)
        fsp_l56 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l57 = self.geom.add_line(fsp_p18, fsp_p23)
        fsp_l58 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l59 = self.geom.add_line(fsp_p21, fsp_p23)
        fsp_l60 = self.geom.add_line(fsp_p22, fsp_p24)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l30, -fsp_l26, fsp_l25, fsp_l28, fsp_l34, -fsp_l32])  # face 1 (0 0 1)
        fsp_ll2 = self.geom.add_line_loop([-fsp_l46, -fsp_l40, -fsp_l37, fsp_l38, fsp_l42, fsp_l44])  # face 2 (0 0 -1)
        fsp_ll3 = self.geom.add_line_loop([-fsp_l51, -fsp_l29, -fsp_l25, fsp_l27, fsp_l49, -fsp_l53])  # face 3 (0 1 0)
        fsp_ll4 = self.geom.add_line_loop([-fsp_l54, -fsp_l39, fsp_l37, fsp_l41, fsp_l56, fsp_l58])  # face 4 (0 -1 0)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l50, -fsp_l27, fsp_l26, fsp_l31])  # face 5 (1 1 0)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l38, fsp_l39, fsp_l55, -fsp_l43])  # face 6 (-1 -1 0)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l59, -fsp_l31, fsp_l30, fsp_l33, -fsp_l56, fsp_l57])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([fsp_l51, -fsp_l45, -fsp_l42, fsp_l43, fsp_l60, -fsp_l52])  # face 8 (-1 1 0)
        fsp_ll9 = self.geom.add_line_loop([fsp_l54, -fsp_l36, -fsp_l34, fsp_l35, -fsp_l60, -fsp_l55])  # face 9 (-1 -1 1)
        fsp_ll10 = self.geom.add_line_loop([fsp_l59, -fsp_l47, fsp_l46, fsp_l48, -fsp_l49, fsp_l50])  # face 10 (1 1 -1)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l35, -fsp_l28, fsp_l29, fsp_l52])  # face 11 (-1 1 1)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l41, fsp_l40, fsp_l47, -fsp_l57])  # face 12 (1 -1 -1)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l33, fsp_l32, fsp_l36, -fsp_l58])  # face 13 (0 -2 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l48, -fsp_l44, fsp_l45, fsp_l53])  # face 14 (0 2 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class anorthite_habit(object):
    # Anorthite, feldspar habit analogue of albite_733 (l indices doubled for c = 14.17 A), point group -1
    # Source: Anorthite_habit.SHP (title AnoHabit)
    # a, b, c = 8.18, 12.88, 14.17 A; alpha, beta, gamma = 93.1667, 115.8500, 91.2167 deg
    # Forms (distance): (0 1 0) 0.8, (0 0 1) 2.1, (1 1 0) 1.2, (1 -1 0) 1.2, (1 3 0) 1.08, (1 -3 0) 1.08, (-1 0 2) 2.03, (-1 0 1) 1.77, (-1 1 2) 2, (-1 -1 2) 2
    # 20 faces, 36 corners, 54 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.382364 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.837620 # Half-width of the crystal along y (for reference only)
        fsp_z = 2.314757 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.8189378138420846*self.sf, 0.7596636819708266*self.sf, 1.8979173755589041*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.2866842854367149*self.sf, 0.7865420810587177*self.sf, 2.1561206542642113*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.7247741511521831*self.sf, 0.8376199641887762*self.sf, -2.0322303890577902*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.3444830086984120*self.sf, 0.8184155502188276*self.sf, -2.2167146894969711*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.8189378138420846*self.sf, 0.7596636819708266*self.sf, -1.6599560670111011*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.7247741511521831*self.sf, 0.8376199641887760*self.sf, 1.6720841002104048*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.7247741511521831*self.sf, -0.8376199641887760*self.sf, 2.0322303890577897*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([0.3444830086984119*self.sf, -0.8184155502188276*self.sf, 2.2167146894969716*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.8189378138420846*self.sf, -0.7596636819708266*self.sf, -1.8979173755589041*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.2866842854367149*self.sf, -0.7865420810587177*self.sf, -2.1561206542642113*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-0.8189378138420847*self.sf, -0.7596636819708265*self.sf, 1.6599560670111011*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.7247741511521832*self.sf, -0.8376199641887760*self.sf, -1.6720841002104048*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([1.3823642380178776*self.sf, -0.0172036078317533*self.sf, 1.6659978302173204*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([1.0222851932836550*self.sf, 0.6308298196425203*self.sf, 1.8058304718312062*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([0.9791461526211647*self.sf, -0.7021498530644618*self.sf, 1.9006255599250048*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([0.0920180838304595*self.sf, 0.4558627571570295*self.sf, 2.2693945638130772*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.0920180838304595*self.sf, -0.3640550064110840*self.sf, 2.3147565996703179*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-1.3823642380178776*self.sf, 0.0172036078317533*self.sf, -1.6659978302173204*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-1.0222851932836550*self.sf, -0.6308298196425203*self.sf, -1.8058304718312062*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-0.9791461526211647*self.sf, 0.7021498530644618*self.sf, -1.9006255599250048*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.0920180838304595*self.sf, -0.4558627571570295*self.sf, -2.2693945638130772*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-0.0920180838304595*self.sf, 0.3640550064110839*self.sf, -2.3147565996703174*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([1.3823642380178776*self.sf, -0.0172036078317533*self.sf, -1.1130562233113157*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([1.0222851932836550*self.sf, 0.6308298196425204*self.sf, -1.6214036115215689*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([1.0460548968901202*self.sf, 0.5880515388083358*self.sf, -1.6306343082096375*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([-1.3823642380178776*self.sf, 0.0172036078317533*self.sf, 1.1130562233113157*self.sf], self.m_r)
        fsp_p27 = self.geom.add_point([-1.0222851932836550*self.sf, -0.6308298196425203*self.sf, 1.6214036115215691*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([-1.0460548968901202*self.sf, -0.5880515388083359*self.sf, 1.6306343082096375*self.sf], self.m_r)
        fsp_p29 = self.geom.add_point([0.9791461526211646*self.sf, -0.7021498530644620*self.sf, -1.6402690939562248*self.sf], self.m_r)
        fsp_p30 = self.geom.add_point([0.9870819112567916*self.sf, -0.6886693860978133*self.sf, -1.6448868159771837*self.sf], self.m_r)
        fsp_p31 = self.geom.add_point([-0.9791461526211644*self.sf, 0.7021498530644618*self.sf, 1.6402690939562246*self.sf], self.m_r)
        fsp_p32 = self.geom.add_point([-0.9870819112567916*self.sf, 0.6886693860978133*self.sf, 1.6448868159771837*self.sf], self.m_r)
        fsp_p33 = self.geom.add_point([-0.8851421354424900*self.sf, 0.5052086063524643*self.sf, 1.8017715274276673*self.sf], self.m_r)
        fsp_p34 = self.geom.add_point([-0.8851421354424900*self.sf, -0.3147091572156488*self.sf, 1.8471335632849077*self.sf], self.m_r)
        fsp_p35 = self.geom.add_point([0.8851421354424900*self.sf, -0.5052086063524643*self.sf, -1.8017715274276673*self.sf], self.m_r)
        fsp_p36 = self.geom.add_point([0.8851421354424901*self.sf, 0.3147091572156489*self.sf, -1.8471335632849080*self.sf], self.m_r)

        fsp_l37 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l38 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l39 = self.geom.add_line(fsp_p1, fsp_p14)
        fsp_l40 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l41 = self.geom.add_line(fsp_p2, fsp_p16)
        fsp_l42 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l43 = self.geom.add_line(fsp_p3, fsp_p6)
        fsp_l44 = self.geom.add_line(fsp_p3, fsp_p20)
        fsp_l45 = self.geom.add_line(fsp_p4, fsp_p5)
        fsp_l46 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l47 = self.geom.add_line(fsp_p5, fsp_p24)
        fsp_l48 = self.geom.add_line(fsp_p6, fsp_p31)
        fsp_l49 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l50 = self.geom.add_line(fsp_p7, fsp_p12)
        fsp_l51 = self.geom.add_line(fsp_p7, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l53 = self.geom.add_line(fsp_p8, fsp_p17)
        fsp_l54 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l55 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l56 = self.geom.add_line(fsp_p9, fsp_p19)
        fsp_l57 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l58 = self.geom.add_line(fsp_p10, fsp_p21)
        fsp_l59 = self.geom.add_line(fsp_p11, fsp_p27)
        fsp_l60 = self.geom.add_line(fsp_p12, fsp_p29)
        fsp_l61 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l62 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l63 = self.geom.add_line(fsp_p13, fsp_p23)
        fsp_l64 = self.geom.add_line(fsp_p14, fsp_p24)
        fsp_l65 = self.geom.add_line(fsp_p15, fsp_p29)
        fsp_l66 = self.geom.add_line(fsp_p16, fsp_p17)
        fsp_l67 = self.geom.add_line(fsp_p16, fsp_p33)
        fsp_l68 = self.geom.add_line(fsp_p17, fsp_p34)
        fsp_l69 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l70 = self.geom.add_line(fsp_p18, fsp_p20)
        fsp_l71 = self.geom.add_line(fsp_p18, fsp_p26)
        fsp_l72 = self.geom.add_line(fsp_p19, fsp_p27)
        fsp_l73 = self.geom.add_line(fsp_p20, fsp_p31)
        fsp_l74 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l75 = self.geom.add_line(fsp_p21, fsp_p35)
        fsp_l76 = self.geom.add_line(fsp_p22, fsp_p36)
        fsp_l77 = self.geom.add_line(fsp_p23, fsp_p25)
        fsp_l78 = self.geom.add_line(fsp_p23, fsp_p30)
        fsp_l79 = self.geom.add_line(fsp_p24, fsp_p25)
        fsp_l80 = self.geom.add_line(fsp_p25, fsp_p36)
        fsp_l81 = self.geom.add_line(fsp_p26, fsp_p28)
        fsp_l82 = self.geom.add_line(fsp_p26, fsp_p32)
        fsp_l83 = self.geom.add_line(fsp_p27, fsp_p28)
        fsp_l84 = self.geom.add_line(fsp_p28, fsp_p34)
        fsp_l85 = self.geom.add_line(fsp_p29, fsp_p30)
        fsp_l86 = self.geom.add_line(fsp_p30, fsp_p35)
        fsp_l87 = self.geom.add_line(fsp_p31, fsp_p32)
        fsp_l88 = self.geom.add_line(fsp_p32, fsp_p33)
        fsp_l89 = self.geom.add_line(fsp_p33, fsp_p34)
        fsp_l90 = self.geom.add_line(fsp_p35, fsp_p36)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l40, -fsp_l37, fsp_l38, -fsp_l45, -fsp_l42, fsp_l43])  # face 1 (0 1 0)
        fsp_ll2 = self.geom.add_line_loop([fsp_l57, -fsp_l50, fsp_l49, fsp_l52, -fsp_l55, fsp_l54])  # face 2 (0 -1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l51, -fsp_l62, fsp_l61, -fsp_l39, fsp_l37, fsp_l41, fsp_l66, -fsp_l53, -fsp_l49])  # face 3 (0 0 1)
        fsp_ll4 = self.geom.add_line_loop([fsp_l56, -fsp_l69, fsp_l70, -fsp_l44, fsp_l42, fsp_l46, -fsp_l74, -fsp_l58, -fsp_l54])  # face 4 (0 0 -1)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l61, fsp_l63, fsp_l77, -fsp_l79, -fsp_l64])  # face 5 (1 1 0)
        fsp_ll6 = self.geom.add_line_loop([fsp_l83, -fsp_l81, -fsp_l71, fsp_l69, fsp_l72])  # face 6 (-1 -1 0)
        fsp_ll7 = self.geom.add_line_loop([fsp_l85, -fsp_l78, -fsp_l63, fsp_l62, fsp_l65])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l70, fsp_l71, fsp_l82, -fsp_l87, -fsp_l73])  # face 8 (-1 1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l38, fsp_l39, fsp_l64, -fsp_l47])  # face 9 (1 3 0)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l72, -fsp_l56, fsp_l55, fsp_l59])  # face 10 (-1 -3 0)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l65, -fsp_l51, fsp_l50, fsp_l60])  # face 11 (1 -3 0)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l43, fsp_l44, fsp_l73, -fsp_l48])  # face 12 (-1 3 0)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l66, fsp_l67, fsp_l89, -fsp_l68])  # face 13 (-1 0 2)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l75, fsp_l74, fsp_l76, -fsp_l90])  # face 14 (1 0 -2)
        fsp_ll15 = self.geom.add_line_loop([-fsp_l88, -fsp_l82, fsp_l81, fsp_l84, -fsp_l89])  # face 15 (-1 0 1)
        fsp_ll16 = self.geom.add_line_loop([-fsp_l80, -fsp_l77, fsp_l78, fsp_l86, fsp_l90])  # face 16 (1 0 -1)
        fsp_ll17 = self.geom.add_line_loop([-fsp_l67, -fsp_l41, fsp_l40, fsp_l48, fsp_l87, fsp_l88])  # face 17 (-1 1 2)
        fsp_ll18 = self.geom.add_line_loop([-fsp_l85, -fsp_l60, -fsp_l57, fsp_l58, fsp_l75, -fsp_l86])  # face 18 (1 -1 -2)
        fsp_ll19 = self.geom.add_line_loop([-fsp_l83, -fsp_l59, -fsp_l52, fsp_l53, fsp_l68, -fsp_l84])  # face 19 (-1 -1 2)
        fsp_ll20 = self.geom.add_line_loop([-fsp_l76, -fsp_l46, fsp_l45, fsp_l47, fsp_l79, fsp_l80])  # face 20 (1 1 -2)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16, fsp_s17, fsp_s18, fsp_s19, fsp_s20]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class anorthite_bfdh_I1bar(object):
    # Anorthite, BFDH morphology on the I-1 lattice, point group -1
    # Source: Anorthite_BFDH_I1bar.SHP (title AnoBFDH-I)
    # a, b, c = 8.18, 12.88, 14.17 A; alpha, beta, gamma = 93.1667, 115.8500, 91.2167 deg
    # Forms (distance): (0 -1 1) 1, (0 1 1) 1.07445, (-1 0 1) 1.16081, (1 -1 0) 1.43775, (0 1 0) 1.46063, (1 1 0) 1.5017, (1 0 1) 1.73083
    # 14 faces, 24 corners, 36 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.693881 # Half-width of the crystal along x (for reference only)
        #fsp_y = 1.502089 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.724199 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([-0.2118300728871147*self.sf, 0.1165671342184663*self.sf, 1.7241986140513237*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([1.5037553930487684*self.sf, 0.0299313706924597*self.sf, 0.8919452450018910*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-1.0199515740673095*self.sf, -1.2561882256523702*self.sf, 0.6369156574863046*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.7841209583986929*self.sf, -1.5020887570399557*self.sf, -0.4002601807904866*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([1.5922466895250802*self.sf, -0.1293262117628414*self.sf, 0.6870284669340870*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.9314560476448048*self.sf, -1.4154534207320919*self.sf, 0.4319890842524252*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.2118300728871147*self.sf, -0.1165671342184663*self.sf, -1.7241986140513237*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-1.5037553930487684*self.sf, -0.0299313706924597*self.sf, -0.8919452450018910*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([1.0199515740673095*self.sf, 1.2561882256523702*self.sf, -0.6369156574863046*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.7841209583986929*self.sf, 1.5020887570399557*self.sf, 0.4002601807904866*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-1.5922466895250804*self.sf, 0.1293262117628414*self.sf, -0.6870284669340873*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.9314560476448048*self.sf, 1.4154534207320919*self.sf, -0.4319890842524252*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([-0.8857592353600590*self.sf, 1.3294358957276560*self.sf, 0.6870256229586675*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.7841209583986929*self.sf, 1.5020887570399557*self.sf, 0.4319927515671504*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([0.9314560476448048*self.sf, 1.4154534207320921*self.sf, -0.4002565134757612*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([1.6053894400639417*self.sf, 0.2025770465984814*self.sf, 0.6369229874790427*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.8857592353600590*self.sf, -1.3294358957276560*self.sf, -0.6870256229586673*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([0.7841209583986929*self.sf, -1.5020887570399557*self.sf, -0.4319927515671504*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.9314560476448048*self.sf, -1.4154534207320921*self.sf, 0.4002565134757614*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-1.6053894400639417*self.sf, -0.2025770465984814*self.sf, -0.6369229874790427*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-1.6938807365402537*self.sf, -0.0433194641431803*self.sf, -0.4002573336063520*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([1.6938807365402537*self.sf, 0.0433194641431803*self.sf, 0.4002573336063520*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([1.6938807365402537*self.sf, 0.0433194641431803*self.sf, 0.4320062094112390*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-1.6938807365402537*self.sf, -0.0433194641431803*self.sf, -0.4320062094112390*self.sf], self.m_r)

        fsp_l25 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l26 = self.geom.add_line(fsp_p1, fsp_p3)
        fsp_l27 = self.geom.add_line(fsp_p1, fsp_p13)
        fsp_l28 = self.geom.add_line(fsp_p2, fsp_p5)
        fsp_l29 = self.geom.add_line(fsp_p2, fsp_p16)
        fsp_l30 = self.geom.add_line(fsp_p3, fsp_p6)
        fsp_l31 = self.geom.add_line(fsp_p3, fsp_p21)
        fsp_l32 = self.geom.add_line(fsp_p4, fsp_p5)
        fsp_l33 = self.geom.add_line(fsp_p4, fsp_p6)
        fsp_l34 = self.geom.add_line(fsp_p4, fsp_p18)
        fsp_l35 = self.geom.add_line(fsp_p5, fsp_p23)
        fsp_l36 = self.geom.add_line(fsp_p6, fsp_p19)
        fsp_l37 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l38 = self.geom.add_line(fsp_p7, fsp_p9)
        fsp_l39 = self.geom.add_line(fsp_p7, fsp_p17)
        fsp_l40 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l41 = self.geom.add_line(fsp_p8, fsp_p20)
        fsp_l42 = self.geom.add_line(fsp_p9, fsp_p12)
        fsp_l43 = self.geom.add_line(fsp_p9, fsp_p22)
        fsp_l44 = self.geom.add_line(fsp_p10, fsp_p11)
        fsp_l45 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l46 = self.geom.add_line(fsp_p10, fsp_p14)
        fsp_l47 = self.geom.add_line(fsp_p11, fsp_p24)
        fsp_l48 = self.geom.add_line(fsp_p12, fsp_p15)
        fsp_l49 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l50 = self.geom.add_line(fsp_p13, fsp_p21)
        fsp_l51 = self.geom.add_line(fsp_p14, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l53 = self.geom.add_line(fsp_p16, fsp_p23)
        fsp_l54 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l55 = self.geom.add_line(fsp_p17, fsp_p22)
        fsp_l56 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l57 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l58 = self.geom.add_line(fsp_p20, fsp_p24)
        fsp_l59 = self.geom.add_line(fsp_p21, fsp_p24)
        fsp_l60 = self.geom.add_line(fsp_p22, fsp_p23)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l28, -fsp_l25, fsp_l26, fsp_l30, -fsp_l33, fsp_l32])  # face 1 (0 -1 1)
        fsp_ll2 = self.geom.add_line_loop([fsp_l45, -fsp_l42, -fsp_l38, fsp_l37, fsp_l40, -fsp_l44])  # face 2 (0 1 -1)
        fsp_ll3 = self.geom.add_line_loop([-fsp_l51, -fsp_l49, -fsp_l27, fsp_l25, fsp_l29, -fsp_l52])  # face 3 (0 1 1)
        fsp_ll4 = self.geom.add_line_loop([-fsp_l41, -fsp_l37, fsp_l39, fsp_l54, fsp_l56, fsp_l57])  # face 4 (0 -1 -1)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l26, fsp_l27, fsp_l50, -fsp_l31])  # face 5 (-1 0 1)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l55, -fsp_l39, fsp_l38, fsp_l43])  # face 6 (1 0 -1)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l35, -fsp_l32, fsp_l34, -fsp_l54, fsp_l55, fsp_l60])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l50, fsp_l49, -fsp_l46, fsp_l44, fsp_l47, -fsp_l59])  # face 8 (-1 1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l45, fsp_l46, fsp_l51, -fsp_l48])  # face 9 (0 1 0)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l34, fsp_l33, fsp_l36, -fsp_l56])  # face 10 (0 -1 0)
        fsp_ll11 = self.geom.add_line_loop([fsp_l53, -fsp_l60, -fsp_l43, fsp_l42, fsp_l48, fsp_l52])  # face 11 (1 1 0)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l36, -fsp_l30, fsp_l31, fsp_l59, -fsp_l58, -fsp_l57])  # face 12 (-1 -1 0)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l53, -fsp_l29, fsp_l28, fsp_l35])  # face 13 (1 0 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l40, fsp_l41, fsp_l58, -fsp_l47])  # face 14 (-1 0 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class anorthite_bfdh_P1bar(object):
    # Anorthite, BFDH morphology on the P-1 lattice, point group -1
    # Source: Anorthite_BFDH_P1bar.SHP (title AnoBFDH-P)
    # a, b, c = 8.18, 12.88, 14.17 A; alpha, beta, gamma = 93.1667, 115.8500, 91.2167 deg
    # Forms (distance): (0 1 0) 1, (0 0 1) 1.0098, (0 -1 1) 1.36928, (-1 0 1) 1.58947, (1 0 0) 1.74697, (-1 1 1) 1.86949, (1 -1 0) 1.96868
    # 14 faces, 24 corners, 36 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.746970 # Half-width of the crystal along x (for reference only)
        #fsp_y = 1.084401 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.590745 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([1.7469700000000001*self.sf, 0.9130536161031366*self.sf, 0.2220967890638680*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([-0.8444831045492915*self.sf, 1.0439200304044447*self.sf, 1.4792452202078772*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([0.9453454192626120*self.sf, 0.9535350488223820*self.sf, -1.4792465711213374*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-1.6461039943116607*self.sf, 1.0844012767322537*self.sf, -0.2220999305184906*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([1.7469699999999997*self.sf, 0.9130536161031365*self.sf, -0.3231279751988761*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-1.6461039943116607*self.sf, 1.0844012767322539*self.sf, 0.3231319474815340*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-1.7469700000000001*self.sf, -0.9130536161031366*self.sf, -0.2220967890638680*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([0.8444831045492915*self.sf, -1.0439200304044447*self.sf, -1.4792452202078772*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.9453454192626120*self.sf, -0.9535350488223820*self.sf, 1.4792465711213374*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([1.6461039943116607*self.sf, -1.0844012767322537*self.sf, 0.2220999305184906*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-1.7469699999999997*self.sf, -0.9130536161031365*self.sf, 0.3231279751988761*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([1.6461039943116607*self.sf, -1.0844012767322539*self.sf, -0.3231319474815340*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([-0.8624746726862006*self.sf, -0.8127625766469387*self.sf, 1.5907445948946952*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([1.7289747408880716*self.sf, -0.9436288045568104*self.sf, 0.3335979542918484*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.8624746726862007*self.sf, 1.0133577677293963*self.sf, 1.4897142962209013*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([1.7469699999999997*self.sf, -0.9130602720253943*self.sf, 0.3231267305450785*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.8624746726862006*self.sf, 0.8127625766469387*self.sf, -1.5907445948946952*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-1.7289747408880716*self.sf, 0.9436288045568104*self.sf, -0.3335979542918484*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([0.8624746726862007*self.sf, -1.0133577677293963*self.sf, -1.4897142962209013*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-1.7469699999999999*self.sf, 0.9130602720253944*self.sf, -0.3231267305450786*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-1.7469699999999999*self.sf, 0.9130602720253945*self.sf, 0.2220980337176657*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([-1.6640955624485700*self.sf, 1.0538390140572058*self.sf, 0.3336010234945579*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([1.7469699999999999*self.sf, -0.9130602720253945*self.sf, -0.2220980337176658*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([1.6640955624485700*self.sf, -1.0538390140572058*self.sf, -0.3336010234945579*self.sf], self.m_r)

        fsp_l25 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l26 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l27 = self.geom.add_line(fsp_p1, fsp_p16)
        fsp_l28 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l29 = self.geom.add_line(fsp_p2, fsp_p15)
        fsp_l30 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l31 = self.geom.add_line(fsp_p3, fsp_p5)
        fsp_l32 = self.geom.add_line(fsp_p3, fsp_p17)
        fsp_l33 = self.geom.add_line(fsp_p4, fsp_p6)
        fsp_l34 = self.geom.add_line(fsp_p4, fsp_p18)
        fsp_l35 = self.geom.add_line(fsp_p5, fsp_p23)
        fsp_l36 = self.geom.add_line(fsp_p6, fsp_p22)
        fsp_l37 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l38 = self.geom.add_line(fsp_p7, fsp_p11)
        fsp_l39 = self.geom.add_line(fsp_p7, fsp_p20)
        fsp_l40 = self.geom.add_line(fsp_p8, fsp_p12)
        fsp_l41 = self.geom.add_line(fsp_p8, fsp_p19)
        fsp_l42 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l43 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l44 = self.geom.add_line(fsp_p9, fsp_p13)
        fsp_l45 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l46 = self.geom.add_line(fsp_p10, fsp_p14)
        fsp_l47 = self.geom.add_line(fsp_p11, fsp_p21)
        fsp_l48 = self.geom.add_line(fsp_p12, fsp_p24)
        fsp_l49 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l50 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l51 = self.geom.add_line(fsp_p14, fsp_p16)
        fsp_l52 = self.geom.add_line(fsp_p15, fsp_p22)
        fsp_l53 = self.geom.add_line(fsp_p16, fsp_p23)
        fsp_l54 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l55 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l56 = self.geom.add_line(fsp_p18, fsp_p20)
        fsp_l57 = self.geom.add_line(fsp_p19, fsp_p24)
        fsp_l58 = self.geom.add_line(fsp_p20, fsp_p21)
        fsp_l59 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l60 = self.geom.add_line(fsp_p23, fsp_p24)

        fsp_ll1 = self.geom.add_line_loop([-fsp_l28, -fsp_l25, fsp_l26, -fsp_l31, fsp_l30, fsp_l33])  # face 1 (0 1 0)
        fsp_ll2 = self.geom.add_line_loop([-fsp_l42, fsp_l43, -fsp_l38, fsp_l37, fsp_l40, -fsp_l45])  # face 2 (0 -1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l49, fsp_l51, -fsp_l27, fsp_l25, fsp_l29, -fsp_l50])  # face 3 (0 0 1)
        fsp_ll4 = self.geom.add_line_loop([-fsp_l41, -fsp_l37, fsp_l39, -fsp_l56, -fsp_l54, fsp_l55])  # face 4 (0 0 -1)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l44, fsp_l42, fsp_l46, -fsp_l49])  # face 5 (0 -1 1)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l34, -fsp_l30, fsp_l32, fsp_l54])  # face 6 (0 1 -1)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l47, -fsp_l43, fsp_l44, fsp_l50, fsp_l52, -fsp_l59])  # face 7 (-1 0 1)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l57, -fsp_l55, -fsp_l32, fsp_l31, fsp_l35, fsp_l60])  # face 8 (1 0 -1)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l35, -fsp_l26, fsp_l27, fsp_l53])  # face 9 (1 0 0)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l39, fsp_l38, fsp_l47, -fsp_l58])  # face 10 (-1 0 0)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l52, -fsp_l29, fsp_l28, fsp_l36])  # face 11 (-1 1 1)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l40, fsp_l41, fsp_l57, -fsp_l48])  # face 12 (1 -1 -1)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l53, -fsp_l51, -fsp_l46, fsp_l45, fsp_l48, -fsp_l60])  # face 13 (1 -1 0)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l36, -fsp_l33, fsp_l34, fsp_l56, fsp_l58, fsp_l59])  # face 14 (-1 1 0)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class anorthite_bfdh_C1bar_subcell(object):
    # Anorthite, BFDH morphology on the 7 A C-1 subcell (indices given on the 14.17 A cell), point group -1
    # Source: Anorthite_BFDH_C1bar_subcell.SHP (title AnoBFDH-sub)
    # a, b, c = 8.18, 12.88, 14.17 A; alpha, beta, gamma = 93.1667, 115.8500, 91.2167 deg
    # Forms (distance): (1 -1 0) 1, (0 1 0) 1.01591, (0 0 1) 1.02587, (1 1 0) 1.04448, (-1 -1 2) 1.1275, (-1 1 2) 1.15314, (0 -1 1) 1.39106
    # 14 faces, 24 corners, 36 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.178147 # Half-width of the crystal along x (for reference only)
        #fsp_y = 1.044746 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.275309 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5453827494555822*self.sf, -1.0447459632662712*self.sf, -0.5663740283595546*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.5453827494555820*self.sf, -1.0447459632662710*self.sf, 0.7723050424383265*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([1.1781473800576754*self.sf, 0.0301308181671510*self.sf, 0.5663800703093128*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.6295784399945483*self.sf, -0.9017228098214031*self.sf, 0.8855857053158115*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([1.1781473800576752*self.sf, 0.0301308181671510*self.sf, -0.7723100672094197*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([1.0939574315061609*self.sf, -0.1128825813681383*self.sf, -0.8855830045591413*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.5453827494555822*self.sf, 1.0447459632662712*self.sf, 0.5663740283595546*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.5453827494555820*self.sf, 1.0447459632662710*self.sf, -0.7723050424383265*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-1.1781473800576754*self.sf, -0.0301308181671510*self.sf, -0.5663800703093128*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.6295784399945483*self.sf, 0.9017228098214031*self.sf, -0.8855857053158115*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-1.1781473800576752*self.sf, -0.0301308181671510*self.sf, 0.7723100672094197*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-1.0939574315061613*self.sf, 0.1128825813681383*self.sf, 0.8855830045591412*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.6478598195452708*self.sf, 0.9844881194632357*self.sf, 0.7723112720435450*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([0.2690018824891375*self.sf, 1.0036201576487189*self.sf, 0.9561003054474155*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([0.6478598195452707*self.sf, 0.9844881194632357*self.sf, -0.5663788654751876*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.1665133284245447*self.sf, 1.0256133451487741*self.sf, -0.9560996468714428*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([-0.6478598195452708*self.sf, -0.9844881194632357*self.sf, -0.7723112720435450*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-0.2690018824891375*self.sf, -1.0036201576487189*self.sf, -0.9561003054474155*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.6478598195452707*self.sf, -0.9844881194632357*self.sf, 0.5663788654751876*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.1665133284245447*self.sf, -1.0256133451487741*self.sf, 0.9560996468714428*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.2795727995614414*self.sf, 0.0717567757505863*self.sf, 1.2753092816470022*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([0.2507090189635109*self.sf, -0.8825901917039055*self.sf, 1.0693803097489276*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([0.2795727995614414*self.sf, -0.0717567757505863*self.sf, -1.2753092816470024*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-0.2507090189635109*self.sf, 0.8825901917039055*self.sf, -1.0693803097489276*self.sf], self.m_r)

        fsp_l25 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l26 = self.geom.add_line(fsp_p1, fsp_p6)
        fsp_l27 = self.geom.add_line(fsp_p1, fsp_p18)
        fsp_l28 = self.geom.add_line(fsp_p2, fsp_p4)
        fsp_l29 = self.geom.add_line(fsp_p2, fsp_p20)
        fsp_l30 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l31 = self.geom.add_line(fsp_p3, fsp_p5)
        fsp_l32 = self.geom.add_line(fsp_p3, fsp_p13)
        fsp_l33 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l34 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l35 = self.geom.add_line(fsp_p5, fsp_p15)
        fsp_l36 = self.geom.add_line(fsp_p6, fsp_p23)
        fsp_l37 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l38 = self.geom.add_line(fsp_p7, fsp_p12)
        fsp_l39 = self.geom.add_line(fsp_p7, fsp_p14)
        fsp_l40 = self.geom.add_line(fsp_p8, fsp_p10)
        fsp_l41 = self.geom.add_line(fsp_p8, fsp_p16)
        fsp_l42 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l43 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l44 = self.geom.add_line(fsp_p9, fsp_p17)
        fsp_l45 = self.geom.add_line(fsp_p10, fsp_p24)
        fsp_l46 = self.geom.add_line(fsp_p11, fsp_p12)
        fsp_l47 = self.geom.add_line(fsp_p11, fsp_p19)
        fsp_l48 = self.geom.add_line(fsp_p12, fsp_p21)
        fsp_l49 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l50 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l51 = self.geom.add_line(fsp_p14, fsp_p21)
        fsp_l52 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l53 = self.geom.add_line(fsp_p16, fsp_p24)
        fsp_l54 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l55 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l56 = self.geom.add_line(fsp_p18, fsp_p23)
        fsp_l57 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l58 = self.geom.add_line(fsp_p20, fsp_p22)
        fsp_l59 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l60 = self.geom.add_line(fsp_p23, fsp_p24)

        fsp_ll1 = self.geom.add_line_loop([fsp_l30, -fsp_l28, -fsp_l25, fsp_l26, -fsp_l34, -fsp_l31])  # face 1 (1 -1 0)
        fsp_ll2 = self.geom.add_line_loop([fsp_l46, -fsp_l38, fsp_l37, fsp_l40, -fsp_l42, fsp_l43])  # face 2 (-1 1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l52, -fsp_l41, -fsp_l37, fsp_l39, -fsp_l49, fsp_l50])  # face 3 (0 1 0)
        fsp_ll4 = self.geom.add_line_loop([fsp_l54, -fsp_l27, fsp_l25, fsp_l29, -fsp_l57, -fsp_l55])  # face 4 (0 -1 0)
        fsp_ll5 = self.geom.add_line_loop([fsp_l59, -fsp_l33, -fsp_l30, fsp_l32, fsp_l49, fsp_l51])  # face 5 (0 0 1)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l54, -fsp_l44, fsp_l42, fsp_l45, -fsp_l60, -fsp_l56])  # face 6 (0 0 -1)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l32, fsp_l31, fsp_l35, -fsp_l50])  # face 7 (1 1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l47, -fsp_l43, fsp_l44, fsp_l55])  # face 8 (-1 -1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l59, -fsp_l48, -fsp_l46, fsp_l47, fsp_l57, fsp_l58])  # face 9 (-1 -1 2)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l52, -fsp_l35, fsp_l34, fsp_l36, fsp_l60, -fsp_l53])  # face 10 (1 1 -2)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l51, -fsp_l39, fsp_l38, fsp_l48])  # face 11 (-1 1 2)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l26, fsp_l27, fsp_l56, -fsp_l36])  # face 12 (1 -1 -2)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l29, fsp_l28, fsp_l33, -fsp_l58])  # face 13 (0 -1 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l45, -fsp_l40, fsp_l41, fsp_l53])  # face 14 (0 1 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class plagioclase_An86_mints(object):
    # Plagioclase An86 (bytownite), habit of the MINTS plagioclase drawings: {001}, {010}, {110}, {1-10}
    # Cell (Wikipedia bytownite): a, b, c = 8.178, 12.870, 14.187 A; alpha, beta, gamma = 93.5, 115.9, 90.63 deg; P-1/I-1
    # Forms (distance): (0 0 1) 1, (0 1 0) 0.42857, (1 1 0) 0.71429, (1 -1 0) 0.71429  (proportions matched by eye to the drawing)
    # 8 faces, 12 corners, 18 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 0.822992 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.452569 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.516215 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5924519066361870*self.sf, 0.4040752929702711*self.sf, 0.8004801733816792*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([-0.5626687302290471*self.sf, 0.4525693657301654*self.sf, 1.3618706830254386*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.5924519066361874*self.sf, -0.4040752929702711*self.sf, 1.4288164619108645*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.5626687302290471*self.sf, -0.4525693657301653*self.sf, 0.8674259522671049*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.8229917570247307*self.sf, -0.0085242439277198*self.sf, 0.7130811381817721*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.8229917570247306*self.sf, 0.0085242439277195*self.sf, 1.5162154971107713*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.5924519066361871*self.sf, 0.4040752929702712*self.sf, -1.4288164619108643*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.5626687302290471*self.sf, 0.4525693657301653*self.sf, -0.8674259522671051*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.5924519066361870*self.sf, -0.4040752929702711*self.sf, -0.8004801733816792*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([0.5626687302290471*self.sf, -0.4525693657301654*self.sf, -1.3618706830254386*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.8229917570247306*self.sf, -0.0085242439277198*self.sf, -1.5162154971107717*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-0.8229917570247306*self.sf, 0.0085242439277198*self.sf, -0.7130811381817718*self.sf], self.m_r)

        fsp_l13 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l14 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l15 = self.geom.add_line(fsp_p1, fsp_p7)
        fsp_l16 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l17 = self.geom.add_line(fsp_p2, fsp_p8)
        fsp_l18 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l19 = self.geom.add_line(fsp_p3, fsp_p6)
        fsp_l20 = self.geom.add_line(fsp_p3, fsp_p9)
        fsp_l21 = self.geom.add_line(fsp_p4, fsp_p5)
        fsp_l22 = self.geom.add_line(fsp_p4, fsp_p10)
        fsp_l23 = self.geom.add_line(fsp_p5, fsp_p11)
        fsp_l24 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l25 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l26 = self.geom.add_line(fsp_p7, fsp_p11)
        fsp_l27 = self.geom.add_line(fsp_p8, fsp_p12)
        fsp_l28 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l29 = self.geom.add_line(fsp_p9, fsp_p12)
        fsp_l30 = self.geom.add_line(fsp_p10, fsp_p11)

        fsp_ll1 = self.geom.add_line_loop([fsp_l21, -fsp_l14, fsp_l13, fsp_l16, -fsp_l19, fsp_l18])  # face 1 (0 0 1)
        fsp_ll2 = self.geom.add_line_loop([fsp_l29, -fsp_l27, -fsp_l25, fsp_l26, -fsp_l30, -fsp_l28])  # face 2 (0 0 -1)
        fsp_ll3 = self.geom.add_line_loop([-fsp_l13, fsp_l15, fsp_l25, -fsp_l17])  # face 3 (0 1 0)
        fsp_ll4 = self.geom.add_line_loop([-fsp_l18, fsp_l20, fsp_l28, -fsp_l22])  # face 4 (0 -1 0)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l15, fsp_l14, fsp_l23, -fsp_l26])  # face 5 (1 1 0)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l20, fsp_l19, fsp_l24, -fsp_l29])  # face 6 (-1 -1 0)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l23, -fsp_l21, fsp_l22, fsp_l30])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l16, fsp_l17, fsp_l27, -fsp_l24])  # face 8 (-1 1 0)

        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class bytownite_dhz_fig210e(object):
    # Bytownite (An86 cell), habit fitted to Deer, Howie and Zussman (1992) Fig. 210e
    # Faces: {001}, {010}, {110}, {1-10} and {-1 0 2} on the 14 A cell (= {-1 0 1}, x, on the 7 A morphological axes)
    # Cell (Wikipedia bytownite): a, b, c = 8.178, 12.870, 14.187 A; alpha, beta, gamma = 93.5, 115.9, 90.63 deg
    # Forms (distance): (0 0 1) 1, (0 1 0) 0.493, (1 1 0) 0.887, (1 -1 0) 0.87, (-1 0 2) 1.144  (distances fitted to the drawing, rms 16-22 px)
    # 10 faces, 16 corners, 24 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.012369 # Half-width of the crystal along x (for reference only)
        #fsp_y = 0.522903 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.222645 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.7581483892839178*self.sf, 0.4616058052787491*self.sf, 0.7160072343355126*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([-0.1601171235111266*self.sf, 0.5001562676020880*self.sf, 1.1622857681649115*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([0.7019311941067589*self.sf, -0.5229026112507509*self.sf, 0.8036883469367573*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.1601171235111266*self.sf, -0.4867122493933483*self.sf, 1.2226452324083750*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([1.0123686976521529*self.sf, 0.0066251096375471*self.sf, 0.6196307959793665*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.7019311941067589*self.sf, 0.5229026112507509*self.sf, -0.8036883469367573*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([0.1601171235111266*self.sf, 0.4867122493933483*self.sf, -1.2226452324083750*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.7581483892839178*self.sf, -0.4616058052787490*self.sf, -0.7160072343355126*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([0.1601171235111266*self.sf, -0.5001562676020880*self.sf, -1.1622857681649115*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-1.0123686976521529*self.sf, -0.0066251096375471*self.sf, -0.6196307959793667*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.7581483892839177*self.sf, 0.4616058052787489*self.sf, -0.9361367792235372*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-0.7019311941067586*self.sf, 0.5229026112507509*self.sf, 0.9027101904908289*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([-0.7581483892839179*self.sf, -0.4616058052787491*self.sf, 0.9361367792235374*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([0.7019311941067586*self.sf, -0.5229026112507509*self.sf, -0.9027101904908289*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([1.0123686976521529*self.sf, 0.0066251096375470*self.sf, -0.7871683236416991*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-1.0123686976521529*self.sf, -0.0066251096375471*self.sf, 0.7871683236416992*self.sf], self.m_r)

        fsp_l17 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l18 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l19 = self.geom.add_line(fsp_p1, fsp_p11)
        fsp_l20 = self.geom.add_line(fsp_p2, fsp_p4)
        fsp_l21 = self.geom.add_line(fsp_p2, fsp_p12)
        fsp_l22 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l23 = self.geom.add_line(fsp_p3, fsp_p5)
        fsp_l24 = self.geom.add_line(fsp_p3, fsp_p14)
        fsp_l25 = self.geom.add_line(fsp_p4, fsp_p13)
        fsp_l26 = self.geom.add_line(fsp_p5, fsp_p15)
        fsp_l27 = self.geom.add_line(fsp_p6, fsp_p7)
        fsp_l28 = self.geom.add_line(fsp_p6, fsp_p10)
        fsp_l29 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l30 = self.geom.add_line(fsp_p7, fsp_p9)
        fsp_l31 = self.geom.add_line(fsp_p7, fsp_p11)
        fsp_l32 = self.geom.add_line(fsp_p8, fsp_p9)
        fsp_l33 = self.geom.add_line(fsp_p8, fsp_p10)
        fsp_l34 = self.geom.add_line(fsp_p8, fsp_p13)
        fsp_l35 = self.geom.add_line(fsp_p9, fsp_p14)
        fsp_l36 = self.geom.add_line(fsp_p10, fsp_p16)
        fsp_l37 = self.geom.add_line(fsp_p11, fsp_p15)
        fsp_l38 = self.geom.add_line(fsp_p12, fsp_p16)
        fsp_l39 = self.geom.add_line(fsp_p13, fsp_p16)
        fsp_l40 = self.geom.add_line(fsp_p14, fsp_p15)

        fsp_ll1 = self.geom.add_line_loop([fsp_l23, -fsp_l18, fsp_l17, fsp_l20, -fsp_l22])  # face 1 (0 0 1)
        fsp_ll2 = self.geom.add_line_loop([fsp_l33, -fsp_l28, fsp_l27, fsp_l30, -fsp_l32])  # face 2 (0 0 -1)
        fsp_ll3 = self.geom.add_line_loop([-fsp_l21, -fsp_l17, fsp_l19, -fsp_l31, -fsp_l27, fsp_l29])  # face 3 (0 1 0)
        fsp_ll4 = self.geom.add_line_loop([fsp_l35, -fsp_l24, fsp_l22, fsp_l25, -fsp_l34, fsp_l32])  # face 4 (0 -1 0)
        fsp_ll5 = self.geom.add_line_loop([-fsp_l19, fsp_l18, fsp_l26, -fsp_l37])  # face 5 (1 1 0)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l36, -fsp_l33, fsp_l34, fsp_l39])  # face 6 (-1 -1 0)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l26, -fsp_l23, fsp_l24, fsp_l40])  # face 7 (1 -1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l29, fsp_l28, fsp_l36, -fsp_l38])  # face 8 (-1 1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l25, -fsp_l20, fsp_l21, fsp_l38, -fsp_l39])  # face 9 (-1 0 2)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l35, -fsp_l30, fsp_l31, fsp_l37, -fsp_l40])  # face 10 (1 0 -2)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())


class high_albite_bfdh(object):
    # High albite, BFDH morphology on the C-1 cell (distance = 1/d_eff; h+k odd planes halved by C-centring)
    # Source: HighAlbite_BFDH_C1bar.SHP (title HiAlbBFDH-C1bar)
    # a, b, c = 8.16, 12.87, 7.11 A; alpha, beta, gamma = 93.5, 116.4, 90.3 deg (beta given as 12.87; 116.4 used)
    # Forms (distance): (1 -1 0) 1, (0 1 0) 1.00482, (0 0 1) 1.0153, (1 1 0) 1.03163, (-1 -1 1) 1.10616, (-1 1 1) 1.14179, (0 -2 1) 1.377
    # 14 faces, 24 corners, 36 edges

    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Half-length of the crystal along z (= c axis); the crystal spans -z..+z and the rest is scaled accordingly.
        self.m_r = m_r # Mesh resolution

        #fsp_x = 1.168871 # Half-width of the crystal along x (for reference only)
        #fsp_y = 1.025803 # Half-width of the crystal along y (for reference only)
        fsp_z = 1.262443 # Half-length of the crystal along z, which is parallel to the c axis

        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5607675702054730*self.sf, -1.0258028371697276*self.sf, -0.5518278821593898*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.5607675702054731*self.sf, -1.0258028371697276*self.sf, 0.7601588174596381*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([1.1688714828938054*self.sf, 0.0216371406833128*self.sf, 0.5518320060394046*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.6425156246080107*self.sf, -0.8849943692314384*self.sf, 0.8699075587628318*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([1.1688714828938054*self.sf, 0.0216371406833128*self.sf, -0.7601620833288377*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([1.0871272242929308*self.sf, -0.1191647891052832*self.sf, -0.8699057286764418*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.5607675702054730*self.sf, 1.0258028371697276*self.sf, 0.5518278821593898*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.5607675702054731*self.sf, 1.0258028371697276*self.sf, -0.7601588174596381*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-1.1688714828938054*self.sf, -0.0216371406833128*self.sf, -0.5518320060394046*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.6425156246080107*self.sf, 0.8849943692314384*self.sf, -0.8699075587628318*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-1.1688714828938052*self.sf, -0.0216371406833128*self.sf, 0.7601620833288376*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-1.0871272242929313*self.sf, 0.1191647891052832*self.sf, 0.8699057286764418*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.6335470407614586*self.sf, 0.9825188697028532*self.sf, 0.7601604025050789*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([0.2627127173247049*self.sf, 0.9959585284358579*self.sf, 0.9443650129149759*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([0.6335470407614585*self.sf, 0.9825188697028531*self.sf, -0.5518336868631635*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.1899256551653928*self.sf, 1.0123629033042678*self.sf, -0.9443671988482304*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([-0.6335470407614586*self.sf, -0.9825188697028532*self.sf, -0.7601604025050789*self.sf], self.m_r)
        fsp_p18 = self.geom.add_point([-0.2627127173247049*self.sf, -0.9959585284358579*self.sf, -0.9443650129149759*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.6335470407614585*self.sf, -0.9825188697028531*self.sf, 0.5518336868631635*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.1899256551653928*self.sf, -1.0123629033042678*self.sf, 0.9443671988482304*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([-0.2636469367627530*self.sf, 0.0893204803714134*self.sf, 1.2624428594320281*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([0.2716737095679305*self.sf, -0.8715544353659783*self.sf, 1.0541159401514244*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([0.2636469367627530*self.sf, -0.0893204803714134*self.sf, -1.2624428594320281*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([-0.2716737095679305*self.sf, 0.8715544353659783*self.sf, -1.0541159401514244*self.sf], self.m_r)

        fsp_l25 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l26 = self.geom.add_line(fsp_p1, fsp_p6)
        fsp_l27 = self.geom.add_line(fsp_p1, fsp_p18)
        fsp_l28 = self.geom.add_line(fsp_p2, fsp_p4)
        fsp_l29 = self.geom.add_line(fsp_p2, fsp_p20)
        fsp_l30 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l31 = self.geom.add_line(fsp_p3, fsp_p5)
        fsp_l32 = self.geom.add_line(fsp_p3, fsp_p13)
        fsp_l33 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l34 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l35 = self.geom.add_line(fsp_p5, fsp_p15)
        fsp_l36 = self.geom.add_line(fsp_p6, fsp_p23)
        fsp_l37 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l38 = self.geom.add_line(fsp_p7, fsp_p12)
        fsp_l39 = self.geom.add_line(fsp_p7, fsp_p14)
        fsp_l40 = self.geom.add_line(fsp_p8, fsp_p10)
        fsp_l41 = self.geom.add_line(fsp_p8, fsp_p16)
        fsp_l42 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l43 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l44 = self.geom.add_line(fsp_p9, fsp_p17)
        fsp_l45 = self.geom.add_line(fsp_p10, fsp_p24)
        fsp_l46 = self.geom.add_line(fsp_p11, fsp_p12)
        fsp_l47 = self.geom.add_line(fsp_p11, fsp_p19)
        fsp_l48 = self.geom.add_line(fsp_p12, fsp_p21)
        fsp_l49 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l50 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l51 = self.geom.add_line(fsp_p14, fsp_p21)
        fsp_l52 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l53 = self.geom.add_line(fsp_p16, fsp_p24)
        fsp_l54 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l55 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l56 = self.geom.add_line(fsp_p18, fsp_p23)
        fsp_l57 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l58 = self.geom.add_line(fsp_p20, fsp_p22)
        fsp_l59 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l60 = self.geom.add_line(fsp_p23, fsp_p24)

        fsp_ll1 = self.geom.add_line_loop([fsp_l30, -fsp_l28, -fsp_l25, fsp_l26, -fsp_l34, -fsp_l31])  # face 1 (1 -1 0)
        fsp_ll2 = self.geom.add_line_loop([fsp_l46, -fsp_l38, fsp_l37, fsp_l40, -fsp_l42, fsp_l43])  # face 2 (-1 1 0)
        fsp_ll3 = self.geom.add_line_loop([fsp_l52, -fsp_l41, -fsp_l37, fsp_l39, -fsp_l49, fsp_l50])  # face 3 (0 1 0)
        fsp_ll4 = self.geom.add_line_loop([fsp_l54, -fsp_l27, fsp_l25, fsp_l29, -fsp_l57, -fsp_l55])  # face 4 (0 -1 0)
        fsp_ll5 = self.geom.add_line_loop([fsp_l59, -fsp_l33, -fsp_l30, fsp_l32, fsp_l49, fsp_l51])  # face 5 (0 0 1)
        fsp_ll6 = self.geom.add_line_loop([-fsp_l54, -fsp_l44, fsp_l42, fsp_l45, -fsp_l60, -fsp_l56])  # face 6 (0 0 -1)
        fsp_ll7 = self.geom.add_line_loop([-fsp_l50, -fsp_l32, fsp_l31, fsp_l35])  # face 7 (1 1 0)
        fsp_ll8 = self.geom.add_line_loop([-fsp_l43, fsp_l44, fsp_l55, -fsp_l47])  # face 8 (-1 -1 0)
        fsp_ll9 = self.geom.add_line_loop([-fsp_l59, -fsp_l48, -fsp_l46, fsp_l47, fsp_l57, fsp_l58])  # face 9 (-1 -1 1)
        fsp_ll10 = self.geom.add_line_loop([-fsp_l52, -fsp_l35, fsp_l34, fsp_l36, fsp_l60, -fsp_l53])  # face 10 (1 1 -1)
        fsp_ll11 = self.geom.add_line_loop([-fsp_l51, -fsp_l39, fsp_l38, fsp_l48])  # face 11 (-1 1 1)
        fsp_ll12 = self.geom.add_line_loop([-fsp_l26, fsp_l27, fsp_l56, -fsp_l36])  # face 12 (1 -1 -1)
        fsp_ll13 = self.geom.add_line_loop([-fsp_l29, fsp_l28, fsp_l33, -fsp_l58])  # face 13 (0 -2 1)
        fsp_ll14 = self.geom.add_line_loop([-fsp_l45, -fsp_l40, fsp_l41, fsp_l53])  # face 14 (0 2 -1)

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

        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14]

        return surfaces

    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())

    def phys_surface(self):
        return self.geom.add_physical(self.surfaces())

    def volume(self):
        return self.geom.add_volume(self.surfaceloop())

    def phys_volume(self):
        return self.geom.add_physical(self.volume())

