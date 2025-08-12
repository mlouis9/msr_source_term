# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..', 'modules'))

from input_templating import lowercase_isotope_name
from math import pi, sqrt
import openmc
import openmc.deplete
import numpy as np
import openmc.mgxs as mgxs
import pandas as pd

###############################################################################
#               Materials Parameters (Taken from MCRE-ENG-PRSNT-0029Rev1A)
###############################################################################

def fuel_density(T): # Temp in K
    return 4.2126E+03 - 1.0686E+00 * T

def create_model(element_position):
    T_fuel = 900 # K
    kg_per_m3_to_cm_per_cm3 = 0.001
    density = fuel_density(T_fuel) * kg_per_m3_to_cm_per_cm3

    mgo_th_fraction = 0.8
    mgo_th_density = 3.58 # g/cm3
    mgo_density = mgo_th_fraction * mgo_th_density

    ###############################################################################

    # Materials definitions

    # HEU chloride fuel salt
    heu = openmc.Material(name='HEUF')
    heu.temperature = T_fuel
    heu.set_density('g/cm3', density)
    components = { # Eutectic NaCl-UCl3 (0.67-0.33 mole%)
        'U': {'percent': 0.33*1/4,
            'enrichment': 93.2},
        'Cl': 0.67*1/2+0.33*1/3,
        'Na': 0.67,
    }
    heu.add_components(components)

    # MgO Reflector
    mgo = openmc.Material(name='reflector')
    mgo.temperature = T_fuel
    mgo.set_density('g/cm3', mgo_density )
    components = {
        'Mg': 0.5,
        'O': 0.5,
    }
    mgo.add_components(components)

    # Incoloy 625 (Material composition from https://en.wikipedia.org/wiki/Inconel_625, material from MCRE-ENG-PRSNT-0029Rev1A)
    inc = openmc.Material(name='vessel')
    inc.temperature = T_fuel
    components = {
        'Ni': 0.58,
        'Cr': 0.215,
        'Mo': 0.09,
        'Fe': 0.05,
        'Nb': 0.01775,
        'Ta': 0.01775,
        'Co': 0.01,
        'Mn': 0.005,
        'Si': 0.005,
        'Al': 0.004,
        'Ti': 0.004,
        'C': 0.001,
        'P': 0.00015,
        'S': 0.00015,
    }   
    inc.add_components(components, percent_type = 'wo')

    # Air Material (from https://www.pnnl.gov/main/publications/external/technical_reports/pnnl-15870rev1.pdf Page 18)
    air = openmc.Material(name='air')
    air.temperature = T_fuel
    air.set_density('g/cm3', 0.001205)
    components = {
        'C': 0.000150,
        'N': 0.784431,
        'O': 0.210748,
        'Ar': 0.004671,
    }
    air.add_components(components)

    # B4C neutron absorber material for control rods and poison rods
    # From from https://www.pnnl.gov/main/publications/external/technical_reports/pnnl-15870rev1.pdf Pages 51,52
    b4c = openmc.Material(name='b4c')
    b4c.temperature = T_fuel  #temperature in Kelvin
    b4c.set_density('g/cm3', 2.52)
    components = {
        'B': {
            'percent': 0.799981,
            'enrichment': 80, # specified in MCRE-ENG-PRSNT-0029Rev1A
            'enrichment_type': 'wo',
            'enrichment_target': 'B10',
        },
        'C': 0.200019,
    }
    b4c.add_components(components)

    # SS 316 Material (from https://www.pnnl.gov/main/publications/external/technical_reports/pnnl-15870rev1.pdf Page 288)
    ss316 = openmc.Material(name='ss316')
    ss316.temperature = T_fuel
    components = {
        'C': 0.001900,
        'Si': 0.010048,
        'P': 0.000413,
        'S': 0.000260,
        'Cr': 0.181986,
        'Mn': 0.010274,
        'Fe': 0.666811,
        'Ni': 0.113803,
        'Mo': 0.014504,
    }
    ss316.add_components(components)

    ###############################################################################
    #                   Add cross section tallying nuclides
    ###############################################################################

    # Read isotopes to generate cross sections for
    isotopes = list(pd.read_csv('../source_term/data/simulation_parameters/isotopes.csv').keys())
    isotopes = [ lowercase_isotope_name(isotope) for isotope in isotopes ]

    data_library = openmc.data.DataLibrary.from_xml()
    data_library_nuclides = [ entry['materials'][0] for entry in data_library ]
    nuclide_list = list(set(isotopes) & set(data_library_nuclides))

    eps = 1e-12
    tally_mat = heu
    for isotope in isotopes:
        if isotope in tally_mat.get_nuclides() or isotope not in data_library_nuclides:
            continue
            
        tally_mat.add_nuclide(isotope, eps)

    ###############################################################################

    # Instantiate a Materials collection and export to xml
    materials = openmc.Materials([heu,mgo,inc,air,b4c,ss316])
    materials.export_to_xml()

    ###############################################################################
    #                           Geometry Definitions
    ###############################################################################
    def sphere_joining_cylinders(l_1, l_2, h):
        """Calculates the radius and offset from the center of a sphere that joins two parallel pipes
        of specified diameters a specified distance from each other. In 2D, this is equivalent to finding a circle with two parallel chords 
        a given distance from each other and with prescribed lengths. Solves the equations
        
        (l_1/2)^2 + x_1^2 = r^2
        (l_2/2)^2 + (x_1 + h)^2 = r^2
        """

        x_1 = 1/(2*h)*((l_1/2)**2 - (l_2/2)**2 - h**2)
        r = sqrt((l_1/2)**2 + x_1**2)
        return x_1, r

    # -------------------------
    # Capsule cylinder geometry
    # -------------------------
    capsule_cyl_ir = 40          # chosen to preserve total fuel volume as in MCRE-ENG-PRSNT-0029Rev1A
    capsule_cyl_thickness = 0.5  # from MCRE-ENG-PRSNT-0029Rev1A
    capsule_sphere_thickness = 1 # from MCRE-ENG-PRSNT-0029Rev1A
    capsule_cyl_height = 91.5    # from MCRE-ENG-PRSNT-0029Rev1A
    capsule_tot_height = 151     # from MCRE-ENG-PRSNT-0029Rev1A
    capsule_cyl_or = capsule_cyl_ir + capsule_cyl_thickness
    capsule_cyl_to_piping_height = (capsule_tot_height - capsule_cyl_height)/2 
    inner_capsule_cyl = openmc.ZCylinder(r=capsule_cyl_ir)
    outer_capsule_cyl = openmc.ZCylinder(r=capsule_cyl_or)
    outer_capsule_cyl_sphere_thickness = openmc.ZCylinder(r=capsule_cyl_ir + capsule_sphere_thickness)
    capsule_upper_sphere_plane = openmc.ZPlane(z0=capsule_cyl_height/2)
    capsule_lower_sphere_plane = openmc.ZPlane(z0=-capsule_cyl_height/2)
    capsule_top = openmc.ZPlane(capsule_tot_height/2)
    outer_sphere_capsule_top = openmc.ZPlane(capsule_tot_height/2 + capsule_sphere_thickness)
    capsule_bot = openmc.ZPlane(-capsule_tot_height/2)
    outer_sphere_capsule_bot = openmc.ZPlane(-capsule_tot_height/2 - capsule_sphere_thickness)

    # ---------------
    # Piping geometry
    # ---------------
    # Pump column: 6" NPS SCH.160 (https://pipestd.com/schedule-160-pipe-6-inch-dn150-mm/)
    pump_col_piping_ir = 6.589
    pump_col_piping_thickness = 1.826
    pump_col_piping_or = pump_col_piping_ir + pump_col_piping_thickness

    # Loop: 4" NPS SCH.40 (https://pipestd.com/schedule-40-pipe-4-inch-dn100-mm/)
    loop_piping_ir = 5.113
    loop_piping_thickness = 0.602
    loop_piping_or = loop_piping_ir + loop_piping_thickness

    # ---------------------------
    # Capsule top sphere geometry
    # ---------------------------
    # First get parameters for inner joining sphere
    x_1, r = sphere_joining_cylinders(capsule_cyl_ir*2, pump_col_piping_ir*2, capsule_cyl_to_piping_height)
    inner_capsule_top_sphere = openmc.Sphere(z0=capsule_cyl_height/2 - x_1, r=r)

    outer_capsule_top_sphere = openmc.Sphere(z0=capsule_cyl_height/2 - x_1, r=r + capsule_sphere_thickness)

    # ------------------------------
    # Capsule bottom sphere geometry
    # ------------------------------
    # First get parameters for inner joining sphere
    x_1, r = sphere_joining_cylinders(capsule_cyl_ir*2, loop_piping_ir*2, capsule_cyl_to_piping_height)
    inner_capsule_bot_sphere = openmc.Sphere(z0=-capsule_cyl_height/2 + x_1, r=r)

    outer_capsule_bot_sphere = openmc.Sphere(z0=-capsule_cyl_height/2 + x_1, r=r + capsule_sphere_thickness)

    # --------------------
    # Loop piping geometry
    # --------------------
    loop_height = 185.5 # centerline to centerline (from MCRE-ENG-PRSNT-0029Rev1A)
    loop_x_length = 175 # centerline to centerline (from MCRE-ENG-PRSNT-0029Rev1A)
    transverse_pipe_length = loop_x_length - 2*loop_piping_ir
    lower_descending_pipe_length = (loop_height - capsule_tot_height)/2 - loop_piping_ir # Last term is to account for elbow

    # Descending pipe
    inner_descending_pipe = openmc.ZCylinder(r=loop_piping_ir)
    outer_descending_pipe = openmc.ZCylinder(r=loop_piping_ir + loop_piping_thickness)

    # Lower transverse pipe
    z0 = -capsule_tot_height/2 -lower_descending_pipe_length -loop_piping_ir
    inner_lower_transverse_pipe = openmc.XCylinder(r=loop_piping_ir,z0=z0)
    outer_lower_transverse_pipe = openmc.XCylinder(r=loop_piping_ir + loop_piping_thickness,z0=z0)

    # Ascending pipe
    x0 = 2*loop_piping_ir + transverse_pipe_length
    inner_ascending_pipe = openmc.ZCylinder(r=loop_piping_ir,x0=x0)
    outer_ascending_pipe = openmc.ZCylinder(r=loop_piping_ir + loop_piping_thickness,x0=x0)

    # Upper transverse pipe
    z0 = capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir
    inner_upper_transverse_pipe = openmc.XCylinder(r=loop_piping_ir,z0=z0)
    outer_upper_transverse_pipe = openmc.XCylinder(r=loop_piping_ir + loop_piping_thickness,z0=z0)

    # Descending to lower transverse elbow
    major_radius = loop_piping_or
    z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
    x0 = loop_piping_ir+loop_piping_thickness
    inner_descending_to_lower_transverse_elbow = openmc.YTorus(a=major_radius,b=loop_piping_ir,c=loop_piping_ir,z0=z0,x0=x0)
    end_descending_pipe = openmc.ZPlane(z0=z0)
    begin_lower_transverse_pipe = openmc.XPlane(x0=x0)

    z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
    x0 = loop_piping_ir+loop_piping_thickness
    outer_descending_to_lower_transverse_elbow = openmc.YTorus(a=major_radius,b=loop_piping_or,c=loop_piping_or,z0=z0,x0=x0)

    # Lower transverse to ascending elbow 
    z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
    x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
    inner_lower_transverse_to_ascending_elbow = openmc.YTorus(a=major_radius,b=loop_piping_ir,c=loop_piping_ir,z0=z0,x0=x0)
    end_lower_transverse_pipe = openmc.XPlane(x0=x0)
    begin_ascending_pipe = openmc.ZPlane(z0=z0)

    major_radius = loop_piping_ir + loop_piping_thickness
    z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
    x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
    outer_lower_transverse_to_ascending_elbow = openmc.YTorus(a=major_radius,b=loop_piping_or,c=loop_piping_or,z0=z0,x0=x0)

    # Ascending to upper transverse elbow
    z0 = capsule_tot_height/2 + lower_descending_pipe_length - loop_piping_thickness
    x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
    inner_ascending_to_upper_transverse_elbow = openmc.YTorus(a=major_radius,b=loop_piping_ir,c=loop_piping_ir,z0=z0,x0=x0)
    begin_upper_transverse_pipe = openmc.XPlane(x0=x0)
    end_ascending_pipe = openmc.ZPlane(z0=z0)

    major_radius = loop_piping_ir + loop_piping_thickness
    z0 = capsule_tot_height/2 + lower_descending_pipe_length - loop_piping_thickness
    x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
    outer_ascending_to_upper_transverse_elbow = openmc.YTorus(a=major_radius,b=loop_piping_or,c=loop_piping_or,z0=z0,x0=x0)

    # --------------------
    # Pump column geometry
    # --------------------
    pump_column_pipe_length = lower_descending_pipe_length + loop_piping_ir - pump_col_piping_ir # Such that elbow centerline aligns with loop piping centerline

    # Pump column
    inner_pump_column_pipe = openmc.ZCylinder(r=pump_col_piping_ir)
    outer_pump_column_pipe = openmc.ZCylinder(r=pump_col_piping_ir + pump_col_piping_thickness)

    # Pump column elbow
    major_radius = pump_col_piping_or
    z0 = capsule_tot_height/2 + pump_column_pipe_length - pump_col_piping_thickness
    x0 = pump_col_piping_ir + pump_col_piping_thickness
    inner_column_to_connector_elbow = openmc.YTorus(a=major_radius,b=pump_col_piping_ir,c=pump_col_piping_ir,z0=z0,x0=x0)
    end_column = openmc.ZPlane(z0=z0)
    begin_connector = openmc.XPlane(x0=x0)

    outer_column_to_connector_elbow = openmc.YTorus(a=major_radius,b=pump_col_piping_or,c=pump_col_piping_or,z0=z0,x0=x0)

    # Pump column to loop piping conical connector
    conical_connector_length = 25 # Arbitrary, but change to match th model
    inner_cone_slope = ((pump_col_piping_ir - loop_piping_ir)/(conical_connector_length - pump_col_piping_thickness))**2
    inner_zero_point_x = pump_col_piping_ir/sqrt(inner_cone_slope) + pump_col_piping_ir + pump_col_piping_thickness
    inner_zero_point_z = capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir
    inner_column_elbow_to_piping_connector = openmc.model.XConeOneSided(r2=inner_cone_slope,up=False,x0=inner_zero_point_x,z0=inner_zero_point_z)

    outer_cone_slope = ((pump_col_piping_ir + pump_col_piping_thickness - (loop_piping_ir + loop_piping_thickness))/(conical_connector_length - pump_col_piping_thickness))**2
    outer_zero_point_x = (pump_col_piping_ir + pump_col_piping_thickness)/sqrt(outer_cone_slope) + pump_col_piping_ir + pump_col_piping_thickness
    outer_zero_point_z = capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir
    outer_column_elbow_to_piping_connector = openmc.model.XConeOneSided(r2=outer_cone_slope,up=False,x0=outer_zero_point_x,z0=outer_zero_point_z)
    end_upper_transverse_pipe = openmc.XPlane(x0=pump_col_piping_ir + conical_connector_length)

    # ------------------
    # Reflector geometry
    # ------------------
    reflector_height = 226            # taken from MCRE-ENG-PRSNT-0029Rev1A
    reflector_long_side_length = 165  # taken from MCRE-ENG-PRSNT-0029Rev1A
    reflector_short_side_length = 115 # taken from MCRE-ENG-PRSNT-0029Rev1A
    cavity_height = 200               # Assumed, not given in model description
    cavity_width = 100                # Assumed, not given in model description
    reflector = openmc.model.CruciformPrism(
        distances=[reflector_short_side_length/2,reflector_long_side_length/2],
        boundary_type='transmission',
    )
    reactor_cavity = openmc.model.RectangularParallelepiped(
        xmin=-cavity_width/2,xmax=cavity_width/2,
        ymin=-cavity_width/2,ymax=cavity_width/2,
        zmin=-cavity_height/2,zmax=cavity_height/2,
    )

    # ------------
    # KCS geometry
    # ------------
    element_ir = 2.851       # 2" NPS SCH.5 (from MCRE-ENG-PRSNT-0029Rev1A)
    element_or = 3.016       # 2" NPS SCH.5 (from MCRE-ENG-PRSNT-0029Rev1A)
    element_height = 95.3    # from MCRE-ENG-PRSNT-0029Rev1A
    gt_ir = 3.49             # 2.5" NPS SCH.5 (from MCRE-ENG-PRSNT-0029Rev1A)
    gt_or = 3.65             # 2.5" NPS SCH.5 (from MCRE-ENG-PRSNT-0029Rev1A)
    element_clearance = 2.54 # from MCRE-ENG-PRSNT-0029Rev1A
    element_rectangular_offset = (capsule_cyl_or + element_clearance + element_or)/sqrt(2)

    NE_element_inner = openmc.ZCylinder(r=element_ir,x0=element_rectangular_offset,y0=element_rectangular_offset)
    NE_element_outer = openmc.ZCylinder(r=element_or,x0=element_rectangular_offset,y0=element_rectangular_offset)
    NE_gt_inner = openmc.ZCylinder(r=gt_ir,x0=element_rectangular_offset,y0=element_rectangular_offset)
    NE_gt_outer = openmc.ZCylinder(r=gt_or,x0=element_rectangular_offset,y0=element_rectangular_offset)

    SE_element_inner = openmc.ZCylinder(r=element_ir,x0=element_rectangular_offset,y0=-element_rectangular_offset)
    SE_element_outer = openmc.ZCylinder(r=element_or,x0=element_rectangular_offset,y0=-element_rectangular_offset)
    SE_gt_inner = openmc.ZCylinder(r=gt_ir,x0=element_rectangular_offset,y0=-element_rectangular_offset)
    SE_gt_outer = openmc.ZCylinder(r=gt_or,x0=element_rectangular_offset,y0=-element_rectangular_offset)

    SW_element_inner = openmc.ZCylinder(r=element_ir,x0=-element_rectangular_offset,y0=-element_rectangular_offset)
    SW_element_outer = openmc.ZCylinder(r=element_or,x0=-element_rectangular_offset,y0=-element_rectangular_offset)
    SW_gt_inner = openmc.ZCylinder(r=gt_ir,x0=-element_rectangular_offset,y0=-element_rectangular_offset)
    SW_gt_outer = openmc.ZCylinder(r=gt_or,x0=-element_rectangular_offset,y0=-element_rectangular_offset)

    NW_element_inner = openmc.ZCylinder(r=element_ir,x0=-element_rectangular_offset,y0=element_rectangular_offset)
    NW_element_outer = openmc.ZCylinder(r=element_or,x0=-element_rectangular_offset,y0=element_rectangular_offset)
    NW_gt_inner = openmc.ZCylinder(r=gt_ir,x0=-element_rectangular_offset,y0=element_rectangular_offset)
    NW_gt_outer = openmc.ZCylinder(r=gt_or,x0=-element_rectangular_offset,y0=element_rectangular_offset)

    element_bot = openmc.ZPlane(z0=-element_height/2 + element_position)
    element_top = openmc.ZPlane(z0=element_height/2 + element_position)
    gt_bot = openmc.ZPlane(z0=-element_height/2)
    gt_top = openmc.ZPlane(z0=element_height/2)

    # Universe boundary
    eps = 1 # 1 cm of extra room to avoid issues with reflector surfaces that would otherwise be coincident with some of these boundaries
    xmin = -reflector_long_side_length/2
    xmax = reflector_long_side_length/2 + loop_x_length
    ymin = -reflector_long_side_length/2
    ymax = reflector_long_side_length/2
    zmin = -reflector_height/2
    zmax = max(reflector_height/2, element_height/2 + element_position) # To allow extra room for fully removed control element
    x1 = openmc.XPlane(xmin - eps, boundary_type='vacuum')
    x2 = openmc.XPlane(xmax + eps, boundary_type='vacuum')
    y1 = openmc.YPlane(ymin - eps, boundary_type='vacuum')
    y2 = openmc.YPlane(ymax + eps, boundary_type='vacuum')
    z1 = openmc.ZPlane(zmin - eps, boundary_type='vacuum')
    z2 = openmc.ZPlane(zmax + eps, boundary_type='vacuum')

    ###############################################################################
    #                           Region Definitions
    ###############################################################################

    # ------------
    # Loop regions
    # ------------
    inner_descending_pipe_region = -inner_descending_pipe & -capsule_bot & +end_descending_pipe
    inner_descending_to_lower_transverse_elbow_region = -inner_descending_to_lower_transverse_elbow & -end_descending_pipe & -begin_lower_transverse_pipe
    inner_lower_transverse_pipe_region = -inner_lower_transverse_pipe & +begin_lower_transverse_pipe & -end_lower_transverse_pipe
    inner_lower_transverse_to_ascending_elbow_region = -inner_lower_transverse_to_ascending_elbow & +end_lower_transverse_pipe & -begin_ascending_pipe
    inner_ascending_pipe_region = -inner_ascending_pipe & +begin_ascending_pipe & -end_ascending_pipe
    inner_ascending_to_upper_transverse_elbow_region = -inner_ascending_to_upper_transverse_elbow & +end_ascending_pipe & +begin_upper_transverse_pipe
    inner_upper_transverse_pipe_region = -inner_upper_transverse_pipe & -begin_upper_transverse_pipe & +end_upper_transverse_pipe
    inner_loop_region = inner_descending_pipe_region | \
                        inner_descending_to_lower_transverse_elbow_region | \
                        inner_lower_transverse_pipe_region | \
                        inner_lower_transverse_to_ascending_elbow_region | \
                        inner_ascending_pipe_region | \
                        inner_ascending_to_upper_transverse_elbow_region | \
                        inner_upper_transverse_pipe_region

    outer_descending_pipe_region = -outer_descending_pipe & +inner_descending_pipe & -capsule_bot & +end_descending_pipe
    outer_descending_to_lower_transverse_elbow_region = (-outer_descending_to_lower_transverse_elbow & +inner_descending_to_lower_transverse_elbow & -end_descending_pipe & -begin_lower_transverse_pipe)
    outer_lower_transverse_pipe_region = -outer_lower_transverse_pipe & +inner_lower_transverse_pipe & +begin_lower_transverse_pipe & -end_lower_transverse_pipe
    outer_lower_transverse_to_ascending_elbow_region = (-outer_lower_transverse_to_ascending_elbow & +inner_lower_transverse_to_ascending_elbow & +end_lower_transverse_pipe & -begin_ascending_pipe)
    outer_ascending_pipe_region = -outer_ascending_pipe & +inner_ascending_pipe & +begin_ascending_pipe & -end_ascending_pipe
    outer_ascending_to_upper_transverse_elbow_region = (-outer_ascending_to_upper_transverse_elbow & +inner_ascending_to_upper_transverse_elbow & +end_ascending_pipe & +begin_upper_transverse_pipe)
    outer_upper_transverse_pipe_region = -outer_upper_transverse_pipe & +inner_upper_transverse_pipe & -begin_upper_transverse_pipe & +end_upper_transverse_pipe
    outer_loop_region = outer_descending_pipe_region | \
                        outer_descending_to_lower_transverse_elbow_region | \
                        outer_lower_transverse_pipe_region | \
                        outer_lower_transverse_to_ascending_elbow_region | \
                        outer_ascending_pipe_region | \
                        outer_ascending_to_upper_transverse_elbow_region | \
                        outer_upper_transverse_pipe_region

    # ---------------------------------
    # Pump column and connector regions
    # ---------------------------------
    inner_pump_column_pipe_region = -inner_pump_column_pipe & +capsule_top & -end_column
    inner_column_to_connector_elbow_region = -inner_column_to_connector_elbow & +end_column & -begin_connector # & ~column_to_connector_elbow_corner_region
    inner_column_elbow_to_piping_connector_region = -inner_column_elbow_to_piping_connector & +begin_connector & -end_upper_transverse_pipe
    inner_pump_col_and_connector_region = inner_pump_column_pipe_region | \
                                        inner_column_to_connector_elbow_region | \
                                        inner_column_elbow_to_piping_connector_region

    outer_pump_column_pipe_region = -outer_pump_column_pipe & +inner_pump_column_pipe & +capsule_top & -end_column
    outer_column_to_connector_elbow_region = (-outer_column_to_connector_elbow & +inner_column_to_connector_elbow & +end_column & -begin_connector) # | column_to_connector_elbow_corner_region
    outer_column_elbow_to_piping_connector_region = -outer_column_elbow_to_piping_connector & +inner_column_elbow_to_piping_connector & +begin_connector & -end_upper_transverse_pipe
    outer_pump_col_and_connector_region = outer_pump_column_pipe_region | \
                                        outer_column_to_connector_elbow_region | \
                                        outer_column_elbow_to_piping_connector_region

    # ---------------
    # Capsule regions
    # ---------------
    inner_capsule_cyl_region = -inner_capsule_cyl & -capsule_upper_sphere_plane & +capsule_lower_sphere_plane
    inner_capsule_region = (inner_capsule_cyl_region | -inner_capsule_top_sphere | -inner_capsule_bot_sphere) & -capsule_top & +capsule_bot & -inner_capsule_cyl

    outer_capsule_cyl_region = -outer_capsule_cyl & +inner_capsule_cyl & -capsule_upper_sphere_plane & +capsule_lower_sphere_plane & -outer_capsule_cyl
    outer_capsule_top_sphere_region = -outer_capsule_top_sphere & +inner_capsule_top_sphere & -outer_capsule_cyl_sphere_thickness & +capsule_upper_sphere_plane
    outer_capsule_bot_sphere_region = -outer_capsule_bot_sphere & +inner_capsule_bot_sphere & -outer_capsule_cyl_sphere_thickness & -capsule_lower_sphere_plane
    outer_capsule_region = ((outer_capsule_cyl_region) | \
                            (outer_capsule_top_sphere_region) | \
                            (outer_capsule_bot_sphere_region)) & -outer_sphere_capsule_top & +outer_sphere_capsule_bot & ~inner_pump_column_pipe_region & ~inner_descending_pipe_region

    # ------------------
    # Guide tube/clad regions
    # ------------------
    NE_gt_region           = -NE_gt_outer      & +NE_gt_inner      & +gt_bot & -gt_top
    NE_element_clad_region = -NE_element_outer & +NE_element_inner & +element_bot & -element_top
    SE_gt_region           = -SE_gt_outer      & +SE_gt_inner      & +gt_bot & -gt_top
    SE_element_clad_region = -SE_element_outer & +SE_element_inner & +element_bot & -element_top
    SW_gt_region           = -SW_gt_outer      & +SW_gt_inner      & +gt_bot & -gt_top
    SW_element_clad_region = -SW_element_outer & +SW_element_inner & +element_bot & -element_top
    NW_gt_region           = -NW_gt_outer      & +NW_gt_inner      & +gt_bot & -gt_top
    NW_element_clad_region = -NW_element_outer & +NW_element_inner & +element_bot & -element_top

    # -----------------------
    # Control element regions
    # -----------------------
    NE_element_region = -NE_element_inner & +element_bot & -element_top
    SE_element_region = -SE_element_inner & +element_bot & -element_top
    SW_element_region = -SW_element_inner & +element_bot & -element_top
    NW_element_region = -NW_element_inner & +element_bot & -element_top

    # -----------------
    # Composite regions
    # -----------------
    fuel_region = inner_loop_region | \
                inner_capsule_region | \
                inner_pump_col_and_connector_region
    vessel_region = (outer_loop_region | \
                    outer_capsule_region | \
                    outer_pump_col_and_connector_region)
    reflector_region = -reflector & ~-reactor_cavity & +z1 & -z2 & ~fuel_region & ~vessel_region
    control_element_region = (NE_element_region | \
                            SE_element_region | \
                            SW_element_region | \
                            NW_element_region) & ~reflector_region & ~fuel_region & ~vessel_region
    gt_and_clad_region = (NE_element_clad_region | NE_gt_region | \
                        SE_element_clad_region | SE_gt_region | \
                        SW_element_clad_region | SW_gt_region | \
                        NW_element_clad_region | NW_gt_region) & ~reflector_region & ~fuel_region & ~vessel_region
    air_region = +x1 & -x2 & +y1 & -y2 & +z1 & -z2 & ~fuel_region & ~control_element_region & \
                ~gt_and_clad_region & ~vessel_region & ~reflector_region


    ###############################################################################
    #                           Cell Definitions
    ###############################################################################

    heu_cell = openmc.Cell(name='reactor_cell')
    heu_cell.fill = heu
    heu_cell.region = fuel_region

    vessel_cell = openmc.Cell(name='vessel_cell')
    vessel_cell.fill = inc
    vessel_cell.region = vessel_region

    air_cell = openmc.Cell(name='air_cell')
    air_cell.fill = air
    air_cell.region = air_region

    mgo_cell = openmc.Cell(name='reflector_cell')
    mgo_cell.fill = mgo
    mgo_cell.region = reflector_region

    ss316_cell = openmc.Cell(name='ss316_cell')
    ss316_cell.fill = ss316
    ss316_cell.region = gt_and_clad_region

    b4c_cell = openmc.Cell(name='b4c_cell')
    b4c_cell.fill = b4c
    b4c_cell.region = control_element_region

    root_universe = openmc.Universe()
    root_universe.add_cell(heu_cell)
    root_universe.add_cell(vessel_cell)
    root_universe.add_cell(air_cell)
    root_universe.add_cell(mgo_cell)
    root_universe.add_cell(ss316_cell)
    root_universe.add_cell(b4c_cell)

    geometry = openmc.Geometry(root_universe)

    geometry.export_to_xml()

    ###############################################################################
    #                           Simulation Parameters
    ###############################################################################

    batches =  500
    inactive = 50
    particles = 2000

    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles

    # Instantiate Source
    point = openmc.stats.Point((0, 0, 0))
    source = openmc.Source(space=point)
    settings_file.source = source


    settings_file.export_to_xml()

    ###############################################################################
    #                           Generate Plots
    ###############################################################################
    xy_plot = openmc.Plot(name='xy_plot')
    xy_plot.basis = 'xy'
    xy_plot.color_by = 'material'
    xy_plot.pixels = (int(5*(xmax - xmin)), int(5*(ymax - ymin)))
    xy_plot.width = (xmax - xmin, ymax - ymin)
    xy_plot.origin = ((xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2)
    xy_plot.colors = {
        heu: 'red',
        mgo: 'blue',
        inc: 'black',
        air: 'green',
      ss316: 'gray',
        b4c: 'cyan'
    }

    yz_plot = openmc.Plot(name='yz_plot')
    yz_plot.basis = 'yz'
    yz_plot.color_by = 'material'
    yz_plot.pixels = (int(5*(ymax - ymin)), int(5*(zmax - zmin)))
    yz_plot.width = (ymax - ymin, zmax - zmin)
    yz_plot.origin = ((xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2)
    yz_plot.colors = {
        heu: 'red',
        mgo: 'blue',
        inc: 'black',
        air: 'green',
      ss316: 'gray',
        b4c: 'cyan'
    }

    xz_plot = openmc.Plot(name='xz_plot')
    xz_plot.basis = 'xz'
    xz_plot.color_by = 'material'
    xz_plot.pixels = (int(5*(xmax - xmin)), int(5*(zmax - zmin)))
    xz_plot.origin = ((xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2)
    xz_plot.width = (xmax - xmin, zmax - zmin)
    xz_plot.colors = {
        heu: 'red',
        mgo: 'blue',
        inc: 'black',
        air: 'green',
      ss316: 'gray',
        b4c: 'cyan'
    }

    xz_plot_kcs = openmc.Plot(name='xz_plot_kcs')
    xz_plot_kcs.basis = 'xz'
    xz_plot_kcs.color_by = 'material'
    xz_plot_kcs.pixels = (int(5*(xmax - xmin)), int(5*(zmax - zmin)))
    xz_plot_kcs.origin = ((xmin + xmax)/2, element_rectangular_offset, (zmin + zmax)/2)
    xz_plot_kcs.width = (xmax - xmin, zmax - zmin)
    xz_plot_kcs.colors = {
        heu: 'red',
        mgo: 'blue',
        inc: 'black',
        air: 'green',
      ss316: 'gray',
        b4c: 'cyan'
    }

    plots = openmc.Plots()
    plots += [xy_plot, yz_plot, xz_plot, xz_plot_kcs]

    plots.export_to_xml()
    openmc.plot_geometry()

    ###############################################################################
    #                   Initialize cross section tallies
    ###############################################################################

    # Instantiate a 2-group EnergyGroups object
    group_edges = np.array([0.0, 20.0E+06])
    groups = mgxs.EnergyGroups(group_edges)

    # Instantiate a few different cross sections
    reaction_list = [
        'transport', 'absorption', 'capture', 
        '(n,a)', '(n,2a)', '(n,2n)', '(n,3n)', '(n,4n)',
        '(n,np)', '(n,p)', '(n,d)', '(n,t)', '(n,3He)',
        'scatter', 'nu-scatter matrix',
        'fission', 'nu-fission', 'kappa-fission',
        'chi', 'chi-prompt', 'chi-delayed', 'beta', 'decay-rate', 'inverse-velocity'
    ]
    library = mgxs.Library(geometry, name='mgxs', mgxs_types=reaction_list)
    library.correction = 'P0'
    library.energy_groups = groups
    library.num_delayed_groups = 6 # Must be set or .add_to_tallies_file() will fail
    library.domain_type = 'cell'
    library.by_nuclide = True
    library.nuclides = nuclide_list
    fuel = heu_cell
    library.domains = [fuel]
    tallies_file = openmc.Tallies()
    library.build_library()
    library.add_to_tallies_file(tallies_file, merge=True)

    # Fission rate tally
    fuel_filter = openmc.CellFilter(fuel.id)
    energy_filter = openmc.EnergyFilter(group_edges)
    fission_rate_tally = openmc.Tally(name='fission_rate')
    fission_rate_tally.filters = [fuel_filter]
    fission_rate_tally.scores = ['fission']
    tallies_file.append(fission_rate_tally)

    # Absorption rate tally
    abs_rate_tally = openmc.Tally(name='absorption_rate')
    abs_rate_tally.filters = [fuel_filter]
    abs_rate_tally.scores = ['absorption']
    tallies_file.append(abs_rate_tally)

    # Energy deposition tally for the fuel cell
    heating_tally = openmc.Tally(name='heating')
    heating_tally.filters = [fuel_filter]
    heating_tally.scores = ['fission-q-recoverable'] 
    tallies_file.append(heating_tally)

    # Create energy filter using SHEM-361 group structure
    shem_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES['SHEM-361'])

    tally_shem = openmc.Tally(name='flux')
    tally_shem.filters = [shem_filter]
    tally_shem.scores = ['flux']

    tallies_file.append(tally_shem)

    # Export to "tallies.xml"
    tallies_file.export_to_xml()

    ###############################################################################
    #                   Initialize and run depletion calculation
    ###############################################################################

    model = openmc.Model(geometry=geometry, settings=settings_file)
    return model

model = create_model(0)
model.run(threads=112)
# crit_height, guesses, keffs = openmc.search_for_keff(create_model, bracket=[0,121], tol=1e-2, print_iterations=True)
# print(f"Critical height: {crit_height} cm")

# # Create depletion "operator"
# # openmc.config['chain_file'] = 'chain_simple.xml'
# openmc.config['chain_file'] = '/projects/MCRE_studies/louime/data/depletion/chain_endfb80_sfr.xml'
# op = openmc.deplete.CoupledOperator(model)

# # Scale power to the fraction of the total fuel volume that was explicitly modeled
# integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='s')
# integrator.integrate()
