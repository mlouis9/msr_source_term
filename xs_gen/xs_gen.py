# Add python modules directory to path
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..', 'modules'))

from input_templating import lowercase_isotope_name
from math import pi, sqrt
import openmc
import openmc.deplete
import matplotlib.pyplot as plt
import numpy as np
import openmc.mgxs as mgxs
import pandas as pd
import yaml
import pdb
import matplotlib.pyplot as plt

###############################################################################
#               Materials Parameters (Taken from MCRE-ENG-PRSNT-0029Rev1A)
###############################################################################

def fuel_density(T): # Temp in K
    return 4.2126E+03 - 1.0686E+00 * T

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
    'B': 0.799981,
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

x1 = openmc.XPlane(-76.2, boundary_type='vacuum')
x2 = openmc.XPlane(76.2, boundary_type='vacuum')
y1 = openmc.YPlane(-76.2, boundary_type='vacuum')
y2 = openmc.YPlane(76.2, boundary_type='vacuum')
z1 = openmc.ZPlane(-96.52, boundary_type='vacuum')
z2 = openmc.ZPlane(96.52, boundary_type='vacuum')

# -------------------------
# Capsule cylinder geometry
# -------------------------
capsule_cyl_ir = 0.4
capsule_cyl_thickness = 0.005
capsule_sphere_thickness = 0.010
capsule_cyl_or = capsule_cyl_ir + capsule_cyl_thickness
capsule_cyl_height = 0.915
capsule_tot_height = 1.51
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
pump_col_piping_ir = 0.06589
pump_col_piping_thickness = 0.01826

# Loop: 4" NPS SCH.40 (https://pipestd.com/schedule-40-pipe-4-inch-dn100-mm/)
loop_piping_ir = 0.05113
loop_piping_thickness = 0.00602

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
loop_height = 1.855 # centerline to centerline
loop_x_length = 1.75 # centerline to centerline
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
elbow_r = loop_piping_ir + loop_piping_thickness/2
z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
x0 = loop_piping_ir+loop_piping_thickness
inner_descending_to_lower_transverse_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
inner_end_descending_pipe = openmc.ZPlane(z0=z0)
inner_begin_lower_transverse_pipe = openmc.XPlane(x0=x0)
inner_descending_to_lower_transverse_elbow_corner_region = +openmc.XPlane(x0=x0-loop_piping_thickness) & +openmc.ZPlane(z0=z0-loop_piping_thickness)

elbow_r = loop_piping_ir + loop_piping_thickness
z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
x0 = loop_piping_ir+loop_piping_thickness
outer_descending_to_lower_transverse_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
cone_z0 = z0 - loop_piping_ir - loop_piping_thickness
cone_x0 = 0
outer_end_descending_pipe = openmc.model.ZConeOneSided(r2=1,up=True,z0=cone_z0,x0=cone_x0)
outer_begin_lower_transverse_pipe = openmc.model.XConeOneSided(r2=1,up=True,z0=cone_z0,x0=cone_x0)

# Lower transverse to ascending elbow
elbow_r = loop_piping_ir + loop_piping_thickness/2
z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
inner_lower_transverse_to_ascending_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
inner_end_lower_transverse_pipe = openmc.XPlane(x0=x0)
inner_begin_ascending_pipe = openmc.ZPlane(z0=z0)
inner_lower_transverse_to_ascending_elbow_corner_region = -openmc.XPlane(x0=x0+loop_piping_thickness) & +openmc.ZPlane(z0=z0-loop_piping_thickness)

elbow_r = loop_piping_ir + loop_piping_thickness
z0 = -capsule_tot_height/2 - lower_descending_pipe_length + loop_piping_thickness
x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
outer_lower_transverse_to_ascending_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
cone_z0 = z0 - loop_piping_ir - loop_piping_thickness
cone_x0 = x0 + loop_piping_ir + loop_piping_thickness
outer_end_lower_transverse_pipe = openmc.model.XConeOneSided(r2=1,up=False,x0=cone_x0,z0=cone_z0)
outer_begin_ascending_pipe = openmc.model.ZConeOneSided(r2=1,up=True,x0=cone_x0,z0=cone_z0)

# Ascending to upper transverse elbow
elbow_r = loop_piping_ir + loop_piping_thickness/2
z0 = capsule_tot_height/2 + lower_descending_pipe_length - loop_piping_thickness
x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
inner_ascending_to_upper_transverse_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
inner_begin_upper_transverse_pipe = openmc.XPlane(x0=x0)
inner_end_ascending_pipe = openmc.ZPlane(z0=z0)
inner_ascending_to_upper_transverse_elbow_corner_region = -openmc.XPlane(x0=x0+loop_piping_thickness) & -openmc.ZPlane(z0=z0+loop_piping_thickness)

elbow_r = loop_piping_ir + loop_piping_thickness
z0 = capsule_tot_height/2 + lower_descending_pipe_length - loop_piping_thickness
x0 = loop_piping_ir + transverse_pipe_length - loop_piping_thickness
outer_ascending_to_upper_transverse_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
cone_z0 = z0 + loop_piping_ir + loop_piping_thickness
cone_x0 = x0 + loop_piping_ir + loop_piping_thickness
outer_begin_upper_transverse_pipe = openmc.model.XConeOneSided(r2=1,up=False,x0=cone_x0,z0=cone_z0)
outer_end_ascending_pipe = openmc.model.ZConeOneSided(r2=1,up=False,x0=cone_x0,z0=cone_z0)

# --------------------
# Pump column geometry
# --------------------
pump_column_pipe_length = lower_descending_pipe_length + loop_piping_ir - pump_col_piping_ir # Such that elbow centerline aligns with loop piping centerline

# Pump column
inner_pump_column_pipe = openmc.ZCylinder(r=pump_col_piping_ir)
outer_pump_column_pipe = openmc.ZCylinder(r=pump_col_piping_ir + pump_col_piping_thickness)

# Pump column elbow
elbow_r = pump_col_piping_ir + pump_col_piping_thickness/2
z0 = capsule_tot_height/2 + pump_column_pipe_length - pump_col_piping_thickness
x0 = pump_col_piping_ir + pump_col_piping_thickness
inner_column_to_connector_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
inner_end_column = openmc.ZPlane(z0=z0)
inner_begin_connector = openmc.XPlane(x0=x0)
inner_column_to_connector_elbow_corner_region = -openmc.ZPlane(z0=z0+pump_col_piping_thickness) & +openmc.XPlane(x0=x0-pump_col_piping_thickness)

elbow_r = pump_col_piping_ir + pump_col_piping_thickness
z0 = capsule_tot_height/2 + pump_column_pipe_length - pump_col_piping_thickness
x0 = pump_col_piping_ir + pump_col_piping_thickness
outer_column_to_connector_elbow = openmc.YTorus(a=elbow_r,b=elbow_r,c=elbow_r,z0=z0,x0=x0)
outer_end_column = openmc.ZPlane(z0=z0)
outer_begin_connector = openmc.XPlane(x0=x0)

# Pump column to loop piping conical connector
conical_connector_length = 0.25 # Arbitrary, but change to match th model
inner_cone_slope = ((pump_col_piping_ir - loop_piping_ir)/(conical_connector_length - pump_col_piping_thickness))**2
inner_zero_point_x = pump_col_piping_ir/sqrt(inner_cone_slope) + pump_col_piping_ir + pump_col_piping_thickness
inner_zero_point_z = capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir
inner_column_elbow_to_piping_connector = openmc.model.XConeOneSided(r2=inner_cone_slope,up=False,x0=inner_zero_point_x,z0=inner_zero_point_z)

outer_cone_slope = ((pump_col_piping_ir + pump_col_piping_thickness - (loop_piping_ir + loop_piping_thickness))/(conical_connector_length - pump_col_piping_thickness))**2
outer_zero_point_x = (pump_col_piping_ir + pump_col_piping_thickness)/sqrt(outer_cone_slope) + pump_col_piping_ir + pump_col_piping_thickness
outer_zero_point_z = capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir
outer_column_elbow_to_piping_connector = openmc.model.XConeOneSided(r2=outer_cone_slope,up=False,x0=outer_zero_point_x,z0=outer_zero_point_z)
end_upper_transverse_pipe = openmc.XPlane(x0=pump_col_piping_ir + conical_connector_length)

###############################################################################
#                           Region Definitions
###############################################################################

# ------------
# Loop regions
# ------------
inner_descending_pipe_region = -inner_descending_pipe & -capsule_bot & +inner_end_descending_pipe
inner_descending_to_lower_transverse_elbow_region = -inner_descending_to_lower_transverse_elbow & -inner_end_descending_pipe & -inner_begin_lower_transverse_pipe & ~inner_descending_to_lower_transverse_elbow_corner_region
inner_lower_transverse_pipe_region = -inner_lower_transverse_pipe & +inner_begin_lower_transverse_pipe & -inner_end_lower_transverse_pipe
inner_lower_transverse_to_ascending_elbow_region = -inner_lower_transverse_to_ascending_elbow & +inner_end_lower_transverse_pipe & -inner_begin_ascending_pipe & ~inner_lower_transverse_to_ascending_elbow_corner_region
inner_ascending_pipe_region = -inner_ascending_pipe & +inner_begin_ascending_pipe & -inner_end_ascending_pipe
inner_ascending_to_upper_transverse_elbow_region = -inner_ascending_to_upper_transverse_elbow & +inner_end_ascending_pipe & +inner_begin_upper_transverse_pipe & ~inner_ascending_to_upper_transverse_elbow_corner_region
inner_upper_transverse_pipe_region = -inner_upper_transverse_pipe & -inner_begin_upper_transverse_pipe & +end_upper_transverse_pipe
inner_loop_region = inner_descending_pipe_region | \
                    inner_descending_to_lower_transverse_elbow_region | \
                    inner_lower_transverse_pipe_region | \
                    inner_lower_transverse_to_ascending_elbow_region | \
                    inner_ascending_pipe_region | \
                    inner_ascending_to_upper_transverse_elbow_region | \
                    inner_upper_transverse_pipe_region

outer_descending_pipe_region = -outer_descending_pipe & +inner_descending_pipe & -capsule_bot & -outer_end_descending_pipe
outer_descending_to_lower_transverse_elbow_region = -outer_descending_to_lower_transverse_elbow & +inner_descending_to_lower_transverse_elbow & +outer_end_descending_pipe & +outer_begin_lower_transverse_pipe
outer_lower_transverse_pipe_region = -outer_lower_transverse_pipe & +inner_lower_transverse_pipe & -outer_begin_lower_transverse_pipe & -outer_end_lower_transverse_pipe
outer_lower_transverse_to_ascending_elbow_region = -outer_lower_transverse_to_ascending_elbow & +inner_lower_transverse_to_ascending_elbow & +outer_end_lower_transverse_pipe & +outer_begin_ascending_pipe
outer_ascending_pipe_region = -outer_ascending_pipe & +inner_ascending_pipe & -outer_begin_ascending_pipe & -outer_end_ascending_pipe
outer_ascending_to_upper_transverse_elbow_region = -outer_ascending_to_upper_transverse_elbow & +inner_ascending_to_upper_transverse_elbow & +outer_end_ascending_pipe & +outer_begin_upper_transverse_pipe
outer_upper_transverse_pipe_region = -outer_upper_transverse_pipe & +inner_upper_transverse_pipe & -outer_begin_upper_transverse_pipe & +end_upper_transverse_pipe
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
inner_pump_column_pipe_region = -inner_pump_column_pipe & +capsule_top & -inner_end_column
inner_column_to_connector_elbow_region = -inner_column_to_connector_elbow & +inner_end_column & -inner_begin_connector & ~inner_column_to_connector_elbow_corner_region
inner_column_elbow_to_piping_connector_region = -inner_column_elbow_to_piping_connector & +inner_begin_connector & -end_upper_transverse_pipe
inner_pump_col_and_connector_region = inner_pump_column_pipe_region | \
                                      inner_column_to_connector_elbow_region | \
                                      inner_column_elbow_to_piping_connector_region

outer_pump_column_pipe_region = -outer_pump_column_pipe & +inner_pump_column_pipe & +capsule_top & -outer_end_column
corner_chunk = +inner_end_column & -inner_begin_connector & inner_column_to_connector_elbow_corner_region
outer_column_to_connector_elbow_region = -outer_column_to_connector_elbow & +inner_column_to_connector_elbow & +outer_end_column & -outer_begin_connector | corner_chunk
outer_column_elbow_to_piping_connector_region = -outer_column_elbow_to_piping_connector & +inner_column_elbow_to_piping_connector & +outer_begin_connector & -end_upper_transverse_pipe
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
                        (outer_capsule_bot_sphere_region)) & -outer_sphere_capsule_top & +outer_sphere_capsule_bot

# -----------------
# Composite regions
# -----------------
fuel_region = inner_loop_region | \
              inner_capsule_region | \
              inner_pump_col_and_connector_region
vessel_region = outer_loop_region | \
                outer_capsule_region | \
                outer_pump_col_and_connector_region
air_region = +x1 & -x2 & +y1 & -y2 & +z1 & -z2 & ~fuel_region & ~vessel_region

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

# test_cell = openmc.Cell(name='test')
# test_cell.fill = ss316
# test_cell.region = outer_lower_transverse_pipe_region

root_universe = openmc.Universe()
root_universe.add_cell(heu_cell)
root_universe.add_cell(vessel_cell)
root_universe.add_cell(air_cell)
# root_universe.add_cell(test_cell)

geometry = openmc.Geometry(root_universe)

geometry.export_to_xml()

# OpenMC simulation parameters

batches =  500
inactive = 50
particles = 20000

settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles

# Instantiate Source
point = openmc.stats.Point((0, 0, 0))
source = openmc.Source(space=point)
settings_file.source = source


settings_file.export_to_xml()

a = 22.9176

# Plot generation
plot1 = openmc.Plot()
plot1.basis = 'xy'
plot1.color_by = 'material'
plot1.pixels = (1000, 1000)
plot1.colors = {
    heu: 'red',
    # mgo: 'blue',
    inc: 'black',
    air: 'green',
  ss316: 'black',
    # b4c: 'gray'
}

plot2 = openmc.Plot()
plot2.basis = 'yz'
plot2.color_by = 'material'
plot2.pixels = (1000, 1000)
plot2.colors = {
    heu: 'red',
    # mgo: 'blue',
    inc: 'black',
    air: 'green',
  ss316: 'black',
    # b4c: 'gray'
}

plot3 = openmc.Plot()
plot3.basis = 'xz'
plot3.color_by = 'material'
plot3.pixels = (5000, 5000)
# plot3.width = (0.75, 0.75)
# plot3.origin = (pump_col_piping_ir, 0, capsule_tot_height/2 + lower_descending_pipe_length + loop_piping_ir)
plot3.colors = {
    heu: 'red',
    # mgo: 'blue',
    inc: 'black',
    air: 'green',
  ss316: 'black',
    # b4c: 'gray'
}

plots = openmc.Plots()
plots.append(plot1)
plots += [plot2,plot3]

# vox_plot = openmc.Plot()
# vox_plot.type = 'voxel'
# vox_plot.width = (100., 100., 200.)
# vox_plot.pixels = (400, 400,  800)
# plots.append(vox_plot)

plots.export_to_xml()
openmc.plot_geometry()

# openmc.voxel_to_vtk('plot_4.h5')

###############################################################################
#                   Initialize cross section tallies
###############################################################################

# Load params
with open('params.yaml', 'r') as f:
    params = yaml.safe_load(f)

# Instantiate a 2-group EnergyGroups object
group_edges = np.array(params['group_edges'])
groups = mgxs.EnergyGroups(group_edges)

# Instantiate a few different cross sections
namespace = {'params': params}
exec(params['reaction_list'], namespace)
reaction_list = [
    'transport', 'absorption', 'capture', 
    '(n,a)', '(n,2a)', '(n,2n)', '(n,3n)', '(n,4n)',
    '(n,np)', '(n,p)', '(n,d)', '(n,t)', '(n,3He)',
    'scatter', 'nu-scatter matrix',
    'fission', 'nu-fission', 'kappa-fission',
    'chi', 'chi-prompt', 'chi-delayed', 'beta', 'decay-rate', 'inverse-velocity'
]
library = mgxs.Library(geometry, name=params['library_name'], mgxs_types=reaction_list)
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
energies_shem = openmc.mgxs.GROUP_STRUCTURES['SHEM-361']
shem_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES['SHEM-361'])

tally_shem = openmc.Tally(name='flux')
tally_shem.filters = [shem_filter]
tally_shem.scores = ['flux']

tallies_file.append(tally_shem)

# Export to "tallies.xml"
tallies_file.export_to_xml()

# Export library object as a pickle, to be loaded later by read_output.py
library.dump_to_file()

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

model = openmc.Model(geometry=geometry, settings=settings_file)

# Create depletion "operator"
# openmc.config['chain_file'] = 'chain_simple.xml'
openmc.config['chain_file'] = '/projects/MCRE_studies/louime/data/depletion/chain_endfb80_sfr.xml'
op = openmc.deplete.CoupledOperator(model)

# Perform simulation using the predictor algorithm
time_steps = params['time_steps']  # seconds

# Scale power to the fraction of the total fuel volume that was explicitly modeled
tp_vol = 289116.4440803906
vol_fraction = heu.volume / tp_vol
power = params['power']  # W/cm
power = [ power_val * vol_fraction for power_val in power ]
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power, timestep_units='s')
integrator.integrate()
