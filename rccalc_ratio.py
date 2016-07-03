#
#   ribbed plate python automated FEM
#
from solid import *
from scad import *
import math

import subprocess
import triangle.plot
import matplotlib.pyplot as plt

import numpy.ma as ma
import numpy as np
import sys
import tvtk
import os

import re
import shutil
import time
import shlex

from triangle import *
from scipy.interpolate import griddata
from threading import Timer
from easyprocess import EasyProcess
from getfem import *
#
#    this is the library of modules
#
import rccalc_lib
#       
#      Main program
#
number_of_ribs = 32
number_of_ribs = 16
#number_of_ribs = 8
plate_radius =    17.0    # mm
plate_thickness = 0.7     # mm
rib_width = 0.7           # mm
rib_width = 0.95           # mm
rib_thickness = 1.4     # mm

rib_length = 0.585*plate_radius
#rib_length = 0.9*plate_radius
rib_radius = math.sqrt(1.0/5.0)*plate_radius
#rib_radius = 1.2*rib_radius
rib_boss_length = 0.65*rib_length
#rib_boss_length = 0.1*rib_length

rib_boss_width = 1.0*rib_width
flare_len = plate_radius - rib_length/2.0-rib_radius
flare_center_r = plate_radius - flare_len/2.0
flare_angle = 360.0/number_of_ribs/2.0*3.0

solfolder = "solution_folder"
fileprefix = "rplate"

run_aborted = False
#
#   make the plate here
#
rib_plus_boss = []
plate = cylinder(r=plate_radius, h=plate_thickness, center=True, segments=60)
rib = cube([rib_length,rib_width,rib_thickness],center=True)
rib_boss = cube([rib_boss_length,rib_boss_width,rib_thickness],center=True)
rib_boss2 = cube([rib_boss_length/2.0,rib_boss_width*2.0,rib_thickness],center=True)
flare = cube([flare_len,rib_width,rib_thickness])
flare = Translate(y=-rib_width/2.0, z =-rib_thickness/2.0) (flare)
flare_up = Rotate(z=flare_angle) (flare)
flare_down = Rotate(z=-flare_angle) (flare)
flare_up = Translate(x=flare_center_r-rib_radius-flare_len/2.0) ( flare_up )
flare_down = Translate(x=flare_center_r-rib_radius-flare_len/2.0) ( flare_down )
#rib_boss2 = Translate(x=rib_length/6.0) (rib_boss2)
rib_plus_boss.append(rib)
#rib_plus_boss.append(rib_boss)
rib_plus_boss.append(rib_boss2)
#rib_plus_boss.append(flare_up)
#rib_plus_boss.append(flare_down)

rib_plus_boss = Union()(*rib_plus_boss)
rib_plus_boss = Translate(x=(rib_radius), z=(rib_thickness+plate_thickness)/2.0-0.05) (rib_plus_boss)
ribbed_plate = []
ribbed_plate.append(plate)
for i in range(0,number_of_ribs):
   _rib = Rotate(z=(360.0 / number_of_ribs) * i) (rib_plus_boss)
   ribbed_plate.append(_rib)
                 
total_shape = Union()(*ribbed_plate)

#
#   now make a folder 
#
#   clean it of it exists
#
if os.path.exists(solfolder):
    shutil.rmtree(solfolder, ignore_errors=True)
    
# now make it again
    
if not os.path.exists(solfolder):
    os.mkdir(solfolder)

#oscad_file = solfolder + "\\" + fileprefix + ".scad"
#stl_file = solfolder + "\\" + fileprefix+".stl"

oscad_file = fileprefix + ".scad"
stl_file = fileprefix+".stl"

sys.stdout.write('printing to oscad file ... ')
total_shape.render(oscad_file)
#  add the cylinder segment number option
add_segment_number_to_oscad(oscad_file,32)
print'done\n'

## if file exists, delete it ##
if os.path.isfile(stl_file):
        os.remove(stl_file)

sys.stdout.write('printing stl file ... ')
sys.stdout.flush()

# openscad to STL conversion

oscad_to_stl_cmd = "openscad.exe -o " + stl_file+ " " + oscad_file

status = subprocess.call(oscad_to_stl_cmd,shell=True)
if (status != 0 ):
    aborted_run = True
    print("oscad to stl failed !")
    sys.exit("Stopping here")
print'done\n'

tetgen_initial_time_tick = time.time()
print "starting tetgen at ", time.asctime(time.localtime(tetgen_initial_time_tick))
sys.stdout.write('meshing with tetgen ... ')
sys.stdout.flush()
#stl_to_mesh_tetgen_meshing_cmd = "tetgen.exe -pQgqa2.0 " + stl_file 
stl_to_mesh_tetgen_meshing_cmd = "tetgen.exe -pgqa2.0 " + stl_file 

pr = EasyProcess(stl_to_mesh_tetgen_meshing_cmd).call(timeout=10)
returncode = pr.return_code
stdoutdata =  pr.stdout
stderrdata =  pr.stderr

tetgen_end_time_tick = time.time()

status = returncode
print stdoutdata

if (status != 0 ):
    aborted_run = True
#  
    print("stl to tetgen mesh failed !")
    print
    sys.exit("Stopping here")

if (status != 0 ):
    aborted_run = True
    print("stl to tetgen mesh timeout !")
    print
    sys.exit("Stopping here")

print'done\n'
sys.stdout.flush()
print "ending tetgen at ", time.asctime(time.localtime(tetgen_end_time_tick))

# converting .mesh to .msh file 

fileposttet = fileprefix + ".1.mesh"
sys.stdout.write('converting mesh to msh file with gmsh ... ')
sys.stdout.flush()
mesh_to_msh_cmd = "gmsh " + fileposttet + " -0 -o " + "mesh3d.msh" 
p = subprocess.call(mesh_to_msh_cmd,shell=True)
print'done\n'

# here strip triangles
sys.stdout.write('stripping triangles from mesh ... ')
#strip_2D_elements_from_mesh("mesh3dproc.msh","meshsolid.msh")
strip_2D_elements_from_mesh("mesh3d.msh","meshsolid.msh")
print'done\n'

# here convert back to .mesh format using gmsh
filepoststrip = "meshsolid.msh"
meshfile      = "meshsolid.mesh"

sys.stdout.write('converting msh to mesh file with gmsh ... ')
sys.stdout.flush()
msh_to_mesh_cmd = "gmsh " + filepoststrip + " -0 -o " + "meshsolid.mesh" 
p = subprocess.call(mesh_to_msh_cmd,shell=True)
print'done\n'

# here convert mesh to abaqus to .inp format using gmsh
filepoststrip = "meshsolid.msh"
abaqusmeshfile      = "meshsolid.inp"

sys.stdout.write('converting msh to inp file with gmsh ... ')
sys.stdout.flush()
# set tolerance to 1e-4 mm so it prints a float which calculix needs
#msh_to_inp_cmd = "gmsh " + filepoststrip + " -0 -tol 0.00001 -o " + "meshsolid_t.inp" 
msh_to_inp_cmd = "gmsh " + filepoststrip + " -0 -o " + "meshsolid_t.inp" 
p = subprocess.call(msh_to_inp_cmd,shell=True)
print'done\n'

sys.stdout.write('cleaning inp file ... ')
sys.stdout.flush()
clean_inp_file("meshsolid_t.inp","meshsolid.inp")
print'done\n'

##################################################################
#
#   here read mesh with getfem interface to find boundary pts
#
##################################################################

sys.stdout.write('importing mesh ... ')
m=Mesh('import','gmsh','meshsolid.msh')
print 'done!'

# first collect the mesh points 

P=m.pts()
num_el = m.nbcvs()  # this is the number of tetrahedra
print "mesh has ", num_el, "tetrahedra"

# find the centroid coordinates for all of the mesh points

centroids = find_solid_centroids(m)
#
#       boundary selection
#
# P[2] contains the z coordinate of the points
#
# anything z >= plate_thickness/2 belongs to top
ctop=((P[2,:] - plate_thickness/2.0) > -1.0e-5*plate_thickness); 
# anything at z=-plate_thickness/2 is part of bottom
cbot=(abs(P[2,:] + plate_thickness/2.0) < 1.0e-5*plate_thickness); 
# anything at x^2+y^2 >=r^2 is part of side
# all points from the faces must be recognized 
# hence it must be on a band 
R = (P[0,:]*P[0,:] +P[1,:]*P[1,:])**0.5
#cside=(abs(R-plate_radius) < 0.015*plate_radius); 
cside=(abs(R[:]-plate_radius) < 0.05*plate_radius); 
# 
# now find bottom faces centroids
#
border = m.outer_faces()
#
pidtop=compress(ctop, range(0, m.nbpts()))
pidbot=compress(cbot, range(0, m.nbpts()))
pidside=compress(cside, range(0, m.nbpts()))
#
fside=m.faces_from_pid(pidside)

ftop=m.faces_from_pid(pidtop)
fbot=m.faces_from_pid(pidbot)
#
fnor = m.normal_of_faces(fside)
fnor1 = m.normal_of_faces(fbot)
#
# find the bottom facet centroid coordinates for all of the solids
#  return(ct,ta,fac_pts)
#
bottom_facet_centroids,area_facets,bot_el_pts,pcp = find_solid_bottom_facet_centroids1(m,fbot,-plate_thickness/2.0,1.0e-5)
#
#   Here identify and refine the edge BC 
#   to make sure that thet are outer faces.
#
fside2 = []
fside1 = fside.tolist()
borderlist = border.tolist()
#
#           correct the edge boundary
#
for index in range(0,len(fside1[0])):
    if (abs(fnor[2,index]) < 0.1 and 
        inlist(fside1[0],borderlist[0],index)==True):
        str1 =[ fside[0,index],fside[1,index] ]
        fside2.append(str1)
 
fside3 = array(fside2)
fside4 = fside3.transpose()
#
#  here are the point IDs for the side boundary
#
side_pts_ids = m.pid_in_faces(fside4)
#
fnor2 =  m.normal_of_faces(fside4)
fnor4 =  m.normal_of_faces(fbot)
#
#     Set the boundaries and multiple forces
#  
#  now get the points for each set 
#
fix_side_point_ids = m.pid_in_faces(fside4)
#
#  here we get the bottom elements concentrated forces
#  first we get the points from every bottom face element
#  and assign the corresponding 1/3 force 
#  from the centroid pressure force and element area 
#  and keep adding force values to account for pressure 
#  of adjacent elements
#
bot_pts = m.pid_in_faces(fbot)
##################################################################
#
#     first find the coma deflection
#
commonfactor = 1.0
magfactor2 = 1.0
#
weight1 = 0.0  # actuator force in gr
#
# piston pressure (for weight1 gr) in mN/mm^2
ppiston = -weight1*commonfactor*1.0e-3*9.8/(math.pi*plate_radius*plate_radius)*1000.0   
#
#  calculate the density pressure due to radius of glycerol
#  need pressure gradient parameter p/r (see Timoshenko,
#  Theory of plates and shells, pg. 285)
#
rhog = 1260.0        # glycerol density in kg/m^3
pres_gly = -rhog*9.8*plate_radius*1.0e-3*1.0e-6*1000.0  # in mN/mm^2
pgmax   =   commonfactor*magfactor2*pres_gly #    maximum pressure difference to center
a0 = pgmax/plate_radius  # coefficient for pressure gradient

print "coma loading ..."
print "ppiston =",ppiston,"mN/mm^2, a0*rad =",a0*plate_radius, "mN/mm^2"
# 
#      find the pressure distribution and concentrated loads
#
press_centroid,pmsh = pressure_at_centroids(bottom_facet_centroids,ppiston,a0)
force_centroid = force_at_centroids(press_centroid,bottom_facet_centroids,area_facets)
plot_pressure(pmsh,20)

bottom_pts_ids,bottom_pts_cload = find_solid_bottom_facet_cloads(m,fbot,force_centroid)
cmsh = force_at_bottom_points(m,bottom_pts_ids,bottom_pts_cload)
plot_cload(cmsh,20)
# 
#      now assemble the calculix input deck and files 
#
matname = "PDMS"
E=1.0e3  # mN/mm^2 PDMS Young modulus (Sylgard 184)
Nu=0.3

calculix_input_deck = "calculix_coma_run.inp"
sys.stdout.write('assembling coma input deck ... ')
sys.stdout.flush()
assemble_input_deck(num_el,matname,E,Nu,calculix_input_deck,abaqusmeshfile,side_pts_ids,bottom_pts_ids,bottom_pts_cload)
print'done\n'
#
# here run the solver
#
calculix_jobname = "calculix_coma_run"
sys.stdout.write('running calculix ... ')
sys.stdout.flush()
calculix_solve_cmd = "ccx -i " + calculix_jobname
#p = subprocess.call(calculix_solve_cmd,shell=True)
p = subprocess.Popen(calculix_solve_cmd,shell=True,stdout=subprocess.PIPE)
result = p.communicate()[0]
print'done\n'
print result
#
# find the extreme dz displacement for uniform pressure
#
dzmin_coma,dzmax_coma = calculix_extreme_z("calculix_coma_run.dat")
#
#     second find the uniform deflection
#
commonfactor = 1.0
magfactor2 = 0.0
#
weight1 = 25.0  # actuator force in gr
#
# piston pressure (for weight1 gr) in mN/mm^2
ppiston = -weight1*commonfactor*1.0e-3*9.8/(math.pi*plate_radius*plate_radius)*1000.0   
#
#  calculate the density pressure due to radius of glycerol
#  need pressure gradient parameter p/r (see Timoshenko,
#  Theory of plates and shells, pg. 285)
#
rhog = 1260.0        # glycerol density in kg/m^3
pres_gly = -rhog*9.8*plate_radius*1.0e-3*1.0e-6*1000.0  # in mN/mm^2
pgmax   =   commonfactor*magfactor2*pres_gly #    maximum pressure difference to center
a0 = pgmax/plate_radius  # coefficient for pressure gradient

print "uniform loading ..."
print "ppiston =",ppiston,"mN/mm^2, a0*rad =",a0*plate_radius, "mN/mm^2"
# 
#      find the pressure distribution and concentrated loads
#
press_centroid,pmsh = pressure_at_centroids(bottom_facet_centroids,ppiston,a0)
force_centroid = force_at_centroids(press_centroid,bottom_facet_centroids,area_facets)
plot_pressure(pmsh,20)
#
bottom_pts_ids,bottom_pts_cload = find_solid_bottom_facet_cloads(m,fbot,force_centroid)
cmsh = force_at_bottom_points(m,bottom_pts_ids,bottom_pts_cload)
plot_cload(cmsh,20)
# 
#      now assemble the calculix input deck and files 
#
matname = "PDMS"
E=1.0e3  # mN/mm^2 PDMS Young modulus (Sylgard 184)
Nu=0.3
#
calculix_input_deck = "calculix_uni_run.inp"
sys.stdout.write('assembling uniform pressure input deck ... ')
sys.stdout.flush()
assemble_input_deck(num_el,matname,E,Nu,calculix_input_deck,abaqusmeshfile,side_pts_ids,bottom_pts_ids,bottom_pts_cload)
print'done\n'
#
# here run the solver
#
calculix_jobname = "calculix_uni_run"
sys.stdout.write('running calculix ... ')
sys.stdout.flush()
calculix_solve_cmd = "ccx -i " + calculix_jobname
#
p = subprocess.Popen(calculix_solve_cmd,shell=True,stdout=subprocess.PIPE)
result = p.communicate()[0]
print'done\n'

print result

# find the extreme dz displacement for uniform pressure

dzmin_p,dzmax_p = calculix_extreme_z("calculix_uni_run.dat")

print "*********************************************************"
print "extreme coma dz = ", dzmin_coma,dzmax_coma
print "extreme uniform dz = ", dzmin_p,dzmax_p

avgc = (abs(dzmin_coma)+abs(dzmax_coma))/1.0
maxdzp = abs(dzmin_p)
print
ratio_coma = avgc/maxdzp
print "coma cont. ratio = ", ratio_coma
print "coma decrease factor = ",ratio_coma/0.0827*100.0/2.0, "%"
print "coma power = ",3.0*ratio_coma, "diopters"
print "*********************************************************"


sys.exit("Stopping here")

print 'You can view the tripod with (for example) mayavi:'
print 'mayavi -d ./rplate.vtk -f WarpVector -m BandedSurfaceMap'
print 'or'
print 'mayavi2 -d rplate.vtk -f WarpScalar -m Surface'
print 'or'
print 'gmsh rplate.pos'


