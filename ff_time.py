#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog='ff_time', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/ff_time
**********************************************
	
[ DESCRITPION ]

This script estimates when the flipflopping lipids identited by ff_detect actually
flip-flops.
	
[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib

[ NOTES ]

1. It's a good idea to trjconv the xtc first and only outputs the relevant bead of each
   phospholipid as the script will run MUCH faster.	

2. Identification of the bilayer leaflets is controlled via 3 options:
   (a) selection of particles
    By default, the PO4 group is used to detect lipids and assess their flip-flop status.
    This bead name can be changed via the --bead flag, the selection string being:
    -> "name " + beadname
   
    Note that only lipids which contain the bead mentioned in the selection string
    will be taken into account to identify leaflets.
        
   (b) leaflet finding method: reference file
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use the default 15 Angstrom cutoff
    directly (without optimising):
     -> '--leaflets 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflets large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

3. The same approach as in ff_detect is used. Flip-flopping lipids are considered to have
   flip-flopped when more of their neighbours within the distance --neighbours belong to 
   the opposite leaflet than to their starting leaflet.
   Obviously this is a simple approach which does not take into account "false starts", i.e.
   flip-flopping lipids which almost flip-flop at a given time but don't actually do until
   quite later in the simulation.
   The trajectory can be browsed backwards using --reverse to check for consistency.


[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: reference structure [.gro] (required)
-x			: trajectory file [.xtc] (required, must be in current folder unless -o is used)
-o			: name of output folder
-b			: beginning time (ns)
-e			: ending time (ns)	
-t 		[10]	: process every t-frames

Leaflets identification (see note 2)
-----------------------------------------------------
--bead		[PO4]	: lipids bead name
--leaflets		: leaflet identification

Flip-flops identification
-----------------------------------------------------
--flipflops		: input file with flipflopping lipids (output of ff_detect)
--neighbours	[15]	: flip-flops detection method, see note 3
--reverse		: browse xtc backwards
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)

#leaflets identification
parser.add_argument('--bead', nargs=1, dest='beadname', default=['PO4'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#flip-flops identification
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', help=argparse.SUPPRESS, required=True)
parser.add_argument('--neighbours', nargs=1, dest='neighbours', default=[15], help=argparse.SUPPRESS)
parser.add_argument('--reverse', dest='reverse', action='store_true', help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start=args.t_start[0]
args.t_end=args.t_end[0]
args.frames_dt=args.frames_dt[0]
#leaflets identification
args.beadname = args.beadname[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
#flip-flops identification
args.selection_file_ff = args.selection_file_ff[0]
args.neighbours = float(args.neighbours[0])

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.selections
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)
if not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.neighbours < 0:
	print "Error: --neighbours should be greater than 0, see note 3"
	sys.exit(1)
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder=="no":
	args.output_folder="ff_times_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folders
	os.mkdir(args.output_folder)

	#create log
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/ff_time.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[ff_time v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python ff_time.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():

	global leaflet_sele_string

	#set default beads
	leaflet_sele_string = "name " + str(args.beadname)

	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global frames_to_write
	global nb_frames_to_process
	global f_start
	global f_end	
	f_start = 0
	
	print "\nLoading trajectory..."
	U = Universe(args.grofilename, args.xtcfilename)
	U_timestep = U.trajectory.dt
	all_atoms = U.selectAtoms("all")
	nb_atoms = all_atoms.numberOfAtoms()
	nb_frames_xtc = U.trajectory.numframes		
	#sanity check
	if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
		print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
		sys.exit(1)
	if U.trajectory.numframes < args.frames_dt:
		print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

	#rewind traj (very important to make sure that later the 1st frame of the xtc will be used for leaflet identification)
	U.trajectory.rewind()
	
	#create list of index of frames to process
	if args.t_end != -1:
		f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
		if f_end < 0:
			print "Error: the starting time specified is before the beginning of the xtc."
			sys.exit(1)
	else:
		f_end = nb_frames_xtc - 1		
	if args.t_start != -1:
		f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
		if f_start > f_end:
			print "Error: the starting time specified is after the end of the xtc."
			sys.exit(1)
	if (f_end - f_start)%args.frames_dt == 0:
		tmp_offset = 0
	else:
		tmp_offset = 1
	frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
	nb_frames_to_process = len(frames_to_process)
			
	return
def identify_ff():
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 4:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_leaflets():
	print "\nIdentifying leaflets..."
	
	#declare variables
	global U_lip
	global leaflet_sele
	global leaflet_sele_atoms
	global upper_resnums
	global lower_resnums
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	U_lip = U.selectAtoms(leaflet_sele_string)
	if U_lip.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"] = L.group(0)
			leaflet_sele["lower"] = L.group(1)
		else:
			leaflet_sele["upper"] = L.group(1)
			leaflet_sele["lower"] = L.group(0)
		leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog:
	else:
		leaflet_sele["both"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"].centerOfGeometry()[2]
		leaflet_sele["upper"] = leaflet_sele["both"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"] = leaflet_sele["both"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'

	#store resnums (avoid repetitive access)
	upper_resnums =leaflet_sele["upper"].resnums()
	lower_resnums =leaflet_sele["lower"].resnums()
		
	return

#=========================================================================================
# core functions
#=========================================================================================

def check_ff(f_nb, t):
	
	global ff_times
	global ff_nb_u2l
	global ff_nb_l2u
	
	#upper to lower
	if f_nb > 0:
		ff_nb_u2l[f_nb] = ff_nb_u2l[f_nb-1]
	for l_index in lipids_ff_u2l_index:
		if ff_times[l_index] == 0:
			tmp_neighbours = U_lip.selectAtoms("around " + str(args.neighbours) + " (resid " + str(lipids_ff_info[l_index][1]) + " and resname " + str(lipids_ff_info[l_index][0]) + ")").resnums()
			tmp_neighbours = np.in1d(tmp_neighbours, lower_resnums)
			if len(tmp_neighbours)> 0:
				tmp_ratio = len(tmp_neighbours[tmp_neighbours==True]) / len(tmp_neighbours)
				if tmp_ratio > 0.5:
					ff_times[l_index] = t
					ff_nb_u2l[f_nb] += 1
	
	#lower to upper
	if f_nb > 0:
		ff_nb_l2u[f_nb] = ff_nb_l2u[f_nb-1]
	for l_index in lipids_ff_l2u_index:
		if ff_times[l_index] == 0:
			tmp_neighbours = U_lip.selectAtoms("around " + str(args.neighbours) + " (resid " + str(lipids_ff_info[l_index][1]) + " and resname " + str(lipids_ff_info[l_index][0]) + ")").resnums()
			tmp_neighbours = np.in1d(tmp_neighbours, upper_resnums)
			if len(tmp_neighbours) > 0:
				tmp_ratio = len(tmp_neighbours[tmp_neighbours==True]) / len(tmp_neighbours)
				if tmp_ratio > 0.5:
					ff_times[l_index] = t
					ff_nb_l2u[f_nb] += 1

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg_times():
	
	print " -writing flip-flopping times..."	
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/ff_times.txt'
	output_txt = open(filename_txt, 'w')
	output_txt.write("[times of lipid flip-flops - written by ff_times v" + str(version_nb) + "]\n")
	#upper to lower
	output_txt.write("")
	output_txt.write("upper to lower\n")
	output_txt.write("--------------\n")
	for l_index in lipids_ff_u2l_index:
		output_txt.write(str(lipids_ff_info[l_index][0]) + "," + str(lipids_ff_info[l_index][1]) + "," + str(ff_times[l_index]) + "\n")
	
	#lower to upper
	output_txt.write("\n")
	output_txt.write("lower to upper\n")
	output_txt.write("--------------\n")
	for l_index in lipids_ff_l2u_index:
		output_txt.write(str(lipids_ff_info[l_index][0]) + "," + str(lipids_ff_info[l_index][1]) + "," + str(ff_times[l_index]) + "\n")
	output_txt.close()
	return

def write_xvg_evolution():
	
	print " -writing evolution of number of flip-flops..."	
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/ff_evolution.xvg'
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of number of flip-flops\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"number of flip-flops\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 2\n")
	output_xvg.write("@ s0 legend \"upper to lower\"\n")
	output_xvg.write("@ s1 legend \"lower to upper\"\n")
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index]) + "	" + str(ff_nb_u2l[f_index]) + "	" + str(ff_nb_l2u[f_index])
		output_xvg.write(results + "\n")
	output_xvg.close()

	return

################################################################################################################################################
# ALGORITHM
################################################################################################################################################

#load ata
set_lipids_beads()
load_MDA_universe()
identify_ff()
identify_leaflets()

#create data structures
global ff_times
global ff_nb_u2l
global ff_nb_l2u
global frames_time
ff_times = np.zeros(lipids_ff_nb)
ff_nb_u2l = np.zeros(nb_frames_to_process)
ff_nb_l2u = np.zeros(nb_frames_to_process)
frames_time = np.zeros(nb_frames_to_process)

#browse trajectory
print "\nChecking for flip-flopping status..."
for f_index in range(0,nb_frames_to_process):
	ts = U.trajectory[frames_to_process[f_index]]
	progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
	sys.stdout.flush()
	sys.stdout.write(progress)
		
	frames_time[f_index] = ts.time/float(1000)
	check_ff(f_index, ts.time/float(1000))

#create outputs
print "\n\nWriting results files..."
write_xvg_times()
write_xvg_evolution()

#exit
print "\nFinished successfully!" "Check results in", str(args.output_folder)
sys.exit(0)
