ff_time
=======

Python utility to estimate at what time flip-flopping lipids flip-flop.
To print the help below: ```python ff_times.py --help```

```
[ DESCRITPION ]

This script estimates when the flipflopping lipids identified by ff_detect actually
flip-flops. Three times are calculated for each flip-flopping lipid:
 - an estimated starting time of the flip-flop
 - an estimated ending time of the flip-flop
 - an estimated mid-point time
	
A file listing the flip-flopping lipids must be supplied with the --flipflops option.
Each line of this file should follow the format (time in ns):

 -> 'resname,resid,starting_leaflet,z_bead'

where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
'z_bead' particle is used to track the position of the lipid. The script will then
produce an updated file with the following format:

 -> 'resname,resid,starting_leaflet,z_bead,t_start,t_end'

The script also outputs the estimated time at which each flip-flop occurs (mid-point time)
and, based on these, the evolution with time of the number of lipids having flip-flopped.


[ REQUIREMENTS ]

The following python modules are needed:
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
   
    In very large systems (more than ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflets large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

3. 
   (a) The same approach as in ff_detect is used. Flip-flopping lipids are considered to have
   flip-flopped when more of their neighbours within the distance --ngh_dist belong to 
   the opposite leaflet than to their starting leaflet. The first time at which this happens
   is the mid-point time.
   
   (b) To detect when flip-flops start/end the distance of the flip-flopping lipids --bead
   particle to the center of geometry (COG) of the --ngh_nb nearest neighbours in each
   leaflet is calculated and a local core membrane layer is defined relative to this distance
   using a fraction of the distance, --ngh_frac,  between those inter-leaflets COGs: t_start
   corresponds to when the flip-flopping lipid enters this core layer (false starts, ie if it
   comes out again, are taken into account) and t_end corresponds to when it leaves it.
   
   NB: --ngh_dist and --ngh_nb are unrelated.
    

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
--ngh_dist	[15]	: distance (Angstrom) within which to consider neighbours for ff time, see note 3(a)
--ngh_nb	[5]	: nb of closest neighbours in each leaflet, see note 3(b)
--ngh_frac 	[0.1]	: fraction of distance between interleaflet neighbours COGs, see note 3(b)
--reverse		: browse xtc backwards [TO DO]
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
```
