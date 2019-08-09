**** very much a work in progress... ****

** Save each frame of the simulation as a pdb file.
If you're using VMD and get a huge file containing all of the frames use split_MDfile.py to split it into individual files.
Check the last one to make sure it's not blank - VMD likes to do that.

** If the frames need to be aligned on a specific protein and aren't already use my align_pdbs_on_body.py script at https://github.com/attamatti/movement_analysis. 

************* Then use lipid_thickness_v3.py *************

** edit the section of the script after:
#########  CAs to draw on the final map  ############

designate which chains are the MSP proteins
MSP_chains = ['F','G']

update the barreldraw variable with a dictionary of other things to draw in this format:
{name1:[chain[AAs]],name2:[chain2,[AAs2]]}

for example:
barreldraw = {  'barrel':['A',range(436,794)],
                'lateralgate':['A',[423,424,425,426,427,428,429,430,801,802,803,804,805,806,807,808,809,810]],
                'membrane_contact1':['A',[310,311,312]],
                'membrane_contact2':['B',[215,216,345,350]],
                }

Finally specify the order to draw the additional elements in
EX:
draworder = ['barrel', 'membrane_contact1', 'lateralgate', 'membrane_contact2']

- I will update this later to read this info from a file like a civilized piece of software should, someday...

** Then run the script:
./lipid_thickness_v3.py frames*.mrc
  - or whatever search string that will find all of your MD frames...

** If you want to compare multiple simulations on the same scale the script replot_thickness_data.py will read the data written by lipid_thickness_v3.py and replot it with a designated scale.  Use it with:

replot_thickness_data.py <thickness min> <thickness max> <thickness increment> <std min> <std max> <std increment>

if run without any arguments it will return the stats of that particular simulation, so you can look at each and then pick the appropriate scale to use to compare them.

MGI
