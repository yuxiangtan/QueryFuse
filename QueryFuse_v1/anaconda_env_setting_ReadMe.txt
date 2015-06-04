#First, you need to download and install the dependencies from channels.
#The command line is:
conda create -c dan_blanchard -c file://restricted/projectnb/montilab-p/conda_channel --file QF_v1_anaconda_requirements.txt -n QF_v1

#In "QF_v1_anaconda_requirements.txt" file, all the dependencies are listed.

#Second, before you run QueryFuse, you need to activate the environment if you did not.
#The command line is:
source activate QF_v1