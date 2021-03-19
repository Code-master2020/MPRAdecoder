# MPRAdecoder
Script used to analyze raw MPRA data in the Laboratory of Cell Division, IMCB SB RAS, Novosibirsk, Russia.

Current version (1.0) of this script was designed by Anna Letiagina, Evgeniya Omelina and Anton Ivankin from Alexey Pindyurin's research group. 

Dependencies  
To succesfully implement this script you need to create Spyder environment for python with conda.

Installing conda on linux: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Spyder installation guide: https://docs.spyder-ide.org/current/installation.html

To activate spyder environment:
conda activate spyder-env

You need to have the next packages in spyder environment:
Bio, scipy, pandas, matplotlib, seaborn, openpyxl. 

To install them please input in command line: 

pip install package, e.g.: pip install Bio


To start analysis, you need to create directory containing files MPRAdecoder.py and parameters.py. 
File MPRAdecoder.py contains script for analysis of raw MPRA data and you should not change something there.
File parameters.py is configuration file and you can input your settings here. 

To test this script, download files from directories MPRAdecoder/example/input_data and MPRAdecoder/example/code. 




