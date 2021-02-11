# MPRAdecoder
Script used to analyze raw MPRA data in the Laboratory of Cell Division, IMCB SB RAS, Novosibirsk, Russia.

Current version (2.0) of this script was designed by Anna Letiagina and Evgeniya Omelina from Alexey Pindyurin's research group. 
Version 1.0 of this script was developed by Anton Ivankin.

Dependencies  
To succesfully implement this script you need to create Spyder environment for python with conda https://docs.spyder-ide.org/current/installation.html and activate it:
conda activate spyder-env

You need to have the next packages in spyder environment:
bio 0.3.0,
biopython 1.78,
matplotlib 3.3.4,
numpy 1.20.0,
openpyxl 3.0.6,
pandas 1.2.1,
scipy 1.6.0,
seaborn 0.11.1.

To start analysis, you need to create directory containing files MPRAdecoder.py and parameters.py. 
File MPRAdecoder.py contains script for analysis of raw MPRA data and you should not change something there.
File parameters.py is configuration file and you can input your settings here. 

To test this script, download files from directories MPRAdecoder/example/input_data and MPRAdecoder/example/code. 




