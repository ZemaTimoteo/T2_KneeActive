Github - Correr MSE Dictionary Mapping

1 - Create a structure of folder with: 'dicom', 'segmentation', 'results'
2 - Copy files from scanner in DICOM for the 'Data' Folder
3 - Copy mask for the 'Segmentation' Folder & change mask of segmentation name to : 'Knee_Segment.nii.gz', if do not have, you can do it manually in the T2_MSE_EPG.py
4 - Open 'T2_MSE_EPG.py'
5 - Select Directories for dictionaries (if you don't have dicitonary files you will have to create them - so DictTest = 'True') and for RF pulse shape (where the subset of folders for all subjects will be).
6 - Set '0.5 - Settings' - see options
7 - Run function T2_MSE_EPG.py
8 - Select folder of subject to analyse
9 - The code, read the mask/ask you so select a mask, run the dictionaries, run template matching, depicts a figure with maps, save the maps in file '.mat' and in '.png'
