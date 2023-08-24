T2_KneeAtive - Reconstruction Pipeline

-----------------

Get T2 maps (for Phantom data, but works similarly for Knee Images) - Get T2 maps & values within each vial mask for T2me MSE data:
	1. Create a folder for subject and create subfolders
		1.Subfolder: 'imagResults_orig'
		2.Subfolder: 'imagResults_preproc'
		3.Subfolder: 'results'
	2. Copy the '.seq' file originally created when running the function 'test_MSE_SEQfile.py' and the 'sequence_info.mat' file
	3. Copy data (from scanner '.dat' file) to the respective folder of data - create a data structure with the following files:
	4. Get .mat file from raw data: Run python code - 'Github\Reconstruction\T2_EPG\T2_MSE_EPG_Phantom_all.py' OR 'Github\Reconstruction\T2_EPG\T2_MSE_EPG.py'
	5. Recon GRAPPA: Run matlab code - 'Github\Reconstruction\MSE_preproc_GRAPPA_recon.m'
	6. Get mask for the 14 vials or for knee
	7. Run python to generate maps for each vial - 'Github\Reconstruction\T2_EPG\T2_MSE_EPG_Phantom_vials.py' or 'Github\Reconstruction\T2_EPG\T2_MSE_EPG.py'
	8. Save '.mat' files with maps for each vial
	9. Copy all maps to a single folder of T2me MSE
	10. Do subsequente analysis on the maps (examples - ISMRM23)



