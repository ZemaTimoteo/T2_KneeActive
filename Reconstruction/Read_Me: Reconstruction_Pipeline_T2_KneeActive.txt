T2_KneeActive - Reconstruction Pipeline

-----------------

Get T2 maps (for Phantom data, but works similarly for Knee Images) - Get T2 maps & values within each vial mask/Knee Articular Cartilage mask for MSE '.dat' data:
	1. Create a folder for the subject and create subfolders
		1.Subfolder: 'imagResults_orig'
		2.Subfolder: 'imagResults_preproc'
		3.Subfolder: 'results'
	2. Copy the '.seq' file originally created when running the function 'test_MSE_cFA_SEQfile.py' and the 'sequence_info.mat' file
	3. Copy data (from scanner '.dat' file) to the respective folder of data - create a data structure with the following files:
	4. Get .mat file from raw data: Run python code - 'Github\Reconstruction\T2_EPG\T2_MSE_EPG_Phantom.py' OR 'Github\Reconstruction\T2_EPG\T2_MSE_EPG.py'
	5. Recon GRAPPA or LORAKS: Run matlab code - 'Github\Reconstruction\MSE_preproc_recon_Phantom.m' OR 'Github\Reconstruction\MSE_preproc_recon_Knee.m'
	6. Get mask for the 14 vials or for knee (in the examples, segmentation was already performed)
	7. Run python to generate maps for each vial - 'Github\Reconstruction\T2_EPG\T2_MSE_EPG_Phantom.py' or 'Github\Reconstruction\T2_EPG\T2_MSE_EPG.py'
		 (change directories within the functions for calling the respective resources, inside the following functions: 'template_match.py' - changing 'GitHub/
Reconstruction' by where you added the repository).
	8. The code automaticly saves '.mat' files with maps for each vial (Phantom) or for each slice (Knee)



