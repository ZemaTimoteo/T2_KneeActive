% Compile .c files from RF Tools Toolbox

cd([toolbox_Path filesep 'rf_tools' filesep 'mex5'])
mex abrx.c
cd(DIR_sequenceOptimization)