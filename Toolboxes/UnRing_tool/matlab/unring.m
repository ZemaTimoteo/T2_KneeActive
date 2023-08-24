%%%%%%%%
%
% unring - tool for removal of the Gibbs ringing artefact
% Usage: outvol = unring(invol,params)
% Options: invol - input volume 
%          params - 3x1 array with [minW maxW nsh]
%                     nsh discretization of subpixel spaceing (default 20)
%                     minW  left border of window used for TV computation (default 1)
%                     maxW  right border of window used for TV computation (default 3)


function v = unring(v,params)

if nargin == 1,
    params = [1 3 20];
end;

v = double(v);
cd('D:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\UnRing_tool\matlab')
mex ringRm.cpp -LD:\Tiago\Trabalho\2021_2025_PhD\Projects\qMRI_Joint\Code\matlabCode\Toolboxes\UnRing_tool\matlab\fftw-3.3.5-dll64\fftw3 -compatibleArrayDims
v = ringRm(v,params); % mex ringRm.cpp -lfftw3 (sually links to MATLAB's libfftw3 library)