----------------  T2_KneeAtive - (T2 Knee Accessible Quantitative)  ------------------------



INTRODUCTION

Toolbox for Quantitative Knee Imaging - T2 mapping - based on Multi (Turbo) Spin-Echo Sequence.


T2 is a promising MRI biomarker, with applications including detecting early cartilage degeneration. Fitting a mono-exponential model to multi-spin-echo (MSE) measurements is typical for clinical T2 mapping. However, as the signal deviates from this model due to echoes caused by unwanted pathways, dictionary-based estimation approaches were proposed. A dictionary of echo modulation curves (EMC) is predicted using the Extended Phase Graphs (EPG) formalism, considering all acquisition parameters and a range of expected tissue parameters (T1 and T2); mapping consists of finding the best match to the measured signal. We proposed an MSE pulse sequence optimization for dictionary-based cartilage T2-mapping, following up on promising results in the knee and hip cartilages. 



GOAL

This toolbox is divided into two accessible tools: 

 - Sequence Optimization:
Code that helps determine the best combination of refocusing flip angle (FA), inter-echo spacing (TE), and total number of echoes (ETL - echo train length).

 - Image Reconstruction:
Helps reconstruct the T2 maps from Knee images - from the '.dat' file from the Siemens/GE scanners.


The present implementation aims to contribute towards the development of open-source tools for MRI in order to make MRI more accessible, particularly on quantitative T2 maps.



DISCLAIMER 

We use the following toolboxes (available in this github):

    GRAPPA_tutorial - by Mark Griswold and  Felix Breuer (breuer@mr-bavaria.de)

    UnRing_tool

    ...

    https://github.com/mriphysics/verse-mb

    https://github.com/imr-framework/pypulseq

    openjournals/joss-reviews#2478


** To run the Matlab codes, download and add the following toolboxes to the Tools folder:

    BART - https://mrirecon.github.io/bart/
    Gridding - https://github.com/ndwork/griddingCodes

** To run the Python codes create a Toolbox folder with the following toolboxes:

    code from pypulseq - https://github.com/imr-framework/pypulseq
    code from pulseqDiffusion - https://github.com/ritagnunes/PulseqDiffusion
