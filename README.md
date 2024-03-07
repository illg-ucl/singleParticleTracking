# Single Particle Tracking code

Matlab code for the detection and tracking in two dimensions (2D) of bright spots in microscope images (on a dark background). The code can be used to analyse fluorescence microscopy videos of living cells with fluorescently labelled proteins that appear on the images as bright spots.

**Related publication**:  
*Single-molecule in vivo imaging of bacterial respiratory complexes indicates delocalized oxidative phosphorylation*, I. Llorente-Garcia, T. Lenn, H. Erhardt, O. L. Harriman, L.-N. Liu, A. Robson, S.-W. Chiu, S. Matthews, N. J. Willis, C. D. Bray, S.-H. Lee, J. Yen Shin, C. Bustamante, J. Liphardt, T. Friedrich, C. W. Mullineaux and M. C. Leake. Biochimica et Biophysica Acta (BBA): Bioenergetics 1837, 811-824 (2014). https://www.sciencedirect.com/science/article/pii/S0005272814000309?via%3Dihub. 

Example:

<p align="center">
  <img src="https://github.com/illg-ucl/singleParticleTracking/blob/master/single_particle_tracking_flu.png" width=55% height=55%>
</p>

# Copyright and License

Copyright (c) 2017. Isabel Llorente-Garcia, Dept. of Physics and Astronomy, University College London, United Kingdom.

Released and licensed under a BSD 2-Clause License:

https://github.com/illg-ucl/singleParticleTracking/blob/master/LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the BSD 2-Clause License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the BSD 2-Clause License for more details. You should have received a copy of the BSD 2-Clause License along with this program.

Citation: If you use this software for your data analysis please acknowledge it in your publications and cite as follows.

          -Citation example 1: 
           singleParticleTracking_illgLab software. (Version). 2017. Isabel Llorente-Garcia, 
           Dept. of Physics and Astronomy, University College London, United Kingdom.
           https://github.com/illg-ucl/singleParticleTracking/. (Download date).
           
          -Citation example 2:
           @Manual{... ,
            title  = {singleParticleTracking_illgLab software. (Version).},
            author       = {{Isabel Llorente-Garcia}},
            organization = { Dept. of Physics and Astronomy, University College London, United Kingdom.},
            address      = {Gower Place, London, UK.},
            year         = 2017,
            url          = {https://github.com/illg-ucl/singleParticleTracking/
            }

# Guide to run code and basic steps of analysis

* Help folder **"GuideToRunCode/"** contains guide **"GuideToRunCode/Guide to run Isabel's code.pdf"** with detailed instructions on how to run the code.
  
* Help folder **"scriptsThatHelpRunCode/"** contains various Matlab scripts that help run the code and do analysis of multiple image files, analysis of output excel files, etc.

* **Basic steps** needed to track bright spots on videos (more details in the above mentioned guide):
    - Set-up Matlab and correct data folders and current directory;
    - Set values of relevant parameters;
    - Find spot trajectories in image sequence and output them to an Excel file in the current directory using function **"FindTrajects.m"**.   
    - Link trajectory segments found into longer trajectories using **"linkTrajSegments.m"**.
    - Inspect tracks visually (on a video) and manually to validate good tracks by deciding which to accept as good using function **"goThroughTracksVideo.m"**. 
    - Analyse each track separatedly and produce one analysis Excel file and graph per track. This is based on functions **"showTrajAnalysis2.m"** and **"showManyTrajAnalysis2.m**. The **analysis includes**: trajectory number, mean x and y positions, number of points in track, track duration, spot intensity versus time with exponential fit, particle trajectory on x-y plane, first frame with trajectory overlayed, calculation of particle mean square displacement (MSD) versus delta time with error bars, fit of MSD curve (mobility and diffusion analysis), analysis of (bacteria) cell sizes and centre (useful when imaging multiple bacteria), etc.


# Matlab function folders

- **"openImageSequences/"**: functions to open different formats of image sequences (.sif, .dv, .tif or .mat data) and return image data in a useful form.

- **"operationsOnImages/"**: various useful auxiliary functions for operations on images, for example, for calculating frame averages, for circular masking, for getting bright spot masks, for getting cell masks and cell boundaries, for creating a moving average video, for calculating the gradient of an image, to calculate 2D Gaussian fits (unused), etc.

- **"saveDisplayOthers/"**: various useful auxiliary functions for saving or displaying images, image sequences (videos), trajectories, etc.

- **"spotDetection/"**: functions to detect bright spots on a single image frame. **"findCandidateSpots.m"** finds candidate bright spots on a single image frame (image of fluorescently labelled proteins in living cells). **"findSpotCentre1frame.m"** iteratively finds the centre of a bright spot on an image frame given initial candidate spot positions (x and y). **"eliminateCoincidentSpots.m"** eliminates coincident bright spot positions (closer than a chosen number of pixels). **"calculateSpotIntensity.m"** calculates the fluorescent spot intensity at a fixed position within a given frame using a square subarray and a circular mask centred on fixed position.

- **"Tracking/" - Find trajectories and ouput them to Excel file**: functions for tracking moving bright spots in an image sequence, finding their trajectories (x, y coordinates and time) and outputting them to an Excel files. 
The two functions used are **"FindTrajects.m"** and **"linkTrajSegments.m"**. The first one finds bright spots, finds their centres and joins spots in subsequent frames into trajectory segments. The second one joins segments into longer trajectories and outputs the resulting tracks to an Excel file. The **parameter files** used by these functions are within folder "Eg_ParameterScripts/" (see below) and need to be set to the correct values for the image sequences under analysis.
Details and explanations can be found in the pdf guide document "GuideToRunCode/Guide to run Isabel's code.pdf" within folder "GuideToRunCode/".

- **"Tracking/" - Inspect and validate tracks**: functions for visually inspecting and validating all bright-spot tracks found in an image sequence.
Function **"goThroughTracksVideo.m"** allows visual inspection and validation of the found tracks by showing a video of the found and accepted tracks (as circles overlayed on the original bright spots) so that you can label each track as good/valid (entering 1) or bad/invalid (entering 0). This process generates files such as "good_track_nums_1030_4_TIRF.mat"  in the current directory. This visual inspection is useful in case there are fluorescent objects (e.g. beads) on the image that you want to exclude, or tracks found outside the cell region, bad tracking or other anomalies. 
Details and explanations can be found in the pdf guide document "GuideToRunCode/Guide to run Isabel's code.pdf" within folder "GuideToRunCode/".

- **"Tracking/" - Analyse each good track**: functions to produce one analysis Excel file and one graph per analysed bright-spot track. A folder is generated which contains the Excel and graph files for the analysis of each track in the image sequence. The function used is **"showManyTrajAnalysis2.m"** (which calls **"showTrajAnalysis2.m"**). The **parameter file** used by this function is "paramsForShowTrajAnalysis2.m" and can be found within folder "Eg_ParameterScripts/" (see below). A copy of this file set to the desired parameter values needs to be placed inside the current Matlab directory for the code to run correctly.

- **"Tracking/" - Other**. There are various other functions inside directory "Tracking/". For example: **"saveTrackVideos.m"** and **"saveVideoManyTracks.m"** (plot and save video (.avi file) of a given image sequence with tracks overlaid on top), **"fullCellIntensity.m"** (get integrated fluorescence intensity over time for a full cell, see also guiding script "scriptsThatHelpRunCode/runToAnalyse_wholeCell_I_2.m"), **"analyseSetOfTracks2.m"** (open the corresponding Excel files from a set of tracks and use to back-fit initial intensity).

- **"Eg_ParameterScripts/"**: folder containing **parameter files**, i.e., lists of parameters needed to run different functions: **"paramsForFindTrajects.m"** (parameters for "FindTrajects.m"), **"paramsForLinkTrajSegments.m"** (parameters for "linkTrajSegments.m") and **"paramsForShowTrajAnalysis2.m"** (parameters for function "ShowTrajAnalysis2.m"). These parameter files need to be copied and placed into the analysis folder where the image sequences you want to analyse are saved ("current directory" for Matlab). The parameter values in these files should be modified to find the optimum ones for the image sequences you are looking at, saving the changes before running other functions.

- **"Diffusion-MSDcalculation/"**. Contains main two main functions: **"analyseTraj.m"** and **"getDisplacement.m"**. The functions analyse all trajectory data and calculate the mean square displacements (MSD) as a function of time-interval delta-t for each bright-spot trajectory (track) in a given image sequence (video).

- **"pair-wise-difference-calculation/"**: functions related to the calculation of pair-wise differences (PwD) for an input vector or series of values. Used for the calculation of the Mean Square Displacement (MSD).
   
- **"colocalisationAnalysis/"**. See guiding script "scriptsThatHelpRunCode/runToAnalyseColocalis.m". Analysis of whether tracks in the top half of the image sequence (channel for one of the fluorophores that labels one type of protein) colocalise with tracks in the bottom half of the images (channel for the different colour fluorophore that labels a different protein). Results and explanations of the colocalisation analysis can be found in this paper:
*Single-molecule in vivo imaging of bacterial respiratory complexes indicates delocalized oxidative phosphorylation*, I. Llorente-Garcia, T. Lenn, H. Erhardt, O. L. Harriman, L.-N. Liu, A. Robson, S.-W. Chiu, S. Matthews, N. J. Willis, C. D. Bray, S.-H. Lee, J. Yen Shin, C. Bustamante, J. Liphardt, T. Friedrich, C. W. Mullineaux and M. C. Leake. Biochimica et Biophysica Acta (BBA): Bioenergetics 1837, 811-824 (2014). https://www.sciencedirect.com/science/article/pii/S0005272814000309?via%3Dihub. 
  
- **"simulatedImages/"**: functions to generate synthetic (simulated) image sequences of bright spots (of known intensity and Gaussian width) on top of a noisy background and functions to validate the tracking code and the error it makes compared to truth. See also "scriptsThatHelpRunCode/ForSimulatedImages/". Used to quantify the mean uncertainty of the bright-spot tracking process for various signal to noise ratios.

- **"InVitroDataAnalysis/"**: functions to analyse fluorescent labels which are fixed (glued) onto a glass surface when imaged (this means tracking of moving bright-spots is not necessary). A frame average is obtained to find the positions of the fixed bright spots. The spots' fluorescence intensity level is measured over time. The aim is to extract the typical single-molecule intensity for a given fluorophore. See "runToAnalyseScript.m" inside this folder for the steps needed to complete the analysis. See also guiding script **"scriptsThatHelpRunCode/runToAnalyseInVitro.m"**.

- "FourierSpectrumCalculation/". Contains two functions: (1) "FourierAndFindPeaks.m", that calculates spectrum (Fourier transform) of histogram data (e.g., of pair-wise difference distribution) and find peaks in the power spectrum (to get the intensity step size); and (2) "trySpectrum.m", to try and identify the spacing between periodic peaks in the histogram of intensity pair wise differences (generates sum of periodic Gaussian peaks with noise and extract the spacing by calculating Fourier power spectrum).

- "Alex-codeMSD/". Folder with preliminary functions to extract Mean Squared Displacement (MSD) from bright-spot trajectories (tracks). Created by Alex Robson; modified and commented by Isabel Llorente Garcia. The functions in this folder are, in general, not used in the final code operations, but functions in folder "Diffusion-MSDcalculation" are based on it.

- "myChungKennedyFilter/": function by I. Llorente Garcia for Chung Kennedy filter. Mostly unused in the previous functions but there is an option to use it. 

- "FromMarkChungKennedyFilter/": original code from Mark Leake for Chung Kennedy Filter.
  
- "checksQuanResults/": ignore folder.

- "some_pdfs/": ignore folder (contains pdf files of a two Matlab functions).
