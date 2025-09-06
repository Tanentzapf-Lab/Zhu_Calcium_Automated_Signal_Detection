# Zhu_Calcium_Automated_Signal_Detection
Code associated with automated signal detection of calcium signals

Run all python code from a main directory containing an 'Analysis' folder with subfolders named: 'CalciumSignal', 'Mask', 'TIF', 'TIF_avrg5frames', 'analyzeddata'

1) Photobleach correction
Open red and green channels separately in Fiji
For each channel, [Image > Adjust > Bleach Correction > Exponential Fit]
For green channel, save resulting .tif in 'CalciumSignal' with suffix '_greenphotobleachexpcorrected' 
For red channel, save resulting .tif in 'TIF' with suffix '_redphotobleachexpcorrected' 

2) Average every 5 frames of red channel
Open 'tif_avrg_every5frames_runfromcurrentfolder.ijm' by dragging into Fiji, click 'Run'
Select the 'Analysis' folder 
This will automatically go through every .tif stack in 'TIF' and create new files in 'TIF_avrg5frames' with the suffix '_AVG' appended

3) Create cell label mask in TrackMate 
Open 'redphotobleachexpcorrected_AVG' file in Fiji
[Plugins > Tracking > Trackmate] Select LoG detector: Estimated size 3.5 micron, Quality threshold 15
Set Quality score 15 for spots
Then use Overlap Tracker (Precise, Min IoU 0.3, scale factor 1)
At this point it is possible to manually edit tracks with TrackScheme. Undesired cells such as those belonging to secondary lymph glands can be deleted, and missed spots or connections can be added manually.
Export Label Image > Export only spots in tracks
Save label image in 'Mask' as .tif, prefix 'LblImg_' will be automatically added

4) Determining background fluorescence value - this step is necessary when changing imaging settings, and should be run once for every new set of imaging conditions
An imaging run should be taken with a sample that does not have green fluorophore.
Run this file in 'backgrounddetermination.py' (name file in first section of code)

5) Analyze experimental data
Open 'calciumsignaling_exportbasalsignal.py' and input filename and background fluorescence value (as determined from stop 4) in first section of code

6) To repeat exact graphing and analysis from publication, run 'downstreamgraphing.py'
