# Sequence-Learning
Code for sequence learning experiments.  See Gavornik &amp; Bear, 2014. 

Use this code with a Plexon electrophysiology system, Plexon's MATLAB Offline SDK, the Psychtoolbox, the Retinotopic-Mapping repository, and Jeff Gavornik's StimulusSuite. The StimulusSuite is used to send event times through a DAQ to the Plexon system. SequenceAnalysis calls a function named MyReadall.m, which is identical to Plexon's readall.m function, except that MyReadall.m saves the converted .plx file to a .mat file with the same name.  See the readme from Retinotopic-Mapping repository for more information.  

In our lab, we implant a single 3mm microelectrode into binocular primary visual cortex, both the left and right hemisphere. In essence, you run the SequenceStim.m function and record the local field potential (LFP) from the two microelectrodes.  The Plexon system outputs a .plx file that the SequenceAnalysis.m function converts to .mat format, via Plexon's MATLAB Offline SDK.  If a mouse has a unique identifier (e.g. 26881), then the function SequenceStim.m can be called as SequenceStim(26881).  The function will search for that mouse's retinotopic map created by the Retinotopy.m and MapRetinotopy.m functions.  It will then use the RetinoMap26881.mat file to generate a stimulus of flashing circles that target the retinotopic location of the implanted LFP-recording electrode. Over the course of 3 days, the mouse will see the same stimulus of flashing circles, each stimulus with 4 elements in the same locations and with the same relative timing, repeated 200 times per day. The number of elements, timing, and size of the circles can be altered by changing the SequenceVars.mat file. 

#Steps:
1) Run Retinotopy(26881) to mouse #26881 while recording the LFP using Plexon system. See the Retinotopic-Mapping repository for details on that process.
#
4) Soon after, run SequenceStim(26881) while recording the LFP with the Plexon system.  This will automatically load the retinotopic map created the previous day and target the map's center of mass. Save the .plx file as SeqData20160706_26881.plx (or whatever date is applicable) to CloudStation. The SequenceStim.m code will automatically save the parameters of the stimulus to CloudStation in a file named SeqStim20160706_26881.mat .
#
5) Do this for 3 consecutive days.
#
6) You can run SequenceAnalysis(26881,20160706) for each of the days, SequenceAnalysis(26881,20160707), etc. or use something like DailySequences.m to automatically do that for you. At the end of the experiment, with all of the data on the MATLAB path, run SeqAnWrapper(26881) .  This will output a series of figures that show the VEPs in response to the full sequence and to each element of the sequence, and how those changed over the 3 days.
