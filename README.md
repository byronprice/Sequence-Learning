# Sequence-Learning
Code for sequence learning experiments.  See Gavornik &amp; Bear, 2014. 

Use this code with a Plexon electrophysiology system, Plexon's MATLAB Offline SDK, the Psychtoolbox, and the Retinotopic-Mapping repository. See the readme for that repository for first steps.  

In our lab, we implant a single 3mm microelectrode into binocular primary visual cortex, both the left and right hemisphere. In essence, you run the SequenceStim.m function and record the local field potential (LFP) from the two microelectrodes.  The Plexon system outputs a .plx file that the SequenceAnalysis.m function converts to .mat format, via Plexon's MATLAB Offline SDK.  If a mouse has a unique identifier (e.g. 26881), then the function SequenceStim.m can be called as SequenceStim(26881).  The function will search for that mouse's retinotopic map created by the Retinotopy.m and MapRetinotopy.m functions.  It will then use the RetinoMap26881.mat file to generate a stimulus of flashing circles that target the retinotopic location of the implanted LFP-recording electrode. Over the course of 3 days, the mouse will see the same stimulus of flashing circles, each stimulus with 4 elements in the same locations and with the same relative timing, repeated 200 times per day.  As each circle flashes, the LFP shows a visually-evoked potential (VEP) that we measure over the course of the three days.  The VEPs potentiate to the training sequence, but not to other sequences.

#Steps:
1) Run Retinotopy(26881) to mouse #26881 while recording the LFP using Plexon system. Name the .plx file as RetinoDataDate_AnimalName.plx (e.g. July 5, 2016 is 20160705 so RetinoData20160705_26881.plx).  The Retinotopy.m file will output a file with the stimulus parameters named RetinoStim20160705_26881.mat .
#
2) Run MapRetinotopy(26881,20160705) as long as the RetinoData file and the RetinoStim file are both on the MATLAB path.

3) This will output a figure and a file named RetinoMap20160705_26881.mat . Clear the workspace, open that file, view the figure, and select the channel with the best retinotopy.  It should be a heat map with a clear center of greatest activity with a relatively circular decay around that center (like a 2D Gaussian). Write Channel = 1 if channel 1 is the best channel, or Channel = 2 if that is best (use the channel name at the top of the figure, ignoring the fact that some systems start at Channel 6 or what have you). Save all of the variables in the workspace as RetinoMap26881.mat .
#
4) The following day, run SequenceStim(26881) while recording the LFP with the Plexon system.  This will automatically load the retinotopic map created the previous day and target the center of mass of the VEPs. Save the .plx file as SeqData20160706_26881.plx (or whatever date is applicable).
#
5) Do this for 3 consecutive days.
#
6) If the first day of the 3-day sequence experiment was July 6, then run SeqAnWrapper(26881,20160706,3) .  This will output a series of figures that show the VEP in response to the sequence and how it changed over the 3 days.
