+-------------------------------------------------------------------------+
|	NOTES                                                   19 Oct 2016   |
+-------------------------------------------------------------------------+

     The synthetic data set was generated using SimTB. It is formed by 
300 images of 100x100 taking in intervals of 2 seconds.

These data are grouped in the next files with some extra information:

File "Data.mat" -----------------------------------------------------------
	
     - Data  >> [double 300x100x100] contains all the images (100x100 voxels)
                 along 300 instances with 2 seconds between pictures.
     - X     >> [double 300x10000] is the same matrix as Data but each images was
			     unfolded into the same row.

  --- Extras ---
     - Clean >> [double 300x10000] matrix which only contains the mixture of sources
			     without adding any artifact.
     - Art   >> [double 300x10000] the implemented artifacts data set selected 
                 form the data set in [1].
     - m     >> [bool 300x10000] is a specific mask implemented to restrict the analysis
                 to the point of interest. 


File "Sources.mat" --------------------------------------------------------
     > This file contains the specific 15 implemented sources with their respective
time courses, they appear grouped in two different matrices:


     - TC    >> [double 300 x 15] contains the 15 different activation patters by columns 

     - SC    >> [double 15 x 100 x 100] It contains grouped the 15 different patters
                 formed by squared images of 100x100.


   References
---------------------------------------------------------------------------
[1] N. Correa, T. Adalı, and Y. Li, “Comparison of blind source separation algorithms 
for fMRI sing a new mat-lab toolbox: GIFT,” in Acoustics, Speech, and Signal Processing. 
Proceedings. ICASSP 2005.

