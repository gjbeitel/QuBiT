# QuBiT - Quantification Tool for Biological Tubes


Overview:
We developed QuBiT, or Quantification tool for Biological Tubes, to measure thin biological tubes accurately and reproducibly in 3D. QuBiT is a series of Matlab scripts to be used with segmented confocal stacks, such as the output of Ilastik (http://ilastik.org). 

Minimum system requirements:

Hardware:
RAM: 8GB
CPU: 2.5GHz
OS: MacOSX 10.5 (Leopard)

Software:
Matlab 2014b
  Matlab code base for multi-stencil fast marching (MSFM) (see: Hassouna and Farag 2007)




This is a list of all Matlab script used in QuBiT, in roughly sequential order of the described image analysis method.

Module 1: Import tubes and skeletonize
s1_import_tubes
 
 This module imports the segmentation data for the tube (e.g. the *.h5 segmentation file exported from Ilastik) and skeletonizes them using MSFM. You can select multiple files to process overnight. The resulting datasets are split and stored in individual structures which will contain all post-processing data for each dataset.

 
 
Module 2: Display tubes
s2_display_tubes
 
 This module displays the segmented and skeletonized tubes. This is a check for proper skeletonization by visualizing the result in 3D. Segments of the skeleton are individually colored and indexed. Use the indices to specify the ROI in the next module.
 
 
Module 3: Specify tube ROI
s3_tube_ROI
 
 This module is required to specify the tube regions of interest. Enter the segment numbers sequentially as they appear in the figure created in module 2.
 
 
Module 4: Calculate tube and segment length
s4_calc_tube_length
 
 This module calculates length of the ROI using the skeleton and gathers information on coordinates of individual tube segments, which will be used in later calculations for e.g. normalizing segment data.
 
 
Module 5: Calculate tube cross-sections
s5_calc_cross_section
 
 This module calculates individual cross-sections, defined as the plane of the tube that is orthogonal to the local skeleton tangent vector. A figure containing the shapes of the successive cross-sections as well as their areas and circumferences can be displayed.
 
 
 
Module 6: Import surface data and mask them on the tube surface
s6_handle_junctions
 
 This module allows you to input an image stack of cell junctions or any tube surface data from the specified folder. It will then invert and project this image onto a mask of the tube with parameters specified by dilation_radius and dilation_value. (Default values result in 1 pixel-wide mask directly above the surface of the tube.)

This module can display the resulting masked data, overlaying the tube image and skeleton. This is done to ensure that the masked data appear acceptable for analysis and that the three components (tube skeleton, tube image, cell image) are all correctly aligned.
 
 
Module 7: Calculate surface object data
s7_calc_cells

 
 This module calculates properties of interest of cells or surface objects defined in the previous step. Properties include apical cell size, aspect ratio, and orientation relative to the tube axis.
 
 
Module 8: Normalize (and filter) surface object data
s8_normalize_cell_segments

 This module normalizes the data to cell segments specified in module 3 and applies user-created filters. 
 

Module 9: Unroll tube (create 2D surface plot)
s9_unroll_tube
 
 This module indexes the tube surface as if it were a thin layer and plots the data as a flat 2D image. This is particularly useful to visualize patterns on surface data. Since unrolling a non-cylindrical tube will deform it, direct calculations on the 2D image should be avoided.


Module: Reroll tube

 This module plots the 2D surface plot data from module 9 on a thin conical surface in 3D.



Supporting scripts

These scripts were used for debugging purposes and/or generating figures specifically for the QuBiT paper and may not be directly applicable to other projects. They are included for reproducibility and as potential starting points for other analyses using QuBiT.

s7x_plot_cells
s7y_bulk_reapply_filters
s8_view_segments
s9a_refine_tube (resulting data not used in paper)
s10a_set_colors
s10b1_aggregate_tube_data
s10b2_plot_tube_data
s10c1_aggregate_cell_data
s10c2_plot_cell_data
s10d2_plot_individual_celldata
s10e1_aggregate_2d_data
s10e2_plot_2d_data
s10f1_aggregate_reroll_data
s10f2_plot_reroll_data
