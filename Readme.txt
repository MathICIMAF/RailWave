
%%%%%%%%%%%%%%%%%%%%%%%%%%
RailWave, 1.0, description, and/or features of the program.
Qt Creator.
Install, uninstall, configuration, and operating instructions.
Files list.
Credit (authors),ICIMAF, CUBA. ahmedlp9@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%


RailWave v1.0
-----

The RailWave program simulates the propagation of an ultrasonic wave through a rail. Given a lane specified in an external file with the manufacturing parameters, a mesh (triangular) is generated that discretizes the cross section. From this mesh are the dispersion curves and choosing a point on one of them the propagation of the corresponding wave on the rail is visualized.


1. Load rail
--------------

1.1 In the "File" menu a rail is loaded from an external file that contains the rail manufacturing parameters (standard).

1.2 In the "Mesh Refinement" box the density of the mesh that discretizes the cross section of the rail is specified with two possibilities: Medium or High. The mesh is displayed in the lower left region "Section".

1.3 In the lower right region "Rail 3D" the rail is displayed as the extension (as a prism) of a cross section. In the attached panel "Modify Section" you select the number of sections (Sections Number) and the spacing ("Spacing") between them.

1.3.1 With the "Zoom In" and "Zoom Out" buttons you can scale the rail size, as well as the mouse scroll. With the arrow keys on the keyboard you can move the rail to different positions.


2. Find dispersion curves
------------------------------

2.1 In the upper right panel "Graph Properties" the maximum number of waves and the number of dispersion curves are selected.

2.2 The "Compute" button shows the dispersion curves in the upper region.

2.3 Enabling the "Show Curves" option in the "Display" menu (by default) the curves are displayed continuously. With this option disabled, the points on the curves corresponding to frequency values / phase speeds are shown.


3. Animation of the propagation of a wave
------------------------------------------

3.1 Selecting a point on a dispersion curve in the upper region, a wave is chosen whose propagation is visualized with the "Animate" button in the "Rail 3D" box.

3.2 When the propagation in the rail is animated, the button changes its name to "Stop Animation" and then allows to stop the animation.
