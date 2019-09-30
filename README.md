# RailWave #
**RailWave** is a software that implements a calculation technique with a semi-analytical finite element method 
for guided waves and its application to simulation and modal analysis of wave propagation in a rail with an arbitrary cross-section. Dispersion curves for any kinds 
of structures can be calculated by the **SAFEM**. Visualization results of guided wave propagation are shown. Also, this software allows to design the cross section of a rail with the standard parameters.


## 1. Load rail


### 1.1

In the **File** menu a rail is loaded from an external file that contains the rail manufacturing parameters (standard).

### 1.2

In the **Mesh Refinement** box the density of the mesh that discretizes the cross section of the rail is specified with two possibilities: **Medium** or **High**. The mesh is displayed in the lower left region **Section**.
![](/images/init.jpg)

### 1.3
In the lower right region **Rail 3D** the rail is displayed as the extension (as a prism) of a cross section. In the attached panel **Modify Section** you select the number of sections (**Sections Number**) and the spacing (**Spacing**) between them.

#### 1.3.1
With the **Zoom In** and **Zoom Out** buttons you can scale the rail size, as well as the mouse scroll. With the arrow keys on the keyboard you can move the rail to different positions.

## 2.Design a rail

### 2.1
In the **File menu** of the main window, specifically in the **Create Rail** submenu opens a new window that allows you to create a new rail by entering the dimensions that it will have. Once the dimensions are entered, the rail can be visualized, as well as the mesh that approximates them.

### 2.2
Once the rail is created, it can be added to the main window by pressing the **Add to Workspace button** or saved to a file through the **File -> Save Rail** menu.

![](/images/designrail.jpg)

## 3. Calculation of dispersion curves

### 3.1
In the upper right panel **Graph Properties** the maximum number of waves and the number of dispersion curves are selected.
### 3.2
The **Compute** button shows the dispersion curves in the upper region.
### 3.3
Enabling the **Show Curves** option in the **Display** menu (by default) the curves are displayed continuously. With this option disabled, the points on the curves corresponding to frequency values / phase speeds are shown.
![](/images/curves.jpg)
## 4.Animation of the propagation of a wave
### 4.1
Selecting a point on a dispersion curve in the upper region, a wave is chosen whose propagation is visualized with the **Animate** button in the **Rail 3D** box.
![](/images/moving.jpg)
### 4.2
When the propagation in the rail is animated, the button changes its name to **Stop Animation** and then allows to stop the animation.
## 5. About RailWave

RailWave was developed using the **C++** programming language and the **Qt Creator** development platform. To visualize the rail, the dispersion curves and the displacements, the **OpenGL** graphic library was used. The **ARPACK** library was used to solve the problems of Linear Algebra that are obtained from the formulation of the problem.
RailWave was conceived by a group of researchers from the **Instituto de Cibernética, Matemática y Física (ICIMAF)** of **Cuba** as a tool for the analysis of ultrasonic wave dispersion curves.















