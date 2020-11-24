# MWSegEval –An image analysis toolbox for microwave breast images

***

## What is MWSegEval?

MWSegEval is an image analysis toolbox that employs methods to automatically segment medical microwave breast images into regions of interest corresponding to various tissue types. The segmentation leads to the decomposition of the breast interior into disjoint tissue masks. An array of region and distance-based metrics are applied to compare masks extracted from reconstructed images and ground truth models. The quantitative results reveal the accuracy with which the geometric and dielectric properties are reconstructed, and are supplemented with qualitative information. Consequently, the toolbox provides a framework that effectively furnishes quantitative and qualitative assessment of regions that contain a specific tissue. The resulting information facilitates comparisons that provide valuable insight into complex issues that impact image quality. This document describes how to apply the techniques furnished by the toolbox to microwave breast images.

## Software Requirements

The toolbox requires Matlab R2017a or above and Matlab’s image processing and statistics toolboxes. The software can run on Windows (32 and 64 bit), Mac OS X (64 bit), Linx (64 bit) operating systems. The toolbox includes the export_fig toolbox[1] for exporting figures from Matlab to standard image and document formats.

![](https://github.com/djkurran/MWSegEval/blob/main/figures/figureReadMe_1.png)

Figure 1. Directory structure of toolbox

## Installation

Installation of the MWSegEval toolbox is achieved by unpacking the MWSegEval libraries at the user’s choice of directory. The libraries can be obtained as a *.zip* archive or by cloning the tool repository (http://www.github.com/djkurran/MWSegEval). The toolbox expects the MWSegEval libraries to be present in the directory structure shown in figure 1. A test data set is also provided in the *testData* folder and is comprised of models represented with square mesh elements and triangular mesh elements that are stored in the *SquareMesh_ScatteredDensityBreast* and *TriangularMesh_HeterogenouslyDenseBreast* folders, respectively.

![](https://github.com/djkurran/MWSegEval/blob/main/figures/figure2.png)
Figure 2. GUI allows user to set parameters for tasks implemented with functions and classes contained in MWSegEval libraries.

## How to use MWSegEval toolbox

The graphical user interface (GUI) shown in figure 2 may be started by right clicking on *MWSegEval_Toolbox.mlapp* (figure 1), and executing the *‘Run’* command. The GUI allows the user to set parameters and to run the tasks that are implemented with functions and classes contained in the *MWSegEval* libraries. When the user selects a task to run, the GUI reads all of the user parameters and calls function *taskManager()*. The *taskManager()* reads the user parameters gathered by the GUI, validates this information, and implements a state machine to run the requested task. The state machine constructs the required objects that call object methods to complete the task and the data structure of the objects hold the results. The state machine carries out error handling and passes error flags and the status of various tasks to the GUI so that the user can monitor the progress of the image analysis. The user requests a sequence of tasks to process and analyze an image. The collective set of tasks comprise a workflow is described in the [on-line manual](https://github.com/djkurran/MWSegEval/wiki) hosted in the project wiki page.

## On-line manual and videos

A comprehensive [on-line manual](https://github.com/djkurran/MWSegEval/wiki) providing a detailed description of how to input images and apply the toolbox methods to the images is hosted by the [repository wiki page](https://github.com/djkurran/MWSegEval/wiki). The on-line manual is supplemented with a series of videos of examples on the complete workflow to process a number of medical microwave breast images. 

## Licensing

Apache License, 2.0

## Citing MWSegEval

D. Kurrant, N. Abdollahi, M. Omer, P. Mojabi, E. Fear, J. LoVetri, “MWSegEval –An image analysis toolbox for microwave breast images”, SoftwareX, 2020 (submitted for review) 




