# Data and configuration files for *Heterogenously Dense* breast 

***

These are the data and configuration files related to an *Heterogenously dense* breast. The files are used to provide an example of microwave inversion results produced from the finite element method contrast source (FEM-CS) inversion technique. The forward and inverse models are represented with various sized triangular mesh elements. The models are input to the MWSegEval toolbox by selecting the *Triangular mesh* option. The procedure for inputing the forward and inverse models is described in the [Input forward and inverse models with triangular mesh elements](https://github.com/djkurran/MWSegEval/wiki/2.1-Workflow-input-models) section of the on-line manual.


## Configuration folder

Mesh information and the value of the complex contrast at nodes of the model are provided for the forward model with the *FwdMesh.msh* and *FwdMaterials.txt* files, respectively. Likewise, mesh information for the inverse models is provided by the *InvMesh.msh* file. The complex permittivity of the background is provided by the *PriorMaterials.txt* file. More general information about how the data were collected is provided by the *input_datacollect.xml* (forward model) and *input_csi_Img_Rec.xml* (inverse model). A frequency file (*frequency.txt*) indicating the incident field frequency is provided. Information about the rectangular mesh for which the triangular elements are transformed to is furnished with the *xNodes.txt*, and *yNodes.txt* files.  

## data folder

Ccontrast permittivity at nodes on the reconstructed model mesh at specific iterations of the inversion algorithm is provided in contrast files in the *data* folder in the *data* folder. The user must load the iteration number that identifies the contrast file by pressing the *Press to get contrast iteration file list* button. Once the list of available contrast iterations is displayed, the user selects the desired contrast iteration from the list.

## Licensing

Apache License, 2.0

## Citing MWSegEval

D. Kurrant, N. Abdollahi, M. Omer, P. Mojabi, E. Fear, J. LoVetri, "MWSegEval - An image analysis toolbox for microwave breast images", *SoftwareX*, 2020 (submitted for review)

## Citing electromagnetic model

Kurrant D, Fear E. "Regional estimation of the dielectric properties of inhomogeneous objects using near-field reflection data", *Inverse Probl*, 2012;28:075001. 

