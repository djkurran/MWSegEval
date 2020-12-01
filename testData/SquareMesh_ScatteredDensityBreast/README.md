# Data and configuration files for *Scattered Density* breast 

***
These are the data and configuration files related to a *Scattered Density* breast. The files are used to provide an example of microwave inversion results produced from various microwave breast imaging algorithms for which the forward and inverse models use or have already been transformed by the user to rectangular (or square) mesh elements. The models are input to the MWSegEval toolbox by selecting the 'Rectangular mesh (Matlab)' option. Refer to [A. Rectangular mesh elements (Matlab)](https://github.com/djkurran/MWSegEval/wiki/2.1-Workflow-input-models) for a detailed description on how to input the forward and inverse models into the MWSegEval toolbox.


## Configuration folder

The frequency and axes files are stored in the *configuration* folder. The frequency file (*frequency.txt*) stores the incident field frequency measured in Hz used to collect the data. X and Y axes node files (*xNodes.txt*, and *yNodes.txt*) used to mesh the imaging domain into rectangular elements and to identify coordinates of each element (or node) of the model are also provided. These files represent distance along an axis in meters. 

## data folder

Values of the electrical properties (complex permittivity) for each element (or node) of the meshed 2D space used to represent the forward and reconstruction models are stored in matrices with variable names *input_forwardModel*, and *input_reconstructedModel*, and are saved in matlab formatted files *input_forwardModel.mat*, and *input_reconstructedModel.mat*, respectively, in the *data* folder. 

## Licensing

Apache License, 2.0

## Citing MWSegEval

D. Kurrant, N. Abdollahi, M. Omer, P. Mojabi, E. Fear, J. LoVetri, 'MWSegEval - An image analysis toolbox for microwave breast images', *SoftwareX*, 2020 (submitted for review)

## Citing electromagnetic model

Kurrant D, Fear E. "Regional estimation of the dielectric properties of inhomogeneous objects using near-field reflection data", *Inverse Probl*, 2012;28:075001. 

