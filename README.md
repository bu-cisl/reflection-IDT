# Inverse Scattering for reflection intensity phase microscopy

This is the MATLAB implementation of the reflection intensity phase microscopy reconstruction pipeline from "Inverse scattering for reflection intensity phase microscopy" published in Biomedical Optics Express in 2020. This github provides example code for reconstructing the complex refractive index contrast of an object evaluated in a reflection imaging geometry. 

## Citation

If using this code or elements of it in your own research, please cite our corresponding publication:
[Alex Matlock, Anne Sentenac, Patrick C. Chaumet, Ji Yi, and Lei Tian, "Inverse scattering for reflection intensity phase microscopy," Biomed. Opt. Express 11, 911-926 (2020)](https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-11-2-911&id=425954)


## Implementation Requirements

This code was implemented entirely in MATLAB and has been coded and implemented successfully with MATLAB R2018b.

## Demo overview

The github repository includes the functions and script necessary for running the reconstruction pipeline recovering Henrietta Lacks (HeLa) cell phase and absorption features from experimental measurements with a reflection imaging setup. The reconstruction recovers the reflection images provided in Figure 5 of our publication.

To run the demo:
1. Download this repository to your computer and add them to your path.
2. Download the raw image data (.m file) from our [Google Drive]().
3. Open the 'main' script in the 'Scripts' folder.
4. In section 1. of the file, modify the 'file.dPath' to match the location of the raw image data.
5. Run the code. 

The file will save the reconstructions in the same location as the raw image data, and a Figure showing the object's recovered complex refractive index contrast will appear.
If issues arise in the code, you cannot access the data, or some other problem occurs, feel free to contact us or start a new Issue on this repository. 

## Abstract

Reflection phase imaging provides label-free, high-resolution characterization of biological samples, typically using interferometric-based techniques. Here, we investigate reflection phase microscopy from intensity-only measurements under diverse illumination. We evaluate the forward and inverse scattering model based on the first Born approximation for imaging scattering objects above a glass slide. Under this design, the measured field combines linear forward-scattering and height-dependent nonlinear back-scattering from the object that complicates object phase recovery. Using only the forward-scattering, we derive a linear inverse scattering model and evaluate this model’s validity range in simulation and experiment using a standard reflection microscope modified with a programmable light source. Our method provides enhanced contrast of thin, weakly scattering samples that complement transmission techniques. This model provides a promising development for creating simplified intensity-based reflection quantitative phase imaging systems easily adoptable for biological research.

## System Design


