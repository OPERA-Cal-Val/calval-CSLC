# OPERA Coregistered Single Look Complex (CSLC) validation tools
Repository of tools intended to address the Cal/Val validation plan for the Level 2 OPERA CSLC product suite. Namely, verify that products meet the established geolocation accuracy requirements.

## Contents
  
1. [Setup](#setup)
    -   [Installing With Conda](#installing-with-conda)
    -   [Running notebooks](#running-notebooks)
2. [Overview of notebooks](#overview-of-notebooks)
3. [Absolute Geolocation Error](#absolute-geolocation-error)
    -   [ALE of COMPASS derived products](#ale-of-compass-derived-products)
    -   [ALE of provisional OPERA CSLC products](#ale-of-provisional-opera-cslc-products)
4. [Relative Geolocation Error](#relative-geolocation-error)
5. [References](#references)
6. [Contributors](#contributors)

## Setup

Refer to `environment.yml` for an explicit list of software dependencies.

### Installing With Conda

To install the dependencies we recommend creating a software environment using __[Conda](https://docs.conda.io/en/latest/index.html)__, a cross-platform way to use Python that allows you to setup and use "virtual environments." These can help to keep dependencies for different sets of code separate. Specifically, we recommend using __[Miniforge](https://github.com/conda-forge/miniforge)__, a conda environment manager that uses conda-forge as its default code repo. Alternatively, see __[here](https://docs.anaconda.com/anaconda/install/)__ for help installing Anaconda and __[here](https://docs.conda.io/en/latest/miniconda.html)__ for installing Miniconda.

Using conda to install dependencies:
```
conda install mamba
mamba env create -f environment.yml
conda activate calval_CSLC
```

### Running notebooks

After successful installation, you may run these tools and visualize outputs and statistics through the provided Jupyter notebooks.

But first, you will have to install and initalize the kernel by activating your environment (i.e. `conda activate calval_CSLC`) and then entering the following command:

`python -m ipykernel install --user --name calval_CSLC`


Then after opening the desired notebook in your preferred code editor, make sure to set the kernel to `calval_CSLC`.

In regards to the use of a specific code editor, we recommend __[Visual Studio Code](https://code.visualstudio.com/)__.


## Overview of notebooks

The Los Angeles Basin in Southern California is a sedimentary basin characterized by significant structural relief and complex faulting. Introduced below are a suite of notebooks which visualize select CSLC products over this region and evaluate their geolocation accuracy with respect to corner reflectors.

## Absolute Geolocation Error

The Absolute Geolocation Error (<B>ALE</B>) can be evaluated using corner reflector (CR) arrays as in [Gisinger et al. 2021]. Specifically, we leverage the OPERA CRs shown in Figure 1 and listed in Table 1. We evaluate the ALE for each CSLC image in a given stack to obtain a time series of offset. For each CSLC image, the image coordinates of the CRs are identified at sub-pixel accuracy using point-target analysis and then compared to field-measured CR coordinates (i.e. [Gisinger et al. 2021], [Schubert et al. 2017]). For example, Figure 2 shows the estimated amplitude peak of the CR arrays for Sentinel-1 A/B for a corner reflector in California; this can be compared to the geolocation of the specific corner reflector and determine if the relevant requirement is met.


<p align="center">
  <img width="90%" height="90%" src="https://user-images.githubusercontent.com/13227405/214785618-d7186ebd-67c6-44f6-b771-c8180b1c87ab.png">
</p>
<B>Figure 1.</B> Proposed sites of the corner reflectors (triangles) to be deployed by OPERA across San Andreas Fault’s creeping segment with different spatial baselines (i.e. 9 - 25 km) and in a subsiding region in Central Valley, CA. Inset map shows the Sentinel-1 (solid lines) and planned NISAR (broken lines) swaths for both their ascending (blue) and descending (red) tracks. The fault trace and the inventory of public open space and the Protected Areas Database are openly available from USGS repositories.
<br />
<br />


<p align="center">
  <img width="80%" height="80%" src="https://user-images.githubusercontent.com/13227405/214785522-2f632ec8-375c-46b7-ac2f-ca47e219b943.png">
</p>
<B>Table 1.</B> Summary of Validation validation sites.
<br />
<br />


<p align="center">
  <img width="90%" height="90%" src="https://user-images.githubusercontent.com/13227405/214784945-7dc6fd73-259f-4628-a83c-1ebb82283578.png">
</p>
<B>Figure 2.</B> (a) One of the 2.4-m corner reflectors in Rosamond, CA, as seen from a Sentinel-1 descending track CSLC with the estimated phase center (yellow cross). (b) Example of an absolute location error analysis for the Rosamond CR array using a stack of Sentinel-1A/B CSLCs spanning Jan 6, 2021 to Dec 14, 2021. 
<br />
<br />


For CR arrays, position data can be accessed through CEOS or on their public repositories (e.g. Rosamond CR array) or -- where applicable -- using co-located cGNSS stations.

To characterize ALE, the difference between predicted target locations in Sentinel-1 images and its measured location are estimated in range and azimuth after multiple corrections. All validations are made on corner reflectors in Rosamond, CA.
<br />
<br />
<br />
<br />

### ALE of COMPASS derived products

ALE assessment performed on CSLC products derived from the COregistered Multi-temPorAl Sar Slc (COMPASS) __[software](https://github.com/opera-adt/COMPASS)__, through the (__[_ALE_.ipynb](https://github.com/OPERA-Cal-Val/calval-CSLC/blob/dev/_ALE_.ipynb)__) notebook. COMPASS is a high level package developed around the InSAR Scientific Computing Environment 3 (__[ISCE3](https://github.com/isce-framework/isce3)__) framework to generate coregistered multi-temporal SAR SLCs in either geographic or radar coordinates.
<br />
<br />
<br />
<br />

### ALE of provisional OPERA CSLC products

ALE assessment performed directly on OPERA CSLC provisional products, through the (__[ALE.ipynb](https://github.com/OPERA-Cal-Val/calval-CSLC/blob/dev/ALE.ipynb)__) notebook.
<br />
<br />
<br />
<br />


## Relative Geolocation Error

The Relative Geolocation Error (<B>RLE</B>) can be evaluated from each stack of CSLCs following state-of-the-art methods in geolocation accuracy assessment for stack coregistration (e.g. [Yunjun et al., 2022], [Fattahi et al., 2016], [Lei et al., 2021]). The accuracy of determining the misregistration between two SAR images is known to decrease as temporal decorrelation increases due to the increased time difference between the pairs. The issue can be circumvented by using a fully connected network of small temporal baseline image pairs that overlap in time to obtain an unbiased time series of misregistration/offset [Fattahi et al., 2016]. For each image in the CSLC stack, we pair each CSLC image with its N-neighbors in time as illustrated in Figure 3. Then, we densify grids of sub-pixel differential offsets by cross-correlating each of the image pairs and apply a simple least squares method to obtain the time series of RLEs between each reference image and the secondary image within the stack. The resulting RLEs in the range and azimuth directions for each track are then statistically analyzed (e.g. L2-norm residual analysis). Unreliable pixels (i.e. pixels that have large residuals or have low coherence values over time) are discarded based on a predefined threshold value (TBD) prior to comparing against the requirements. This comparison provides metric offsets in each coordinate direction to determine if the requirement is satisfied.


<p align="center">
  <img width="90%" height="90%" src="https://user-images.githubusercontent.com/13227405/214784837-7aa19c05-ac42-4a55-b8a2-072068542e9b.png">
</p>
<B>Figure 3.</B> (a) Examples of the sequential network of CSLCs, with circles showing the SAR acquisitions and arcs representing the pairs of CSLCs [Fattahi et al., 2016]. (b) Temporal matching using normalized cross-correlation, from  [Lei et al., 2021]. 
<br />
<br />
<br />
<br />

<I>To be added</I>. To characterize RLE, range and azimuth offsets are estimated from neighboring pairs of coregistered Sentinel-1 SLCs and the offsets are measured by ISCE-2's AmpCor module.

## References

Fattahi, H., Agram, P., and Simons, M. (2016). A network-based enhanced spectral diversity approach for TOPS time-series analysis. IEEE Transactions on Geoscience and Remote Sensing, 55(2), 777-786, https://doi.org/10.1109/TGRS.2016.2614925

Gisinger, C., Schubert, A., Breit, H., Garthwaite, M., Balss, U., Willberg, M., Small, D., Eineder, M. and Miranda, N. (2020). In-depth verification of Sentinel-1 and TerraSAR-X geolocation accuracy using the Australian corner reflector array. IEEE Transactions on Geoscience and Remote Sensing, 59(2), https://doi.org/10.1109/TGRS.2019.2961248.

Kellndorfer, J., Cartus, O., Lavalle, M., Magnard, C., Milillo, P., Oveisgharan, S., Osmanoglu, B., Rosen, P.A. and Wegmüller, U. (2022). Global seasonal Sentinel-1 interferometric coherence and backscatter data set. Scientific Data, 9(1), 73, https://doi.org/10.1038/s41597-022-01189-6

Lei, Y., Gardner, A., and Agram, P. (2021). Autonomous repeat image feature tracking (autoRIFT) and its application for tracking ice displacement. Remote Sensing, 13(4), 749, https://doi.org/10.3390/rs13040749

Schubert, A., Miranda, N., Geudtner, D., and Small, D. (2017). Sentinel-1A/B combined product geolocation accuracy. Remote Sensing, 9(6), 607, https://doi.org/10.3390/rs9060607

Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., and Aoki, Y. (2022). Range Geolocation Accuracy of C-/L-Band SAR and its Implications for Operational Stack Coregistration. IEEE Transactions on Geoscience and Remote Sensing, 60, 1-19, https://doi.org/10.1109/TGRS.2022.3168509


## Contributors

Grace Bato
<br />
Zhong Lu
<br />
Jin Woo Kim
<br />
Simran Sangha
