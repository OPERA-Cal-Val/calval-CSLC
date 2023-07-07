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
There are two approaches to validate the OPERA Coregistered Single-Look Complex (<B>CSLC</B>) products and ensure the products meet the specified requirements which include the: (1) Absolute Location Error (<B>ALE</B>) estimation through point-target analysis using corner reflectors network in North America and (2) Relative Location Error (<B>RLE</B>) estimation by cross-correlating pairs of OPERA CSLC products.

### 0. Building the environments using conda
To run the ALE and RLE notebooks we need to create two separate environments and install relevant python dependencies. To install the dependencies, we recommend using __[Conda](https://docs.conda.io/en/latest/index.html)__, a cross-platform way to use Python that allows you to setup and use "virtual environments." These can help to keep dependencies for different sets of code separate. Specifically, we recommend using __[Miniforge](https://github.com/conda-forge/miniforge)__, a conda environment manager that uses conda-forge as its default code repo. Alternatively, see __[here](https://docs.anaconda.com/anaconda/install/)__ for help installing Anaconda and __[here](https://docs.conda.io/en/latest/miniconda.html)__ for installing Miniconda.

### 1. Install the python dependencies to run the ALE notebooks

Refer to `environment_ALE.yml` for an explicit list of software dependencies.

Use conda to install dependencies via your regular terminal:
```
conda install mamba
mamba env create -f environment_ALE.yml
conda activate calval_CSLC_ALE
python -m ipykernel install --user --name calval_CSLC_ALE
```

Whenever you are running the <B>ALE</B> notebooks make sure to activate the `calval_CSLC_ALE` environment by typing the following in your terminal:
```
conda activate calval_CSLC_ALE
```

Don't forget to deactivate the environment before running the RLE notebooks using:
```
conda deactivate
```

### 2. Install PyCuAmpcor and the python dependencies to run the RLE notebooks
PyCuAmpcor is included in ISCE2, and can be compiled/installed by CMake or SCons, together with ISCE2 or as a standalone version. An installation guide can be found at [isce-framework](https://github.com/isce-framework/isce2#building-isce).

Refer to `environment_RLE.yml` for an explicit list of software dependencies.

Again, use conda to install dependencies via your regular terminal:
```
conda install mamba
mamba env create -f environment_RLE.yml
conda activate calval_CSLC_RLE
python -m ipykernel install --user --name calval_CSLC_RLE
```

Whenever you are running the <B>RLE</B> notebooks make sure to activate the `calval_CSLC_RLE` environment by typing the following in your terminal:
```
conda activate calval_CSLC_RLE
```

Don't forget to deactivate the environment before running the ALE notebooks using:
```
conda deactivate
```

### Overview of the Validation Approaches

## Absolute Location Error

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

To characterize ALE, the difference between predicted target locations in Sentinel-1 images and its measured location are estimated in range and azimuth after multiple corrections. All validations are made on corner reflectors in North America (i.e., Rosamond, CA and Oklahoma)
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


## Key Contributors (Alphabetical)

M Grace Bato (JPL)
<br />
Liang Kang (SMU)
<br /> 
Jinwoo Kim (SMU)
<br />
Zhong Lu (SMU)
<br />
Simran Sangha (JPL)
