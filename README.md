# scripts_for_2021_paper
> The scripts are used for the paper 2021ApJS..254....3M (see https://ui.adsabs.harvard.edu/abs/2021ApJS..254....3M/abstract), not all, but some examples.
## Description 
### 1. clusters.py
- Used for identifying molecular clouds in a 3D 12CO or 13CO datacube.
#### Outputs
- dendro.fits: dendrogram 
- cluster.fits: masks of the identified clouds
- dendro_catalog.cat: catalog of the identified hierarchical dendrogram structures
- clusters.cat: catalog of the identified clouds, containing information on the PPV centroid positions, angular size, velocity dispersion, etc.
- parameters.txt: text file recording the input parameters for structure identification


### 2. plt_ave.py
- Make plots of the average spectrum for a series of molecular clouds. 


### 3. maskphy.pro 
- Used to calculate the physical maps such as excitation temperature, optical depth, and column density of the molecular gas in different layers of velocity, i.e., with different distances, and make histograms of the derived properties. The three histograms of column density distributions of 12CO, 13CO, and C18O are given as examples (in the directory './example/').

### 4. pdf_panels.pro
- Plot the N-PDFs of selected molecular clouds (make panel plots), fit each N-PDF with both LN and PL functions using mpfit, and select a better model with a smaller reduced chi-squared. The figure 'pdf_panel_loc.eps' in the example folder is given as an example.

### 5. plot_clouds.pro 
- Overlay the boundaries of the identified molecular clouds on an m0 map and a p-v map of the 12CO or 13CO data, with each cloud corresponding to one random color. (see 'layer1_m0.eps' in the directory './example/')
