# Seurat to RNA-Velocity
#### *By Basil Khuder*

- [Introduction](#introduction)
- [Generating Loom files](#generating-loom-files)
  * [Kallisto Bustools](#kallisto-bustools)
    + [Installation](#installation)
    + [Usage](#usage)
  * [Velocyto](#velocyto)
    + [Installation](#installation-1)
    + [Usage](#usage-1)
- [Extracting Meta-data From a Seurat Object](#extracting-meta-data)
- [Integrating Loom File and Meta-data](#integrating-loom-file-and-meta-data)
  * [Multiple-Sample Integration](#multiple-sample-integration)
- [Running RNA Velocity](#running-rna-velocity)
- [FAQ](#faq)

## Introduction
This guide will demonstrate how to use a processed/normalized Seurat object in conjunction with an RNA Velocity analysis. Keep in mind that although Seurat is R-based, all of the available RNA Velocity software/packages are Python, so we will be moving back and forth between the two. We will be using the following programs: <br>

- [scVelo](https://github.com/theislab/scvelo) (For RNA Velocity)
- [Velocyto](http://velocyto.org/) or [Kallisto Bustools](https://www.kallistobus.tools/introduction) (To produce our initial RNA Velocity Object) 
- [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/) (For manipulation of our RNA Velocity object)
- Seurat
- Samtools 
 -- *optional* (Velocyto will run Samtools sort on unsorted .bam)
 

 ## Generating Loom files

To start, we will be generating loom files (a file format designed for genomics datasets such as single-cell) for every single-cell sample you used in your Seurat analysis. The two methods I will discuss for doing this will be Velocyto's ```run``` command and Kallisto Bustools (KB.) 

### Kallisto Bustools 
#### Installation
As of 8/17, to get the most up-to-date Kallisto Bustools (KB), you should install it from Github as shown below. Some of these parameters won't work if you install KB from pypi.org

```
pip install git+https://github.com/pachterlab/kb_python@devel
```
#### Usage
Using KB is a two-step process with you first creating a reference and then generating the counts table. You'll need FASTA and GTF files for your species (download them from Ensembl if you need them: https://useast.ensembl.org/index.html)

The ```-g```, ```-f1```, ```-f2```, ```c1```, ```-c2``` are all files that will be generated from the kb ref step (you don't need to provide them, just give them names. ) 
```
kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno -n 4 \
fasta.fa \
gtf.gtf
```
We'll then use ```kb count``` to generate our LOOM file. -x specifies the single-cell technology (use ```kb --list``` for all technology parameters) and --lamanno specifies we want to calculate RNA velocity. 

```
kb count -i transcriptome.idx -g t2g.txt -x 10xv2 --workflow lamanno --loom -c1 cdna_t2c.txt -c2 intron_t2c.txt read_1.fastq.gz read_2.fastq.gz  
```
### Velocyto
#### Installation
Velocyto can be installed through pip. 
```
#Download dependencies first
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install velocyto
```
#### Usage
Note that I have found it easier to use ```velocyto run``` 
for whichever scRNA-seq chemistry you are working with rather than Velocyto's "ready-to-use subcommands." 
``` 
velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf bam_file.bam annotation.gtf
```
*If time or memory is limited, you can disregard the -m parameter. Just note that this could contribute a downstream confounding factor
later on in the analysis.*

Once this step has finished and your loom file is generated, we can go ahead and use anndata to import our loom file and make the necessary adjustments/additions. 

But, before we do that, we will export all of the necessary meta-data from our Seurat object that is needed for our RNA Velocity object. This includes:

- Filtered Cell Ids
- UMAP or TSNE coordinates
- Clusters (Optional)
- Cluster Colors (Optional)

## Extracting Meta-data
<i> I put the meta-extraction steps in this [Google Collab](https://colab.research.google.com/drive/1LDnUGrALuR6ulXUWtH9K1mPe6ESEKk3K?usp=sharing) document if you're interested in playing around with the code. </i><br>
One way we can access our filtered cell id's is through Seurat's ```Cells``` function:

``` r
%%R
write.csv(Cells(seurat_object), file = "cellID_obs.csv", row.names = FALSE)
```

*If you have a Seurat object that is composed of multiple single-cell samples, you either can use the code above, and then later
some type of pattern to extract each sample (for example, if you added unique cell prefixes to each sample then you
could use that pattern.) Likewise, you can also create a cell ID observation file for every sample, and use each one 
individually to filter each RNA Velocity object.*

To get UMAP or TSNE coordinates, we use the ```Embeddings``` function:

``` r
%%R
write.csv(Embeddings(seurat_object, reduction = "umap"), file = "cell_embeddings.csv")
```
And finally we can extract our clusters with: 

```r
%%R
write.csv(seurat_object@meta.data$seurat_clusters, file = "clusters.csv")
```
## Integrating Loom File and Meta-data
We now can import our loom file(s) and all of our Seurat meta-data using anndata

```python
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
%load_ext rpy2.ipython

sample_one = anndata.read_loom("sample_one.loom")
.... 
sample_n = anndata.read_loom("sample_n.loom")

sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters_obs.csv")
```

With our extracted Cell IDs from Seurat, we'll need to filter our uploaded loom (now as an anndata object) based upon them.  

``` python
sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["x"])]
```
### Multiple-Sample Integration
-------------------------------------------------
If you have individual observation files for every sample, you'll do the filtering above one by one.  If you have a combined observation file, you'll want to filter it based upon the cell pattern and then use that to filter the RNA Velocity sample. For example, if these
were your Cell IDs:
<br><br>

| Sample Cell IDs | 
| ------------- | 
| sample1_ACTCACT |
| sample1_ACTCCAC |
|  .....          |
| sample2_CACACTG |

<br>

You could use the pattern sample1_, sample2_, to filter as such:

``` python
cellID_obs_sample_one = cellID_obs[cellID_obs_sample_one[0].str.contrains("sample1_")]
cellID_obs_sample_two = cellID_obs[cellID_obs_sample_two[0].str.contrains("sample2_")]
sample_one = sample_one[np.isin(sample_one.obs.index,cellID_obs_sample_one)]
sample_two = sample_one[np.isin(sample_two.obs.index,cellID_obs_sample_two)]
```
Once all the samples have been properly filtered, we can merge them into one. 
``` python
sample_one = sample_one.concatenate(sample_two,sample_three,sample_four)
```
-------------------------------------------------
Now that we have our Velocity file filtered based upon our Seurat object, we can go ahead and add UMAP coordinates. We'll first upload them:
``` python
umap = pd.read_csv("umap.csv")
```

With the coordinates, we will need to make sure we add them so they match the order of the Cell IDs in our anndata object. Our Cell IDs are rownames in the observation layer of our object, so we can view them by using the following:

``` python
sample_one.obs.index
```
Let's cast our index as a data frame and change the column name

``` python
sample_one_index = pd.DataFrame(sample_one.obs.index]
sample_one = sample_one_index.rename(columns = {0:'Cell ID'})
```
Let's also change the first column of our UMAP data frame to the same name:

``` python
umap = umap.rename(Columns = {'Unnamed: 0':'Cell ID'})
```
Now if we merge our index dataframe with our UMAP, the order will match our anndata object.
``` python
umap_ordered = sample_one_index.merge(umap,on="Cell ID")
```
Since we're certain the orders are the same, we can remove the first column of the data frame and add the UMAP coordinates to our anndata object.
``` python
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
```
Clusters and their cluster colors can be added in the same fashion (and again, they must match the order of the Cell IDs.) Instead of adding them as an multidimensional observation ('obsm'), we'd add them under the unstructured annotation 'uns.'

``` python
sample_one.uns['Cluster_colors']
```

## Running RNA Velocity

At this point, we can now run the scVelo commands and generate our RNA Velocity plot based upon our Seurat UMAP coordinates.

``` python
scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis='umap')
```
If you want to incorporate your clusters and cluster colors in the embedding plot, the parameters you would add would be:
```
color = sample_one.uns['Cluster_colors']
```

## FAQ 

*For Kallisto Bustools, why should I generate a loom file rather than an h5ad file?*
<br>
Either one is fine. I recommended generating a loom file just in case the user wanted to use Velocyto for their RNA velocity, but the h5ad format makes more sense if you're using scVelo (the Anndata format is an extension of h5ad.)

