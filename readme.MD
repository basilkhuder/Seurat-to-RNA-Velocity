# Seurat to RNA-Velocity
#### *By Basil Khuder*

This guide will demonstrate how to use a processed/normalized Seurat object in conjunction with an RNA Velocity analysis. Keep in mind that although Seurat is R-based, all of the available RNA Velocity software/packages are Python, so we will be moving back and forth between the two. We will be using the following programs: <br>

- [scVelo](https://github.com/theislab/scvelo) (For RNA Velocity)
- [Velocyto](http://velocyto.org/) (To produce our initial RNA Velocity Object)
- [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/) (For manipulation of our RNA Velocity object)
- Seurat
- Samtools 
 -- *optional* (Velocyto will run Samtools sort on unsorted .bam)

To start, we will be generating loom files (a file format designed for genomics datasets such as single-cell) for every single-cell sample you used in your Seurat analysis using Velocyto's ```velocyto run```

*Although there is a plethora of more up-to-date and less memory intensive software to produce the spliced/unspliced matrix required
for RNA Velocity, Velocyto's run command is still a more straightforward approach.*

Note that I have found it easier to use ```velocyto run``` 
for whichever scRNA-seq chemistry you are working with rather than Velocyto's "ready-to-use subcommands." 
```
velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf bam_file.bam annotation.gtf
```
*If time or memory is limited, you can disregard the -m parameter. Just note that this could contribute a downstream confounding factor
later on in the analysis.*

Once this step has finished and your loom file is generated, we can go ahead and use anndata to import our loom file and make the necessary adjustments/additions. 

But, before we do that, we will export all of the necessary meta-data from our Seurat object that is needed for our RNA Velocity object. This includes:

- Filtered Cell IDS
- UMAP or TSNE coordinates
- Clusters (Optional)
- Cluster Colors (Optional)

One way we can access our filtered Cell ID's is through Seurat's ```GetAssayData``` function and extract the column names:

```
%%R
write.csv(colnames(GetAssayData(seurat_object)), file = "cellID_obs.csv")
```
*If you have a Seurat object that is composed of multiple single-cell samples, you either can use the code above, and then later
on use some type of pattern to extract each sample (for example, if you added unique cell prefixes to each sample then you
could use that pattern.)  Likewise, you can also create a cell ID observation file for every sample, and use each one 
individually to filter each RNA Velocity object.*

To get UMAP or TSNE coordinates, we use the ```Embeddings``` function:
```
%%R
write.csv(Embeddings(seurat_object, reduction = "umap"), file = "cell_embeddings.csv")
```
And finally we can extract our clusters with: 

```
%%R
write.csv(seurat_object@meta.data$seurat_clusters, file = "clusters.csv")
```

We now can import our loom file(s) and all of our Seurat meta-data using anndata

```
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

```
sample_one = sample_one[sample_one[np.isin(sample_one.obs.index,cellID_obs[0])]]
```
#### *If your Seurat object has multiple samples*<br>
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

```
cellID_obs_sample_one = cellID_obs[cellID_obs_sample_one[0].str.contrains("sample1_")]
cellID_obs_sample_two = cellID_obs[cellID_obs_sample_two[0].str.contrains("sample2_")]
sample_one = sample_one[sample_one[np.isin(sample_one.obs.index,cellID_obs_sample_one)]]
sample_two = sample_two[sample_one[np.isin(sample_two.obs.index,cellID_obs_sample_two)]]
```
Once all the samples have been properly filtered, we can merge them into one. 
```
sample_one = sample_one.concatenate(sample_two,sample_three,sample_four)
```
-------------------------------------------------
Now that we have our Velocity file filtered based upon our Seurat object, we can go ahead and add UMAP coordinates. We'll first upload them:
```
umap = pd.read_csv("umap.csv")
```

With the coordinates, we will need to make sure we add them so they match the order of the Cell IDs in our anndata object. Our Cell IDs are rownames in the observation layer of our object, so we can view them by using the following:

```
sample_one.obs.index
```
Let's cast our index as a data frame and change the column name

```
sample_one_index = pd.DataFrame(sample_one.obs.index]
sample_one = sample_one_index.rename(columns = {0:'Cell ID'})
```
Let's also change the first column of our UMAP data frame to the same name:

```
umap = umap.rename(Columns = {'Unnamed: 0':'Cell ID'})
```
Now if we merge our index dataframe with our UMAP, the order will match our anndata object.
```
umap_ordered = sample_one_index.merge(umap,on="Cell ID")
```
Since we're certain the orders are the same, we can remove the first column of the data frame and, finally, add the UMAP coordinates to our anndata object.
```
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
```
At this point, we can now run the scVelo commands and generate our RNA Velocity plot based upon our Seurat UMAP coordinates.

```
scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis='umap')
```
