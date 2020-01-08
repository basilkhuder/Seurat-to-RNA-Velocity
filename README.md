# Seurat to RNA-Velocity

This guide will show you how you can use a processed/normalized Seurat object in conjunction with your RNA Velocity analysis. 
We will be using the following programs: <br>

- [scVelo](https://github.com/theislab/scvelo) (For RNA Velocity)
- [Velocyto](http://velocyto.org/) (To produce our initial RNA Velocity Object)
- [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/) (For manipulation of our RNA Velocity object)
- Seurat (To extract meta data from processed Seurat object)
- Samtools *optional* (Velocyto will run Samtools sort on unsorted .bam)

To start, we'll be generating loom files for every single-cell sample you used in your Seurat analysis using Velocyto's <br>```velocyto run```

*Although there are a plethora of more up-to-date and less memory intensive software to produce the spliced/unspliced matrix required
for RNA Velocity, Velocyto's run command is still a more straightforward approach.*

Note that I've found it better to use ```velocyto run``` for whichever scRNA-seq chemistry you are working with rather than Velocyto's "ready-to-use subcommands." 
The only caveat is you'll have to use a [filtered barcode file](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-) in conjunction with the run command, but you can easily find those from
your sequencing company's website. I've uploaded 10Xv2 and 10Xv3 filtered barcodes in this repo that you can use if that is the
chemistry for your samples. The command is then as follows:

```
velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf bam_file.bam annotation.gtf
```
Once this step has finished and your loom file is generated, we can go ahead and use anndata to import our loom file and make the necessary adjustments/additions. If you'd like a more interactive visualization for the rest of this guide,
an interactive jupyter notebook is provided here.

But, before we do that, we will export all of the necessary meta data from our Seurat object that we will be incorporating into
our RNA Velocity object. This includes:

- Filtered Cell IDS
- UMAP or TSNE Embeddings
- Clusters

One way we can access our filtered Cell ID's is through Seurat's ```GetAssayData``` method and extract the column names:

```
write.table(colnames(GetAssayData(seurat_object)), file = "cellID_obs.csv")
```
To get UMAP coordinates, we use the ```Embeddings``` method:
```
write.table(Embeddings(seurat_object, reduction = "umap"), file = "cell_embeddings.csv")

```

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
```

  

```
sample_obs = pd.read_csv("sample_obs.csv")
sample_one = sample_one[sample_one[np.isin(sample_one.obs.index,sample_obs[0])]]
```

