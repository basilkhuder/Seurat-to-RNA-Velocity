# Seurat to RNA-Velocity

This guide will show you how you can use a processed/normalized Seurat object in conjunction with your RNA Velocity analysis. 
We will be using the following programs: <br>

- [scVelo](https://github.com/theislab/scvelo)
- [Anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/)
- [Velocyto](http://velocyto.org/)
- Seurat
- Samtools

To start, we'll be generating loom files for every single-cell sample you used in your Seurat analysis using Velocyto's <br>```velocyto run```

Note that I've found it better to use ```velocyto run``` for which ever scRNA-seq chemistry you are working with rather than Velocyto's "ready-to-use subcommands." 
The only caveat is you'll have to use a filtered barcode file in conjunction with the run command, but you can easily find those from
your sequencing company's website. I've uploaded 10XV2 and 10XV3 filtered barcodes in this repo that you can use if that is the
chemistry for your samples. The command is then as follows:

```
velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf bam_file.bam annotation.gtf
```
Once this step has finished and your loom file is generated, we will go ahead and use anndata to import our loom file and make the necessary adjustments/additions. If you'd like a more interactive visualization for the rest of this guide,
an interactive jupyter notebook is provided here.

