# Seurat to RNA-Velocity

This guide will show you how you can use a processed/normalized Seurat object in conjunction with your RNA Velocity analysis. 
This guide will use the following programs: <br>

- [scVelo](https://github.com/theislab/scvelo)
- [anndata](https://icb-anndata.readthedocs-hosted.com/en/stable/)
- [velocyto](http://velocyto.org/)
- Seurat
- Samtools

To start, we'll be generating loom files for every single-cell sample you used in your Seurat analysis using velocyto's <br>```velocyto run```

Note that I've found it better to use ```velocyto run``` for which ever scRNA-seq chemistry you are working with rather than velocyto's "ready-to-use subcommands." 
The only caveat is you'll have to use a filtered barcode file in conjunction with the run command, but you can easily find those from
your sequencing company's website. I've uploaded 10XV2 and 10XV2 filtered barcodes in this directory that you can use if that is the
chemistry for your samples. The command is then as follows:

```
velocyto run_smartseq2 -o OUTPUT -m repeat_msk.gtf -e MyTissue plateX/*/*.bam mm10_annotation.gtf
```
