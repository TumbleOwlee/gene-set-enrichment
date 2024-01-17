# GeneSetEnrichment

[![Singularity Build](https://github.com/TumbleOwlee/gene-set-enrichment/actions/workflows/singularity.yml/badge.svg?branch=main)](https://github.com/TumbleOwlee/gene-set-enrichment/actions/workflows/singularity.yml)

This repository contains the R markdown file and environment to perform GO and KEGG enrichment analysis and generate plots based on DESeq2 input data. The concepts are based on the [GeneSetEnrichment Guide](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/) from the [New York University](https://www.nyu.edu/) with various modification.

# Quickstart

The repository provides the R markdown file `./src/gene_set_enrichment.Rmd`. Just upload it to your R studio instance and execute the chunks. But keep in mind to update the configuration since the defaults are for my use case.

# Configuration

The configuration consists of multiple parameters. Each one has to be carefully chosen to work with your provided data and to generate the correct results.

**Parameters:**
* `quiet`: If true, this option just disables excessive logging to prevent terminal flood.
* `organism.annotation`: This is the database to use for annotation retrieval. It's used by all functions provided by `clusterProfiler`.
* `deseq.file`: Path to the DESeq2 file.
* `annotated.file`: Path to the output file of eggNOG. Since we have to know the GO terms for the locus tags, we utilize the eggNOG service to provide them.
* `locus.tag.map`: Path to the locus tag map. E.g. downloaded from [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/299/455/GCF_000299455.1_ASM29945v1/) for `esl`. Check out [this](https://www.genome.jp/kegg/catalog/org_list.html) for other links.
* `key.type`: Key type to use for annotation retrieval.
* `ont`: The ontology to use. Can be either "BP", "MF", "CC" or "ALL".
* `min.gs.size`: See `?gseGO` and `?gseKEGG`.
* `max.gs.size`: See `?gseGO` and `?gseKEGG`.
* `pvalue.cutoff`: See `?gseGO` and `?gseKEGG`.
* `padjust.method`: See `?gseGO` and `?gseKEGG`.
* `kegg.code`: The KEGG code as listed [here](https://www.genome.jp/kegg/catalog/org_list.html)
* `show.categories`: Used to limit the number of elements in plots.
* `pathway.id`: The pathway id for the pathview (see `?pathview`).
* `pathway.gene.idtype`: The gene id type (see `?pathview`).

# Requirements

To simplify the scripts and minimize code duplication, the package `ifutil` is used heavily. You can find it [here](https://github.com/tumbleowlee/ifutil). You can install it easily by executing the following command in your `R` shell.

```R
> remotes::install_github("tumbleowlee/ifhutil")
```

Afterwards you should be able to load the library with `library(ifhutil)`. Documentation for each provided function is also added, thus you can check out `?ifh.<function>` easily.

# Environment

You can execute the scripts in your local environment or create a singularity or docker based environment providing all required `R` libraries.

## Singularity

If you have `singularity` installed, this repository provides the recipe to create the necessaary `r-lang` environment that contains all necessary packages. Please refer to [the Singularity guide](https://docs.sylabs.io/guides/2.6/user-guide/installation.html) for installation instruction if `singularity` is currently missing. You can create the `renv.sif` file with `singularity` as shown below.

```bash
sudo singularity build renv.sif singularity/renv.def
```

The build will take some time since multiple `r-lang` packages have to compile. By default the compilation is executed by `make` using a single core. But you can also assign more cores and thus speed up the build time. To achieve this, just edit the `sed` line in `renv.def` and replace the `-j 1` with the core count (e.g. `-j 4` to assign 4 cores or `-j $(nproc --all)` to assign all available cores).

Afterwards you will have a `renv.sif` file that contains the R environment. You can execute any command in this environment by executing `./renv.sif <command>`. To get to the interactive R shell, just run `./renv.sif R` and use it as if it was installed locally.

Alternatively you can check out the `Actions` tab since the configured workflow will perform a build whenever the file `./singularity/renv.def` changes. This way you can simply download the attached file and unpack the `renv.sif` image. Make sure to look for any workflow that took more than 5 minutes since a workflow is started by any commit, but the image is only created if the definition file changes as part of the commit.

## Docker

If you prefer to use a docker container, you can simply use the `./singularity/singularity2docker.sh` script. Just execute the command below. The script is provided by [SingularityHub](https://github.com/singularityhub/singularity2docker).

```bash
./singularity/singularity2docker.sh -n renv:latest renv.sif
```

Afterwards you will have an image in your docker image list. To get the image id execute the command below.

```bash
bash> docker image ls
REPOSITORY                  TAG         IMAGE ID       CREATED          SIZE
renv                        latest      e62ceaad5db5   11 minutes ago   2.3GB
```

Using the `IMAGE ID` (e.g. `e62ceaad5db5`) you can run a container and e.g. start `R` in interactive mode.

```bash
      docker run -it -v .:/home e62ceaad5db5 bash -c 'cd /home/ && R'
             |-| |-| |--------| |----------| |-----| |--------------|
              |   |       |          |          |           |
  start container |       |       image id      |   change directory into /home/ and run R
                  |       |                     |
       run interactively  |            execute given command
                          |          (used to combine cd and R)
                          |
        mount . directory as /home/ in container
```

## Resources

* List of available KEGG codes: [https://www.genome.jp/kegg/catalog/org_list.html](https://www.genome.jp/kegg/catalog/org_list.html)
* List of available organism annotation databases: [https://bioconductor.org/packages/3.18/BiocViews.html#___OrgDb](https://bioconductor.org/packages/3.18/BiocViews.html#___OrgDb)

