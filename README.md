# GeneSetEnrichment

This repository contains the script and environment to perform GO and KEGG enrichment analysis and generate plots based on DESeq2 input data. The concepts are based on the [GeneSetEnrichment Guide](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/) from the [New York University](https://www.nyu.edu/).

# Environment

You can execute the scripts in your local environment or create a singularity or docker based environment providing all required `R` libraries.

## Singularity

If you have `singularity` installed, this repository provides the recipe to create the necessaary `r-lang` environment that contains all necessary packages. Please refer to [the Singularity guide](https://docs.sylabs.io/guides/2.6/user-guide/installation.html) for installation instruction if `singularity` is currently missing. You can create the `renv.sif` file with `singularity` as shown below.

```bash
cd /path/to/repo/
sudo singularity build renv.sif singularity/renv.def
```

The build will take some time since multiple `r-lang` packages have to compile. By default the compilation is executed by `make` using a single core. But you can also assign more cores and thus speed up the build time. To achieve this, just edit the `sed` line in `renv.def` and replace the `-j 1` with the core count (e.g. `-j 4` to assign 4 cores or `-j $(nproc --all)` to assign all available cores).

Afterwards you will have a `renv.sif` file that contains the R environment. You can execute any command in this environment by executing `./renv.sif <command>`. To get to the interactive R shell, just run `./renv.sif R` and use it as if it was installed locally.

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
