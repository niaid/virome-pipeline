# Using the DiscoVir docker image

The DiscoVir pipeline is available in NIAID's web app for microbiome analysis, [Nephele](https://nephele.niaid.nih.gov).  Nephele publishes docker images for all pipelines in [our AWS ECR gallery](https://gallery.ecr.aws/niaid_nephele/).  Older versions of the pipeline images can also be found in the archive there.

- The [DiscoVir image page](https://gallery.ecr.aws/niaid_nephele/pipeline/nephele_virome) has full instructions for how to download the image and where to mount the external databases as well as some example commands for running the pipeline under the Usage tab.  
- **About the external databases**: While the images contain all the tools for DiscoVir, they do not contain the databases, which can be quite large (hundreds of GBs to TB).  We do not have the resources to host these third party databases for download, so they need to be downloaded and mounted separately in order to use the image to run the pipeline.  We link to instructions for downloading the databases from each of the respective tools below.
- Using the image requires knowledge of using and troubleshooting docker on the command line.   Other users should try out [DiscoVir on the web](https://nephele.niaid.nih.gov/user-guide/pipeline-descriptions/discovir).

- For datasets smaller than 150 gb gzipped in size, we use a machine with 48 cpus and 192 GB memory.
    - The significant memory usage is for iPHoP.  If you are not running that step, you may need less.

- To report any bugs or issues you have with the docker image, use the [Nephele contact form](https://forms.office.com/g/xfMLrmVXDE).

## Databases

All databases should be downloaded inside of a single folder (in our [documentation under usage](https://gallery.ecr.aws/niaid_nephele/pipeline/nephele_virome), we use a folder called *dbs/virome*), and then mounted to the container using `-v /path/to/dbs:/dbs`.  Where the download instructions use the tool itself, you can use the tool in the docker image itself to download the database.

<u>Instructions for downloading each database</u>

- **[CheckV](https://bitbucket.org/berkeleylab/CheckV/src/master/#markdown-header-checkv-database)**
- **[geNomad](https://github.com/apcamargo/genomad/?tab=readme-ov-file#downloading-the-database)**

- [**VirSorter2**](https://github.com/jiarong/VirSorter2?tab=readme-ov-file#download-database-and-dependencies)

- **[DRAM-v](https://github.com/WrightonLabCSU/DRAM?tab=readme-ov-file#i-dont-have-access-to-kegg)**: This database is quite large.  We suggest not downloading the UniRef files which are the largest part, as those are not needed for DRAM-v.
- [**iPHoP**](https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-downloading-iphop-host-database)
- [**Diamond nr (NCBI)**](https://scienceparkstudygroup.github.io/ibed-bioinformatics-page/source/core_tools/ncbi_nr.html#get-the-ncbi-nr-database)