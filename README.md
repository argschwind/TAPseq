# TAPseq
An R-package to design PCR primers for TAP-seq.

## Installation

This package requires local installations of Primer3 and BLASTn. TAPseq was developed and tested
using Primer3 v.2.5.0 and blastn v.2.6.0. It's strongly suggested to use Primer3 >= 2.5.0! Earlier
versions require a primer3_config directory, which needs to be provided whenever calling functions
interacting with Primer3. Source code and installation instructions can be found under:

Primer3: <https://github.com/primer3-org/primer3/releases>  
BLASTn: <https://www.ncbi.nlm.nih.gov/books/NBK279690/>

Please install these tools first and add them to your PATH. If you don't want to add the tools to
your "global" PATH, you can add the following code to an ~/.Rprofile file. This should add the tools
to your PATH in R whenever you start a new session.
```
Sys.setenv(PATH = paste("/full/path/to/primer3-x.x.x/src", Sys.getenv("PATH"), sep = ":"))

Sys.setenv(PATH = paste("/full/path/to/blast+/ncbi-blast-x.x.x+/bin", Sys.getenv("PATH"), 
                        sep = ":"))
```

The R-package itself can be installed from source by using the devtools package. This also allows
building vignettes.
```
# minimal install without vignettes
devtools::install_github("argschwind/TAPseq")

# install with building vignettes and also installing suggested dependencies
devtools::install_github("argschwind/TAPseq", build_vignettes = TRUE, dependencies = TRUE)
```

Note: If missing Bioconductor dependencies fail to install from GitHub, please install them from
Bioconductor (see: <http://bioconductor.org/install/>).

## Examples
An example of the TAPseq primer design workflow can be found in a vignette. To view the vignette,
run the following command (assuming vignettes were built when the package was installed).
```
vignette("tapseq_primer_design", package = "TAPseq")
```

Examples of how to select and evaluate target genes to identify cell populations can be found in
a separate vignette. This requires that the additional dependencies are installedl, which should be
the case if the package was installed with building vignettes and suggested dependencies.
```
vignette("tapseq_target_genes", package = "TAPseq")
```
