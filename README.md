# TAPseq
An R-package to design PCR primers for TAP-seq.

## Installation
The package can be installed from source by using the devtools package. This also builds all
vignettes.
```
devtools::install_github("argschwind/TAPseq", build_vignettes = TRUE)
```

This package requires local installations of Primer3 and BLASTn. TAPseq was developed and tested
using Primer3 v.2.5.0 and blastn v.2.6.0. It's strongly suggested to use Primer3 >= 2.5.0! Earlier
versions require a primer3_config directory, which needs to be specified whenever calling functions
interacting with Primer3. Source code and installation instructions can be found under:

Primer3: <https://github.com/primer3-org/primer3/releases>  
BLASTn: <https://www.ncbi.nlm.nih.gov/books/NBK279690/>

Please install these tools first and add them to your PATH. If you don't want to add the tools to
your "global" PATH, you can add the following code to your .Rprofile file. This should add the tools
to your PATH in R whenever you start a new session.
```
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/full/path/to/primer3-x.x.x/src", sep = ":"))

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/full/path/to/blast+/ncbi-blast-x.x.x+/bin",
                        sep = ":"))
```

Alternatively you can specify the paths to 3rd party software as arguments when calling TAPseq
functions (TAPseqInput(), designPrimers(), checkPrimers()).

## Examples
An example of the TAPseq primer design workflow can be found in a vignette. To view the vignette,
run the following command (assuming vignettes were built when the package was installed).
```
vignette("tapseq_primer_design", package = "TAPseq")
```

Examples of how to select and evaluate target genes to identify cell populations can be found in
the examples of following functions. This requires that the suggested dependencies are installed as
well, which should be the case if vignettes were built when installing the package.
```
?selectTargetGenes()
?plotTargetGenes()
```
