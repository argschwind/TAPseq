## Check that all required 3rd party programs are installed

#' Check if required tools are installed
#'
#' Check if requried software is installed and return paths to executables (if found).
#'
#' @return A \code{list} containing paths to required tools, if found, else NA.
#' @keywords internal
check_tool_installation <- function() {

  tools <- list()

  # get paths to required executables
  tools[["primer3_core"]] <- Sys.which("primer3_core")
  tools[["makeblastdb"]] <- Sys.which("makeblastdb")
  tools[["blastn"]] <- Sys.which("blastn")

  # NA if a tool is not found
  tools[tools == ""] <- NA

  # return list with path to tools
  return(tools)

}

#' Actions to perform when package is attached
#'
#' Check for required tools and set tool paths as package options.
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {

  # check for required tools
  tools <- check_tool_installation()

  # raise warning if some tools are not found, else print message with paths to used tools
  na_tools <- names(tools[is.na(tools)])
  if (length(na_tools) > 0) {

    warning("Following required software is not installed or not in PATH:",
            "\n\n\t", paste0(na_tools, collapse = "\n\t"),
            "\n\nPlease install these tools before trying to use this package!",
            call. = FALSE)

  }else{

    # create printable string listing tool names and path to tools
    tools_print <- lapply(1:length(tools), FUN = function(x) {
      paste(names(tools)[x], tools[[x]] , sep = ": ")
    })
    tools_print <- paste(unlist(tools_print), collapse = "\n")

    # print start up message providing name and paths of used tools
    packageStartupMessage("\nUsing the following tools:\n", tools_print, "\n")

  }

  # create tool paths for package options
  tool_opts <- tools
  names(tool_opts) <- paste0("TAPseq.", names(tool_opts))

  ### TO DO: implement better! ---------------------------------------------------------------------
  # add default path to primer3 thermodynamic parameters
  therm_params <- paste0(dirname(tools[["primer3_core"]]), "/primer3_config/")
  tool_opts[["TAPseq.thermodynamic_params_path"]] <- therm_params
  ### ----------------------------------------------------------------------------------------------

  # create package options
  op <- options()
  op.primer3 <- tool_opts

  # set options
  toset <- !(names(op.primer3) %in% names(op))
  if(any(toset)) options(op.primer3[toset])

  invisible()

}
