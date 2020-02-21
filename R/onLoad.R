## Check that all required 3rd party programs are installed

#' Check required tools
#'
#' Check if required software is installed and return paths to executables (if found).
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

# actions to perform when package is LOADED. check for required tools and set tool paths as package
# options
.onLoad <- function(libname, pkgname) {

  # check for required tools
  tools <- check_tool_installation()

  # raise warning if some tools are not found
  na_tools <- names(tools[is.na(tools)])
  if (length(na_tools) > 0) {

    warning("Following software required by TAPseq is not installed or not in PATH:",
            "\n\n\t", paste0(na_tools, collapse = "\n\t"),
            "\n\nPlease install these tools before trying to use this package!",
            call. = FALSE)

  }

  # create tool paths for package options
  tool_opts <- tools
  names(tool_opts) <- paste0("TAPseq.", names(tool_opts))

  # create package options
  op <- options()
  op.tapseq <- tool_opts

  # set options
  toset <- !(names(op.tapseq) %in% names(op))
  if(any(toset)) options(op.tapseq[toset])

  invisible()

}

# print startup message if package is ATTACHED
.onAttach <- function(libname, pkgname) {

  # get all tools used by TAPseq
  tool_ops <- c("TAPseq.primer3_core", "TAPseq.makeblastdb", "TAPseq.blastn")
  tools <- unlist(lapply(tool_ops, FUN = getOption))

  # create printable strings listing tool names and path to tools
  tools_print <- lapply(seq_along(tools), FUN = function(x) {
    paste(names(tools)[x], tools[[x]] , sep = ": ")
  })
  tools_print <- paste(unlist(tools_print), collapse = "\n")

  # print start up message providing name and paths of used tools
  packageStartupMessage("\nTAPseq is using the following tools:\n", tools_print, "\n")

}
