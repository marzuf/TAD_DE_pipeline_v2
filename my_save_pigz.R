pigz_exec_path <- file.path(setDir, "/mnt/ed4/marie/pigz-2.4/pigz")

# function code taken from fastSave package !

my_save.pigz <- function (..., list = character(), file = stop("'file' must be specified"), 
          pigz_exec_path = "",
          envir = parent.frame(), n.cores = 4, eval.promises = TRUE, 
          precheck = TRUE) 
{
  stopifnot(file.exists(pigz_exec_path))
  if (.Platform$OS.type == "unix") {
    if (system(paste0("command -v ", pigz_exec_path), wait = T, ignore.stdout = TRUE, 
               ignore.stderr = TRUE) > 0) {
      stop("The pigz command is not available on this system!")
    }
  }
  else {
    stop("Platform is not a unix system!")
  }
  if (!is.numeric(n.cores)) 
    stop("'n.cores' mut be numeric")
  names <- as.character(substitute(list(...)))[-1L]
  if (missing(list) && !length(names)) 
    warning("nothing specified to be save()d")
  list <- c(list, names)
  if (precheck) {
    ok <- vapply(list, exists, NA, envir = envir)
    if (!all(ok)) {
      n <- sum(!ok)
      stop(sprintf(ngettext(n, "object %s not found", "objects %s not found"), 
                   paste(sQuote(list[!ok]), collapse = ", ")), domain = NA)
    }
  }
  on.exit(close(con))
  n.cores.arg <- paste0(" --processes ", n.cores)
  con <- pipe(paste0(pigz_exec_path, " ", n.cores.arg, " > ", file))
  save(..., list = list, file = con, ascii = FALSE, version = 2, 
       envir = envir, eval.promises = eval.promises, precheck = precheck)
}
# x = 2
# save.pigz(x, file="x.Rdata")
# myx=5
# my_save.pigz(myx, pigz_exec_path = pigz_exec_path, file="myx.Rdata")
# x <- get(load("x.Rdata"))
# x
# myx <- get(load("myx.Rdata"))
# myx