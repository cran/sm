.First.lib <- function(library, pkg)
{
    if(version$major == 0 && version$minor < 64)
        stop("This package requires R 0.64 or later")
    library.dynam("sm", pkg, library)
    assign(".sm.home", file.path(library, pkg), envir = .GlobalEnv)
    cat("Library `sm'; Copyright (C) 1997 A.W.Bowman & A.Azzalini\n")
    cat("See `Licence' for details and conditions for use\n")
    invisible()
}

isMatrix <- function(x) length(dim(x)) == 2


sm.script <- function(name, path)
{
    if (missing(path)) path <- file.path(.sm.home, "scripts")
    else path <- as.character(substitute(path))
    if (missing(name)) {
      file.show(file.path(path, Index.doc))
    } else {
        name <- as.character(substitute(name))
        if(length(name) == 3 && name[1] == "<-")
            name <- paste(name[2:3], collapse="_")
        source(file.path(path, paste(name, ".q", sep = "")))
    }
    invisible()
}
