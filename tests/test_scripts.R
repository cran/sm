scripts <- list.files("../sm/scripts", ".*\.q$")
## these need gam
omit <- match(c("trwlgam1.q", "trwlgam2.q", "mackgam.q", "trwlgam3.q",
                "smackgam.q"), scripts)
## these are interactive
omit2 <- match(c("bissell3.q", "dogs.q"), scripts)
scripts <- scripts[-c(omit, omit2)]
library(sm)
if(.Platform$OS.type == "unix") options(pager="cat")
postscript(file="test_scripts.ps")
for(z in scripts) {
    cat("\n============ running script `", z, "' ============\n", sep="")
    source(file.path("../sm/scripts", z), echo=TRUE)
    rm(list = ls(all = TRUE))
}
