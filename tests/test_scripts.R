## Note: R CMD check may run these scripts from an installed package
scripts <- list.files(system.file("scripts", package = "sm"), ".*\.q$")
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
    set.seed(123)
    source(system.file("scripts", z, package = "sm"), echo=TRUE)
    rm(list = ls(all = TRUE))
}
