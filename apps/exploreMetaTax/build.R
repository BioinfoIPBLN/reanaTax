#!/usr/bin/env Rscript
#
# Rebuild the self-contained (shinylive) exploreMetaTax bundle in `site/`.
#
#   Rscript apps/exploreMetaTax/build.R [destdir]
#
# Run it with an R that has the `shinylive` package and every package the app
# imports, e.g. at the FGCZ:
#
#   /misc/ngseq12/miniforge3/envs/ps_jlruiz/bin/Rscript apps/exploreMetaTax/build.R
#
# shinylive does not run the app: it reads the source, works out the package
# list, and downloads a WebAssembly build of each from https://repo.r-wasm.org.
# The local installs are consulted only for their metadata (name and version),
# which is what the LifemapR workaround below is about.
#
# Expect ~180 MB and a few minutes on a cold cache.

args <- commandArgs(trailingOnly = TRUE)
app_dir <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "app")
dest_dir <- if (length(args) >= 1) args[[1]] else file.path(dirname(app_dir), "site")

stopifnot(dir.exists(app_dir))
if (!requireNamespace("shinylive", quietly = TRUE)) {
  stop("The 'shinylive' package is required to build the bundle.", call. = FALSE)
}

# ── LifemapR workaround ──────────────────────────────────────────────────────
# LifemapR is normally installed from GitHub, so its metadata carries
# RemoteType=github and RemoteRef=HEAD. shinylive takes that as "fetch a
# WebAssembly build from this repo's GitHub releases", finds no release tagged
# HEAD, and aborts - even though LifemapR is in the webR CRAN mirror and works
# fine from there.
#
# So: copy the installed package into a throwaway library, drop the Remote*
# fields, and put that library first. Note the fields have to be removed from
# Meta/package.rds, NOT from DESCRIPTION - packageDescription() reads the
# former for an installed package, so editing only the latter silently changes
# nothing. Nothing outside tempdir() is touched.
shim_lib <- file.path(tempdir(), "shinylive_shim_lib")
if (requireNamespace("LifemapR", quietly = TRUE)) {
  dir.create(shim_lib, recursive = TRUE, showWarnings = FALSE)
  src <- find.package("LifemapR")
  if (!identical(normalizePath(dirname(src)), normalizePath(shim_lib))) {
    file.copy(src, shim_lib, recursive = TRUE)
    meta_file <- file.path(shim_lib, "LifemapR", "Meta", "package.rds")
    if (file.exists(meta_file)) {
      meta <- readRDS(meta_file)
      meta$DESCRIPTION <- meta$DESCRIPTION[!grepl("^Remote", names(meta$DESCRIPTION))]
      meta$DESCRIPTION["Repository"] <- "CRAN"
      saveRDS(meta, meta_file)
    }
    desc_file <- file.path(shim_lib, "LifemapR", "DESCRIPTION")
    if (file.exists(desc_file)) {
      writeLines(grep("^Remote", readLines(desc_file), value = TRUE, invert = TRUE), desc_file)
    }
    .libPaths(c(shim_lib, .libPaths()))
    message("Using a de-GitHub-ified copy of LifemapR from ", shim_lib)
  }
}

message("Exporting ", app_dir, " -> ", dest_dir)
unlink(dest_dir, recursive = TRUE)
shinylive::export(appdir = app_dir, destdir = dest_dir)

message("\nDone. Serve it with any static web server, e.g.:")
message("  Rscript -e 'httpuv::runStaticServer(\"", dest_dir, "\", port = 8080)'")
message("Opening index.html straight from disk (file://) does NOT work - the")
message("service worker webR relies on is only served over http(s).")
