# exploreMetaTax — self-contained (shinylive) build

An offline, serverless build of the **exploreMetaTax** Shiny app, for exploring
the taxonomic tables this pipeline produces. It runs entirely in the browser:
[shinylive](https://posit-dev.github.io/r-shinylive/) ships R itself as
WebAssembly, so there is no R server, no ShinyProxy, and no upload of your data
to anything. Everything you load stays on your machine.

> [!IMPORTANT]
> The hosted app at <https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax> is
> likely to be **more up to date and considerably faster** than this bundle.
> Use it whenever you can, and fall back to this copy when the hosted instance
> is unreachable, when your data must not leave your machine, or when you want
> a build pinned alongside a particular pipeline run. If this bundle feels slow
> or misbehaves, try the hosted app before reporting a problem — R in
> WebAssembly is several times slower than native R, and large datasets feel it.

## Layout

| Path       | What it is                                                                 |
| ---------- | -------------------------------------------------------------------------- |
| `app/`     | The app source the bundle is built from: `app.R` plus its data file.        |
| `site/`    | The generated static site (~180 MB). Not tracked by git — see below.        |
| `build.R`  | Regenerates `site/` from `app/`.                                            |

## Running it

Any static web server will do. From the pipeline root:

```bash
Rscript -e 'httpuv::runStaticServer("apps/exploreMetaTax/site", port = 8080)'
# or
python3 -m http.server 8080 --directory apps/exploreMetaTax/site
```

then open <http://localhost:8080>.

Opening `site/index.html` directly as a `file://` URL **does not work**: webR
loads through a service worker, which browsers only allow over `http(s)`. The
same directory can be published as-is to GitHub Pages or any static host.

The first load fetches the R runtime and ~90 packages from `site/` and is
therefore slow (tens of seconds); afterwards the browser caches them.

## Feeding it pipeline output

The app is upload-driven, which is what makes a serverless build possible. Point
it at the files in your `--outdir`:

- `kraken2/*.kraken2.report.txt` — per-sample Kraken2 reports
- `bracken/*.bracken_S.tsv` — per-sample Bracken tables
- `bracken/bracken_combined_S.txt`, `kraken2/kraken2_combined_report.txt` — the combined tables
- a metadata TSV/CSV of your own, to group and colour samples

## Rebuilding

```bash
/misc/ngseq12/miniforge3/envs/ps_jlruiz/bin/Rscript apps/exploreMetaTax/build.R
```

Needs an R with `shinylive` and every package the app imports. shinylive never
runs the app — it reads the source, resolves the package list, and downloads a
WebAssembly build of each from <https://repo.r-wasm.org>; the local installs are
consulted only for their name and version.

## How this build differs from the server one

`app/app.R` is a copy of `shinyproxy_apps/exploreMetaTax/app.R` with three
changes, all forced by the browser sandbox:

1. **No telemetry.** The `shiny.telemetry` wiring is removed outright — it wrote
   to a database on the FGCZ ShinyProxy host, which does not exist here.
2. **Pure-R Kraken2 lineage parser.** The original compiles a C++ helper with
   `Rcpp::cppFunction()` at startup. There is no C++ compiler in WebAssembly, so
   the build falls back to an R implementation of the same function, verified to
   return identical output. It is slower on very large reports but correct, and
   the C++ path is still used automatically wherever a compiler exists.
3. **No Krona/sunburst tab.** That tab needs
   [`taxplore`](https://github.com/markschl/taxplore), which is GitHub-only and
   has no WebAssembly build, so it cannot ship in a shinylive bundle. The tab
   shows its usual "package not installed" notice. This costs you little here:
   the pipeline already writes standalone interactive Krona charts to `krona/`.

Everything else — upload, filtering, alpha/beta diversity, rarefaction, PCA,
the LifemapR tree, plots and exports — behaves as it does on the server.

The `?data=` gstore URL loader is inert in this build: it requires a ShinyProxy
username and an LDAP lookup, neither of which exists in a browser. Use the
upload tab.
