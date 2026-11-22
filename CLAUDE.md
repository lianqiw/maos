# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

MAOS (Multithreaded Adaptive Optics Simulator) is an end-to-end adaptive optics simulator written in C11 (optional CUDA for GPUs). Manual: https://lianqiw.github.io/maos (sources in `docs/*.md`, built with doxygen).

## Build

Autotools, always built **out of tree**. This checkout (`~/work/maos/src`) already has configured build dirs as siblings: `../optim` (default optimized), `../debug` (`--enable-debug`), `../release` (`--enable-release`), `../cpp` (`CC=clang++`, compiles the C as C++; note it is currently configured against `../src_master`, not this tree).

```
./autogen.sh                         # only needed after editing configure.ac / Makefile.am (or fresh clone)
mkdir ../optim && cd ../optim && ../src/configure [--enable-debug|--enable-release|--disable-double|--with-cuda=DIR]
make -j8                             # executables land in <builddir>/bin (maos, skyc, drawdaemon, monitor, aolib.so, ...)
```

Useful configure switches: `--enable-debug`, `--enable-release` (portable build, disables scheduler), `--disable-double` (single precision CPU+GPU), `--disable-scheduler`, `--enable-gcov`, `--with-matlab=DIR`, `--disable-gtk-3/4`. Missing deps (FFTW, BLAS/LAPACK, CHOLMOD, libwebsockets) are auto-downloaded by configure unless `--disable-download`.

## Tests

- Unit tests (`test/cmocka_*.c`, CMocka; include `test/unittest.h`): `make check` in the build dir. Tests must not read or write temporary files; generate inputs in memory. Each run updates the status report `docs/unittest_status.md` in the build dir; `make doc` reruns the tests and puts it (with gcov line coverage of `sys`/`math`/`lib` when configured with `--enable-gcov`) in the manual in place of the placeholder `docs/96_unittest.md`.
- Single test: `make -C test check TESTS=cmocka_mat` (from the build dir), or run `test/cmocka_mat` directly.
- End-to-end regression: from `<builddir>/bin`, run `../../src/scripts/maos_check.sh 5` (D=5 m; also `10`, `30`). It runs many configurations and compares wavefront error against hard-coded `REF_GPU`/`REF_CPU` arrays in the script; update those arrays when a change intentionally alters results. Add `-g-1` to force CPU.
- CI (`.github/workflows/`) builds with gcc and clang (`--enable-release`), runs `make check`, `maos_check.sh 5`, and a macOS build with `CC=clang++` — code must also compile as C++.

## Running

`maos [-c base.conf] [key=value ...] [extra.conf ...] -o outdir`. Defaults come from `config/maos/default.conf` (found via `MAOS_CONFIG_PATH`, `~/.aos/`, or the source tree). `-g-1` disables GPU, `-n N` threads, `-d` detaches. Env vars like `MAOS_LOG_LEVEL`, `MAOS_PARALLEL`, `MAOS_MEM_DEBUG=1` (leak accounting) are documented in `docs/30_run.md`. Sky coverage post-processing is the separate `skyc` program (`config/skyc`).

## Architecture

Layered libraries, each depending only on the ones above it (`SUBDIRS` order in `Makefile.am`):

1. `sys/` – OS layer: memory accounting (`mem.c`), threads/thread pool, sockets, the `.bin` file format and zstd/gzip IO (`bin.c`, `zstdio.c`), scheduler client, per-OS process info.
2. `math/` – numeric containers (`dmat`, `cmat`, `smat`, `zmat`, `lmat`, sparse `dsp`/`csp`, `cell` arrays) with BLAS/FFT/CHOLMOD wrappers. **Type-generic via macros**: `defs.h` defines `X(name)` → `d`/`s`/`c`/`z`/`l` prefixes; `dmath.c`, `smath.c`, `cmath.c`, `zmath.c`, `lmath.c` each set `COMP_*` flags and `#include` the same template sources (`mat.c`, `matmath.c`, `blas.c`, `sp.c`, ...). Edit the template files, not the per-type wrappers; write new functions as `X(foo)` with `T`/`R` types. With `--disable-double`, `d*`-prefixed functions operate on `float` (`real` typedef). Index with `P(a,i,j)`; free with the matching `dfree`/`cellfree`.
3. `lib/` → `libaos`: physics/AO building blocks reusable outside maos (`loc_t`/`pts_t`/`map_t` in `loc.c`/`map.c`, derived from `dmat` and declared in `math/type.h` so generic math routines accept them; turbulence, ray tracing `accphi`, `mkg/mkh/mkz/mkw` operators, PCG, servo, PSD, pyramid WFS, `readcfg` config parser, `draw` plotting client). Functions here take basic types, not maos structs. `lib/aolib.c` builds `bin/aolib.so` for Python/Matlab.
4. `cuda/` (`math`, `recon`, `sim`) – GPU versions of simulation and reconstruction, mirroring the CPU code paths; selected at runtime (`parms->gpu.*`).
5. `maos/` – the simulator. Flow (see `docs/93_maos.md`): `maos.c:main` → `setup_parms()` (`parms.c`, reads every key into `parms_t`) → `maos_setup()` (`setup.c`: `aper`, `powfs`, `recon_setup*`, `surf`, GPU init) → `maos_sim()` (`sim.c`) which per step runs `perfevl`, `wfsgrad`, `reconstruct`, `filter_dm`, serial or parallel depending on `PARALLEL`. Runtime state lives in `sim_t` (`types.h`); static setup in `parms_t` (`parms.h`), `powfs_t`, `recon_t`, `aper_t`.
6. `skyc/` – sky coverage (separate executable, own `parms`).
7. `tools/` – `scheduler` (job queue + websocket server), `monitor` (GTK and `monitor.html/js` browser UI), `drawdaemon` (GTK plotting), `drawres`, `bin2fits`.
8. `mex/`, `scripts/` – Matlab/Python bindings. `scripts/lib2mex.py` and `scripts/lib2py.py` parse headers (`lib/*.h`, `mex/aolib.h`, `maos/parms.h`, `maos/types.h`) to generate `mex/aolib.c` and `scripts/libaos.py`; struct/doc-comment formatting in those headers matters to the parser. `scripts/readbin.py`, `maos_result.py` read outputs.

### Configuration parameters

Every runtime option is a `key=value` in `config/maos/*.conf`, mapped 1:1 onto `parms_t` fields (`wfs.thetax` → `parms->wfs.thetax`) via `READ_INT/READ_DBL/READ_DMAT/...` macros in `maos/parms.c`. A key given on the command line or in a user `.conf` that maos does not read is an error (a warning for `dbg.*`, or with `MAOS_PERMISSIVE=1`); a key read but never set only warns and is taken as 0 (`lib/readcfg.c`). So a new parameter needs both a `READ_*` in `parms.c` (plus any derived values/validation in the `setup_parms_postproc_*` functions) and a documented default in the appropriate `.conf` (e.g. `sim.conf`, `recon.conf`, `powfs_common.conf`). Conf files `include=` each other; later values override, `key+=value` appends; command-line `key=value` beats `.conf` files on the command line, which beat the main `-c` conf, regardless of order.

## Conventions (from `docs/91_basics.md`)

- Allocate structs zeroed (`mycalloc`/`calloc`); free everything (memory accounting reports leaks).
- In/out arguments first, except when there is a main object; prefer a single return.
- `lib/` utilities take basic types; `maos/` functions may take `sim_t`/`parms_t`.
- No non-const static variables; no hard-coded tunable values (make them parameters); no function-pointer casting.
- Mark temporary debugging changes with a `warning(...)` so they are not left in silently.
- Logging uses `dbg`, `info`, `warning`, `error` macros (from `sys/`), filtered by `MAOS_LOG_LEVEL`.
- Each source file carries the GPL copyright header (`scripts/addcopy.sh`).
