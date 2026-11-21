# Changes by Claude, 2026-11-20

Two commits on top of `dev` (93bc69583):

1. "Fix bugs found by code review and unit tests": bug fixes from code reviews
   of sys, math, lib, maos, skyc, tools, cuda and mex, from the new unit tests,
   from a documentation audit, and from macOS, UndefinedBehaviorSanitizer and
   valgrind runs; the code builds as C++; the cmocka unit tests and their
   report in `make doc`; a CUDA build in the Linux CI workflow.
2. "Update source code documentation": comment-only changes, the user manual
   in docs/ checked against the code, the doxygen page layout, and CLAUDE.md.
   The code tokens are unchanged.

## Impact on simulation results {#claude_results}

`scripts/maos_check.sh 5` on GPU and CPU is within run-to-run noise (MCAO
cases vary by about 5% between runs), so the default configurations are not
measurably changed. The fixes below change results only in the listed cases.

### maos (CPU and GPU)
| Change | When it matters |
|---|---|
| Reconstructor updates from Cn2 estimation or dithering (`postproc_isim`) were applied only with PARALLEL=2 | `cn2.tomo=1` or dithering with PARALLEL=0 or 1 (default). Example: 40 steps with `cn2.tomo=1` went from 155.4 nm (same as `cn2.tomo=0`) to 164.2 nm |
| `mkh` bilinear interpolation now renormalizes edge points and applies alpha, as `prop()` does; the cubic version accumulates its weight sum | Edge points of HXW/HA/HXF and the PWFS/mkdtf `si` |
| `Tomo_iprop_do` traced forward instead of the transpose | CPU tomography with `tomo.square=1` and `dbg.tomo_hxw=0` |
| The periodic Laplacian regularization used the x size of the grid for both dimensions | `tomo.square=1` with a non-square reconstruction grid |
| `tomo.cxxalg=1/2`: invpsd and fractal take the full `cxxscale*TOMOSCALE`, not its square root | `tomo.cxxalg=1` or 2 |
| `clipdm_ia` clip was lost | `sim.dmclipia` with `fit.square=0` |
| Automatic `sim.eplo` (`servo_optim_margin`) was discarded, always 0.5 | Automatic low order gain with `dtrat_lo>0` |
| `readcfg_siglev` skipped the wavelength index for an empty powfs | Later powfs got wrong wvlwts, throughput, mag and magbkgrnd when an empty powfs came first |
| `readcfg_fov` left (0,0) directions for clipped points | Rectangular or square FoV with clip=1 |
| A single string for all powfs/WFS (e.g. `powfs.qe=["f"]`) was ignored | String arrays given with one entry |
| Transposed non-square config arrays (`` [...]' ``) were parsed wrong | No shipped config uses one |
| Matched filter: residual background added twice, and its terms not scaled by kpx/kpy | `powfs.bkgrndc!=1` or `bkgrnd2!=bkgrnd2c` |
| CPU `addnoise` multiplied by QE instead of dividing (matches GPU now) | CPU physical optics with `powfs.qe` set |
| `i0save`: CPU doubled the images (i0save=2) or saved nothing (i0save=1) | `powfs.i0save` |
| `pixoff` x/y pixel sizes swapped | `powfs.pixoffx/y` with non-square pixels |
| CORR/CORRS gradients used WFS 0's reference image | Correlation centroiding with several WFS |
| atmscale applied in CPU evaluation with `atm.dtrat>0` | `atm.dtrat>0` on CPU |
| Dithering drift buffer overflowed and only the first dithered powfs was used | Several dithered powfs |
| Dither amplitude reduction never happened (integer division) | PWFS multi-mode dithering (`powfs.dither<-1`, `recon.modal=1`) |
| `ahst_wt=4` used the wrong temporary and column bound | `recon.ahst_wt=4` with 2 or more DMs |
| Focus tracking NEA used the physical WFS index; other GLAO index fixes (petal, twfs, dither, ploc, GA, save.gcov) | GLAO |
| SH subaperture selection kept the whole padded grid | `powfs.saat=0` only (e.g. 60x60 goes from 3844 to 2904 subapertures); all default powfs use saat>0 |
| `sa_reduce` summed `saasum` wrong | SH WFS with `sigmatch=2` when subapertures are dropped |
| `pistatout` accumulated before `pistatstart` (and crashed with `pistatstc=0`) | `pistatout` |
| Threaded MVM assembly shared FFT buffers (race) | MVM with `tomo.cxxalg=1/2` and several threads |

### GPU only
| Change | When it matters |
|---|---|
| MOAO DM wiring for WFS | MOAO on WFS |
| PWFS interaction matrix passed a displacement as a direction | Altitude DMs with an off-axis PWFS |
| `svd_pow` truncation used host memory and was never applied | `svd_pow` with NX>500 |
| Uplink turbulence traced at `llt->misreg`, not the launch position (GPU and CPU now agree to 2e-12) | Side-launched or off-axis LLT |
| Uplink tip/tilt added only when both x and y were nonzero | `llt_ws` (all on x) |
| `pistatout` overwrote the FSM command; idealfsm also added wind shake | `pistatout` or idealfsm |
| Ideal NGS modes applied to `powfs.lo` instead of `powfs.skip` | `tomo.ahst_idealngs` |
| Background column test inverted | `powfs.bkgrnd`/`bkgrndc` from files |
| PWFS with several wavelengths gave 2*w0+w1 | PWFS with nwvl>1 |
| DM misregistration and MOAO embedding with ghost points; missing `__syncthreads()` in tip/tilt evaluation | WFS misregistration; MOAO with `fit.square=0` |
| Partial atmosphere screen path (and OOM fallback) was dead; it now matches the full screen to 6e-7 | When the full screen does not fit in GPU memory |

### lib and math (shared by maos and skyc)
| Change | When it matters |
|---|---|
| CPU PWFS modulation shift dropped PSF rows/columns; `pywfs.psx/psy` operator precedence | CPU PWFS intensities, optical gain (also at setup with GPU), CPU interaction matrix |
| `cn2est` y frequency grid in the wrong order | Cn2 estimation with baselines not along x, y or 45 degrees |
| `X(interp1)` on non-uniform grids; `interp1linear`/`interp1log` at the last grid point | PSD interpolation and skyc nonlinearity lookups |
| `X(add)` with `ac==0` now scales B by `bc` | servo integrator with gain `ep==0` |
| `loc_calc_mod` used uninitialized memory | Single precision builds (perfevl mode RMS) |
| `X(circle)` margin units, `create_metapupil` and loc/map functions with dx!=dy | dx!=1 or dx!=dy |
| Petal screen edge weights | Petal screens |
| Type II `servo_optim` rejected every gain by round off | skyc with servo type 2 |
| `spatial_psd` did not zero a reused output | Repeated screen generation with the spectrum cache |
| `gensei` i0scale lost the wavelength weighting | `i0scale` with several wavelengths |
| Bessel K and I go from about 1e-9 to 1e-14 relative error | von Karman covariance changes by about 1e-9 |
| `fresnel_prop` angular spectrum method built an aliased kernel | Fresnel propagation method 2 |
| Complex random numbers depended on the compiler's argument order (gcc swapped real and imaginary parts) | Complex `randn`/`randu`; gcc results change |
| Complex Gram-Schmidt now uses the conjugate inner product | Complex `X(gramschmidt)` |
| Complex `X(spmulsp2)` swapped 't' and 'c' | FDPCG is updated to match and unchanged |

### skyc
| Change | When it matters |
|---|---|
| NEA kept only the last wavelength | `skyc.neaaniso=1` |
| Wind shake added twice to phyres and to "ATM NGS" | `psd_ws` with `addws=0` |
| Wavelength weight used the wrong index | Wavelengths not matched in order |
| Nonlinearity read sanea as interleaved x/y | `skyc.neanonlin=1` (delete the `powfs%d_nonlin` cache) |
| Detector half a pixel off for odd `pixpsa` | Odd `skyc.pixpsa` |
| minrad, multirate SNR rate, pistat curve for ig==0 | `skyc.minrad>0`, `skyc.multirate=1`, `skyc.neanonlin=1` |
| `sortstars` read the next star's H band with one band | Star fields with one band |

### Configuration and file format
- `dbg.wfs_iac` is read as a real number.
- The mcao and keck examples that failed the `powfs.astscale` check have the
  base asterism halved and `astscale` doubled (guide stars unchanged). The
  check skips powfs whose WFS are all on axis, which fixes the strata
  examples. skyc_ttf.conf uses current keys.
- Single precision complex FITS files had re and im exchanged for other
  readers (astropy, cfitsio). **zmat FITS files written by earlier versions
  now read back with re and im exchanged.** `.bin` files are not affected.
- FITS: long keywords use CONTINUE cards; int8 data padding; files of 2 GB
  or more; DM grid maps keep iac, dratio and aoi.

## Other changes (no effect on results)
- Sizes, element indices and offsets are `long` instead of `int`
  (clang `-Wshorten-64-to-32` reports nothing outside cuda), so arrays
  larger than 2^31 elements no longer overflow in ray tracing, matrix and
  sparse loops, operator assembly, WFS OPD offsets and MVM transfers. Small
  counts (WFS, DM, wavelengths, threads) stay `int` with explicit casts.
  Header changes: `pts_t.nxsa/nysa`, `cn2est_t.nsa`, `dtf_t` sizes,
  `atm.nx/ny/nxnmax`, skyc `nstep`, and the size or index arguments of
  `gensepsf`, `X(embed_wvf)`, `X(fft2partial)`, `zfarr_push`,
  `mklaplacian_map`, `nea_check`, `maxapriori`, `ngsmod2wvf`,
  `petal_connect` and `wrap_seq`.
- Robustness fixes: crashes, out of bounds access, leaks, NULL
  checks, return after `error()` (which is recoverable in tests and mex),
  zstd seeking, shared memory files, sockets, thread pool waits, and
  scheduler/monitor/drawdaemon issues.
- Security: an empty `~/.aos/htcookie` bypassed the scheduler's web
  authentication; htpasswd/htcookie were world readable.
- The scheduler connection regression from 88beaf2da (non-blocking
  `connect_port`) is fixed.
- C++ build (`CC=g++`): configure library checks, `ADDR()` for compound
  literals, no function pointer casts (UBSan), `const void*` qsort
  comparators.
- macOS: process information (memory, load, CPU, launch time) is implemented;
  `SO_NOSIGPIPE` on sockets; `realloc_maos(p,0)` frees on every platform;
  `mymalloc`/`myrealloc` parenthesize the count.
- configure: `--disable-threadlib` with OpenBLAS, the CUDA and mex checks.
- Python: maos_client.py protocol fixes, maos_result.py with NumPy >= 1.24.
- CI: the Linux gcc job builds with CUDA (compiled, not run).
- `maos -h` and the documentation describe the current options.

## Unit tests
- The old test_* programs are replaced by 34 cmocka programs with 324 test
  cases that assert results and use no temporary files. Line coverage of
  sys, math and lib is 86% (was 54%).
- `make check` writes docs/unittest_status.md in the build folder; `make doc`
  includes it, with coverage when configured with `--enable-gcov`.

## Verification {#claude_verification}
- No warnings in C with CUDA, C++ without CUDA, and C++ without OpenMP.
- `make check` passes in these builds and with clang
  `-fsanitize=undefined` (no reports). valgrind is clean, except
  cmocka_matops, cmocka_recon and cmocka_servo, which did not finish within
  25 minutes.
- `maos_check.sh 5` on GPU and CPU is within noise, and targeted runs were
  done for the GPU fixes, Cn2 tomography, dithering, GLAO and pistatout.
- Not run: skyc, the scheduler and the GTK tools (compiled only), mex (no
  Matlab), macOS and BSD code (not compiled here).

## Not changed {#claude_not_changed}
- `X(minsearch)`: zero components of the initial guess give a degenerate
  simplex.
- `mkw_annular`: the inner and outer normalizations do not cancel exactly;
  fixing it changes the fitting weights of annular pupils.
- `wfspupmask` in maos/powfs.c is unused; its radius is always d/2.
- skyc/powfs.c calls `dinterp1` with a 2-D `xnew`, which it rejects.
- maos/petal_utils.c keeps its static thread variables; tools/memdbg.c does
  not compile as C++.
