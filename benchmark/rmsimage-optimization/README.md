# rmsimage / bstat optimization

Scripts and notes from the optimization of `Op_rmsimage.rms_mean_map()` and
`functions.bstat()`. See `conversation.md` for the full analysis and results.

All scripts compare two versions of PyBDSF without rebuilding: they copy the
installed package from `venv` (which has the compiled extensions) and replace
its Python sources with those of a git revision (`pkg_orig`) and of the
working tree (`pkg_new`).

| File | Purpose |
| --- | --- |
| `run_comparison.sh <work dir> [rev]` | Runs everything below in order |
| `setup_overlays.sh <work dir> [rev]` | Builds `pkg_orig` (rev, default HEAD) and `pkg_new` (working tree) |
| `analyze_profile.py [prof]` | Prints the parts of `../filter_skymodel.prof` used in the analysis |
| `compare_bstat.py <pkg_orig> <pkg_new> <fits>` | Speed and accuracy of the two `bstat` versions |
| `bench_rms_mean_map.py <pkg> <label> <fits> <out.npz> [--tiled]` | Times `rms_mean_map` and saves the maps |
| `compare_maps.py <ref.npz> <new.npz> [--rtol r]` | Compares maps from two runs |

`run_comparison.sh` also runs `test/tbdsf_process_image.py` with both
versions, compares their outputs at rtol 1e-3 (catalogs included), and times
`../run_pybdsf.py` with both.

Run Python with `-I` (as the scripts do) so the repository's `bdsf/`, which
has no compiled extensions, isn't imported by accident.

Example:

```bash
./run_comparison.sh /tmp/rms-bench HEAD
```
