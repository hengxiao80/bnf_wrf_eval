# WRF `p8w` (interface / "level" pressure) formulation

Reference snippet handed to Claude for the task of fixing `crtm.py`'s
`_layer_pressures` (which previously used a geometric mean of adjacent
layer pressures). This is how WRF itself builds `p8w`, the interface
pressure it passes to its physics/radiation packages: a linear-in-eta
interpolation of the layer full pressure `P + PB` onto the staggered
w-levels using WRF's own vertical-stretch weights `FNM`/`FNP`
(`fzm`/`fzp` in the model source), with `psfc` and the model-top pressure
as the two boundary interfaces. Implemented and verified in commit that
followed (matches WRF's `p_hyd_w` / stored `P_HYD` to ~0.1-0.3 Pa mean);
see the `crtm.py` docstring and CLAUDE.md `crtm.py` section for details.

```python
pm = P + PB                       # shape (N, ...)
pw = np.empty((pm.shape[0] + 1,) + pm.shape[1:])

# FNM/FNP index 0 is not used for an interior interface
pw[1:-1] = (
    FNM[1:, None, None] * pm[1:]
    + FNP[1:, None, None] * pm[:-1]
)

pw[0]  = psfc
pw[-1] = ptop_interface
```

Here `psfc` is surface pressure and `ptop_interface` is model top pressure. As you can see, here index 0 is at the surface.