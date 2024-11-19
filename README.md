# ArgonCubeLUTSim
ArgonCube optical simulations based on Look-up-Table (LUT)

The light propagation model is constructed with this repo: https://github.com/Frappa/ArCubeOptSim
The 2x2 related optical properties are listed in [this sheet](https://docs.google.com/spreadsheets/d/1yGrKdGnJMx8poRzJhDruYQg0A4hAS3w4hp9X0tV6RpM/edit?usp=sharing).

Run the translation script by \
```python3 LUT_root2numpy.py /PATH/TO/GEANT4/LUT -o /OUTPUT```

On NERSC (perlmutter), the light LUT's are stored here: `/global/cfs/cdirs/dune/www/data/2x2/simulation/larndsim_data/light_LUT/`
On s3df, the light LUT's are stored here: `/sdf/data/neutrino/2x2/light_lut/`
