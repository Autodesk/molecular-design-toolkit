## 0.3
##### Features:
- Visualizer displays for force field setup errors and warnings
- Added `viewer.color_by` to set colors in 3D representations using colormaps and categorical data
- Made most objects serializable via pickle
- Better markdown displays for molecules, biounits, and atoms


##### Bug fixes / backend:
 - Drive trajectory animations from Python to prevent overloading browser memory
 - Get the `Deploy` bootstrapping script working consistently
 - Workaround problems with pickling numpy index arrays
 - Overhaul and simplify selection UI classes
 - Remove "frame_interval" option from openmm minimizations - it prevented full convergence
 - Improve performance by overriding how arrays with units are sliced
 - Update pint dependency to 0.7.2
 - Workaround for pint comparison bug (see https://github.com/hgrecco/pint/issues/351)
 - Formalize the way we name atoms, biounits, and molecules (see buckyball/README.md)

 
