## The Molecular Design Toolkit Package

This is the main MDT package. Except for a few core functions, most of the functionality here is organized into subpackages. *If you're not a developer* you don't need to worry which function is in which subpackage - after you run `import moldesign`, all user-facing APIs are available at the top level of the `moldesign` namespace.

The complete internal APIs are exposed by the subpackages (generally).

Modules in this package are general-purpose utilities:
 * The physical unit system (in `units.py`);
 * generic molecular and biomolecular data (in `data.py`);
 * general tools for working with molecules (in `helpers.py`);
 * a standardized set of parameters for interacting with molecular modeling techniques (in `parameters.py`)