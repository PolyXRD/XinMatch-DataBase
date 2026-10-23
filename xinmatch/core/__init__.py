"""XinMatch computation kernel.

All modules here are **GUI-free**: they consume and produce plain data
structures (numpy arrays, dataclasses) so they can be unit-tested in
isolation and driven from either the CLI or the UI.

Subpackages
-----------
model       Data classes (Pattern, Peak, Entry, Phase, Restraints, Project)
io          Diffraction data file I/O (33 formats) + CIF + exporters
preprocess  Smoothing, background, Kα2 stripping, zero point, displacement
peak        Peak search, profile fitting, peak shapes, peak editing
db          Reference database (SQLite over COD) + user database
search      FoM search-match, PFSM, restraints, candidate list
index       Powder indexing (Treor / Dicvol)
refine      Le Bail / Rietveld refinement (incl. GSAS-II / FullProf adapters)
quant       Quantitative analysis (RIR, internal standard, Rietveld)
micro       Crystallite size, microstrain, degree of crystallinity
report      Report data model and rendering
"""
