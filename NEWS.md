# SeuratIntegrate (development version)

* Speed up Dijkstra's algorithm-like used in `ExpandNeighbours` with a new c++
implementation

* The most suited matrix format is automatically chosen for corrected counts
output by integration methods (should be dense matrix most of the time)

* Add support for scGraph metric (https://doi.org/10.1101/2024.04.02.587824)

* Improved speed of `CreateIntegrationGroups` for non-SCT assay in unambiguous cases


# SeuratIntegrate 0.4.0

* Initial public release
