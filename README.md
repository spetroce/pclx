https://github.com/spetroce/pclx

# PCLx
PCLx is a header library for performing fundamental operations on point clouds. The motivation was to create something more flexible and efficient than the PCL library. The bulk of PCLx is templated to allow users to use their own point and color types. The 'mio' library is a dependency of PCLx and provdes default types. A main feature of PCLx is it holds only a single source copy of the 3D data (points and their colors) and all functions operate on indices of the source data. 

##core.h
core.h simply includes core-inl.h, which contains the base CPointCloud class used by all the PCLx functions.
