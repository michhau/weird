#**W**ind Field **E**stimation using thermal **IR** **D**ata

- Data from IR-camera must be available as 3D array stored as variable 'irdata' in an .NetCDF4-file. Dimensions:
  - 1: height of image
  - 2: width of image
  - 3: temporal sequence

##

- 'weird_preproc.jl' performs the preprocessing steps as explained in Haugeneder et al. 2022
- 'weird_proc.jl' performs the actual wind field estimation on multiple cores (designed for the use on HYPERION, the high-performance computing facility at the Swiss Federal Institute for Forest, Snow and Landscape Research WSL)
- 'weird_output_stitch.jl' stitches the single output wind field frame into a 3D array and saves the array

- variables can be changed in the files in /args/ that are used for argument parsing
