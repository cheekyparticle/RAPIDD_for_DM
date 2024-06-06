# RAPIDD_for_DM

Things you need for the code to compile 

- gsl libraries
- pkg-config

1) Go to parent directory and type 

$ sh instructions_cmake.sh

Then these lines 
export PYTHONPATH="<path_to_directory_of_RAPIDD.so_file>:$PYTHONPATH"
export LD_LIBRARY_PATH=<path_to_directory_of_RAPIDD.so_file>:$LD_LIBRARY_PATH

or put them in your ~/.bashrc file.
