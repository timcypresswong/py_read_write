# py_read_write
Python scripts for reading and writing files, mainly for the use of quantum chemistry (Gaussian log file)

### Usage

After cloning to local, it is recommended to export the path of directory of README.md to environment variable, e.g. add the following command to the `~/.bashrc`

`export PYREADPATH=/path/to/py_read_write`

then use it to export the xyz format of the Gaussian log files. For instance,

`python ${PYREADPATH}/clever_split.py gaussianfile.log` to output the last structure of the gaussian.

`python ${PYREADPATH}/clever_split.py gaussianfile.log --all_frame True` to output the xyz trajectory of the gaussian (e.g. job type optimization)

`python ${PYREADPATH}/clever_split.py gaussianfile.log --frame_number 3` to output the 3rd structure of the gaussian (0-based)
