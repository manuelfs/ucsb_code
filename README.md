ucsb_code
==========

Utilities to handle cfA ntuples.

Code may be compiled with ./compile.sh. Executables and scripts are stored in scripts directory 
and are all intended to be run from the root ucsb_code directory (i.e. not from within the scripts directory).

Repository is setup with the assumption that (assuming the repository is locally called ucsb_code) 
there is a directory ucsb_code/../data. This directory is intended to store cfA-style .root files 
containing data from cfA needed in the analysis. The directory can be filled (but not created) by running 
./scripts/make_skims.sh.
