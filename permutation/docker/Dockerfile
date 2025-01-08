# Dockerfile for germline:somatic convergence permutations

# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute.  
# Distributed under terms of the GNU GPL v2.0 License

FROM vanallenlab/rlctools:latest
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install necessary python modules
RUN mamba install --update-all \
	-c conda-forge -c bioconda -c r \
	pandas numpy

# Launch bash at runtime
CMD ["/bin/bash"]
