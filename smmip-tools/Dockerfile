# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_19-R-4.4.1

# Update apt-get
RUN apt-get update \
	## Install the python package tensorflow
	# && pip install tensorflow \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install required Bioconductor package
RUN R -e 'BiocManager::install(c("optparse", "data.table", "dplyr", "parallel", "Rsamtools", "IRanges", "ggplot2", "RColorBrewer", "pwalign"))'

RUN git clone https://github.com/jdemeul/smMIP-tools.git
