# Base image using r2u which installs all R packages as native Ubuntu binaries
# This completely bypasses C++ compilation and prevents out-of-memory errors
FROM rocker/r2u:jammy

# Install all required packages including Bioconductor via standard install command
# The r2u engine automatically intercepts this and fetches precompiled deb packages
RUN R -e "install.packages(c('Seurat', 'Matrix', 'EnhancedVolcano', 'patchwork', 'ggplot2', 'GEOquery', 'DESeq2', 'clusterProfiler', 'org.Mm.eg.db', 'ashr'))"

# Strict validation step to guarantee clusterProfiler is fully functional
RUN R -e "if(!requireNamespace('clusterProfiler', quietly = TRUE)) stop('Build failed! clusterProfiler is missing.')"

# Set the working directory inside the container
WORKDIR /workspace

# Define default command to execute the analysis script
CMD ["Rscript", "analisis.R"]
