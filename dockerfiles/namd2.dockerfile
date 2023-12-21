# Use the Debian base image
FROM debian:bullseye

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    python \
    tar


# Set the working directory
WORKDIR /

# Download and extract NAMD
RUN wget https://www.ks.uiuc.edu/Research/namd/2.14/download/946183/NAMD_2.14_Linux-x86_64-multicore.tar.gz && \
    tar -xzf NAMD_2.14_Linux-x86_64-multicore.tar.gz && \
    rm NAMD_2.14_Linux-x86_64-multicore.tar.gz

# Add NAMD to the PATH
ENV PATH="/NAMD_2.14_Linux-x86_64-multicore/:${PATH}"

# Command to run when the container starts
CMD ["bash"]
