# Issue build command from main directory 
FROM condaforge/miniforge3
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/mosdef-gomc.yml env.yml 
RUN mamba env update --file=env.yml -n base && mamba clean --all
RUN mamba update -c conda-forge procps-ng && mamba clean --all --yes 
RUN mamba init

ENV BASH_ENV=/root/.bashrc

# Set the default command to run when the container starts
CMD ["/bin/bash"]
