FROM mambaorg/micromamba:2.0.5
LABEL maintainer="Antonio Camargo (antoniop.camargo@lbl.gov)"
LABEL version="1.11.0"

RUN micromamba install -y -n base -c conda-forge -c bioconda genomad==1.11.0 && \
    micromamba clean --all --yes
WORKDIR /app
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "genomad"]
