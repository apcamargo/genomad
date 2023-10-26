# Installation

## Installing geNomad

You can install geNomad in you computer using either a general-purpose package managers (`mamba` or `conda`) or a Python-specific package manager (`pipx` or `pip`). `pip` is the standard command to install Python libraries and should be available to everyone with a Python installation.

::::{tab-set}

:::{tab-item} mamba
```bash
mamba create -n genomad -c conda-forge -c bioconda genomad
mamba activate genomad
```
:::

:::{tab-item} conda
```bash
conda create -n genomad -c conda-forge -c bioconda genomad
conda activate genomad
```
:::

:::{tab-item} pipx
```bash
pipx install genomad
```
:::

:::{tab-item} pip
```bash
pip install genomad
```
:::

::::

```{admonition} pipx installation
:class: tip
We recommend using [`pipx`](https://pypa.github.io/pipx/) over `pip` if possible. By using `pipx` you will avoid dependency conflicts that might arise if you try to install geNomad in an existing Python environment.
```

Conda and Mamba will install both the Python dependencies and the third-party software required by geNomad. If you install geNomad using `pip` or `pipx`, make sure to add [`MMseqs2`](https://github.com/soedinglab/MMseqs2/) and [`ARAGORN`](http://www.ansikte.se/ARAGORN/) to your `$PATH`.

## Running geNomad using containers

You can also execute geNomad using containerization tools, such as Docker and Podman. To pull the image, execute the command below.

::::{tab-set}

:::{tab-item} Docker
```bash
docker pull antoniopcamargo/genomad
```
:::

:::{tab-item} Podman
```bash
podman pull docker.io/antoniopcamargo/genomad
```
:::

::::

To start a geNomad container you have to mount a folder from the host system into the container with the `-v` argument. The following command mounts the current working directory (`$(pwd)`) under `/app` inside the container and then executes the `genomad download-database` and `genomad end-to-end` commands.

::::{tab-set}

:::{tab-item} Docker
```bash
docker run -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad download-database .
docker run -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad end-to-end input.fna output genomad_db
```
:::

:::{tab-item} Podman
```bash
podman run -u 0 -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad download-database .
podman run -u 0 -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad end-to-end input.fna output genomad_db
```
:::

::::
