# Installation

## Installing geNomad

You can install geNomad in you computer using either a general-purpose package manager (`pixi`, `mamba`, `conda`) or a Python-specific package manager (`uv` or `pip`).

::::{tab-set}

:::{tab-item} Pixi
[Pixi](https://pixi.sh/latest/) is probably the simplest way to install geNomad. It takes care of all the dependencies for you and doesn't require any environment management, meaning the `genomad` command will be available globally.

```bash
pixi global install -c conda-forge -c bioconda genomad
genomad --help
```
:::

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/) is a package manager that handles all your dependencies for you. To install geNomad using Mamba, you need to create a new environment and activate it before running the `genomad` command.

```bash
mamba create -n genomad -c conda-forge -c bioconda genomad
mamba activate genomad
genomad --help
```
:::

:::{tab-item} Conda
[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that handles all your dependencies for you. To install geNomad using Conda, you need to create a new environment and activate it before running the `genomad` command.

```bash
conda create -n genomad -c conda-forge -c bioconda genomad
conda activate genomad
genomad --help
```
:::

:::{tab-item} uv
[`uv`](https://github.com/astral-sh/uv) is a Python-specific package manager that lets you install geNomad as a global command within an isolated environment. However, it won't take care of the non-Python dependencies required by geNomad (see note below).

```bash
uv tool install genomad
genomad --help
```
:::

:::{tab-item} pip
`pip` is a Python-specific package manager that lets you install geNomad as a global command. It is the standard tool for installing Python libraries and should be available to everyone with a Python installation. However, `pip` won't take care of the non-Python dependencies required by geNomad (see note below).

```bash
pip install genomad
genomad --help
```
:::

::::

```{admonition} Non-python dependencies
:class: attention
Pixi, Mamba, and Conda will install both the Python dependencies and the non-Python dependencies required by geNomad. If you install geNomad using `uv` or `pip`, make sure to add [`MMseqs2`](https://github.com/soedinglab/MMseqs2/) and [`ARAGORN`](http://www.ansikte.se/ARAGORN/) to your `$PATH`.
```

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
