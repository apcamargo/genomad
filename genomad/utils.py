import bz2
import gzip
import hashlib
import io
import json
import lzma
import os
import re
import shutil
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from enum import Enum, auto
from importlib import metadata
from pathlib import Path
from typing import List

import numpy as np
from rich import box
from rich._log_render import LogRender
from rich.console import Console, Group
from rich.padding import Padding
from rich.panel import Panel
from rich.rule import Rule
from rich.tree import Tree
from genomad._paths import GenomadOutputs


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


class HybridConsole:
    def __init__(self, output_file=None, verbose=True, width=None):
        self.output_file = output_file
        self.verbose = verbose
        self.width = width
        # Create the regular console
        if self.verbose:
            self.regular_console = Console(
                width=self.width, highlight=False, record=True
            )
        else:
            devnull_file = open(os.devnull, "w")
            self.regular_console = Console(
                file=devnull_file, width=self.width, highlight=False, record=True
            )
        self.regular_console._log_render = LogRender(
            omit_repeated_times=False,
            show_time=True,
            show_path=False,
            time_format="[%X]",
        )
        # Create the error console
        self.error_console = Console(
            stderr=True, style="red", width=self.width, highlight=False, record=True
        )
        self.error_console._log_render = LogRender(
            omit_repeated_times=False,
            show_time=True,
            show_path=False,
            time_format="[%X]",
        )
        if self.output_file and self.output_file.exists():
            self.output_file.unlink()

    def write_print(self, *args, **kwargs):
        with open(self.output_file, "a") as fout:
            self.writer_console = Console(file=fout, width=self.width)
            self.writer_console._log_render = LogRender(
                omit_repeated_times=False,
                show_time=True,
                show_path=False,
                time_format="[%X]",
            )
            self.writer_console.print(*args, **kwargs)

    def write_log(self, *args, **kwargs):
        with open(self.output_file, "a") as fout:
            self.writer_console = Console(file=fout, width=self.width)
            self.writer_console._log_render = LogRender(
                omit_repeated_times=False,
                show_time=True,
                show_path=False,
                time_format="[%X]",
            )
            self.writer_console.log(*args, **kwargs)

    def print(self, *args, **kwargs):
        self.regular_console.print(*args, **kwargs)
        if self.output_file:
            self.write_print(*args, **kwargs)

    def log(self, *args, **kwargs):
        self.regular_console.log(*args, **kwargs)
        if self.output_file:
            self.write_log(*args, **kwargs)

    def warning(self, *args, **kwargs):
        self.regular_console.log(*args, **kwargs, style="#FFA500")
        if self.output_file:
            self.write_log(*args, **kwargs)

    def error(self, *args, **kwargs):
        self.error_console.log(*args, **kwargs)
        if self.output_file:
            self.write_log(*args, **kwargs)

    def status(self, *args, **kwargs):
        return self.regular_console.status(*args, **kwargs)


def is_compressed(filepath: Path) -> Compression:
    """
    Checks if a file is compressed (gzip, bzip2 or xz).

    Parameters
    ----------
    filepath : Path
        Path object pointing to a file.

    Returns
    -------
    Compression
        Returns a `Compression` enum with one of the following attributes:
        `bzip2`, `gzip`, `xz`, or `uncompressed`
    """
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


def read_file(filepath: Path, skip_header: bool = False) -> str:
    with open_file(filepath) as fin:
        if skip_header:
            next(fin)
        yield from fin


@contextmanager
def suppress_stdout(suppress=True):
    std_ref = sys.stdout
    if suppress:
        sys.stdout = open("/dev/null", "w")
        yield
    sys.stdout = std_ref


def natsort(iterable):
    k = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split("(\d+)", s)]
    return sorted(iterable, key=k)


def check_executables(executables: List[str]) -> List[str]:
    """
    Checks if the executables from the input list are available in the PATH.

    Parameters
    ----------
    executables : list of str
        Names of the executables

    Returns
    -------
    list of str
        Returns a list containing the executables that are not in the PATH.
    """
    return [executable for executable in executables if not shutil.which(executable)]


def get_md5(filepath, size=io.DEFAULT_BUFFER_SIZE):
    m = hashlib.md5()
    with open(filepath, "rb") as fin:
        b = fin.read(size)
        while len(b) > 0:
            m.update(b)
            b = fin.read(size)
        return str(m.hexdigest())


def write_execution_info(
    module_name: str, input_file: Path, parameters: dict, output_file: Path
) -> None:
    input_file_hash = get_md5(input_file)
    start_time = datetime.now(timezone.utc).astimezone().isoformat()
    json_dump = json.dumps(
        {
            "module": module_name,
            "input": input_file.name,
            "input_md5": input_file_hash,
            "start_time": start_time,
            "parameters": parameters,
        },
        indent=4,
    )
    with open(output_file, "w") as fout:
        fout.write(f"{json_dump}\n")


def get_execution_info(input_file: Path):
    with open(input_file) as fin:
        execution_info = json.load(fin)
    input_md5 = execution_info["input_md5"]
    module_name = execution_info["module"]
    parameters = execution_info["parameters"]
    return input_md5, module_name, parameters


def compare_executions(
    input_file: Path,
    parameters: dict,
    execution_info_file: Path,
    only_md5: bool = False,
) -> bool:
    input_md5 = get_md5(input_file)
    previous_input_md5, _, previous_parameters = get_execution_info(execution_info_file)
    if only_md5:
        return input_md5 == previous_input_md5
    else:
        return (parameters == previous_parameters) and (input_md5 == previous_input_md5)


def check_provirus_execution(prefix: str, input_file: Path, output_dir: Path) -> bool:
    outputs = GenomadOutputs(prefix, output_dir)
    if not outputs.find_proviruses_execution_info.exists():
        return False
    input_md5 = get_md5(input_file)
    provirus_input_md5, *_ = get_execution_info(outputs.find_proviruses_execution_info)
    if input_md5 != provirus_input_md5:
        return False
    required_outputs = [
        outputs.find_proviruses_output,
        outputs.find_proviruses_nucleotide_output,
        outputs.find_proviruses_proteins_output,
        outputs.find_proviruses_genes_output,
    ]
    if all(i.exists() for i in required_outputs) and sum(
        1 for _ in read_file(outputs.find_proviruses_output, skip_header=True)
    ):
        return True
    return False


def display_header(
    console: Console,
    module_name: str,
    module_description: str,
    output_dir: Path,
    output_files: List[Path],
    output_descriptions: List[str],
) -> None:
    version = metadata.version("genomad")
    tree = Tree(f"{output_dir}", style="green")
    for f, d in zip(output_files, output_descriptions):
        tree.add(f"{f.name} [dim][white]({d})")
    console.print(
        Panel(
            Group(
                f"Executing [cyan]geNomad {module_name}[/cyan] (v{version}). "
                + module_description,
                Rule(style="dim"),
                "Outputs:",
                Padding(tree, (0, 0, 0, 2)),
            ),
            border_style="dim",
            box=box.ROUNDED,
            padding=(0, 2),
        )
    )


def logistic(x, temperature=1.0):
    return 1 / (1 + np.exp(-x / temperature))


def softmax(x, temperature=1.0):
    assert len(x.shape) == 2
    x /= temperature
    s = np.max(x, axis=1)
    s = s[:, np.newaxis]
    e_x = np.exp(x - s)
    div = np.sum(e_x, axis=1)
    div = div[:, np.newaxis]
    return e_x / div


def entropy(x):
    x = np.array(x)
    n = len(x)
    if not np.any(x):
        return np.log2(n)
    p = x / np.sum(x)
    p = p[p!=0]
    return -1 * np.dot(p, np.log2(p))


def specificity(x):
    x = np.array(x)
    if not np.any(x):
        return 0.0
    n = len(x)
    if n == 1:
        return 0.0
    ss = np.log2(n) - entropy(x)
    return ss / np.log2(n)


def rle_encode(array):
    """
    Run-length encoding of an input array
    """
    value_array = []
    count_array = []
    i = 0
    while i <= len(array) - 1:
        c = 1
        v = array[i]
        j = i
        while j < len(array) - 1 and array[j] == array[j + 1]:
            c += 1
            j += 1
        count_array.append(c)
        value_array.append(v)
        i = j + 1
    return count_array, value_array


def rle_decode(count_array, value_array):
    decoded_array = []
    for c, v in zip(count_array, value_array):
        decoded_array += [v] * c
    return decoded_array
