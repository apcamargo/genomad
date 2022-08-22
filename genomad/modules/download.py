import shutil
import sys
import urllib
from functools import partial
from urllib.request import urlopen

import genomad
from genomad import utils
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)


class DatabaseDownloader:
    def __init__(self, destination, console, version=None):
        self.base_url = "https://portal.nersc.gov/genomad/__data__/"
        self.destination = destination
        self.console = console
        self.version = str(version) if version else self.get_version()
        self.filename = f"genomad_db_v{self.version}.tar.gz"
        self.database_url = self.base_url + self.filename
        self.output_file = self.destination / self.filename

    def get_version(self):
        genomad_version = genomad.__version__.rsplit(".", 1)[0]
        file_contents = (
            urllib.request.urlopen(self.base_url + "releases.txt")
            .read()
            .decode("utf-8")
            .strip()
            .split("\n")
        )
        selected_db_version = None
        for line in file_contents[1:]:
            db_version, package_version = line.strip().split("\t")
            if package_version == genomad_version:
                selected_db_version = db_version
        if selected_db_version:
            return selected_db_version
        else:
            self.console.error("A compatible database version was not found.")
            sys.exit(1)

    def _copy_url(self, task_id, progress):
        progress.console.log(
            f"Requesting [blue link=self.database_url]{self.database_url}[/blue link]."
        )
        response = urlopen(self.database_url)
        progress.update(task_id, total=int(response.info()["Content-length"]))
        with open(self.output_file, "wb") as dest_file:
            progress.start_task(task_id)
            for data in iter(partial(response.read, 32768), b""):
                dest_file.write(data)
                progress.update(task_id, advance=len(data))

    def download(self):
        progress = Progress(
            TextColumn("{task.fields[filename]}", justify="right", style="green"),
            BarColumn(bar_width=None),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "|",
            DownloadColumn(),
            "|",
            TransferSpeedColumn(),
            "|",
            TimeRemainingColumn(elapsed_when_finished=True),
            console=self.console.regular_console,
            transient=True,
        )
        with progress:
            task_id = progress.add_task("download", filename=self.filename, start=False)
            self._copy_url(task_id, progress)

    def extract(self):
        shutil.unpack_archive(self.output_file, self.destination, "gztar")


def main(destination, keep, verbose):
    console = utils.HybridConsole(verbose=verbose)
    database_downloader = DatabaseDownloader(destination, console)
    database_downloader.download()
    console.log(
        f"Finished downloading [green]{database_downloader.output_file.name}[/green]."
    )
    database_path = database_downloader.destination / "genomad_db"
    with console.status(f"Extracting the database to [green]{database_path}[/green]."):
        database_downloader.extract()
        console.log(f"Database extracted to [green]{database_path}[/green].")
    if not keep:
        console.log(f"Deleting [green]{database_downloader.output_file.name}[/green].")
        database_downloader.output_file.unlink()
    console.log(
        f"geNomad database (v{database_downloader.version}) is ready to be used!",
        style="yellow",
    )
