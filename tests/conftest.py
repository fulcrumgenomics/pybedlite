import gzip
from pathlib import Path
from typing import List

import pytest
import requests

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add the `--benchmark` option to pytest."""
    parser.addoption(
        "--benchmark",
        action="store_true",
        default=False,
        help="Run benchmarking tests",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Modify tests that are marked as slow with whether we will run them or not."""
    if config.getoption("--benchmark"):
        return None
    skip_slow = pytest.mark.skip(reason="need --benchmark option to run")
    for item in items:
        if "benchmark" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def bed_records() -> List[BedRecord]:
    return [
        BedRecord(
            chrom="1",
            start=100,
            end=150,
            name="test_record1",
            score=100,
            strand=BedStrand.Positive,
            thick_start=100,
            thick_end=100,
            item_rgb=(0, 0, 0),
            block_count=1,
            block_sizes=[50],
            block_starts=[0],
        ),
        BedRecord(
            chrom="1",
            start=200,
            end=300,
            name="test_record2",
            score=100,
            strand=BedStrand.Negative,
            thick_start=210,
            thick_end=290,
            item_rgb=(0, 0, 0),
            block_count=1,
            block_sizes=[100],
            block_starts=[0],
        ),
        BedRecord(
            chrom="2",
            start=200,
            end=300,
            name="test_record3",
            score=None,
            strand=None,
            thick_start=None,
            thick_end=None,
            item_rgb=None,
            block_count=None,
            block_sizes=None,
            block_starts=None,
        ),
    ]


def fetch_known_genes() -> list[list[str]]:
    """Fetch known genes for hg38 from the UCSC Genome Browser."""
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz"
    response = requests.get(url)
    text = gzip.decompress(response.content).decode("utf-8")
    return [line.split("\t") for line in text.splitlines()]


def write_known_genes(output: Path) -> None:
    """Write known genes for hg38 to an uncompressed BED file."""
    with open(output, "w") as writer:
        for line in fetch_known_genes():
            chrom, start, end, name, strand = line[1], line[3], line[4], line[0], line[2]
            writer.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


@pytest.fixture(scope="function")
def known_genes_hg38(cache: pytest.Cache) -> Path:
    """Fixture to provide known genes BED file for hg38 in a cached location."""
    test_data_dir = cache.mkdir("test_data_dir")
    if not (test_data_dir / "known_genes.hg38.bed").exists():
        write_known_genes(test_data_dir / "known_genes.hg38.bed")
    return test_data_dir / "known_genes.hg38.bed"
