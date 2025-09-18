import gzip
from pathlib import Path

import pytest
import requests


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
