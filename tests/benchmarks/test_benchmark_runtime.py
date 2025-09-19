from pathlib import Path

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from pybedlite import BedRecord
from pybedlite import BedSource
from pybedlite.overlap_detector import OverlapDetector


@pytest.mark.benchmark
def test_query_exon_against_exome(benchmark: BenchmarkFixture, known_genes_hg38: Path) -> None:
    """Benchmark loading all known human genes into the overlap detector and querying each gene."""
    detector = OverlapDetector[BedRecord](BedSource(known_genes_hg38))
    assert len(detector) == 412_044
    benchmark(lambda: [detector.get_overlaps(bed) for bed in BedSource(known_genes_hg38)])
