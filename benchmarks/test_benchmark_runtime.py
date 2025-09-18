from pathlib import Path
from pytest_benchmark.fixture import BenchmarkFixture

from pybedlite import BedSource, BedRecord
from pybedlite.overlap_detector import OverlapDetector

def test_query_exon_against_exome(benchmark: BenchmarkFixture, known_genes_hg38: Path) -> None:
    detector = OverlapDetector[BedRecord]()
    detector.add_all(BedSource(known_genes_hg38))
    assert len(detector) == 412_044
    benchmark(lambda: [detector.get_overlaps(bed) for bed in BedSource(known_genes_hg38)])
