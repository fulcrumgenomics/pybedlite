from typing import List

import pytest

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand


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
