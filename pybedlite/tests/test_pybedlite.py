"""Tests for :py:mod:`~pybedlite.__init__.py`"""

import pytest
from pathlib import Path
from py._path.local import LocalPath as TmpDir
from typing import List
import pybedlite as pybed
from pybedlite.bed_writer import MAX_BED_FIELDS
from pybedlite.bed_writer import BedWriter
from pybedlite.bed_source import BedSource
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


SNIPPET_BED = """\
1	100	150	test_record1	100	+	100	100	0,0,0	1	50	0
1	200	300	test_record2	100	-	210	290	0,0,0	1	100	0
2	200	300	test_record3	.	.	.	.	.	.	.	.
"""


def compare_bed_records(
    record1: BedRecord,
    record2: BedRecord,
    record_number: int,
    num_fields: int = MAX_BED_FIELDS,
) -> None:
    assert record1.chrom == record2.chrom, f"Chromosome didn't match in record {record_number}"
    assert (
        record1.start == record2.start
    ), f"Start coordinate didn't match in record {record_number}"
    assert record1.end == record2.end, f"End coordinate didn't match in record {record_number}"
    if num_fields >= 4:
        assert record1.name == record2.name, f"Name didn't match in record {record_number}"
    if num_fields >= 5:
        assert record1.score == record2.score, f"Score didn't match in record {record_number}"
    if num_fields >= 6:
        assert record1.strand == record2.strand, f"Strand didn't match in record {record_number}"
    if num_fields >= 7:
        assert (
            record1.thick_start == record2.thick_start
        ), f"Thick start didn't match in record {record_number}"
        assert (
            record1.thick_end == record2.thick_end
        ), f"Thick end didn't match in record {record_number}"
    if num_fields >= 9:
        assert (
            record1.item_rgb == record2.item_rgb
        ), f"Item RGB didn't match in record {record_number}"
    if num_fields >= 10:
        assert (
            record1.block_count == record2.block_count
        ), f"Block count didn't match in record {record_number}"
        assert (
            record1.block_sizes == record2.block_sizes
        ), f"Block sizes didn't match in record {record_number}"
        assert (
            record1.block_starts == record2.block_starts
        ), f"Block starts didn't match in record {record_number}"

    assert record1.as_bed_line(num_fields) == record2.as_bed_line(
        num_fields
    ), f"Derived bed lines differed from expectation in record {record_number}"


@pytest.mark.parametrize(
    "bed_field_number",
    [
        3,
        4,
        5,
        6,
        8,
        9,
        12,
    ],
)
def test_bed_parsing(tmpdir: TmpDir, bed_field_number: int, bed_records: List[BedRecord]) -> None:
    tmpdir_path = Path(tmpdir)
    test_bed = tmpdir_path / "test.bed"

    with BedWriter(test_bed, num_fields=bed_field_number) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    with BedSource(test_bed) as test_in:
        for i, parsed_record in enumerate(test_in):
            expected_record = bed_records[i]
            compare_bed_records(
                record1=parsed_record,
                record2=expected_record,
                record_number=i,
                num_fields=bed_field_number,
            )


@pytest.mark.parametrize(
    "bed_field_number",
    [
        3,
        4,
        5,
        6,
        8,
        9,
        12,
    ],
)
def test_preopened_bed_parsing(
    tmpdir: TmpDir, bed_field_number: int, bed_records: List[BedRecord]
) -> None:
    tmpdir_path = Path(tmpdir)
    test_bed = tmpdir_path / "test.bed"

    with BedWriter(test_bed, num_fields=bed_field_number) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    test_bed_fh = test_bed.open("r")

    with pybed.reader(path=test_bed_fh) as test_in:
        for i, parsed_record in enumerate(test_in):
            expected_record = bed_records[i]
            compare_bed_records(
                record1=parsed_record,
                record2=expected_record,
                record_number=i,
                num_fields=bed_field_number,
            )


@pytest.mark.parametrize(
    "bed_field_number",
    [
        3,
        4,
        5,
        6,
        8,
        9,
        12,
    ],
)
def test_bed_writing(tmpdir: TmpDir, bed_field_number: int, bed_records: List[BedRecord]) -> None:
    tmpdir_path = Path(tmpdir)
    test_written_bed = tmpdir_path / "test_written.bed"
    test_premade_bed = tmpdir_path / "test_premade.bed"

    with BedWriter(test_written_bed, num_fields=bed_field_number) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    with test_premade_bed.open("w") as premade_out:
        premade_out.write(SNIPPET_BED)

    with BedSource(test_written_bed) as test_written_in, BedSource(
        test_premade_bed
    ) as test_premade_in:
        for i, (parsed_record, expected_record) in enumerate(
            zip(test_written_in, test_premade_in)
        ):
            compare_bed_records(
                record1=parsed_record,
                record2=expected_record,
                record_number=i,
                num_fields=bed_field_number,
            )


@pytest.mark.parametrize(
    "bed_field_number",
    [
        3,
        4,
        5,
        6,
        8,
        9,
        12,
    ],
)
def test_preopened_bed_writing(
    tmpdir: TmpDir, bed_field_number: int, bed_records: List[BedRecord]
) -> None:
    tmpdir_path = Path(tmpdir)
    test_written_bed = tmpdir_path / "test_written.bed"
    test_premade_bed = tmpdir_path / "test_premade.bed"

    test_written_bed_fh = test_written_bed.open("w")

    with pybed.writer(
        path=test_written_bed_fh,
        num_fields=bed_field_number,
    ) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    with test_premade_bed.open("w") as premade_out:
        premade_out.write(SNIPPET_BED)

    with BedSource(test_written_bed) as test_written_in, BedSource(
        test_premade_bed
    ) as test_premade_in:
        for i, (parsed_record, expected_record) in enumerate(
            zip(test_written_in, test_premade_in)
        ):
            compare_bed_records(
                record1=parsed_record,
                record2=expected_record,
                record_number=i,
                num_fields=bed_field_number,
            )
