"""Tests for :py:mod:`~pybedlite.__init__.py`."""

from pathlib import Path
from typing import List

import pytest

import pybedlite as pybed
from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand
from pybedlite.bed_source import BedSource
from pybedlite.bed_writer import MAX_BED_FIELDS
from pybedlite.bed_writer import BedWriter

SNIPPET_BED = """\
# Test header, with a line with whitespace below it.

# More header
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
    assert record1.start == record2.start, (
        f"Start coordinate didn't match in record {record_number}"
    )
    assert record1.end == record2.end, f"End coordinate didn't match in record {record_number}"
    if num_fields >= 4:
        assert record1.name == record2.name, f"Name didn't match in record {record_number}"
    if num_fields >= 5:
        assert record1.score == record2.score, f"Score didn't match in record {record_number}"
    if num_fields >= 6:
        assert record1.strand == record2.strand, f"Strand didn't match in record {record_number}"
    if num_fields >= 7:
        assert record1.thick_start == record2.thick_start, (
            f"Thick start didn't match in record {record_number}"
        )
        assert record1.thick_end == record2.thick_end, (
            f"Thick end didn't match in record {record_number}"
        )
    if num_fields >= 9:
        assert record1.item_rgb == record2.item_rgb, (
            f"Item RGB didn't match in record {record_number}"
        )
    if num_fields >= 10:
        assert record1.block_count == record2.block_count, (
            f"Block count didn't match in record {record_number}"
        )
        assert record1.block_sizes == record2.block_sizes, (
            f"Block sizes didn't match in record {record_number}"
        )
        assert record1.block_starts == record2.block_starts, (
            f"Block starts didn't match in record {record_number}"
        )

    assert record1.as_bed_line(num_fields) == record2.as_bed_line(num_fields), (
        f"Derived bed lines differed from expectation in record {record_number}"
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
def test_bed_parsing(tmp_path: Path, bed_field_number: int, bed_records: List[BedRecord]) -> None:
    test_bed = tmp_path / "test.bed"

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
    tmp_path: Path, bed_field_number: int, bed_records: List[BedRecord]
) -> None:
    test_bed = tmp_path / "test.bed"

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
def test_bed_writing(tmp_path: Path, bed_field_number: int, bed_records: List[BedRecord]) -> None:
    test_written_bed = tmp_path / "test_written.bed"
    test_premade_bed = tmp_path / "test_premade.bed"

    with BedWriter(test_written_bed, num_fields=bed_field_number) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    with test_premade_bed.open("w") as premade_out:
        premade_out.write(SNIPPET_BED)

    with (
        BedSource(test_written_bed) as test_written_in,
        BedSource(test_premade_bed) as test_premade_in,
    ):
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
    tmp_path: Path, bed_field_number: int, bed_records: List[BedRecord]
) -> None:
    test_written_bed = tmp_path / "test_written.bed"
    test_premade_bed = tmp_path / "test_premade.bed"

    test_written_bed_fh = test_written_bed.open("w")

    with pybed.writer(
        path=test_written_bed_fh,
        num_fields=bed_field_number,
    ) as test_out:
        test_out.write_all(bed_records, truncate=True, add_missing=True)

    with test_premade_bed.open("w") as premade_out:
        premade_out.write(SNIPPET_BED)

    with (
        BedSource(test_written_bed) as test_written_in,
        BedSource(test_premade_bed) as test_premade_in,
    ):
        for i, (parsed_record, expected_record) in enumerate(
            zip(test_written_in, test_premade_in)
        ):
            compare_bed_records(
                record1=parsed_record,
                record2=expected_record,
                record_number=i,
                num_fields=bed_field_number,
            )


def test_bedstrand_opposite() -> None:
    """Test that we can reverse a BedStrand."""
    assert BedStrand.Positive.opposite is BedStrand.Negative
    assert BedStrand.Negative.opposite is BedStrand.Positive


def test_bedrecord_refname() -> None:
    """Test that the alternate property for reference sequence name is correct."""
    assert BedRecord(chrom="chr1", start=1, end=2).refname == "chr1"


def test_bedrecord_negative() -> None:
    """Test that the negative property is set correctly."""
    assert not BedRecord(chrom="chr1", start=1, end=2, strand=None).negative
    assert not BedRecord(chrom="chr1", start=1, end=2, strand=BedStrand.Positive).negative
    assert BedRecord(chrom="chr1", start=1, end=2, strand=BedStrand.Negative).negative

    """Test that BedRecord validates start < end."""
    # Test start == end
    with pytest.raises(ValueError, match="End of interval must be greater than start"):
        BedRecord(chrom="chr1", start=100, end=100)

    # Test start > end
    with pytest.raises(ValueError, match="End of interval must be greater than start"):
        BedRecord(chrom="chr1", start=100, end=50)


def test_bedrecord_validation_thick_start_end_mismatch() -> None:
    """Test that BedRecord validates thick_start and thick_end are both defined or both None."""
    # thick_end without thick_start
    with pytest.raises(ValueError, match="Thick end cannot be defined if thick start is not"):
        BedRecord(chrom="chr1", start=100, end=200, thick_start=None, thick_end=150)

    # thick_start without thick_end
    with pytest.raises(ValueError, match="Thick start cannot be defined if thick end is not"):
        BedRecord(chrom="chr1", start=100, end=200, thick_start=120, thick_end=None)


def test_bedrecord_validation_block_count_without_sizes() -> None:
    """Test that BedRecord validates block_sizes is defined when block_count is."""
    with pytest.raises(ValueError, match="Block sizes cannot be undefined if block count"):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=None,
            block_starts=[0, 50],
        )


def test_bedrecord_validation_block_count_without_starts() -> None:
    """Test that BedRecord validates block_starts is defined when block_count is."""
    with pytest.raises(ValueError, match="Block starts cannot be undefined if block count"):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=[30, 40],
            block_starts=None,
        )


def test_bedrecord_validation_block_sizes_without_count() -> None:
    """Test that BedRecord validates block_count is defined when block_sizes is."""
    with pytest.raises(ValueError, match="Block count must be defined if block sizes"):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=None,
            block_sizes=[30, 40],
            block_starts=[0, 60],
        )


def test_bedrecord_validation_block_starts_without_count() -> None:
    """Test that BedRecord validates block_count is defined when block_starts is."""
    with pytest.raises(ValueError, match="Block count must be defined if block starts"):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=None,
            block_sizes=None,
            block_starts=[0, 60],
        )


def test_bedrecord_validation_block_sizes_length_mismatch() -> None:
    """Test that BedRecord validates block_sizes length matches block_count."""
    with pytest.raises(
        ValueError, match="Number of items in block_sizes .* must match block_count"
    ):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=[30, 40, 30],
            block_starts=[0, 60],
        )


def test_bedrecord_validation_block_starts_length_mismatch() -> None:
    """Test that BedRecord validates block_starts length matches block_count."""
    with pytest.raises(
        ValueError, match="Number of items in block_starts .* must match block_count"
    ):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=[30, 40],
            block_starts=[0, 60, 90],
        )


def test_bedrecord_validation_block_starts_first_nonzero() -> None:
    """Test that BedRecord validates first block_start is zero."""
    with pytest.raises(ValueError, match="Block start at first position should be zero"):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=[30, 40],
            block_starts=[10, 60],
        )


def test_bedrecord_validation_block_end_mismatch() -> None:
    """Test that BedRecord validates last block end matches interval end."""
    with pytest.raises(
        ValueError, match="Overall interval end should be equal to the last defined block's end"
    ):
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            block_count=2,
            block_sizes=[30, 40],
            block_starts=[0, 50],
        )


def test_bedrecord_as_bed_line_invalid_field_count() -> None:
    """Test that as_bed_line raises ValueError for invalid field counts."""
    record = BedRecord(chrom="chr1", start=100, end=200)

    # Test field count < 3
    with pytest.raises(ValueError, match="BED records can only contain between 3 and 12 fields"):
        record.as_bed_line(number_of_output_fields=2)

    # Test field count > 12
    with pytest.raises(ValueError, match="BED records can only contain between 3 and 12 fields"):
        record.as_bed_line(number_of_output_fields=13)


def test_bedsource_close_unopened_file(tmp_path: Path) -> None:
    """Test that BedSource raises ValueError when closing an unopened file."""
    test_bed = tmp_path / "test.bed"
    test_bed.write_text("chr1\t100\t200\n")

    source = BedSource(test_bed)
    with pytest.raises(ValueError, match="Cannot close file .* if it is not already open"):
        source.close()


def test_bedsource_invalid_rgb(tmp_path: Path) -> None:
    """Test that BedSource raises ValueError for invalid RGB values."""
    test_bed = tmp_path / "test.bed"
    # RGB should have 3 components, but we provide 2
    test_bed.write_text("chr1\t100\t200\tname\t100\t+\t100\t200\t255,255\n")

    with pytest.raises(ValueError, match="item_rgb.* must contain 3 comma separated integers"):
        with BedSource(test_bed) as source:
            list(source)


def test_bedsource_too_few_fields(tmp_path: Path) -> None:
    """Test that BedSource raises ValueError for BED records with fewer than 3 fields."""
    test_bed = tmp_path / "test.bed"
    # Only 2 fields provided
    test_bed.write_text("chr1\t100\n")

    with pytest.raises(
        ValueError, match="BED records must conform to specifications.* at least 3 input fields"
    ):
        with BedSource(test_bed) as source:
            list(source)


def test_bedwriter_invalid_num_fields_too_small(tmp_path: Path) -> None:
    """Test that BedWriter raises ValueError for num_fields < 3."""
    test_bed = tmp_path / "test.bed"

    with pytest.raises(ValueError, match="BED files must contain between 3 and 12 columns"):
        BedWriter(test_bed, num_fields=2)


def test_bedwriter_invalid_num_fields_too_large(tmp_path: Path) -> None:
    """Test that BedWriter raises ValueError for num_fields > 12."""
    test_bed = tmp_path / "test.bed"

    with pytest.raises(ValueError, match="BED files must contain between 3 and 12 columns"):
        BedWriter(test_bed, num_fields=13)


def test_bedwriter_close_unopened_file(tmp_path: Path) -> None:
    """Test that BedWriter raises ValueError when closing an unopened file."""
    test_bed = tmp_path / "test.bed"

    writer = BedWriter(test_bed)
    with pytest.raises(ValueError, match="Cannot close file .* if it is not already open"):
        writer.close()


def test_bedrecord_bed_field_num_coverage() -> None:
    """Test bed_field_num property for various field counts."""
    # BED3
    assert BedRecord(chrom="chr1", start=100, end=200).bed_field_num == 3

    # BED4
    assert BedRecord(chrom="chr1", start=100, end=200, name="test").bed_field_num == 4

    # BED5
    assert BedRecord(chrom="chr1", start=100, end=200, name="test", score=100).bed_field_num == 5

    # BED6
    assert (
        BedRecord(
            chrom="chr1", start=100, end=200, name="test", score=100, strand=BedStrand.Positive
        ).bed_field_num
        == 6
    )

    # BED8
    assert (
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            name="test",
            score=100,
            strand=BedStrand.Positive,
            thick_start=110,
            thick_end=190,
        ).bed_field_num
        == 8
    )

    # BED9
    assert (
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            name="test",
            score=100,
            strand=BedStrand.Positive,
            thick_start=110,
            thick_end=190,
            item_rgb=(255, 0, 0),
        ).bed_field_num
        == 9
    )

    # BED12
    assert (
        BedRecord(
            chrom="chr1",
            start=100,
            end=200,
            name="test",
            score=100,
            strand=BedStrand.Positive,
            thick_start=110,
            thick_end=190,
            item_rgb=(255, 0, 0),
            block_count=2,
            block_sizes=[50, 50],
            block_starts=[0, 50],
        ).bed_field_num
        == 12
    )
