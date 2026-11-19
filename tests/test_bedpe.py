"""Tests for BEDPE record reading, writing, and validation."""

import logging
from pathlib import Path
from typing import Any

import pytest

import pybedlite as pybed
from pybedlite.bed_record import BedStrand
from pybedlite.bedpe_record import MAX_BEDPE_FIELDS
from pybedlite.bedpe_record import MIN_BEDPE_FIELDS
from pybedlite.bedpe_record import BedPeRecord
from pybedlite.bedpe_source import BedPeSource
from pybedlite.bedpe_writer import BedPeWriter

SNIPPET_BEDPE = """\
# Test header

chr1\t100\t200\tchr2\t300\t400\trecord1\t500\t+\t-
chr3\t1000\t2000\tchr3\t3000\t4000\trecord2\t0\t+\t+
chrX\t50\t100\tchrY\t200\t250\t.\t.\t.\t.
"""


@pytest.fixture
def bedpe_records() -> list[BedPeRecord]:
    """A set of three BedPeRecords covering 10-field, 10-field, and 6-field cases."""
    return [
        BedPeRecord(
            chrom1="chr1",
            start1=100,
            end1=200,
            chrom2="chr2",
            start2=300,
            end2=400,
            name="record1",
            score=500,
            strand1=BedStrand.Positive,
            strand2=BedStrand.Negative,
        ),
        BedPeRecord(
            chrom1="chr3",
            start1=1000,
            end1=2000,
            chrom2="chr3",
            start2=3000,
            end2=4000,
            name="record2",
            score=0,
            strand1=BedStrand.Positive,
            strand2=BedStrand.Positive,
        ),
        BedPeRecord(
            chrom1="chrX",
            start1=50,
            end1=100,
            chrom2="chrY",
            start2=200,
            end2=250,
        ),
    ]


def compare_bedpe_records(
    record1: BedPeRecord,
    record2: BedPeRecord,
    record_number: int,
    num_fields: int = MAX_BEDPE_FIELDS,
) -> None:
    """Assert field-by-field equality between two BedPeRecords up to num_fields."""
    assert record1.chrom1 == record2.chrom1, f"chrom1 didn't match in record {record_number}"
    assert record1.start1 == record2.start1, f"start1 didn't match in record {record_number}"
    assert record1.end1 == record2.end1, f"end1 didn't match in record {record_number}"
    assert record1.chrom2 == record2.chrom2, f"chrom2 didn't match in record {record_number}"
    assert record1.start2 == record2.start2, f"start2 didn't match in record {record_number}"
    assert record1.end2 == record2.end2, f"end2 didn't match in record {record_number}"
    if num_fields >= 7:
        assert record1.name == record2.name, f"name didn't match in record {record_number}"
    if num_fields >= 8:
        assert record1.score == record2.score, f"score didn't match in record {record_number}"
    if num_fields >= 9:
        assert record1.strand1 == record2.strand1, (
            f"strand1 didn't match in record {record_number}"
        )
    if num_fields >= 10:
        assert record1.strand2 == record2.strand2, (
            f"strand2 didn't match in record {record_number}"
        )

    assert record1.as_bedpe_line(num_fields) == record2.as_bedpe_line(num_fields), (
        f"Derived BEDPE lines differed from expectation in record {record_number}"
    )


@pytest.mark.parametrize("num_fields", [6, 7, 8, 9, 10])
def test_bedpe_parsing_round_trip(
    tmp_path: Path, num_fields: int, bedpe_records: list[BedPeRecord]
) -> None:
    """Write records with a given num_fields then parse them back to verify fields."""
    test_bedpe = tmp_path / "test.bedpe"

    with BedPeWriter(test_bedpe, num_fields=num_fields) as out_fh:
        out_fh.write_all(bedpe_records, truncate=True, add_missing=True)

    with BedPeSource(test_bedpe) as in_fh:
        for i, parsed_record in enumerate(in_fh):
            compare_bedpe_records(
                record1=parsed_record,
                record2=bedpe_records[i],
                record_number=i,
                num_fields=num_fields,
            )


@pytest.mark.parametrize("num_fields", [6, 7, 8, 9, 10])
def test_preopened_bedpe_parsing(
    tmp_path: Path, num_fields: int, bedpe_records: list[BedPeRecord]
) -> None:
    """Parse a BEDPE file via a pre-opened file handle using the module-level bedpe_reader."""
    test_bedpe = tmp_path / "test.bedpe"

    with BedPeWriter(test_bedpe, num_fields=num_fields) as out_fh:
        out_fh.write_all(bedpe_records, truncate=True, add_missing=True)

    with pybed.bedpe_reader(path=test_bedpe.open("r")) as in_fh:
        for i, parsed_record in enumerate(in_fh):
            compare_bedpe_records(
                record1=parsed_record,
                record2=bedpe_records[i],
                record_number=i,
                num_fields=num_fields,
            )


@pytest.mark.parametrize("num_fields", [6, 7, 8, 9, 10])
def test_bedpe_writing(tmp_path: Path, num_fields: int, bedpe_records: list[BedPeRecord]) -> None:
    """Write records to a file and compare against a pre-made BEDPE snippet."""
    test_written = tmp_path / "test_written.bedpe"
    test_premade = tmp_path / "test_premade.bedpe"

    with BedPeWriter(test_written, num_fields=num_fields) as out_fh:
        out_fh.write_all(bedpe_records, truncate=True, add_missing=True)

    with test_premade.open("w") as fh:
        fh.write(SNIPPET_BEDPE)

    with (
        BedPeSource(test_written) as written_in,
        BedPeSource(test_premade) as premade_in,
    ):
        for i, (parsed, expected) in enumerate(zip(written_in, premade_in)):
            compare_bedpe_records(
                record1=parsed,
                record2=expected,
                record_number=i,
                num_fields=num_fields,
            )
        assert i + 1 == len(bedpe_records)


@pytest.mark.parametrize("num_fields", [6, 7, 8, 9, 10])
def test_preopened_bedpe_writing(
    tmp_path: Path, num_fields: int, bedpe_records: list[BedPeRecord]
) -> None:
    """Write via a pre-opened file handle and compare against a pre-made snippet."""
    test_written = tmp_path / "test_written.bedpe"
    test_premade = tmp_path / "test_premade.bedpe"

    with pybed.bedpe_writer(path=test_written.open("w"), num_fields=num_fields) as out_fh:
        out_fh.write_all(bedpe_records, truncate=True, add_missing=True)

    with test_premade.open("w") as fh:
        fh.write(SNIPPET_BEDPE)

    with (
        BedPeSource(test_written) as written_in,
        BedPeSource(test_premade) as premade_in,
    ):
        for i, (parsed, expected) in enumerate(zip(written_in, premade_in)):
            compare_bedpe_records(
                record1=parsed,
                record2=expected,
                record_number=i,
                num_fields=num_fields,
            )
        assert i + 1 == len(bedpe_records)


@pytest.mark.parametrize(
    "kwargs, expected_fields",
    [
        (dict(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4), 6),
        (dict(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4, name="sv_001"), 7),
        (
            dict(
                chrom1="chr1",
                start1=1,
                end1=2,
                chrom2="chr2",
                start2=3,
                end2=4,
                name="sv_001",
                score=10,
            ),
            8,
        ),
        (
            dict(
                chrom1="chr1",
                start1=1,
                end1=2,
                chrom2="chr2",
                start2=3,
                end2=4,
                name="sv_001",
                score=10,
                strand1=BedStrand.Positive,
            ),
            9,
        ),
        (
            dict(
                chrom1="chr1",
                start1=1,
                end1=2,
                chrom2="chr2",
                start2=3,
                end2=4,
                name="sv_001",
                score=10,
                strand1=BedStrand.Positive,
                strand2=BedStrand.Negative,
            ),
            10,
        ),
    ],
)
def test_bedpe_record_field_num(kwargs: dict[str, Any], expected_fields: int) -> None:
    """bedpe_field_num returns the count of populated optional fields and 6 required fields."""
    assert BedPeRecord(**kwargs).bedpe_field_num == expected_fields


@pytest.mark.parametrize("num_fields", [5, 11])
def test_bedpe_record_as_bedpe_line_invalid_field_count(num_fields: int) -> None:
    """as_bedpe_line raises ValueError when number_of_output_fields is outside [6, 10]."""
    record = BedPeRecord(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4)
    with pytest.raises(
        ValueError,
        match=f"BEDPE records can only contain between {MIN_BEDPE_FIELDS} and {MAX_BEDPE_FIELDS}",
    ):
        record.as_bedpe_line(number_of_output_fields=num_fields)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (
            dict(chrom1="chr1", start1=100, end1=100, chrom2="chr2", start2=1, end2=2),
            "End of first interval must be greater than start",
        ),
        (
            dict(chrom1="chr1", start1=200, end1=100, chrom2="chr2", start2=1, end2=2),
            "End of first interval must be greater than start",
        ),
        (
            dict(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=100, end2=100),
            "End of second interval must be greater than start",
        ),
        (
            dict(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=200, end2=100),
            "End of second interval must be greater than start",
        ),
    ],
)
def test_bedpe_record_validation_end_le_start(kwargs: dict[str, Any], match: str) -> None:
    """BedPeRecord raises ValueError when end1 <= start1 or end2 <= start2."""
    with pytest.raises(ValueError, match=match):
        BedPeRecord(**kwargs)


def test_bedpe_source_close_unopened_file(tmp_path: Path) -> None:
    """Closing a BedPeSource that was never opened raises ValueError."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\t4\n")

    source = BedPeSource(test_bedpe)
    with pytest.raises(ValueError, match="Cannot close file .* if it is not already open"):
        source.close()


def test_bedpe_source_too_few_fields(tmp_path: Path) -> None:
    """Parsing a BEDPE line with fewer than six fields raises ValueError."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\n")

    with pytest.raises(
        ValueError,
        match=f"BEDPE records must conform to specifications.*{MIN_BEDPE_FIELDS} input fields",  # noqa: E501
    ):
        with BedPeSource(test_bedpe) as source:
            list(source)


def test_bedpe_source_skips_header_lines(tmp_path: Path) -> None:
    """BedPeSource skips comment, browser, and track lines and yields only data records."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("# comment\nbrowser foo\ntrack name=test\nchr1\t1\t2\tchr2\t3\t4\n")

    with BedPeSource(test_bedpe) as source:
        records = list(source)

    assert len(records) == 1
    assert records[0].chrom1 == "chr1"


@pytest.mark.parametrize("num_fields", [5, 11])
def test_bedpe_writer_invalid_num_fields(tmp_path: Path, num_fields: int) -> None:
    """BedPeWriter raises ValueError when num_fields is outside the valid range of [6, 10]."""
    with pytest.raises(
        ValueError,
        match=f"BEDPE files must contain between {MIN_BEDPE_FIELDS} and {MAX_BEDPE_FIELDS}",
    ):
        BedPeWriter(tmp_path / "test.bedpe", num_fields=num_fields)


def test_bedpe_writer_close_unopened_file(tmp_path: Path) -> None:
    """Closing a BedPeWriter that was never opened raises ValueError."""
    writer = BedPeWriter(tmp_path / "test.bedpe")
    with pytest.raises(ValueError, match="Cannot close file .* if it is not already open"):
        writer.close()


def test_bedpe_writer_truncate_raises_without_flag(tmp_path: Path) -> None:
    """Writing a rec with more fields than allowed raises unless truncate=True."""
    record = BedPeRecord(
        chrom1="chr1",
        start1=1,
        end1=2,
        chrom2="chr2",
        start2=3,
        end2=4,
        name="sv_001",
        score=10,
        strand1=BedStrand.Positive,
        strand2=BedStrand.Negative,
    )
    with BedPeWriter(tmp_path / "test.bedpe", num_fields=6) as writer:
        with pytest.raises(ValueError, match="truncate must be set to True"):
            writer.write(record)


def test_bedpe_writer_add_missing_raises_without_flag(tmp_path: Path) -> None:
    """Writing a rec with fewer fields than expected raises unless add_missing=True."""
    record = BedPeRecord(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4)
    with BedPeWriter(tmp_path / "test.bedpe", num_fields=10) as writer:
        with pytest.raises(ValueError, match="add_missing must be set to True"):
            writer.write(record)


def test_bedpe_source_infers_num_fields(tmp_path: Path) -> None:
    """BedPeSource sets num_fields from the field count of the first record read."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\t4\tname\t100\t+\t-\n")

    with BedPeSource(test_bedpe) as source:
        list(source)
        assert source.num_fields == 10


def test_bedpe_iterator_without_context_manager(tmp_path: Path) -> None:
    """BedPeSource can be iterated directly without a context manager."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\t4\n")

    source = BedPeSource(test_bedpe)
    records = list(source)
    assert len(records) == 1
    assert records[0].chrom1 == "chr1"


def test_bedpe_source_reiteration(tmp_path: Path) -> None:
    """A BedPeSource opened via path can be iterated multiple times, yielding identical records."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\t4\n")

    source = BedPeSource(test_bedpe)
    records1 = list(source)
    records2 = list(source)
    assert records1 == records2


def test_bedpe_source_invalid_path_type() -> None:
    """BedPeSource raises TypeError when passed an unsupported path type."""
    with pytest.raises(TypeError, match="Cannot open"):
        BedPeSource(42)  # type: ignore[arg-type]


def test_bedpe_writer_invalid_path_type() -> None:
    """BedPeWriter raises TypeError when passed an unsupported path type."""
    with pytest.raises(TypeError, match="Cannot open"):
        BedPeWriter(42)  # type: ignore[arg-type]


def test_bedpe_source_invalid_strand(tmp_path: Path) -> None:
    """Parsing a BEDPE line with an unrecognised strand character raises ValueError."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text("chr1\t1\t2\tchr2\t3\t4\tname\t100\t?\t+\n")

    with pytest.raises(ValueError, match="is not a valid BedStrand"):
        with BedPeSource(test_bedpe) as source:
            list(source)


def test_bedpe_record_strand2_without_strand1(caplog: pytest.LogCaptureFixture) -> None:
    """A record with strand2 set and strand1=None logs a warning and serializes strand1 as '.'."""
    with caplog.at_level(logging.WARNING):
        record = BedPeRecord(
            chrom1="chr1",
            start1=1,
            end1=2,
            chrom2="chr2",
            start2=3,
            end2=4,
            strand2=BedStrand.Positive,
        )
    assert "strand2 is set but strand1 is not" in caplog.text
    assert record.bedpe_field_num == 10
    fields = record.as_bedpe_line().split("\t")
    assert fields[8] == "."
    assert fields[9] == "+"


def test_bedpe_writer_write_after_close(tmp_path: Path) -> None:
    """Calling write() on a closed BedPeWriter raises ValueError."""
    record = BedPeRecord(chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4)
    writer = BedPeWriter(tmp_path / "test.bedpe")
    with writer:
        pass
    with pytest.raises(ValueError, match="Cannot write to a closed BedPeWriter"):
        writer.write(record)


def test_bedpe_writer_num_fields_auto_inferred(tmp_path: Path) -> None:
    """When num_fields is not set, BedPeWriter infers it from the first record written."""
    record = BedPeRecord(
        chrom1="chr1",
        start1=1,
        end1=2,
        chrom2="chr2",
        start2=3,
        end2=4,
        name="sv_001",
        score=10,
    )
    test_bedpe = tmp_path / "test.bedpe"
    with BedPeWriter(test_bedpe) as writer:
        assert writer.num_fields is None
        writer.write(record)
        assert writer.num_fields == 8

    with BedPeSource(test_bedpe) as source:
        records = list(source)
    assert len(records) == 1
    assert source.num_fields == 8


def test_bedpe_record_zero_start_is_valid() -> None:
    """A BedPeRecord with start1=0 or start2=0 is valid (0-based coordinates allow zero)."""
    record = BedPeRecord(chrom1="chr1", start1=0, end1=1, chrom2="chr2", start2=0, end2=1)
    assert record.start1 == 0
    assert record.start2 == 0


def test_bedpe_record_score_zero_field_num() -> None:
    """score=0 is treated as a defined field, so bedpe_field_num reflects it."""
    record = BedPeRecord(
        chrom1="chr1", start1=1, end1=2, chrom2="chr2", start2=3, end2=4, name="sv_001", score=0
    )
    assert record.bedpe_field_num == 8


@pytest.mark.parametrize(
    "content",
    [
        "",
        "# comment\nbrowser foo\ntrack name=test\n",
    ],
    ids=["empty", "headers_only"],
)
def test_bedpe_source_yields_no_records(tmp_path: Path, content: str) -> None:
    """An empty file or a file containing only comment/browser/track lines yields zero records."""
    test_bedpe = tmp_path / "test.bedpe"
    test_bedpe.write_text(content)

    with BedPeSource(test_bedpe) as source:
        records = list(source)
    assert records == []


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (
            dict(chrom1="chr1", start1=-1, end1=0, chrom2="chr2", start2=0, end2=1),
            "start1 must be >= 0",
        ),
        (
            dict(chrom1="chr1", start1=0, end1=1, chrom2="chr2", start2=-1, end2=0),
            "start2 must be >= 0",
        ),
    ],
)
def test_bedpe_record_validation_negative_start(kwargs: dict[str, Any], match: str) -> None:
    """BedPeRecord raises ValueError when start1 or start2 is negative."""
    with pytest.raises(ValueError, match=match):
        BedPeRecord(**kwargs)
