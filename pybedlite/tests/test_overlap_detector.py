"""Tests for :py:mod:`~pybedlite.overlap_detector`"""

from dataclasses import dataclass
from typing import List
from typing import Union

import pytest
from typing_extensions import TypeAlias

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand
from pybedlite.overlap_detector import Interval
from pybedlite.overlap_detector import OverlapDetector


def run_test(targets: List[Interval], query: Interval, results: List[Interval]) -> None:
    detector: OverlapDetector[Interval] = OverlapDetector()
    # Use add_all() to covert itself and add()
    detector.add_all(intervals=targets)
    # Test overlaps_any()
    assert detector.overlaps_any(query) == (len(results) > 0)
    # Test get_overlaps()
    assert detector.get_overlaps(query) == results


def test_same_interval() -> None:
    interval = Interval("1", 10, 100)
    run_test(targets=[interval], query=interval, results=[interval])


def test_query_wholly_contained_in_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 11, 99)
    run_test(targets=[target], query=query, results=[target])


def test_target_wholly_contained_in_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 9, 101)
    run_test(targets=[target], query=query, results=[target])


def test_target_overlaps_first_base_of_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 99, 100)
    run_test(targets=[target], query=query, results=[target])


def test_target_overlaps_last_base_of_query() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 10, 11)
    run_test(targets=[target], query=query, results=[target])


def test_query_before_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 9, 10)
    run_test(targets=[target], query=query, results=[])


def test_query_after_target() -> None:
    target = Interval("1", 10, 100)
    query = Interval("1", 100, 101)
    run_test(targets=[target], query=query, results=[])


def test_different_references() -> None:
    target = Interval("1", 10, 100)
    query = Interval("2", 10, 100)
    run_test(targets=[target], query=query, results=[])


def test_multiple_overlaps() -> None:
    interval_a = Interval("1", 10, 20)
    interval_b = Interval("1", 15, 25)
    interval_c = Interval("1", 19, 30)
    interval_d = Interval("1", 24, 35)

    # B overlaps both A and C
    run_test(targets=[interval_a, interval_c], query=interval_b, results=[interval_a, interval_c])
    # C overlaps both A and B
    run_test(targets=[interval_a, interval_b], query=interval_c, results=[interval_a, interval_b])
    # D overlaps only B and C (is after A)
    run_test(
        targets=[interval_a, interval_b, interval_c],
        query=interval_d,
        results=[interval_b, interval_c],
    )


def test_multiple_references() -> None:
    target_chr1 = Interval("1", 10, 20)
    target_chr2 = Interval("2", 10, 20)
    run_test(targets=[target_chr1, target_chr2], query=target_chr1, results=[target_chr1])
    run_test(targets=[target_chr1, target_chr2], query=target_chr2, results=[target_chr2])


def test_same_interval_twice() -> None:
    interval = Interval("1", 10, 100)
    run_test(targets=[interval, interval], query=interval, results=[interval])


def test_wholly_contained_target() -> None:
    target_inner = Interval("1", 50, 60)
    target_outer = Interval("1", 40, 80)

    run_test(
        targets=[target_inner, target_outer],
        query=target_inner,
        results=[target_outer, target_inner],
    )


def test_get_enclosing_intervals() -> None:
    a = Interval("1", 1, 250)
    b = Interval("1", 5, 30)
    c = Interval("1", 10, 99)
    d = Interval("1", 15, 19)
    e = Interval("1", 16, 20)

    detector: OverlapDetector[Interval] = OverlapDetector()
    detector.add_all([a, b, c, d, e])

    assert detector.get_enclosing_intervals(Interval("1", 10, 100)) == [a]
    assert detector.get_enclosing_intervals(Interval("1", 15, 20)) == [a, b, c]
    assert detector.get_enclosing_intervals(Interval("1", 18, 19)) == [a, b, c, d, e]
    assert detector.get_enclosing_intervals(Interval("1", 50, 99)) == [a, c]


def test_get_enclosed() -> None:
    a = Interval("1", 10, 100)
    b = Interval("1", 15, 20)
    c = Interval("1", 18, 19)
    d = Interval("1", 50, 99)

    detector: OverlapDetector[Interval] = OverlapDetector()
    detector.add_all([a, b, c, d])

    assert detector.get_enclosed(Interval("1", 1, 250)) == [a, b, c, d]
    assert detector.get_enclosed(Interval("1", 5, 30)) == [b, c]
    assert detector.get_enclosed(Interval("1", 16, 20)) == [c]
    assert detector.get_enclosed(Interval("1", 15, 19)) == [c]
    assert detector.get_enclosed(Interval("1", 10, 99)) == [b, c, d]


def test_iterable() -> None:
    a = Interval("1", 1, 250)
    b = Interval("1", 5, 30)
    c = Interval("1", 10, 99)
    d = Interval("1", 15, 19)
    e = Interval("1", 16, 20)

    detector: OverlapDetector[Interval] = OverlapDetector()
    detector.add_all([a])
    assert list(detector) == [a]
    detector.add_all([a, b, c, d, e])
    assert list(detector) == [a, a, b, c, d, e]


def test_conversion_to_interval(bed_records: List[BedRecord]) -> None:
    """
    Test that we can convert a BedRecord to an Interval.
    """

    # I don't think pytest.mark.parametrize can accept a fixture and expand over its values.
    # For loop it is.
    for record in bed_records:
        interval = Interval.from_bedrecord(record)

        assert interval.refname == record.chrom
        assert interval.start == record.start
        assert interval.end == record.end
        assert interval.negative is (record.strand is BedStrand.Negative)
        assert interval.name == record.name


def test_construction_from_interval(bed_records: List[BedRecord]) -> None:
    """
    Test that we can convert a BedRecord to an Interval and back.
    """

    # I don't think pytest.mark.parametrize can accept a fixture and expand over its values.
    # For loop it is.
    for record in bed_records:
        new_record = BedRecord.from_interval(Interval.from_bedrecord(record))

        assert new_record.chrom == record.chrom
        assert new_record.start == record.start
        assert new_record.end == record.end
        assert new_record.name == record.name

        if record.strand is None:
            assert new_record.strand is BedStrand.Positive
        else:
            assert new_record.strand is record.strand


def test_arbitrary_interval_types() -> None:
    """
    Test that an overlap detector can receive different interval-like objects and query them too.
    """

    @dataclass(eq=True, frozen=True)
    class ZeroBasedOpenEndedProtocol:
        reference_name: str
        zero_based_start: int
        zero_based_open_end: int

        @property
        def is_negative(self) -> bool:
            return False

    @dataclass(eq=True, frozen=True)
    class OneBasedProtocol:
        contig: str
        start: int
        end: int

        @property
        def reference_name(self) -> str:
            return self.contig

        @property
        def zero_based_start(self) -> int:
            """A 0-based start position."""
            return self.start - 1

        @property
        def zero_based_open_end(self) -> int:
            """A 0-based open-ended position."""
            return self.end

        @property
        def is_negative(self) -> bool:
            """True if the interval is on the negative strand, False otherwise"""
            return False

    @dataclass(eq=True, frozen=True)
    class ZeroBasedUnstranded:
        reference_name: str
        zero_based_start: int
        zero_based_open_end: int

    @dataclass(eq=True, frozen=True)
    class ZeroBasedStranded:
        reference_name: str
        zero_based_start: int
        zero_based_open_end: int
        is_negative: bool

    # Create minimal features of all supported structural types
    zero_based_protocol = ZeroBasedOpenEndedProtocol(
        reference_name="chr1", zero_based_start=1, zero_based_open_end=50
    )
    one_based_protocol = OneBasedProtocol(contig="chr1", start=10, end=60)
    zero_based_unstranded = ZeroBasedUnstranded(
        reference_name="chr1", zero_based_start=20, zero_based_open_end=70
    )
    zero_based_stranded = ZeroBasedStranded(
        reference_name="chr1", zero_based_start=30, zero_based_open_end=80, is_negative=True
    )
    # Set up an overlap detector to hold all the features we care about
    AllKinds: TypeAlias = Union[
        ZeroBasedOpenEndedProtocol, OneBasedProtocol, ZeroBasedUnstranded, ZeroBasedStranded
    ]
    features: List[AllKinds] = [
        zero_based_protocol,
        one_based_protocol,
        zero_based_unstranded,
        zero_based_stranded,
    ]
    detector: OverlapDetector[AllKinds] = OverlapDetector(features)

    assert OverlapDetector._is_negative(zero_based_protocol) is False
    assert OverlapDetector._is_negative(one_based_protocol) is False
    assert OverlapDetector._is_negative(zero_based_unstranded) is False
    assert OverlapDetector._is_negative(zero_based_stranded) is True

    # Query the overlap detector with yet another type
    assert detector.get_overlaps(Interval("chr1", 0, 1)) == []
    assert detector.get_overlaps(Interval("chr1", 0, 9)) == [zero_based_protocol]
    assert detector.get_overlaps(Interval("chr1", 11, 12)) == [
        zero_based_protocol,
        one_based_protocol,
    ]
    assert detector.get_overlaps(Interval("chr1", 21, 27)) == [
        zero_based_protocol,
        one_based_protocol,
        zero_based_unstranded,
    ]
    assert detector.get_overlaps(Interval("chr1", 32, 35)) == [
        zero_based_protocol,
        one_based_protocol,
        zero_based_unstranded,
        zero_based_stranded,
    ]
    assert detector.get_overlaps(Interval("chr1", 54, 55)) == [
        one_based_protocol,
        zero_based_unstranded,
        zero_based_stranded,
    ]
    assert detector.get_overlaps(Interval("chr1", 61, 62)) == [
        zero_based_unstranded,
        zero_based_stranded,
    ]
    assert detector.get_overlaps(Interval("chr1", 78, 79)) == [zero_based_stranded]
    assert detector.get_overlaps(Interval("chr1", 80, 81)) == []


def test_the_overlap_detector_wont_accept_a_non_hashable_feature() -> None:
    """
    Test that an overlap detector will not accept a non-hashable feature.
    """

    @dataclass  # A dataclass missing both `eq` and `frozen` does not implement __hash__.
    class ChromFeature:
        reference_name: str
        zero_based_start: int
        zero_based_open_end: int

    with pytest.raises(ValueError):
        OverlapDetector(
            [ChromFeature(reference_name="chr1", zero_based_start=0, zero_based_open_end=30)]
        )
