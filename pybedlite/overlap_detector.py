"""
Utility Classes for Querying Overlaps with Genomic Regions
----------------------------------------------------------

The :class:`~pybedlite.overlap_detector.OverlapDetector` class detects and returns overlaps between
a set of regions and another region on a reference. The overlap detector may contain a collection
of interval-like Python objects that have the following properties:

  * `refname` (str): The reference sequence name
  * `start` (int): A 0-based start position
  * `end` (int): A 0-based half-open end position

This contract is described by the :class:`~pybedlite.overlap_detector.Span` protocol.

Interval-like Python objects may also contain strandedness information which will be used
for sorting them in :func:`~pybedlite.overlap_detector.OverlapDetector.get_overlaps` using
the following property if it is present:

  * `negative (bool)`: True if the interval is negatively stranded, False if the interval is
    unstranded or positively stranded

This contract is described by the :class:`~pybedlite.overlap_detector.StrandedSpan` protocol.

Examples of Detecting Overlaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> from pybedlite.overlap_detector import Interval, OverlapDetector
    >>> detector = OverlapDetector()
    >>> query = Interval("chr1", 2, 20)
    >>> detector.overlaps_any(query)
    False
    >>> detector.add(Interval("chr2", 1, 100))
    >>> detector.add(Interval("chr1", 21, 100))
    >>> detector.overlaps_any(query)
    False
    >>> detector.add(Interval("chr1", 1, 1))
    >>> detector.overlaps_any(query)
    True
    >>> detector.get_overlaps(query)
    [Interval("chr1", 1, 1)]
    >>> detector.add(Interval("chr1", 3, 10))
    >>> detector.overlaps_any(query)
    True
    >>> detector.get_overlaps(query)
    [Interval("chr1", 1, 1), interval("chr1", 3, 10)]

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedlite.overlap_detector.Interval` -- Represents a region mapping to a reference
        sequence that is 0-based and open-ended
    - :class:`~pybedlite.overlap_detector.OverlapDetector` -- Detects and returns overlaps between
        a set of regions and another region on a reference
"""

import itertools
from pathlib import Path
from typing import Dict
from typing import Generic
from typing import Hashable
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Optional
from typing import Protocol
from typing import Set
from typing import Type
from typing import TypeVar
from typing import Union

import attr

import cgranges as cr
from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand
from pybedlite.bed_source import BedSource


class Span(Hashable, Protocol):
    """
    A structural type for a span on a reference sequence with zero-based open-ended coordinates.
    """

    @property
    def refname(self) -> str:
        """A reference sequence name."""

    @property
    def start(self) -> int:
        """A 0-based start position."""

    @property
    def end(self) -> int:
        """A 0-based open-ended end position."""


class StrandedSpan(Span, Protocol):
    """
    A structural type for a stranded span on a reference sequence with zero-based open-ended
    coordinates.
    """

    @property
    def negative(self) -> bool:
        """
        True if the interval is negatively stranded, False if the interval is unstranded or
        positively stranded.
        """


@attr.s(frozen=True, auto_attribs=True)
class Interval:
    """A region mapping to a reference sequence that is 0-based and open-ended

    Attributes:
        refname (str): the refname (or chromosome)
        start (int): the 0-based start position
        end (int): the 0-based half-open end position
        negative (bool): true if the interval is negatively stranded, False if the interval is
            unstranded or positively stranded
        name (Optional[str]): an optional name assigned to the interval
    """

    refname: str = attr.ib()
    start: int = attr.ib()
    end: int = attr.ib()
    negative: bool = False
    name: Optional[str] = None

    def __attrs_post_init__(self) -> None:
        """Performs simple validation.

        Checks:
            - 0 <= start
            - start < end
        """
        if self.start < 0:
            raise ValueError(f"start is out of range: {self.start}")
        if self.end <= self.start:
            raise ValueError(f"end <= start: {self.end} <= {self.start}")

    def overlap(self, other: "Interval") -> int:
        """Returns the overlap between this interval and the other, or zero if there is none.

        Args:
            other (Interval): the other interval to find the overlap with
        """
        if self.refname != other.refname:
            return 0

        overlap = min(self.end, other.end) - max(self.start, other.start)
        return overlap if overlap > 0 else 0

    def length(self) -> int:
        """Returns the length of the interval."""
        return self.end - self.start

    @classmethod
    def from_bedrecord(cls: Type["Interval"], record: BedRecord) -> "Interval":
        """
        Construct an `Interval` from a `BedRecord` instance.

        Note that when the `BedRecord` does not have a specified strand, the `Interval`'s negative
        attribute is set to False. This mimics the behavior of `OverlapDetector.from_bed()` when
        reading a record that does not have a specified strand.

        Args:
            record: The `BedRecord` instance to convert.

        Returns:
            An `Interval` corresponding to the same region specified in the record.
        """
        return cls(
            refname=record.chrom,
            start=record.start,
            end=record.end,
            negative=record.strand is BedStrand.Negative,
            name=record.name,
        )

    @classmethod
    def from_ucsc(
        cls: Type["Interval"],
        string: str,
        name: Optional[str] = None,
    ) -> "Interval":
        """
        Construct an `Interval` from a UCSC "position"-formatted string.

        The "Position" format (referring to the "1-start, fully-closed" system as coordinates are
        "positioned" in the browser):

            * Written as: chr1:127140001-127140001
            * The location may optionally be followed by a parenthetically enclosed strand, e.g.
              chr1:127140001-127140001(+).
            * No spaces.
            * Includes punctuation: a colon after the chromosome, and a dash between the start and
              end coordinates.
            * When in this format, the assumption is that the coordinate is **1-start,
              fully-closed.**

        https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/  # noqa: E501

        Note that when the string does not have a specified strand, the `Interval`'s negative
        attribute is set to `False`. This mimics the behavior of `OverlapDetector.from_bed()` when
        reading a record that does not have a specified strand.

        Args:
            string: The UCSC "position"-formatted string.
            name: An optional name for the interval.

        Returns:
            An `Interval` corresponding to the same region specified in the string.
            Note that the `Interval` is **zero-based open-ended**.

        Raises:
            ValueError: If the string is not a valid UCSC position-formatted string.
        """
        try:
            if string[-1] == ")":
                interval, strand = string.rstrip(")").rsplit("(", 1)
            else:
                interval, strand = string, "+"

            contig, span = interval.rsplit(":", 1)
            start, end = span.split("-")

            return Interval(contig, int(start) - 1, int(end), negative=(strand == "-"), name=name)

        except Exception as exception:
            raise ValueError(
                f"Not a valid UCSC position-formatted string: {string}"
            ) from exception


SpanType = TypeVar("SpanType", bound=Union[Span, StrandedSpan])
"""
A generic reference sequence span. This type variable is used for describing the generic type
contained within the :class:`~pybedlite.overlap_detector.OverlapDetector`.
"""


class OverlapDetector(Generic[SpanType], Iterable[SpanType]):
    """Detects and returns overlaps between a set of regions and another region on a reference.

    The overlap detector may contain a collection of interval-like Python objects that have the
    following properties:

        * `refname` (str): The reference sequence name
        * `start` (int): A 0-based start position
        * `end` (int): A 0-based half-open end position

    Interval-like Python objects may also contain strandedness information which will be used
    for sorting them in :func:`~pybedlite.overlap_detector.OverlapDetector.get_overlaps` using
    the following property if it is present:

        * `negative (bool)`: True if the interval is negatively stranded, False if the interval is
            unstranded or positively stranded

    The same interval may be added multiple times, but only a single instance will be returned
    when querying for overlaps.

    This detector is the most efficient when all intervals are added ahead of time.
    """

    def __init__(self, intervals: Optional[Iterable[SpanType]] = None) -> None:
        # A mapping from the contig/chromosome name to the associated interval tree
        self._refname_to_tree: Dict[str, cr.cgranges] = {}  # type: ignore
        self._refname_to_indexed: Dict[str, bool] = {}
        self._refname_to_intervals: Dict[str, List[SpanType]] = {}
        if intervals is not None:
            self.add_all(intervals)

    def __iter__(self) -> Iterator[SpanType]:
        """Iterates over the intervals in the overlap detector."""
        return itertools.chain(*self._refname_to_intervals.values())

    def add(self, interval: SpanType) -> None:
        """Adds an interval to this detector.

        Args:
            interval: the interval to add to this detector
        """
        if interval.refname not in self._refname_to_tree:
            self._refname_to_tree[interval.refname] = cr.cgranges()  # type: ignore
            self._refname_to_indexed[interval.refname] = False
            self._refname_to_intervals[interval.refname] = []

        # Append the interval to the list of intervals for this tree, keeping the index
        # of where it was inserted
        interval_idx: int = len(self._refname_to_intervals[interval.refname])
        self._refname_to_intervals[interval.refname].append(interval)

        # Add the interval to the tree
        tree = self._refname_to_tree[interval.refname]
        tree.add(interval.refname, interval.start, interval.end, interval_idx)

        # Flag this tree as needing to be indexed after adding a new interval, but defer
        # indexing
        self._refname_to_indexed[interval.refname] = False

    def add_all(self, intervals: Iterable[SpanType]) -> None:
        """Adds one or more intervals to this detector.

        Args:
            intervals: the intervals to add to this detector
        """
        for interval in intervals:
            self.add(interval)

    def overlaps_any(self, interval: Span) -> bool:
        """Determines whether the given interval overlaps any interval in this detector.

        Args:
            interval: the interval to check

        Returns:
            True if and only if the given interval overlaps with any interval in this
            detector.
        """
        tree = self._refname_to_tree.get(interval.refname)
        if tree is None:
            return False
        else:
            if not self._refname_to_indexed[interval.refname]:
                tree.index()
            try:
                next(iter(tree.overlap(interval.refname, interval.start, interval.end)))
            except StopIteration:
                return False
            else:
                return True

    def get_overlaps(self, interval: Span) -> List[SpanType]:
        """Returns any intervals in this detector that overlap the given interval.

        Args:
            interval: the interval to check

        Returns:
            The list of intervals in this detector that overlap the given interval, or the empty
            list if no overlaps exist. The intervals will be returned sorted using the following
            sort keys:

                * The interval's start (ascending)
                * The interval's end (ascending)
                * The interval's strand, positive or negative (assumed to be positive if undefined)
                * The interval's reference sequence name (lexicographically)
        """
        tree = self._refname_to_tree.get(interval.refname)
        if tree is None:
            return []
        else:
            if not self._refname_to_indexed[interval.refname]:
                tree.index()
            ref_intervals: List[SpanType] = self._refname_to_intervals[interval.refname]
            # NB: only return unique instances of intervals
            intervals: Set[SpanType] = {
                ref_intervals[index]
                for _, _, index in tree.overlap(interval.refname, interval.start, interval.end)
            }
            return sorted(
                intervals,
                key=lambda intv: (
                    intv.start,
                    intv.end,
                    self._negative(intv),
                    intv.refname,
                ),
            )

    @staticmethod
    def _negative(interval: Span) -> bool:
        """
        True if the interval is negatively stranded, False if the interval is unstranded or
        positively stranded.
        """
        return getattr(interval, "negative", False)

    def get_enclosing_intervals(self, interval: Span) -> List[SpanType]:
        """Returns the set of intervals in this detector that wholly enclose the query interval.
        i.e. `query.start >= target.start` and `query.end <= target.end`.

        Args:
            interval: the query interval

        Returns:
            The list of intervals in this detector that enclose the query interval. The intervals
            will be returned sorted using the following sort keys:

                * The interval's start (ascending)
                * The interval's end (ascending)
                * The interval's strand, positive or negative (assumed to be positive if undefined)
                * The interval's reference sequence name (lexicographically)
        """
        results = self.get_overlaps(interval)
        return [i for i in results if interval.start >= i.start and interval.end <= i.end]

    def get_enclosed(self, interval: Span) -> List[SpanType]:
        """Returns the set of intervals in this detector that are enclosed by the query
        interval.  I.e. target.start >= query.start and target.end <= query.end.

        Args:
            interval: the query interval

        Returns:
            The list of intervals in this detector that are enclosed within the query interval.
            The intervals will be returned sorted using the following sort keys:

                * The interval's start (ascending)
                * The interval's end (ascending)
                * The interval's strand, positive or negative (assumed to be positive if undefined)
                * The interval's reference sequence name (lexicographically)
        """
        results = self.get_overlaps(interval)
        return [i for i in results if i.start >= interval.start and i.end <= interval.end]

    @classmethod
    def from_bed(cls, path: Path) -> "OverlapDetector[BedRecord]":
        """Builds a :class:`~pybedlite.overlap_detector.OverlapDetector` from a BED file.
        Args:
            path: the path to the BED file
        Returns:
            An overlap detector for the regions in the BED file.
        """
        detector: OverlapDetector[BedRecord] = OverlapDetector()

        for region in BedSource(path):
            detector.add(region)

        return detector
