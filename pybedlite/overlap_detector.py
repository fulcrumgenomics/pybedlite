"""
Utility Classes for Querying Overlaps with Genomic Regions
----------------------------------------------------------


The :class:`~pybedlite.overlap_detector.OverlapDetector` class detects and returns overlaps between
a set of genomic regions and another genomic region.  The overlap detector may contain any
interval-like Python objects that have the following properties:

  * `refname` (str): The reference sequence name
  * `start` (int): A 0-based start position
  * `end` (int): A 0-based exclusive end position

This is encapsulated in the :class:`~pybedlite.overlap_detector.GenomicSpan` protocol.

Interval-like Python objects may also contain strandedness information which will be used
for sorting them in :func:`~pybedlite.overlap_detector.OverlapDetector.get_overlaps` using
the following property if it is present, otherwise assumed to be positive stranded:

  * `negative (bool)`: Whether the feature is negative stranded or not

This is encapsulated in the :class:`~pybedlite.overlap_detector.StrandedGenomicSpan` protocol.

Examples of Detecting Overlaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

. code-block:: python

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

    - :class:`~pybedlite.overlap_detector.Interval` -- Represents a region mapping to the genome
        that is 0-based and open-ended
    - :class:`~pybedlite.overlap_detector.OverlapDetector` -- Detects and returns overlaps between
        a set of genomic regions and another genomic region
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


class GenomicSpan(Protocol):
    """
    A genomic span which has protected methods that must be implemented by all subclasses to give
    a zero-based open-ended genomic span.
    """

    @property
    def refname(self) -> str:
        """A reference sequence name."""

    @property
    def start(self) -> int:
        """A 0-based start position."""

    @property
    def end(self) -> int:
        """A 0-based open-ended position."""


class StrandedGenomicSpan(GenomicSpan, Protocol):
    @property
    def negative(self) -> bool:
        """True if the interval is on the negative strand, False otherwise"""


@attr.s(frozen=True, auto_attribs=True)
class Interval:
    """A region mapping to the genome that is 0-based and open-ended

    Attributes:
        refname (str): the refname (or chromosome)
        start (int): the 0-based start position
        end (int): the 0-based end position (exclusive)
        negative (bool): true if the interval is on the negative strand, false otherwise
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


GenericGenomicSpan = TypeVar("GenericGenomicSpan", bound=Union[GenomicSpan, StrandedGenomicSpan])
"""
A generic genomic feature. This type variable is used for describing the
generic type contained within the :class:`~pybedlite.overlap_detector.OverlapDetector`.
"""


class OverlapDetector(Generic[GenericGenomicSpan], Iterable[GenericGenomicSpan]):
    """Detects and returns overlaps between a set of genomic regions and another genomic region.

    The overlap detector may contain any interval-like Python objects that have the following
    properties:
      * `chrom` or `contig` or `refname`: The reference sequence name
      * `start`: A 0-based start position
      * `end`: A 0-based exclusive end position

    Interval-like Python objects may also contain strandedness information which will be used
    for sorting them in :func:`~pybedlite.overlap_detector.OverlapDetector.get_overlaps` using
    either of the following properties if they are present:
      * `negative (bool)`: Whether or not the feature is negative stranded or not
      * `strand (BedStrand)`: The BED strand of the feature
      * `strand (str)`: The strand of the feature (`"-"` for negative)

    The same interval may be added multiple times, but only a single instance will be returned
    when querying for overlaps.

    This detector is the most efficient when all intervals are added ahead of time.
    """

    def __init__(self, intervals: Optional[Iterable[GenericGenomicSpan]] = None) -> None:
        # A mapping from the contig/chromosome name to the associated interval tree
        self._refname_to_tree: Dict[str, cr.cgranges] = {}  # type: ignore
        self._refname_to_indexed: Dict[str, bool] = {}
        self._refname_to_intervals: Dict[str, List[GenericGenomicSpan]] = {}
        if intervals is not None:
            self.add_all(intervals)

    def __iter__(self) -> Iterator[GenericGenomicSpan]:
        """Iterates over the intervals in the overlap detector."""
        return itertools.chain(*self._refname_to_intervals.values())

    def add(self, interval: GenericGenomicSpan) -> None:
        """Adds an interval to this detector.

        Args:
            interval: the interval to add to this detector
        """
        if not isinstance(interval, Hashable):
            raise ValueError(f"Interval feature is not hashable but should be: {interval}")

        refname = interval.refname
        if refname not in self._refname_to_tree:
            self._refname_to_tree[refname] = cr.cgranges()  # type: ignore
            self._refname_to_indexed[refname] = False
            self._refname_to_intervals[refname] = []

        # Append the interval to the list of intervals for this tree, keeping the index
        # of where it was inserted
        interval_idx: int = len(self._refname_to_intervals[refname])
        self._refname_to_intervals[refname].append(interval)

        # Add the interval to the tree
        tree = self._refname_to_tree[refname]
        tree.add(refname, interval.start, interval.end, interval_idx)

        # Flag this tree as needing to be indexed after adding a new interval, but defer
        # indexing
        self._refname_to_indexed[refname] = False

    def add_all(self, intervals: Iterable[GenericGenomicSpan]) -> None:
        """Adds one or more intervals to this detector.

        Args:
            intervals: the intervals to add to this detector
        """
        for interval in intervals:
            self.add(interval)

    def overlaps_any(self, interval: GenomicSpan) -> bool:
        """Determines whether the given interval overlaps any interval in this detector.

        Args:
            interval: the interval to check

        Returns:
            True if and only if the given interval overlaps with any interval in this
            detector.
        """
        refname = interval.refname
        tree = self._refname_to_tree.get(refname)
        if tree is None:
            return False
        else:
            if not self._refname_to_indexed[refname]:
                tree.index()
            try:
                next(iter(tree.overlap(refname, interval.start, interval.end)))
            except StopIteration:
                return False
            else:
                return True

    def get_overlaps(self, interval: GenomicSpan) -> List[GenericGenomicSpan]:
        """Returns any intervals in this detector that overlap the given interval.

        Args:
            interval: the interval to check

        Returns:
            The list of intervals in this detector that overlap the given interval, or the empty
            list if no overlaps exist.  The intervals will be return in ascending genomic order.
        """
        refname = interval.refname
        tree = self._refname_to_tree.get(refname)
        if tree is None:
            return []
        else:
            if not self._refname_to_indexed[refname]:
                tree.index()
            ref_intervals: List[GenericGenomicSpan] = self._refname_to_intervals[refname]
            # NB: only return unique instances of intervals
            intervals: Set[GenericGenomicSpan] = {
                ref_intervals[index]
                for _, _, index in tree.overlap(refname, interval.start, interval.end)
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
    def _negative(interval: GenomicSpan) -> bool:
        return getattr(interval, "negative", False)

    def get_enclosing_intervals(self, interval: GenomicSpan) -> List[GenericGenomicSpan]:
        """Returns the set of intervals in this detector that wholly enclose the query interval.
        i.e. `query.start >= target.start` and `query.end <= target.end`.

          Args:
              interval: the query interval
          Returns:
              The list of intervals in this  detector that enclose the query interval.
              The intervals will be returned in ascending genomic order.
        """
        results = self.get_overlaps(interval)
        return [i for i in results if interval.start >= i.start and interval.end <= i.end]

    def get_enclosed(self, interval: GenomicSpan) -> List[GenericGenomicSpan]:
        """Returns the set of intervals in this detector that are enclosed by the query
        interval.  I.e. target.start >= query.start and target.end <= query.end.

          Args:
              interval: the query interval

          Returns:
              The list of intervals in this detector that are enclosed within the query interval.
              The intervals will be return in ascending genomic order.
        """
        results = self.get_overlaps(interval)
        return [i for i in results if i.start >= interval.start and i.end <= interval.end]

    @classmethod
    def from_bed(cls, path: Path) -> "OverlapDetector":
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
