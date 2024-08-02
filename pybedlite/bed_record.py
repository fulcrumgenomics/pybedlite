"""
Lightweight interface for reading and writing BED records
---------------------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedtools.bed_record.BedStrand` -- Enumeration of possible strands for a bed record
    - :class:`~pybedtools.bed_record.BedRecord` -- Lightweight class for storing information
        pertaining to a BED record.
"""

from __future__ import annotations

import enum
from typing import TYPE_CHECKING
from typing import ClassVar
from typing import List
from typing import Optional
from typing import Tuple
from typing import Type

import attr

if TYPE_CHECKING:
    from pybedlite.overlap_detector import Interval

"""Maximum BED fields that can be present in a well formed BED file written to specification"""
MAX_BED_FIELDS: int = 12


@enum.unique
class BedStrand(enum.Enum):
    Positive = "+"
    Negative = "-"

    @property
    def opposite(self) -> BedStrand:
        """Return the opposite strand of the current strand."""
        return BedStrand.Positive if self is BedStrand.Negative else BedStrand.Negative


@attr.s(frozen=True, auto_attribs=True, kw_only=True, slots=True)
class BedRecord:
    """Lightweight class for storing BED records. A more comprehensive description of BED format
    can be found at https://genome.ucsc.edu/FAQ/FAQformat.html#format1. Only `chrom`, `start`, and
    `end` are required.

    Attributes:
        chrom: the reference name of the interval described by the record
        start: the start coordinate in 0-based half open coordinates (inclusive)
        end: the end coordinate in 0-based half open coordinates (exclusive)
        name: the name of the bed record
        score: a score for the interval (Officially in the UCSC spec should be between 0 and 1000,
               however we do not enforce this constraint on the score. We only require that if
               defined it stores an integer)
        strand: defines the strand of the interval described by the record
        thick_start: the starting position at which the bed feature should be drawn thickly
        thick_end: the ending position at which the bed feature should be drawn thickly
        item_rgb: an RGB value of the form (R, G, B)
        block_count: the number of blocks in the bed line
        block_sizes: a list of block (exon) sizes. Number of items must correspond to
            block_count
        block_starts: a list of block (exon) starts relative to start. 0-based inclusive.
            The number of items must correspond to block_count
    """

    chrom: str
    start: int = attr.ib()
    end: int
    name: Optional[str] = None
    score: Optional[int] = None
    strand: Optional[BedStrand] = None
    thick_start: Optional[int] = attr.ib(default=None)
    thick_end: Optional[int] = None
    item_rgb: Optional[Tuple[int, int, int]] = None
    block_count: Optional[int] = attr.ib(default=None)
    block_sizes: Optional[List[int]] = None
    block_starts: Optional[List[int]] = attr.ib(default=None)

    MissingValue: ClassVar[str] = "."

    ##############
    # Validators #
    ##############

    @start.validator
    def _validate_start(self, attribute: str, value: Optional[int]) -> None:
        assert self.end > value, (
            "End of interval must be greater than start of interval."
            + f" start: {value}, end: {self.end}"
        )

    @thick_start.validator
    def _validate_thick_definitions(self, attribute: str, value: Optional[int]) -> None:
        if value is None:
            assert self.thick_end is None, "Thick end cannot be defined if thick start is not"
        else:
            assert self.thick_end is not None, "Thick start cannot be defined if thick end is not"

    @block_count.validator
    def _validate_block_counts(self, attribute: str, value: Optional[int]) -> None:
        if value is None:
            assert (
                self.block_sizes is None
            ), "Block count must be defined if block sizes is defined"
            assert (
                self.block_starts is None
            ), "Block count must be defined if block starts is defined"
        else:
            assert (
                self.block_sizes is not None
            ), "Block sizes cannot be undefined if block count is defined"
            assert (
                self.block_starts is not None
            ), "Block starts cannot be undefined if block count is defined"
            assert (
                len(self.block_sizes) == value
            ), "Number of items in block_sizes must match block_counts"
            assert (
                len(self.block_starts) == value
            ), "Number of items in block_starts must match block_counts"

    @block_starts.validator
    def _validate_bounds(self, attribute: str, value: Optional[List[int]]) -> None:
        if value is not None:
            assert 0 == value[0], "Block start at first position should be zero"
            assert self.end == self.start + value[-1] + self.block_sizes[-1], (
                "overall interval end should be equal to the last defined block's end. "
                + f"Last block end: {self.start + value[-1] + self.block_sizes[-1]}, "
                + f"Interval end: {self.end}"
            )

    @property
    def bed_field_num(self) -> int:
        if self.block_starts is not None:
            return 12
        elif self.item_rgb is not None:
            return 9
        elif self.thick_end is not None:
            return 8
        elif self.strand is not None:
            return 6
        elif self.score is not None:
            return 5
        elif self.name is not None:
            return 4
        else:
            return 3

    @property
    def bed_fields(self) -> List[str]:
        """
        Converts a BED record to a list of it's BED field string equivalents
        """
        return [
            self.chrom,
            f"{self.start}",
            f"{self.end}",
            BedRecord.MissingValue if self.name is None else self.name,
            BedRecord.MissingValue if self.score is None else f"{self.score}",
            BedRecord.MissingValue if self.strand is None else self.strand.value,
            BedRecord.MissingValue if self.thick_start is None else f"{self.thick_start}",
            BedRecord.MissingValue if self.thick_end is None else f"{self.thick_end}",
            (
                BedRecord.MissingValue
                if self.item_rgb is None
                else ",".join([f"{x}" for x in self.item_rgb])
            ),
            BedRecord.MissingValue if self.block_count is None else f"{self.block_count}",
            (
                BedRecord.MissingValue
                if self.block_sizes is None
                else ",".join([f"{x}" for x in self.block_sizes])
            ),
            (
                BedRecord.MissingValue
                if self.block_starts is None
                else ",".join([f"{x}" for x in self.block_starts])
            ),
        ]

    @property
    def refname(self) -> str:
        """The reference name of the interval described by the record."""
        return self.chrom

    @property
    def negative(self) -> bool:
        """
        True if the interval is negatively stranded, False if the interval is unstranded or
        positively stranded.
        """
        return self.strand is BedStrand.Negative

    def as_bed_line(self, number_of_output_fields: Optional[int] = None) -> str:
        """
        Converts a BED record to a tab delimited BED line equivalent, including up to the number of
        fields specified in the output line.

        Args:
            number_of_output_fields: the number of fields that should be output in the bed line.
                i.e. if you'd like a BED6 line, this should be set to 6. Etc.
        """

        assert (
            number_of_output_fields is None or 3 <= number_of_output_fields <= MAX_BED_FIELDS
        ), "BED records can only contain between 3 and 12 fields"

        number_of_output_fields = (
            self.bed_field_num if number_of_output_fields is None else number_of_output_fields
        )
        fields = self.bed_fields[:number_of_output_fields]
        return "\t".join(fields)

    @classmethod
    def from_interval(cls: Type["BedRecord"], interval: Interval) -> "BedRecord":
        """
        Construct a `BedRecord` from a `Interval` instance.

        **Note that `Interval` cannot represent a `BedRecord` with a missing strand.**
        Converting a record with no strand to `Interval` and then back to `BedRecord` will result in
        a record with **positive strand**.

        Args:
            interval: The `Interval` instance to convert.

        Returns:
            A `BedRecord` corresponding to the same region specified in the interval.
        """

        return BedRecord(
            chrom=interval.refname,
            start=interval.start,
            end=interval.end,
            strand=BedStrand.Negative if interval.negative else BedStrand.Positive,
            name=interval.name,
        )
