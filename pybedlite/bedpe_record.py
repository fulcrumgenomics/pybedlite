"""
Lightweight interface for reading and writing BEDPE records.
------------------------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedlite.bedpe_record.BedPeRecord` -- Lightweight class for storing information
        pertaining to a BEDPE record.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import ClassVar

from pybedlite.bed_record import BedStrand

_logger = logging.getLogger(__name__)

#: Minimum BEDPE fields that must be present in a well formed BEDPE file.
MIN_BEDPE_FIELDS: int = 6

#: Maximum BEDPE fields that can be present in a well formed BEDPE file written to specification.
MAX_BEDPE_FIELDS: int = 10


@dataclass(frozen=True, kw_only=True, slots=True)
class BedPeRecord:
    """
    Lightweight class for storing BEDPE records.

    BEDPE format describes a pair of genomic intervals, commonly used for structural variants and
    paired-end sequencing data. A more comprehensive description of BEDPE format can be found at
    https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format. Only
    the first six fields are required.

    Attributes:
        chrom1: the reference name of the first interval
        start1: the start coordinate of the first interval in 0-based half open coordinates
            (inclusive)
        end1: the end coordinate of the first interval in 0-based half open coordinates (exclusive)
        chrom2: the reference name of the second interval
        start2: the start coordinate of the second interval in 0-based half open coordinates
            (inclusive)
        end2: the end coordinate of the second interval in 0-based half open coordinates
            (exclusive)
        name: the name of the record
        score: a score for the record (integer)
        strand1: the strand of the first interval
        strand2: the strand of the second interval
    """

    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str | None = None
    score: int | None = None
    strand1: BedStrand | None = None
    strand2: BedStrand | None = None

    MissingValue: ClassVar[str] = "."

    ##############
    # Validators #
    ##############
    def __post_init__(self) -> None:
        """
        Validate BEDPE record constraints.

        Raises:
            ValueError: If any validation checks fail
        """
        if self.start1 < 0:
            raise ValueError(f"start1 must be >= 0, got {self.start1}")
        if self.start2 < 0:
            raise ValueError(f"start2 must be >= 0, got {self.start2}")
        if self.end1 <= self.start1:
            raise ValueError(
                f"End of first interval must be greater than start of first interval. "
                f"start1: {self.start1}, end1: {self.end1}"
            )
        if self.end2 <= self.start2:
            raise ValueError(
                f"End of second interval must be greater than start of second interval. "
                f"start2: {self.start2}, end2: {self.end2}"
            )
        if self.strand2 is not None and self.strand1 is None:
            # Log rather than raise: the BEDPE spec permits strand1="." alongside a defined
            # strand2, so this state is valid and round-trips correctly. We log because the
            # combination is unusual and likely a caller mistake.
            _logger.warning(
                "strand2 is set but strand1 is not. While spec-compliant (strand1 will be "
                "serialized as '.'), this combination is unusual."
            )

    @property
    def bedpe_field_num(self) -> int:
        """The number of BEDPE fields that are defined in this record."""
        if self.strand2 is not None:
            return 10  # strand1=None is valid; it serializes as "." (see bedpe_fields)
        elif self.strand1 is not None:
            return 9
        elif self.score is not None:
            return 8
        elif self.name is not None:
            return 7
        else:
            return 6

    @property
    def bedpe_fields(self) -> list[str]:
        """The list of field string equivalents for this record, always ten elements long."""
        return [
            self.chrom1,
            f"{self.start1}",
            f"{self.end1}",
            self.chrom2,
            f"{self.start2}",
            f"{self.end2}",
            BedPeRecord.MissingValue if self.name is None else self.name,
            BedPeRecord.MissingValue if self.score is None else f"{self.score}",
            BedPeRecord.MissingValue if self.strand1 is None else self.strand1.value,
            BedPeRecord.MissingValue if self.strand2 is None else self.strand2.value,
        ]

    def as_bedpe_line(self, number_of_output_fields: int | None = None) -> str:
        """
        Converts a BEDPE record to a tab-delimited BEDPE line.

        Args:
            number_of_output_fields: the number of fields that should be output. Must be between
                6 and 10. If None, the number of fields is inferred from the record. If less than
                :attr:`bedpe_field_num`, trailing fields are silently truncated.
                No error is raised.

        Returns:
            A tab-delimited string representation of this record with the specified number of
            fields.

        Raises:
            ValueError: If number_of_output_fields is not between 6 and 10
        """
        if number_of_output_fields is not None and (
            number_of_output_fields < MIN_BEDPE_FIELDS
            or number_of_output_fields > MAX_BEDPE_FIELDS
        ):
            raise ValueError(
                f"BEDPE records can only contain between {MIN_BEDPE_FIELDS} and "
                f"{MAX_BEDPE_FIELDS} fields, got {number_of_output_fields}"
            )
        number_of_output_fields = (
            self.bedpe_field_num if number_of_output_fields is None else number_of_output_fields
        )
        fields = self.bedpe_fields[:number_of_output_fields]
        return "\t".join(fields)
