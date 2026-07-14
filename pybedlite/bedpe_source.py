"""
Reader class for BEDPE files producing BedPeRecords.
-----------------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedlite.bedpe_source.BedPeSource` -- Reader class for parsing BEDPE files and
        iterating over their contained records.
"""

from __future__ import annotations

from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from types import TracebackType
from typing import IO
from typing import Any
from typing import ContextManager
from typing import TypeVar

from pybedlite.bed_record import BedStrand
from pybedlite.bed_source import BedPath
from pybedlite.bed_source import _IOClasses
from pybedlite.bedpe_record import MIN_BEDPE_FIELDS
from pybedlite.bedpe_record import BedPeRecord

T = TypeVar("T")


class BedPeSource(ContextManager, Iterable[BedPeRecord]):
    """
    Reader for BEDPE records stored in a BEDPE file.

    Attributes:
        num_fields: the number of BEDPE fields present for records in this file. This will be set
            to the number of fields present in the first record parsed by this class. Note that
            while the BEDPE spec requires all records in a file to use the same number of fields,
            this reader does not enforce that constraint for records after the first.
    """

    def __init__(self, path: BedPath) -> None:
        """
        Initialize a BedPeSource.

        Args:
            path: Path to the BEDPE file or a file handle.
        """
        from pathlib import Path

        self._path: Path | None
        self._in_fh: IO[str] | None
        self._file_is_open: bool

        if isinstance(path, (str, Path)):
            self._path = Path(path)
            self._in_fh = None
            self._file_is_open = False
        elif isinstance(path, _IOClasses):
            self._path = None
            self._in_fh = path
            self._file_is_open = not self._in_fh.closed
        else:
            raise TypeError(f"Cannot open '{type(path)}' for reading.")

        self.num_fields: int | None = None

    def __enter__(self) -> "BedPeSource":
        """Enter context manager and open the file."""
        return self.open()

    def __exit__(
        self,
        __exc_type: type[BaseException] | None,
        __exc_value: BaseException | None,
        __traceback: TracebackType | None,
    ) -> None:
        """Exit context manager and close the file."""
        self.close()

    def open(self) -> "BedPeSource":
        """
        Open the BedPeSource file for reading.

        Must be called before iterating over the file. Make sure to close when done.
        """
        if self._in_fh is None or (not self._file_is_open and self._path is not None):
            assert self._path is not None, "Assertion present to satisfy mypy"
            self._in_fh = self._path.open("r")
            self._file_is_open = True
        else:
            if not self._file_is_open:
                raise ValueError("File must be pre-opened if filehandle specified")
        return self

    def close(self) -> None:
        """Closes the BedPeSource file. Should be called after iterating over the file."""
        if not self._file_is_open:
            raise ValueError(f"Cannot close file {self._path} if it is not already open!")
        self._file_is_open = False
        if self._in_fh is not None:
            self._in_fh.close()

    def __iter__(self) -> Iterator[BedPeRecord]:  # noqa: C901
        """Iterate over BEDPE records in the file."""

        def helper(fields: list[str], index: int, present_fn: Callable[[str], T]) -> T | None:
            if len(fields) <= index or fields[index] == BedPeRecord.MissingValue:
                return None
            return present_fn(fields[index])

        context_managed_by_iterator: bool = False
        if not self._file_is_open:
            context_managed_by_iterator = True
            self.open()

        if self._in_fh is None:
            raise ValueError("File must be opened before iterating over it!")
        for line_number, line in enumerate(self._in_fh, start=1):
            if line.startswith("#") or line.startswith("browser") or line.startswith("track"):
                continue
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            if len(fields) < MIN_BEDPE_FIELDS:
                raise ValueError(
                    f"BEDPE records must conform to specifications, which requires at least "
                    f"{MIN_BEDPE_FIELDS} input fields. On line {line_number} in {self._path} "
                    f"had only {len(fields)} fields"
                )

            num_fields = len(fields)

            if self.num_fields is None:
                self.num_fields = num_fields

            init_args: dict[str, Any] = {}

            init_args["chrom1"] = fields[0]
            init_args["start1"] = int(fields[1])
            init_args["end1"] = int(fields[2])
            init_args["chrom2"] = fields[3]
            init_args["start2"] = int(fields[4])
            init_args["end2"] = int(fields[5])
            if num_fields >= 7:
                init_args["name"] = helper(fields=fields, index=6, present_fn=lambda x: x)
            if num_fields >= 8:
                init_args["score"] = helper(fields=fields, index=7, present_fn=int)
            if num_fields >= 9:
                init_args["strand1"] = helper(fields=fields, index=8, present_fn=BedStrand)
            if num_fields >= 10:
                init_args["strand2"] = helper(fields=fields, index=9, present_fn=BedStrand)

            yield BedPeRecord(**init_args)

        if context_managed_by_iterator:
            self.close()
