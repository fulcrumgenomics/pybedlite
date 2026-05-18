"""
Writer class for outputting BedPeRecords to a file.
----------------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedlite.bedpe_writer.BedPeWriter` -- Writer class for writing BEDPE files.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from types import TracebackType
from typing import IO
from typing import ContextManager

from pybedlite.bed_source import BedPath
from pybedlite.bed_source import _IOClasses
from pybedlite.bedpe_record import MAX_BEDPE_FIELDS
from pybedlite.bedpe_record import MIN_BEDPE_FIELDS
from pybedlite.bedpe_record import BedPeRecord


class BedPeWriter(ContextManager):
    """
    Writer class for writing BEDPE records to a file.

    Attributes:
        num_fields: The number of BEDPE fields to report. Must be between 6 and 10. If not set
            at construction time, this is inferred from the first record passed to :meth:`write`.
            May be changed between calls to :meth:`write`, but doing so will affect all subsequent
            records written.
    """

    def __init__(
        self,
        path: BedPath,
        num_fields: int | None = None,
    ) -> None:
        """
        Instantiate a BedPeWriter.

        Args:
            path: Path specifying where to write the BEDPE file output by this class
            num_fields: Number of BEDPE columns to write in BEDPE file
        """
        if num_fields is not None and (
            num_fields < MIN_BEDPE_FIELDS or num_fields > MAX_BEDPE_FIELDS
        ):
            raise ValueError(
                f"BEDPE files must contain between {MIN_BEDPE_FIELDS} and {MAX_BEDPE_FIELDS} "
                f"columns, got {num_fields}"
            )

        self._path: Path | None
        self._file_handle: IO[str] | None
        self._file_is_open: bool
        if isinstance(path, (str, Path)):
            self._path = Path(path)
            self._file_handle = None
            self._file_is_open = False
        elif isinstance(path, _IOClasses):
            self._path = None
            self._file_handle = path
            self._file_is_open = not self._file_handle.closed
        else:
            raise TypeError(f"Cannot open '{type(path)}' for writing.")
        self.num_fields: int | None = num_fields

    def __enter__(self) -> "BedPeWriter":
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

    def open(self) -> "BedPeWriter":
        """Open the BedPeWriter's file handle."""
        if self._file_handle is None or (not self._file_is_open and self._path is not None):
            assert self._path is not None, "Assertion present to satisfy mypy"
            self._file_handle = self._path.open("w")
            self._file_is_open = True
        else:
            if not self._file_is_open:
                raise ValueError("File must be pre-opened if filehandle specified")
        return self

    def close(self) -> None:
        """
        Close the BedPeWriter file.

        Should be called after all records to write have been added.
        """
        if not self._file_is_open:
            raise ValueError(f"Cannot close file {self._path} if it is not already open!")
        self._file_is_open = False
        if self._file_handle is not None:
            self._file_handle.close()

    def write(
        self, record: BedPeRecord, truncate: bool = False, add_missing: bool = False
    ) -> None:
        """
        Write a single BedPeRecord to the file.

        Args:
            record: the BEDPE record to write to the file.
            truncate: if false and a BEDPE record is passed with more fields than the writer is set
                to output a `ValueError` will be raised. If true such a record will be written in a
                truncated fashion, with only the number of fields written by this writer.
            add_missing: if false and a BEDPE record is passed with fewer fields than the writer is
                set to output a `ValueError` will be raised. If true such a record will be written
                in a padded fashion, with '.' output for the missing fields up to the number of
                fields written by this writer.
        """
        if not self._file_is_open:
            raise ValueError("Cannot write to a closed BedPeWriter!")

        if self.num_fields is None:
            self.num_fields = record.bedpe_field_num

        if self.num_fields < record.bedpe_field_num and not truncate:
            raise ValueError(
                "To write a record to a BEDPE file with fewer BEDPE fields than are present in "
                "the record truncate must be set to True. "
                f"number of fields expected {self.num_fields}, "
                f"number of fields observed: {record.bedpe_field_num}"
            )
        elif self.num_fields > record.bedpe_field_num and not add_missing:
            raise ValueError(
                "To write a record to a BEDPE file with more fields than are present in the "
                "record add_missing must be set to True. "
                f"number of fields expected {self.num_fields}, "
                f"number of fields observed: {record.bedpe_field_num}"
            )

        if self._file_handle is None:
            raise ValueError("File must be opened before writing!")
        self._file_handle.write(
            f"{record.as_bedpe_line(number_of_output_fields=self.num_fields)}\n"
        )

    def write_all(
        self,
        records: Iterable[BedPeRecord],
        truncate: bool = False,
        add_missing: bool = False,
    ) -> None:
        """
        Write multiple BedPeRecords to a file.

        Args:
            records: the BEDPE records to write to the file (must be iterable)
            truncate: if false and a BEDPE record is passed with more fields than the writer is set
                to output a `ValueError` will be raised. If true such records will be written in a
                truncated fashion, with only the number of fields written by this writer
            add_missing: if false and a BEDPE record is passed with fewer fields than the writer is
                set to output a `ValueError` will be raised. If true such records will be written
                in a padded fashion, with '.' output for the missing fields up to the number of
                fields written by this writer.
        """
        for record in records:
            self.write(record, truncate=truncate, add_missing=add_missing)
