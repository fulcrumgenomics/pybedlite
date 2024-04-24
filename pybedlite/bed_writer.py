"""
Writer class for outputting BedRecords to a file
------------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedtools.bed_source.BedWriter` -- Writer class for writing BED files
"""

from pathlib import Path
from types import TracebackType
from typing import IO
from typing import ContextManager
from typing import Iterable
from typing import Optional
from typing import Type

from pybedlite.bed_record import BedRecord
from pybedlite.bed_source import BedPath
from pybedlite.bed_source import _IOClasses

"""Maximum BED fields that can be present in a well formed BED file written to specification"""
MAX_BED_FIELDS: int = 12


class BedWriter(ContextManager):
    """Writer class for writing BED records to a file.

    Attributes:
        num_fields: The number of BED fields to report. Must be between 3 and 12.
    """

    def __init__(
        self,
        path: BedPath,
        num_fields: Optional[int] = None,
    ) -> None:
        """Instantiates a BedWriter.

        Args:
            outfile: Path specifying where to write the BED file output by this class
            num_fields: Number of BED columns to write in BED file
        """
        assert num_fields is None or (
            num_fields >= 3 and num_fields <= MAX_BED_FIELDS
        ), "BED files must contain between 3 and 12 columns"

        self._path: Optional[Path]
        self._file_handle: Optional[IO]
        self._file_is_open: bool
        if isinstance(path, (str, Path)):
            self._path = Path(path)
            self._file_handle = None
            self._file_is_open = False
        elif isinstance(path, _IOClasses):
            self._path = None
            self._file_handle = path
            self._file_is_open = not self._file_handle.closed
        self.num_fields: int = num_fields

    def __enter__(self) -> "BedWriter":
        return self.open()

    def __exit__(
        self,
        __exc_type: Type[BaseException],
        __exc_value: BaseException,
        __traceback: TracebackType,
    ) -> None:
        self.close()

    def open(self) -> "BedWriter":
        """Opens the BedWriter's file handle."""

        if self._file_handle is None or (not self._file_is_open and self._path is not None):
            self._file_handle = self._path.open("w")
            self._file_is_open = True
        else:
            assert self._file_is_open, "File must be pre-opened if filehandle specified"
        return self

    def close(self) -> None:
        """Closes the BedWriter file. Should be called after all records to write have been
        added.
        """

        assert self._file_is_open, f"Cannot close file {self._path} if it is not already open!"
        self._file_is_open = False
        self._file_handle.close()

    def write(self, record: BedRecord, truncate: bool = False, add_missing: bool = False) -> None:
        """Writes a single BedRecord to the file.

        Args:
            record: the BED record to write to the file.
            truncate: if false and a BED record is passed with more fields than the writer is set
                to output a `ValueError` will be raised. If true such a record will be written in a
                truncated fashion, with only the number of fields written by this writer.
            add_missing: if false and a BED record is passed with fewer fields than the writer is
                set to output a `ValueError` will be raised. If true such a record will be written
                in a padded fashion, with '.' output for the missing fields up to the number of
                fields written by this writer.
        """
        if self.num_fields is None:
            self.num_fields = record.bed_field_num

        if self.num_fields < record.bed_field_num and not truncate:
            raise ValueError(
                "To write a record to a BED file with fewer BED fields than are present in the "
                + "record truncate must be set to True. "
                + f"number of fields expected {self.num_fields}, "
                + f"number of fields observed: {record.bed_field_num}"
            )
        elif self.num_fields > record.bed_field_num and not add_missing:
            raise ValueError(
                "To write a record to a BED file with more fields than are present in the "
                + "record add_missing must be set to True. "
                + f"number of fields expected {self.num_fields}, "
                + f"number of fields observed: {record.bed_field_num}"
            )

        self._file_handle.write(f"{record.as_bed_line(number_of_output_fields=self.num_fields)}\n")

    def write_all(
        self, records: Iterable[BedRecord], truncate: bool = False, add_missing: bool = False
    ) -> None:
        """Writes multiple BedRecords to a file

        Arguments:
            records: the BED records to write to the file (must be iterable)
            truncate: if false and a BED record is passed with more fields than the writer is set
                to output a `ValueError` will be raised. If true such records will be written in a
                truncated fashion, with only the number of fields written by this writer
            add_missing: if false and a BED record is passed with fewer fields than the writer is
                set to output a `ValueError` will be raised. If true such records will be written
                in a padded fashion, with '.' output for the missing fields up to the number of
                fields written by this writer.
        """
        for record in records:
            self.write(record, truncate=truncate, add_missing=add_missing)
