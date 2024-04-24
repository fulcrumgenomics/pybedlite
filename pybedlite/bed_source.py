"""
Reader class for BED files producing BedRecords
-----------------------------------------------

Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedtools.bed_source.BedSource` -- Reader class for parsing BED files and iterate
        over their contained records
"""

import io
from pathlib import Path
from types import TracebackType
from typing import IO
from typing import Any
from typing import Callable
from typing import ContextManager
from typing import Dict
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple
from typing import Type
from typing import TypeVar
from typing import Union

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand

"""The classes that should be treated as file-like classes"""
_IOClasses = (io.TextIOBase, io.BufferedIOBase, io.RawIOBase, io.IOBase)

"""The valid base classes for opening a BED file."""
BedPath = Union[Path, str, IO[Any]]

T = TypeVar("T")


class BedSource(ContextManager, Iterable[BedRecord]):
    """Reader for BED records stored in a BED file

    Attributes:
        num_fields: the number of BED fields present for records in this file. This will be set to
            the number of fields present in the first record parsed by this class. Note that while
            the official BED spec indicates that all BED records in the same file should be written
            with the same number of fields, we do not enforce that this is the case in the BED
            files this reader parses.
    """

    def __init__(self, path: BedPath) -> None:
        self._path: Optional[Path]
        self._in_fh: Optional[IO]
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

        self.num_fields: Optional[int] = None

    def __enter__(self) -> "BedSource":
        return self.open()

    def __exit__(
        self,
        __exc_type: Type[BaseException],
        __exc_value: BaseException,
        __traceback: TracebackType,
    ) -> None:
        self.close()

    def open(self) -> "BedSource":
        """Opens the BedSources file for reading. Must be called before iterating over the file.
        Make sure to close when done.
        """
        if self._in_fh is None or (not self._file_is_open and self._path is not None):
            self._in_fh = self._path.open("r")
            self._file_is_open = True
        else:
            assert self._file_is_open, "File must be pre-opened if filehandle specified"
        return self

    def close(self) -> None:
        """Closes the BedSources file. Should be called after iterating over the file."""
        assert self._file_is_open, f"Cannot close file {self._path} if it is not already open!"
        self._file_is_open = False
        self._in_fh.close()

    def __iter__(self) -> Iterator[BedRecord]:
        def helper(fields: List[str], index: int, present_fn: Callable[[str], T]) -> Optional[T]:
            if len(fields) <= index or fields[index] == BedRecord.MissingValue:
                return None
            return present_fn(fields[index])

        def parse_rgb(x: str) -> Tuple[int, int, int]:
            # This is fairly verbose for what it's doing, but it makes mypy happy because
            # if converted to a tuple directly there's no bounds on its size.
            rgb_split = [int(x) for x in fields[8].split(",")]
            assert (
                len(rgb_split) == 3
            ), "item_rgb, if defined, must contain 3 comma separated integers"
            r, g, b = rgb_split
            return (r, g, b)

        context_managed_by_iterator: bool = False
        if not self._file_is_open:
            context_managed_by_iterator = True
            self.open()
            self._file_is_open = True

        for i, line in enumerate(self._in_fh):
            # Skip header lines
            if line.startswith("#") or line.startswith("browser") or line.startswith("track"):
                continue
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            assert len(fields) >= 3, (
                "BED records must conform to specifications, which requires at least 3 input "
                + f"fields. On line {i} in {self._path} had only {len(fields)} fields"
            )

            num_fields = len(fields)

            if self.num_fields is None:
                self.num_fields = num_fields

            init_args: Dict[str, Any] = {}

            init_args["chrom"] = fields[0]
            init_args["start"] = int(fields[1])
            init_args["end"] = int(fields[2])
            if num_fields >= 10:
                init_args["block_count"] = helper(fields=fields, index=9, present_fn=int)
                init_args["block_sizes"] = helper(
                    fields=fields, index=10, present_fn=lambda x: [int(y) for y in x.split(",")]
                )
                init_args["block_starts"] = helper(
                    fields=fields, index=11, present_fn=lambda x: [int(y) for y in x.split(",")]
                )
            if num_fields >= 9:
                init_args["item_rgb"] = helper(fields=fields, index=8, present_fn=parse_rgb)
            if num_fields >= 7:
                init_args["thick_start"] = helper(fields=fields, index=6, present_fn=int)
                init_args["thick_end"] = helper(fields=fields, index=7, present_fn=int)
            if num_fields >= 6:
                init_args["strand"] = helper(fields=fields, index=5, present_fn=BedStrand)
            if num_fields >= 5:
                init_args["score"] = helper(fields=fields, index=4, present_fn=int)
            if num_fields >= 4:
                init_args["name"] = helper(fields=fields, index=3, present_fn=lambda x: x)

            yield BedRecord(**init_args)

        if context_managed_by_iterator:
            self.close()
            self._file_is_open = False
