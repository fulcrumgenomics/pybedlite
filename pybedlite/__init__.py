"""
Lightweight interfaces for reading and writing BED records
----------------------------------------------------------

Examples of Parsing BED files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> import pybedlite as pybed
    >>> from pathlib import Path
    >>> with pybed.reader(path=Path("infile.bed")) as in_fh:
            for record in in_fh:
                # Do work with records
                pass

Examples of Writing BED files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> import pybedlite as pybed
    >>> from pathlib import Path
    >>> # Get records from somewhere
    >>> records = []
    >>> with pybed.reader(path=Path("infile.bed")) as in_fh:
            for record in in_fh:
                records.append(record)
    >>> # Write records to somewhere
    >>> with pybed.writer(path=Path("outfile.bed"), num_fields=6) as out_fh:
            for record in records:
                out_fh.write(record)


Module Contents
~~~~~~~~~~~~~~~

The module contains the following public classes:

    - :class:`~pybedlite.bed_record.BedStrand` -- Enumeration of possible strands for a bed record
    - :class:`~pybedlite.bed_record.BedRecord` -- Lightweight class for storing information
        pertaining to a BED record.
    - :class:`~pybedlite.bed_source.BedSource` -- Reader class for parsing BED files and iterate
        over their contained records
    - :class:`~pybedlite.bed_writer.BedWriter` -- Writer class for writing BED files

The module contains the following methods:

    - :func:`pybedlite.reader` -- opens a BED file for reading.
    - :func:`pybedlite.writer` -- opens a BED file for writing.

"""

from typing import Optional

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand
from pybedlite.bed_source import BedPath
from pybedlite.bed_source import BedSource
from pybedlite.bed_writer import BedWriter


def reader(path: BedPath) -> "BedSource":
    """Returns BED file for reading. File handle will need to be opened with the

    Args:
        path: a file handle or path to the Bed to read.
    """
    return BedSource(path=path)


def writer(
    path: BedPath,
    num_fields: Optional[int] = None,
) -> "BedWriter":
    """Returns a BED file for writing.

    Args:
        path: a file handle or path to the BED to write.
        num_fields: the number of BED fields to write for each record. If none this value will be
            set to the number of fields present in the first BED record written by this object.
    """
    return BedWriter(path=path, num_fields=num_fields)


__all__ = [
    "reader",
    "writer",
    "BedSource",
    "BedWriter",
    "BedRecord",
    "BedStrand",
]
