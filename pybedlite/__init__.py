"""
Lightweight interfaces for reading and writing BED and BEDPE records.
----------------------------------------------------------------------

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

Examples of Parsing BEDPE files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> import pybedlite as pybed
    >>> from pathlib import Path
    >>> with pybed.bedpe_reader(path=Path("infile.bedpe")) as in_fh:
            for record in in_fh:
                pass

Examples of Writing BEDPE files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    >>> import pybedlite as pybed
    >>> from pathlib import Path
    >>> records = []
    >>> with pybed.bedpe_reader(path=Path("infile.bedpe")) as in_fh:
            for record in in_fh:
                records.append(record)
    >>> with pybed.bedpe_writer(path=Path("outfile.bedpe"), num_fields=10) as out_fh:
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
    - :class:`~pybedlite.bedpe_record.BedPeRecord` -- Lightweight class for storing information
        pertaining to a BEDPE record.
    - :class:`~pybedlite.bedpe_source.BedPeSource` -- Reader class for parsing BEDPE files and
        iterating over their contained records
    - :class:`~pybedlite.bedpe_writer.BedPeWriter` -- Writer class for writing BEDPE files

The module contains the following methods:

    - :func:`pybedlite.reader` -- opens a BED file for reading.
    - :func:`pybedlite.writer` -- opens a BED file for writing.
    - :func:`pybedlite.bedpe_reader` -- opens a BEDPE file for reading.
    - :func:`pybedlite.bedpe_writer` -- opens a BEDPE file for writing.

"""

from typing import Optional

from pybedlite.bed_record import BedRecord
from pybedlite.bed_record import BedStrand
from pybedlite.bed_source import BedPath
from pybedlite.bed_source import BedSource
from pybedlite.bed_writer import BedWriter
from pybedlite.bedpe_record import MAX_BEDPE_FIELDS
from pybedlite.bedpe_record import MIN_BEDPE_FIELDS
from pybedlite.bedpe_record import BedPeRecord
from pybedlite.bedpe_source import BedPeSource
from pybedlite.bedpe_writer import BedPeWriter


def reader(path: BedPath) -> "BedSource":
    """
    Return a BedSource for reading the BED file.

    Args:
        path: a file handle or path to the Bed to read.
    """
    return BedSource(path=path)


def writer(
    path: BedPath,
    num_fields: Optional[int] = None,
) -> "BedWriter":
    """
    Return a BedWriter for writing the BED file.

    Args:
        path: a file handle or path to the BED to write.
        num_fields: the number of BED fields to write for each record. If `None` this value will
            be set to the number of fields present in the first BED record written by this object.
    """
    return BedWriter(path=path, num_fields=num_fields)


def bedpe_reader(path: BedPath) -> "BedPeSource":
    """
    Return a BedPeSource for reading the BEDPE file.

    Args:
        path: a file handle or path to the BEDPE file to read.
    """
    return BedPeSource(path=path)


def bedpe_writer(
    path: BedPath,
    num_fields: Optional[int] = None,
) -> "BedPeWriter":
    """
    Return a BedPeWriter for writing the BEDPE file.

    Args:
        path: a file handle or path to the BEDPE file to write.
        num_fields: the number of BEDPE fields to write for each record. If `None` this value will
            be set to the number of fields present in the first BEDPE record written by this
            object.
    """
    return BedPeWriter(path=path, num_fields=num_fields)


__all__ = [
    "reader",
    "writer",
    "bedpe_reader",
    "bedpe_writer",
    "BedSource",
    "BedWriter",
    "BedRecord",
    "BedStrand",
    "BedPeRecord",
    "BedPeSource",
    "BedPeWriter",
    "MIN_BEDPE_FIELDS",
    "MAX_BEDPE_FIELDS",
]
