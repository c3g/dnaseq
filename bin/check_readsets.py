#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path


logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="Sample",
        readset_col="Readset",
        library_col="Library",
        runtype_col="RunType",
        run_col="Run",
        lane_col="Lane",
        fastq1_col="FASTQ1",
        fastq2_col="FASTQ2",
        bam_col="BAM",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "Sample").
            readset_col (str): The name of the column that contains the (unique)
                readset name (default "Readset").
            library_col (str): The name of the column that contains the library
                name (default "Library").
            runtype_col (str): The name of the column that contains the either the
                string PAIRED_END or SINGLE_END (default "RunType").
            run_col (str): The name of the column that contains the run name.
            lane_col (str): The name of the column that contains the lane name.
            lane_col (str): The name of the column that contains the lane name.
            fastq1_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "FASTQ1").
            fastq2_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "FASTQ2").
            bam_col (str): The name of the column that contains the second (if any)
                BAM file path (default "BAM").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._readset_col = readset_col
        self._library_col = library_col
        self._runtype_col = runtype_col
        self._run_col = run_col
        self._lane_col = lane_col
        self._lane_col = lane_col
        self._fastq1_col = fastq1_col
        self._fastq2_col = fastq2_col
        self._bam_col = bam_col
        self._single_col = single_col
        self._seen = set()
        self._has_fastq = False
        self._has_bam = False
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_readset(row)
        self._validate_sample(row)
        if not self._validate_bam(row):
            self._validate_first_fastq(row)
            self._validate_second_fastq(row)
            self._validate_fastq_pair(row)
        self._seen.add(row[self._readset_col])
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_readset(self, row):
        """Assert that the readset name exists and convert any spaces to underscores."""
        assert len(row[self._readset_col]) > 0, "Readset name is required."
        # Sanitize readset names slightly.
        row[self._readset_col] = row[self._readset_col].replace(" ", "_")

    def _validate_first_fastq(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        # assert len(row[self._fastq1_col]) > 0, "At least th   e first FASTQ file is required."
        self._validate_fastq_format(row[self._fastq1_col])

    def _validate_second_fastq(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._fastq2_col]) > 0:
            self._validate_fastq_format(row[self._fastq2_col])

    def _validate_fastq_pair(self, row):
        """Assert that read pairs have the same file extension and match runtype."""
        if row[self._runtype_col] == "PAIRED_END":
            assert self._fastq2_col in row, (
                f"Readsets file says readset {row[self._readset_col]} is RunType SINGLE_END, but ",
                f"a path is provided in the FASTQ2 column: {row[self._fastq2_col]}"
            )
        if row[self._fastq1_col] and row[self._fastq2_col]:
            row[self._runtype_col] = "PAIRED_END"
            row[self._single_col] = False
            assert (
                Path(row[self._fastq1_col]).suffixes == Path(row[self._fastq2_col]).suffixes
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._runtype_col] = "SINGLE_END"
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

    def _validate_bam(self, row):
        """Assert that if a bam (if provided) has a .bam extension"""
        if len(row[self._bam_col]) > 0:
            self._validate_bam_format(row[self._bam_col])
            return True
        else:
            return False

    def _validate_bam_format(self, filename):
        """Assert that the bam file has the .bam extension"""
        assert filename.endswith(".bam"), (
            f"The BAM file has an unrecognized extension: {filename}\n"
            f"The filename should end with *.bam"
        )

    def validate_unique_readsets(self):
        """
        Assert that the combination of readset name is unique.

        In addition to the validation, also rename the readset if more than one exists.

        """
        assert len(self._seen) == len(self.modified), "The pair of sample name and FASTQ must be unique."
        if len({pair[0] for pair in self._seen}) < len(self._seen):
            counts = Counter(pair[0] for pair in self._seen)
            seen = Counter()
            for row in self.modified:
                sample = row[self._sample_col]
                seen[sample] += 1
                if counts[sample] > 1:
                    row[self._sample_col] = f"{sample}_T{seen[sample]}"


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = handle.read(2048)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    handle.seek(0)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,

            Sample,Readset,Library,RunType,Run,Lane,Adapter1,Adapter2,QualityOffset,BED,FASTQ1,FASTQ2,BAM
            sampleA,readset1,lib0001,PAIRED_END,run100,1,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,33,,dnaseq_GM12878_chr19_R1.fastq.gz,dnaseq_GM12878_chr19_R2.fastq.gz,
            sampleA,readset2,lib0001,PAIRED_END,run100,2,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,33,,,,dnaseq_GM12878_chr19.bam

    """
    required_columns = {"Sample", "Readset", "Library", "RunType", "Run", "Lane"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            logger.critical(f"The sample sheet **must** contain the column headers: {', '.join(required_columns)}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_readsets()
    header = list(reader.fieldnames)
    header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
