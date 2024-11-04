#!/usr/bin/env python
# coding: utf-8

try:
    from itertools import batched
except ImportError:
    from itertools import islice

    def batched(iterable, n, *, strict=False):
        # batched('ABCDEFG', 3) → ABC DEF G
        if n < 1:
            raise ValueError('n must be at least one')
        iterator = iter(iterable)
        while batch := tuple(islice(iterator, n)):
            if strict and len(batch) != n:
                raise ValueError('batched(): incomplete batch')
            yield batch
import os
from xml.etree import ElementTree as ET

import pyarrow as pa
from tqdm.auto import tqdm
import duckdb


NS_UNIPROT = "http://uniprot.org/uniprot"
UNIPROT_PROTEIN_NAME_TYPES = ("recommendedName", "alternativeName", "submittedName")


# Parse XML on the fly, and insert into database
def parse_uniprot_xml(fh):
    # https://stackoverflow.com/questions/12160418/why-is-lxml-etree-iterparse-eating-up-all-my-memory/12161185#12161185
    # get the root element:
    it = ET.iterparse(fh, events=('start', 'end'))

    _, root = next(it)

    for (k, s) in filter(
        lambda ks: (ks[0] == "end") and (ks[1].tag == "{ns}entry".format(ns="{%s}" % NS_UNIPROT)),
        it,
    ):
        assert k == "end"

        dataset = s.attrib['dataset']
        created = s.attrib['created']
        modified = s.attrib['modified']

        accessions = [
            x.text
            for x in s.findall("./accession", namespaces={'': NS_UNIPROT})
        ]

        names = [
            x.text
            for x in s.findall("./name", namespaces={'': NS_UNIPROT})
        ]

        protein_names = {
            key: list()
            for key in UNIPROT_PROTEIN_NAME_TYPES
        }
        for e in s.findall("./protein/*/fullName/..", namespaces={'': NS_UNIPROT}):
            key = e.tag.replace("{%s}" % NS_UNIPROT, '')
            value = e.find("./fullName", namespaces={'': NS_UNIPROT}).text
            protein_names[key].append(value)
        # by the XSD[1], "recommendedName" should appear at most once:
        # [1]: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd
        assert len(protein_names["recommendedName"]) <= 1

        gene_names = [
            x.text
            for x in s.findall("./gene/name", namespaces={'': NS_UNIPROT})
        ]

        subcellular_locations = [
            x.text
            for x in s.findall("./comment/subcellularLocation/location", namespaces={'': NS_UNIPROT})
        ]

        hgnc_refs = [
            x.attrib["id"]
            for x in s.findall("./dbReference[@type='HGNC']", namespaces={'': NS_UNIPROT})
        ]

        ensembl_refs = [
            {
                "tid": x.attrib["id"],
                "pid": x.find("./property[@type='protein sequence ID']", namespaces={'': NS_UNIPROT}).attrib["value"],
                "gid": x.find("./property[@type='gene ID']", namespaces={'': NS_UNIPROT}).attrib["value"],
            }
            for x in s.findall("./dbReference[@type='Ensembl']", namespaces={'': NS_UNIPROT})
        ]

        sequence = s.findall("./sequence", namespaces={'': NS_UNIPROT})[0].text
        proteinexistence = s.find("./proteinExistence", namespaces={'': NS_UNIPROT}).attrib["type"]

        # TODO: also extract the Uniprot data on features (keyed by "type")
        yield {
            "dataset": dataset,
            "created": created,
            "modified": modified,
            "accessions": accessions,
            "names": names,
            "protein_names": protein_names,
            "gene_names": gene_names,
            "subcellular_locations": subcellular_locations,
            "sequence": sequence,
            "hgnc_refs": hgnc_refs,
            "ensembl_refs": ensembl_refs,
            "proteinexistence": proteinexistence,
        }

        root.clear()


def have_uniprot_table(conn):
    return conn.sql("""
        SELECT COUNT(*) as n
        FROM duckdb_tables
        WHERE table_name = 'uniprot'
    """).fetchone()[0] == 1


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    args = ap.parse_args()

    INFN = args.infile
    DUCKDB_FN = args.outfile

    if os.path.splitext(INFN)[1].lower() == ".gz":
        import gzip
        fh = gzip.open(INFN, 'rt')
    else:
        fh = open(INFN, 'rt')

    # conn = duckdb.connect(":memory:")
    conn = duckdb.connect(DUCKDB_FN)
    base_schema = None

    it = tqdm()
    data_iter = parse_uniprot_xml(fh)
    item = next(data_iter)
    # infer schema
    if base_schema is None:
        myarrow = pa.Table.from_pylist([item])

        # first pass: use automatic inference
        # NB: arrow uses "first row" to infer types
        base_schema = pa.schema(myarrow.schema)

        # fix types being set to null upon empty lists:
        base_schema = base_schema.set(
            base_schema.get_field_index("protein_names"),
            pa.field("protein_names", pa.struct([pa.field(col, pa.list_(pa.string())) for col in UNIPROT_PROTEIN_NAME_TYPES]))
        )

    # read data, enforce schema
    myarrow = pa.Table.from_pylist([item], schema=base_schema)

    # write to database
    if have_uniprot_table(conn):
        conn.sql("""
        INSERT INTO uniprot SELECT * FROM myarrow
        """)
    else:
        conn.sql("""
        CREATE TABLE uniprot as SELECT * FROM myarrow
        """)
    it.update(1)
    conn.commit()

    # Iterate over further entries
    for items in batched(data_iter, 25):
        items = list(items)
        myarrow = pa.Table.from_pylist(items, schema=base_schema)
        conn.sql("""
        INSERT INTO uniprot SELECT * FROM myarrow
        """)
        it.update(len(items))

    conn.commit()

    # cleanup of database, infer types, add enums
    # Create basic ENUM types and convert columns:
    conn.sql("""
    CREATE TYPE dataset_type AS ENUM (
        SELECT DISTINCT dataset
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TYPE proteinexistence_type AS ENUM (
        SELECT DISTINCT proteinexistence
        FROM uniprot
    )
    """)

    # And does [this](https://github.com/duckdb/duckdb/issues/3075) work now?
    conn.sql("""
    ALTER TABLE uniprot ALTER dataset TYPE dataset_type
    """)

    conn.sql("""
    ALTER TABLE uniprot ALTER proteinexistence TYPE proteinexistence_type
    """)

    # Create date columns:
    conn.sql(r"""
    ALTER TABLE uniprot
    ADD COLUMN IF NOT EXISTS created_date DATE;

    UPDATE uniprot
    SET created_date = strptime(created, '%Y-%m-%d');

    ALTER TABLE uniprot
    DROP COLUMN IF EXISTS created;
    """)

    conn.sql(r"""
    ALTER TABLE uniprot
    ADD COLUMN IF NOT EXISTS modified_date DATE;

    UPDATE uniprot
    SET modified_date = strptime(modified, '%Y-%m-%d');

    ALTER TABLE uniprot
    DROP COLUMN IF EXISTS modified;
    """)

    # Create `primary_accession` column:
    conn.sql("""
    ALTER TABLE uniprot
    ADD COLUMN IF NOT EXISTS primary_accession VARCHAR;

    UPDATE uniprot
    SET primary_accession = accessions[1];
    """)

    # Create the LUTs
    # First, create ENUM for subcellular_location (which is a list in `uniprot` table)
    conn.sql("""
    CREATE TYPE subcellular_location_type AS ENUM (
        WITH u AS (
            SELECT unlist(subcellular_locations) as l
            FROM uniprot
        )
        SELECT DISTINCT l
        FROM u
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_subcellular_location_lut AS (
        SELECT rowid as uniprot_id, CAST(unnest(subcellular_locations) AS subcellular_location_type) as subcellular_location
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_accession_lut AS (
        SELECT rowid as uniprot_id, unnest(accessions) as accession
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_name_lut AS (
        SELECT rowid as uniprot_id, unnest(names) as name
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_gene_name_lut AS (
        SELECT rowid as uniprot_id, unnest(gene_names) as gene_name
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_hgnc_ref_lut AS (
        SELECT rowid as uniprot_id, unnest(hgnc_refs) as hgnc_ref
        FROM uniprot
    )
    """)

    conn.sql("""
    CREATE TABLE uniprot_ensembl_ref_lut AS (
        SELECT uniprot_id, ensembl_ref.* FROM (
            SELECT rowid as uniprot_id, unnest(ensembl_refs) as ensembl_ref
            FROM uniprot
        )
    )""")

    conn.commit()
    conn.close()

    return


if __name__ == "__main__":
    main()
