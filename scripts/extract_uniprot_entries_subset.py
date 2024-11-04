#!/usr/bin/env python
# coding: utf-8

# Take Uniprot (species) XML dump
# Filter entries to retrieve a small subset
# https://stackoverflow.com/a/15457389

import gzip
from xml.etree import ElementTree as ET

try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(iterator, *args, **kwargs):
        return iterator


NS_UNIPROT = "http://uniprot.org/uniprot"
ENTRY_NAME = "{ns}entry".format(ns="{%s}" % NS_UNIPROT)


def filter_nodes(fh):
    it = ET.iterparse(fh, events=('start', 'end'))
    _, root = next(it)

    for event, item in tqdm(it):
        if event == 'end':
            if item.tag == ENTRY_NAME:
                hgnc_refs = [
                    x.attrib["id"]
                    for x in item.findall("./dbReference[@type='HGNC']", namespaces={'': NS_UNIPROT})
                ]
                if 'HGNC:11892' not in hgnc_refs:
                    item.clear()
                    root.remove(item)
    return root


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")

    args = ap.parse_args()

    FN = args.infile
    OUTFN = args.outfile

    with gzip.open(FN, 'rt') as fh, open(OUTFN, 'wb') as outfh:
        root = filter_nodes(fh)
        outfh.write(ET.tostring(root, encoding='utf8'))

    return


if __name__ == "__main__":
    main()
