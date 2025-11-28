#!/bin/bash

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import sys

def reverse_feature_location(feature, sequence_length):
    # Reverse the location
    new_location = FeatureLocation(
        sequence_length - feature.location.end,
        sequence_length - feature.location.start,
        strand=-1 if feature.location.strand == 1 else 1
    )
    new_feature = SeqFeature(location=new_location, type=feature.type, qualifiers=feature.qualifiers)
    return new_feature

def reverse_genbank(infile, outfile):
    for record in SeqIO.parse(infile, "genbank"):
        reversed_seq = record.seq.reverse_complement()
        reversed_features = [reverse_feature_location(f, len(record.seq)) for f in record.features]

        record.seq = reversed_seq
        record.features = reversed_features
        record.annotations["comment"] = "Reversed sequence and features"

        SeqIO.write(record, outfile, "genbank")

# Usage
reverse_genbank('s0.ctg006697l.gbff', 'rev_s0.ctg006697l.gbff')