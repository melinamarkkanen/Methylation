import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"NZ_FJWZ01000025.1.gbff",
"rev_NZ_CP031449.2.gbff",
"s28039.ctg036605l.gbff",
"s5338.ctg005986l.gbff",
"rev_s2137.ctg002908l.gbff",
"s18438.ctg024208l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {
  "aadA1": "#3CB44B",
  "aadA5": "#4363D8",
  "aph(3'')-Ib": "#F58231",
  "aph(6)-Id": "#46F0F0",
  "Aspartate--tRNA ligase": "#f58231",
  "blaOXA-129": "#F54927",
  "blaOXA": "#F54927",
  "Chromosome-partitioning ATPase Soj": "#f032e6",
  "Class 1 integrase": "#F1F505",
  "Cupin type-2 domain-containing protein": "#fabebe",
  "dfrA21": "#000075",
  "dosC": "#666666",
  "DUF502 domain-containing protein": "#9a6324",
  "DUF86 domain-containing protein": "#fffac8",
  "Gamma-butyrobetaine hydroxylase-like N-terminal": "#800000",
  "Glutamine--tRNA ligase": "#aaffc3",
  "Integrase catalytic domain-containing protein": "#808000",
  "Integrase": "#F1F505",
  "Integrating conjugative element protein": "#E8D27B",
  "Integron gene cassette protein": "#DEB821",
  "Integron integrase": "#F2F7D2",
  "intI": "#F5F78F",
  "intI1": "#F1F505",
  "IS6-like element IS6100 family transposase": "#b0e0e6",
  "merA": "#B22222",
  "merD": "#6A5ACD",
  "merE": "#708090",
  "merT": "#2F4F4F",
  "N-acetyltransferase domain-containing protein": "#9932cc",
  "NTP-binding protein": "#ff4500",
  "Nucleotidyl transferase AbiEii/AbiGii toxin family": "#2e8b57",
  "NYN domain-containing protein": "#daa520",
  "orfB": "#89AFF0",
  "Polymerase nucleotidyl transferase domain-containing": "#87ceeb",
  "qacEdelta1": "#C71585",
  "Relaxase": "#00ced1",
  "repB": "#db7093",
  "Resolvase/invertase-type recombinase catalytic": "#556b2f",
  "site-specific DNA-methyltransferase": "#708090",
  "Streptomycin 3''-O-adenylyltransferase partial": "#d2691e",
  "sul1": "#A38DF0",
  "Tn3 transposase DDE domain-containing protein": "#075EA8",
  "tnp": "#2D96FC",
  "traG": "#bc8f8f",
  "Transcriptional regulator AbiEi antitoxin N-terminal": "#6495ed",
  "Transcriptional regulator": "#dc143c",
  "Transposase for transposon Tn21": "#075EA8",
  "Transposase IS4-like domain-containing protein": "#4682B4",
  "Transposase": "#4682B4",
  "Transposon DNA-invertase": "#9ACD32",
  "Transposon Tn501 resolvase": "#C0C0C0",
  "trnG": "#c71585",
  "Type ISP restriction-modification enzyme LLaBIII": "#b8860b",
  "ubiE": "#32cd32",
  "Uncharacterized protein": "#000000"
}

default_color = "#E3E3E3"

# Plot CDS and misc_feature with custom coloring
for gbff in gbff_list:
    track = gv.add_feature_track(gbff.name, gbff.get_seqid2size(), align_label=False)
    
    for ft in ["CDS", "misc_feature"]:
        for seqid, features in gbff.get_seqid2features(ft).items():
            segment = track.get_segment(seqid)
            for feature in features:
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                label = feature.qualifiers.get("gene", [""])[0]

                # Skip unlabeled or generic features like "CDS1"
                if label.startswith("CDS") and label[3:].isdigit():
                    continue
                if label.startswith("hypothethical protein"):
                    continue

                # Determine color based on gene label
                facecolor = manual_colors.get(label, default_color)


                # Only label if it's a misc_feature OR label is in manual_colors
                add_label = (ft == "misc_feature" or label in manual_colors)

                # Define max width per line
                wrap_width = 20

                # Wrap label text if needed
                if add_label and label:
                                wrapped_label = "\n".join(textwrap.wrap(label, width=wrap_width))
                else:
                                wrapped_label = None

                segment.add_feature(
                                start, end, strand,
                                plotstyle="bigarrow", lw=0.7,
                                label=wrapped_label,
                                facecolor= facecolor,
                                text_kws=dict(rotation=0, vpos="top", hpos="center"))

# Run BLAST alignment & filter by user-defined threshold
align_coords = Blast(gbff_list, seqtype="nucleotide").run()
align_coords = AlignCoord.filter(align_coords, length_thr=100, identity_thr=30)

# Plot BLAST alignment links
if len(align_coords) > 0:
    min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
    color, inverted_color = "#83B7C9", "#83B7C9"
    for ac in align_coords:
        gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident, curve=True)
    gv.set_colorbar([color, inverted_color], vmin=min_ident)

gv.savefig("december_genbank_comparison_by_blast_nucl.svg")
gv.savefig("december_genbank_comparison_by_blast_nucl.png")

