import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"CP026207.1.gbff",
"s137.ctg000200l.gbff",
"s399.ctg000560l.gbff",
"s17.ctg041013l.gbff",
"rev_CP055486.1.gbff",
"rev_s5821.ctg007802l.gbff",
"rev_s30.ctg032113lf.gbff",
"CP021775.1.gbff",
"s21430.ctg025057l.gbff",
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {

# ARGs
  "sul1": "red",
  "aadA1": "#0B9948",
  "aadA2": "#0B9948",
  "ant(2'')-Ia": "#9477A6",
  "blaBEL-1": "#E8FFE9",
  "blaBEL": "#E5FCA9",
  "blaOXA-10": "#581845",
  "catB": "#BF5C0F",
  "catB3": "#BF5C0F",
  "qacEdelta1": "#6A2CC7",
  "N-acetyltransferase domain-containing protein": "#FCA41C",
  "dfrA1": "#808000",

# MGEs
  "Integrase": "#FCFFBF",
  "Integron gene cassette protein": "#BD9E04",
  "intI1": "#FADB41",
  "IS6-like element IS6100 family transposase": "#4157FA",
  "ISXO2-like transposase domain-containing protein": "#daa520",
  "Prophage protein": "#000000",
  "Resolvase/invertase-type recombinase catalytic": "#B2BABB",

  "Transposase for transposon Tn21": "#0B1FB5",
  "Transposase Tn21": "#0B1FB5",
  "Transposase": "#279BF5",
  "Transposon Tn501 resolvase": "#4263B8",

  "3-oxoacyl-ACP synthase": "#e6194b",
  "AAA+ ATPase domain-containing protein": "#3cb44b",
  "ABC transporter domain-containing protein": "#f58231",
  "Abi-like protein": "#911eb4",
  "AbiTii domain-containing protein": "#46f0f0",
  "citrate synthase (unknown stereospecificity)": "#9a6324",
  "CobQ/CobB/MinD/ParA nucleotide binding": "#fffac8",
  "Cupin type-2 domain-containing protein": "#800000",
  "Cytochrome c": "#aaffc3",
  "DNA-invertase hin": "#ffd8b1",
  "Dna-invertase": "#000075",
  "DUF1010 domain-containing protein": "#808080",
  "DUF3018 domain-containing protein": "#ffffff",
  "DUF86 domain-containing protein": "#a9a9a9",
  "Essential protein Yae1 N-terminal domain-containing": "#ff69b4",
  "GGDEF domain-containing protein": "#b0e0e6",
  "hrtB": "#ff7f00",
  "HTH tetR-type domain-containing protein": "#8b0000",
  "mobB": "#cd5c5c",
  "mobC": "#87ceeb",
  "Mobilization protein": "#6a5acd",
  "N-acetyltransferase domain-containing protein": "#00ced1",
  "NTP-binding protein": "#db7093",
  "Peptidase M23 domain-containing protein": "#556b2f",
  "Polymerase nucleotidyl transferase domain-containing": "#708090",
  "Secreted protein": "#5f9ea0",
  "Serpin domain-containing protein": "#bc8f8f",
  "SpoVT-AbrB domain-containing protein": "#6495ed",
  "STAS domain-containing protein": "#dc143c",
  "Transcriptional regulator": "#ffdead",
  "Transmembrane protein": "#20b2aa",
  "vapC": "#32cd32",
  "Zinc finger DksA/TraR C4-type domain-containing": "#000000",
  "istB": "#242411",
  "ydiV": "#035C49",
  "xerC": "#F7F705",
  "tnp": "#0535F7",
  "mph(E)": "#F77605",
  "msr(E)": "#05F7A2"
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

gv.savefig("genbank_comparison_by_blast_nucl.svg")
gv.savefig("genbank_comparison_by_blast_nucl.png")

