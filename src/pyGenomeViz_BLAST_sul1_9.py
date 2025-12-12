import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"CP026207.1.gbff",
"AY963803.6.gbff",
"s137.ctg000200l.gbff",
"s399.ctg000560l.gbff",
"s17.ctg041013l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {

# ARGs
"ant(2'')-Ia": "#9467bd",
"catB3": "#8c564b",
"dfrA1": "#FFB30A",
"floR2": "#F2BEFA",
"mph(E)": "#016620",
"msr(E)": "#9BBA77",
"sul1": "#F54927",
"tet(G)": "#ABA28A",
"tetR(G)": "#82775B",

# MGEs
"Chaperonin/integrase fusion protein": "#CACC91",
"DNA-invertase hin": "#17becf",
"Integrase catalytic domain-containing protein": "#843c39",
"Integrase": "#C9B608",
"intI1": "#F2D335",
"IS6-like element IS6100 family transposase": "#35ADF2",
"istB": "#31a354",
"Putative transposase": "#ce6dbd",
"tnp": "#3B27F2",
"Transposase for transposon Tn21": "#a55194",
"Transposase": "#ff9896",
"xerC": "#F2F5A2",


#"AAA domain-containing protein, putative AbiEii": "#1f77b4",
#"AAA+ ATPase domain-containing protein": "#ff7f0e",
#"Abi-like protein": "#2ca02c",
#"AbiTii domain-containing protein": "#d62728",
#"citrate synthase (unknown stereospecificity)": "#7f7f7f",

#"Essential protein Yae1 N-terminal domain-containing": "#393b79",
#"GGDEF domain-containing protein": "#8c6d31",


#"N-acetyltransferase domain-containing protein": "#636363",
#"NTP-binding protein": "#b5cf6b",
#"Peptidase M23 domain-containing protein": "#9c9ede",
#"Peptidase S8/S53 domain-containing protein": "#6b6ecf",
#"Secreted protein": "#bd9e39",
#"Serpin domain-containing protein": "#fd8d3c",
#"Transcriptional regulator": "#969696",
#"Transmembrane protein": "#c7e9c0",
#"ydiV": "#9edae5"

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

