import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"blaOXA-490.gbff",
"blaOXA-464.gbff",
"CP183406.1.gbff",
"rev_CP053835.1.gbff",
"rev_CP042652.1.gbff",
"rev_s1.ctg044793l.gbff",
"rev_s0.ctg006697l.gbff",
"rev_s11272.ctg014836l.gbff",
"rev_s24657.ctg032319l.gbff",
"rev_s23410.ctg030464l.gbff",
"rev_s2037.ctg002770l.gbff",
"rev_s21411.ctg023591l.gbff",
"rev_s18593.ctg020543l.gbff",
"rev_s0.ctg019872l.gbff",
"rev_s1.ctg003290l.gbff",
"rev_s8353.ctg009320l.gbff", 
"rev_s0.ctg011225l.gbff",
"rev_s5709.ctg006388l.gbff",
"rev_s10.ctg001312l.gbff",
"rev_s16264.ctg021428l.gbff",
"rev_s4806.ctg006467l.gbff",
"rev_s46534.ctg051470l.gbff",
"rev_s42218.ctg057264l.gbff",
"rev_s38688.ctg042595l.gbff",
"CP032100.1.gbff",
"rev_s0.ctg025057l.gbff",
"rev_s1.ctg009132l.gbff",
"rev_s0.ctg000707l.gbff",
"rev_s1.ctg046449l.gbff",
"rev_s31643.ctg041501l.gbff",
"rev_s1.ctg035853l.gbff",
"rev_s27745.ctg036219l.gbff",
"rev_s0.ctg027561l.gbff",
"rev_s1.ctg008922l.gbff",
"rev_s1.ctg008129l.gbff", 
"rev_s0.ctg008951l.gbff",
"rev_s25498.ctg028069l.gbff",
"rev_s0.ctg001594l.gbff",
"rev_s0.ctg001044l.gbff",
"rev_s1.ctg029157l.gbff",
"rev_s1.ctg046422l.gbff",
"rev_s1.ctg023029l.gbff",
"rev_s1.ctg008404l.gbff",
"rev_s0.ctg020806l.gbff",
"rev_s52630.ctg058870l.gbff",
"rev_s10.ctg003382l.gbff",
"rev_s1.ctg001528l.gbff", 
"rev_s16739.ctg021805l.gbff",
"rev_s10.ctg008364l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {
    "blaOXA": "#F54927",
    "blaOXA-464": "#F54927",
    "msrA": "#A327F5",
    "rhlE": "#21179A",
    "Aldose 1-epimerase": "#709C6E",
    "Resolvase/invertase-type recombinase catalytic domain-containing protein": "#F9E45E",
    "putative DNA invertase": "#e1ad4c",
    "istB": "#5E98F9",
    "IS21 family transposase": "#7AC1F5",
    "IS5 family transposase": "#0C7ACC",
    "Transposase IS4-like domain-containing protein": "#A2DEF2",
    "Transposase": "#6FFCF3",
    "Phage protein": "black",
    "Site-specific tyrosine recombinase, phage integrase family (INT_XerDC domain)": "#8A8A79",
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

