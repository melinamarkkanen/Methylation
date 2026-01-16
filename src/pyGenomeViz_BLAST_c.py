import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"rev_CP032274.1.gbff",
"s0.ctg032507l.gbff",
"rev_s12384.ctg016378l.gbff",
"rev_s0.ctg060815l.gbff",
"s19014.ctg024941c.gbff",
"s14420.ctg018836c.gbff",
"rev_s15889.ctg017579c.gbff",
"rev_s0.ctg023587l.gbff",
"s23052.ctg030217l.gbff",
"s29743.ctg039014l.gbff",
"s29686.ctg038937l.gbff",
"s42147.ctg057144l.gbff",
"rev_s38098.ctg041930l.gbff",
"rev_s34351.ctg037774l.gbff",
"s0.ctg007657l.gbff",
"rev_s342.ctg000481l.gbff",
"rev_s294.ctg000339l.gbff",
"rev_s41333.ctg056107l.gbff",
"s3705.ctg005001l.gbff",
"rev_s26928.ctg035146l.gbff",
"rev_s1.ctg000565l.gbff",
"rev_s7899.ctg010523l.gbff",
"rev_s1.ctg017122l.gbff",
"rev_s0.ctg024624l.gbff",
"rev_s1.ctg057230l.gbff",
"Acinetobacter_baumannii.gbff",
"s21891.ctg024122c.gbff",
"rev_s26768.ctg034942l.gbff",
"s0.ctg009177l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {
# ARGs
	"blaMCA": "#F54927",
	"mph(E)": "#CC80F2",
	"msr(E)": "#CF27F5",

# toxin-antitoxin
	"relE": "#00FF97",
	"brnT": "#4A8F69",
	"RelB/DinJ family addiction module antitoxin": "#98F5CF",
	"RelE/StbE family addiction module toxin": "#00ED8B",
	"dinJ": "#10C206",
	"yafQ": "#499C5B",
	"Addiction module killer protein": "#15F4B3",
	"putative addiction module antidote protein": "#1FA4C1",
	"Toxin HigB-2": "#08571F",
	"Toxin-antitoxin system, antitoxin component": "#7CA684",
	"Major facilitator superfamily (MFS) profile domain-containing protein": "#738F80",

# MGE
	"repB": "#FFB030",
	"Replication initiation protein-like protein": "#F2C79B",
	"repM": "#C77522",
	"Initiator Rep protein domain-containing protein": "#BD9E6D",
	"Initiator Rep protein WH1 domain-containing protein": "#BD9E6D",
	"RepB family plasmid replication initiator protein": "#FACB7F",
	"plasmid replication DNA-binding protein": "#FAAC7F",

	"Integrase catalytic domain-containing protein": "#F5D927",
	"Resolvase/invertase-type recombinase catalytic domain-containing protein": "#F0E28B",

	"ftsK": "#F59D22",

	"Transposase": "#223BF5",
	"Transposase IS4-like domain-containing protein": "#22C0F5",
	"tnp": "#22F1F5",

	"Insertion element IS402-like domain-containing": "#7E9AF7",
	"IS982 family transposase": "#0897CF",

	"Mobilization protein": "#E97A8F",
	"MobA/MobL family protein": "#BF4C86",
	"MobA/MobL protein domain-containing protein": "#BF4C86",
	"MobA/VirD2-like nuclease domain-containing protein": "#D18A7D",
	"mobC": "#F09BB0",

	"C|D": "#010101",
	"D|C": "#010101",

# other
	"Short-chain dehydrogenase/reductase SDR": "#B7BF4C"
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
                                facecolor=facecolor,
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

