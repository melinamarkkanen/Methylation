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

#"rev_s39922.ctg053620l.gbff",
#"s0.ctg049933l.gbff",
#"s1.ctg045269l.gbff",
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

# RE
#	"site-specific DNA-methyltransferase (adenine-specific)": "#7A0F2C",
#	"Cytosine-specific methyltransferase": "#8F0033",
#	"Type II restriction endonuclease, subtype P (Modular protein)": "#6A1147",
#	"Restriction endonuclease": "#4C0A2B",
#	"Type II methyltransferase M-CfrBI": "#9B2E4F",
#	"site-specific DNA-methyltransferase": "#8A1A22",
#	"Type II restriction endonuclease, subtype P (Modular": "#A1232D",
#	"Type III restriction enzyme, res subunit": "#7C003F",

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




#    "ampC": "#b96d85",
#"3-oxoacyl-[acyl-carrier-protein] synthase 3": "#e0c2b2",
#"AAA+ ATPase domain-containing protein": "#0f95b0",
#"abrB": "#cbfdfd",
#"Adenylosuccinate lyase": "#b869b9",
#"Aldehyde dehydrogenase domain-containing protein": "#f33f94",
#"Amidase": "#92b46d",
#"Arc family DNA-binding protein": "#cd9a4b",
#"Arsenate reductase": "#54e38f",
#"C2H2-type domain-containing protein": "#fb711e",
#"C4": "#b7d0b5",
#"CD-NTase-associated protein 12/Pycsar effector": "#a564ed",
#"Cytoplasmic protein": "#d1f1e9",
#"dacC": "#3cb6c6",
#"DNA-binding protein HU": "#ea9c63",
#"DNA-binding protein": "#786e85",
#"DUF3387 domain-containing protein": "#95a7fb",
#"Fido domain-containing protein": "#fcf7b2",
#"Gamma carbonic anhydrase family protein": "#c05af3",
#"Helix-turn-helix domain-containing protein": "#7acba0",
#"High frequency lysogenization HflD-like protein": "#568e2f",
#"Hsp20/alpha crystallin family protein": "#dc9e3c",
#"HTH cro/C1-type domain-containing protein": "#62dffb",
#"HTH marR-type domain-containing protein": "#f2b08e",
#"Immunity protein 43 domain-containing protein": "#2f5d2e",
#"Knr4/Smi1-like domain-containing protein": "#dc3a1a",
#"Lipoprotein": "#6c87fa",
#"lysM": "#a94a71",
#"Macro domain-containing protein": "#d0d2fa",
#"Major facilitator superfamily (MFS) profile": "#6a3e2f",
#"mnmA": "#b4f775",
#"Nucleotidyl transferase AbiEii/AbiGii toxin family": "#a0c4a8",
#"nudJ": "#f2decc",
#"OmpR/PhoB-type domain-containing protein": "#7d1bb3",
#"Short-chain dehydrogenase/reductase SDR": "#397a9e",
#"SHSP domain-containing protein": "#fe5b87",
#"surE": "#bc93c3",
#"Tsi6 domain-containing protein": "#516060",
#"vbhA": "#dc6889",
#"Very short patch repair endonuclease": "#fdab48",
#"yafQ": "#8bbfed",
#"yflT": "#4f8a35"

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

