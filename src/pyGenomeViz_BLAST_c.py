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
    "blaMCA": "red",
    "ampC": "#b96d85",
    "ftsK": "#f87c03",
    "mph(E)": "#0b5f23",
    "msr(E)": "#8bca9d",
    "brnT": "#8b98ca",
    "relE": "#a3f107",
    "RelB/DinJ family addiction module antitoxin": "#0265fe",
    "Addiction module killer protein": "#4ff8b5",
    "putative addiction module antidote protein": "#53c4ca",
    "Toxin HigB-2": "#02992c",
    "Toxin-antitoxin system, antitoxin component": "#11e34b",
    "repB": "#500299",
    "Replication initiation protein-like protein": "#ccbeda",
    "repM": "#87709c",
    "Integrase catalytic domain-containing protein": "#fde218",
    "Resolvase/invertase-type recombinase catalytic": "#f7faa7",
    "site-specific DNA-methyltransferase (adenine-specific)": "#eca2fe",
    "Cytosine-specific methyltransferase": "#a8128f",
    "Type II restriction endonuclease, subtype P (Modular protein)": "#864fd5",
    "Restriction endonuclease": "#1244a8",
    "Type II methyltransferase M-CfrBI": "#cea3c7",
    "Transposase": "#fbe9a6",
    "Transposase IS4-like domain-containing protein": "#a1b471",
    "tnp": "#fefe02",
    "Insertion element IS402-like domain-containing": "#e1fe7f",
    "IS982 family transposase": "#83a315",
    "plasmid replication DNA-binding protein": "#caaa8b",
    "Mobilization protein": "#f69788",
    "MobA/MobL family protein": "#894602",
    "MobA/MobL protein domain-containing protein": "#93627f",
    "MobA/VirD2-like nuclease domain-containing protein": "#936b53",
    "mobC": "#e69e05",
    "mobQ": "#cd8e05",
    "C|D": "#010101",
    "D|C": "#010101",

"3-oxoacyl-[acyl-carrier-protein] synthase 3": "#e0c2b2",
"AAA+ ATPase domain-containing protein": "#0f95b0",
"abrB": "#cbfdfd",
"Adenylosuccinate lyase": "#b869b9",
"Aldehyde dehydrogenase domain-containing protein": "#f33f94",
"Amidase": "#92b46d",
"Arc family DNA-binding protein": "#cd9a4b",
"Arsenate reductase": "#54e38f",
"C2H2-type domain-containing protein": "#fb711e",
"C4": "#b7d0b5",
"CD-NTase-associated protein 12/Pycsar effector": "#a564ed",
"Cytoplasmic protein": "#d1f1e9",
"dacC": "#3cb6c6",
"dinJ": "#a8e272",
"DNA-binding protein HU": "#ea9c63",
"DNA-binding protein": "#786e85",
"DUF3387 domain-containing protein": "#95a7fb",
"Fido domain-containing protein": "#fcf7b2",
"Gamma carbonic anhydrase family protein": "#c05af3",
"Helix-turn-helix domain-containing protein": "#7acba0",
"High frequency lysogenization HflD-like protein": "#568e2f",
"Hsp20/alpha crystallin family protein": "#dc9e3c",
"HTH cro/C1-type domain-containing protein": "#62dffb",
"HTH marR-type domain-containing protein": "#f2b08e",
"Immunity protein 43 domain-containing protein": "#2f5d2e",
"Initiator Rep protein domain-containing protein": "#e378dd",
"Initiator Rep protein WH1 domain-containing protein": "#c6e00e",
"Knr4/Smi1-like domain-containing protein": "#dc3a1a",
"Lipoprotein": "#6c87fa",
"lysM": "#a94a71",
"Macro domain-containing protein": "#d0d2fa",
"Major facilitator superfamily (MFS) profile": "#6a3e2f",
"mnmA": "#b4f775",
"Nucleotidyl transferase AbiEii/AbiGii toxin family": "#a0c4a8",
"nudJ": "#f2decc",
"OmpR/PhoB-type domain-containing protein": "#7d1bb3",
"pdif": "#42a7cf",
"RelE/StbE family addiction module toxin": "#f2752e",
"Short-chain dehydrogenase/reductase SDR": "#397a9e",
"SHSP domain-containing protein": "#fe5b87",
"site-specific DNA-methyltransferase": "#96e91c",
"surE": "#bc93c3",
"Tsi6 domain-containing protein": "#516060",
"Type II restriction endonuclease, subtype P (Modular": "#fcbb0d",
"Type III restriction enzyme, res subunit": "#3ed5ba",
"vbhA": "#dc6889",
"Very short patch repair endonuclease": "#fdab48",
"yafQ": "#8bbfed",
"yflT": "#4f8a35"

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

