import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

input_dir = "input_dir_mod"

desired_order = ["MRSN15084_A_baumannii_plasmid_pdif_size_1.gbff", 
"rev_INF2_s1.ctg057230l_plasmid_pdif_size_5.gbff", 
"INF2_s1.ctg032844l_plasmid_pdif_size_1.gbff", 
"INF1_s0.ctg009177l_plasmid_pdif_size_3.gbff", 
"INF1_s0.ctg049933l_plasmid_pdif_size_2.gbff", 
"INF3_s21891.ctg024122c_plasmid_pdif_size_6.gbff", 
"rev_INF1_s0.ctg014782l_plasmid_pdif_size_1.gbff", 
"INF2_s1.ctg045269l_plasmid_pdif_size_7.gbff", 
"rev_INF2_s26768.ctg034942l_plasmid_pdif_size_1.gbff", 
"rev_INF2_s7899.ctg010523l_plasmid_pdif_size_1.gbff", 
"rev_INF2_s26928.ctg035146l_pdif_size_1.gbff", 
"rev_INF1_s0.ctg024624l_pdif_size_13.gbff", 
"rev_INF2_s1.ctg000565l_plasmid_pdif_size_1.gbff", 
"rev_INF2_s1.ctg017122l_plasmid_pdif_size_1.gbff", 
"rev_INF3_s34351.ctg037774l_pdif_size_2.gbff", 
"INF1_s0.ctg007657l_Acinetobacter_pdif_size_2.gbff", 
"rev_INF2_s342.ctg000481l_pdif_size_1.gbff", 
"rev_INF3_s294.ctg000339l_pdif_size_1.gbff", 
"rev_INF3_s10.ctg039353l_plasmid_pdif_size_1.gbff", 
"rev_INF2_s39922.ctg053620l_plasmid_pdif_size_1.gbff", 
"INF1_s42147.ctg057144l_plasmid_pdif_size_1.gbff", 
"INF1_s19014.ctg024941c_pdif_size_2.gbff", 
"rev_INF1_s0.ctg060815l_pdif_size_1.gbff", 
"INF2_s14420.ctg018836c_plasmid_pdif_size_2.gbff", 
"INF1_s29743.ctg039014l_plasmid_pdif_size_1.gbff", 
"INF1_s23052.ctg030217l_plasmid_pdif_size_1.gbff", 
"rev_INF1_s0.ctg023587l_plasmid_pdif_size_1.gbff", 
"INF1_s0.ctg032507l_plasmid_pdif_size_3.gbff", 
"INF1_s29686.ctg038937l_pdif_size_1.gbff", 
"rev_INF1_s12384.ctg016378l_plasmid_pdif_size_1.gbff", 
"rev_INF3_s15889.ctg017579c_pdif_size_3.gbff", 
"rev_INF2_s41333.ctg056107l_plasmid_pdif_size_1.gbff", 
"INF1_s3705.ctg005001l_plasmid_pdif_size_1.gbff", 
"EFF1_s36773.ctg044390l_size_1.gbff", 
"EFF3_s5267.ctg006296l_size_1.gbff", 
"EFF2_s12321.ctg014473l_size_1.gbff", 
"EFF2_s1026.ctg001124l_Aquabacterium_A_size_3.gbff", 
"rev_INF3_s8904.ctg009920l_size_1.gbff", 
"INF1_s10685.ctg014148l_size_1.gbff", 
"rev_INF2_s12179.ctg016005l_size_1.gbff", 
"rev_INF2_s7316.ctg009764l_size_1.gbff", 
"rev_INF1_s12481.ctg016498l_size_1.gbff", 
"rev_INF3_s7624.ctg008510l_size_1.gbff", 
"INF2_s2798.ctg003801l_size_1.gbff", 
"INF2_s762.ctg001037l_A_celticus_size_2.gbff", 
"rev_INF2_s8338.ctg011098l_size_1.gbff", 
"rev_INF3_s24314.ctg026770l_size_1.gbff", 
"rev_INF3_s5778.ctg006463l_size_1.gbff", 
"rev_INF1_s0.ctg006759l_size_2.gbff", 
"INF3_s8137.ctg009080l_size_1.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="center", fig_track_height=0.5, feature_track_ratio=0.4)
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
    "Transposase IS4-like domain-containing protein": "#ffe9b8",
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
}
default_color = "#cbcdce"

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

                # Only label if it's a misc_feature AND label is in manual_colors
                add_label = (ft == "misc_feature" or label in manual_colors)

                segment.add_feature(
                    start, end, strand,
                    plotstyle="bigarrow", lw=0.5,
                    label=label if add_label else None,
                    facecolor=facecolor,
                    text_kws=dict(rotation=15, vpos="top", hpos="center")
                )

# Run BLAST alignment & filter by user-defined threshold
align_coords = Blast(gbff_list, seqtype="protein").run()
align_coords = AlignCoord.filter(align_coords, length_thr=100, identity_thr=30)

# Plot BLAST alignment links
if len(align_coords) > 0:
    min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
    color, inverted_color = "#38b1ee", "red"
    for ac in align_coords:
        gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident, curve=True)
    gv.set_colorbar([color, inverted_color], vmin=min_ident)


gv.savefig("genbank_comparison_by_blast.png")