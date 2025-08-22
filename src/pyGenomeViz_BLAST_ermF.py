import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"CP054002.1.gbff",
"s10.ctg005298l.gbff",
"rev_s252.ctg000287l.gbff",
"s11680.ctg016811l.gbff",
"rev_s3603.ctg005124l.gbff",
"rev_s1.ctg011090l.gbff",
"s1.ctg049768l.gbff",
"s1.ctg023398l.gbff",
"rev_s3.ctg008624l.gbff",
"rev_s45902.ctg073735l.gbff",
"rev_s2030.ctg002865l.gbff",
"rev_s1.ctg021279l.gbff",
"s3.ctg005828l.gbff",
"rev_s1.ctg052735l.gbff",
"rev_s1.ctg001235l.gbff",
"s4.ctg000209l.gbff",
"s4.ctg034248l.gbff",
"rev_s4.ctg016419l.gbff",
"s30909.ctg045958l.gbff",
"rev_s854.ctg000934l.gbff",
"rev_s3.ctg021958l.gbff",
"rev_s29774.ctg043727l.gbff",
"rev_s2396.ctg003397l.gbff",
"rev_s498.ctg000705l.gbff",
"rev_s3.ctg007810l.gbff",
"rev_s25641.ctg036837l.gbff",
"s35793.ctg053218l.gbff",
"s40323.ctg061682l.gbff",
"rev_s4.ctg005206l.gbff",
"s3.ctg002540l.gbff",
"rev_s3.ctg034885l.gbff",
"s1.ctg010434l.gbff",
"s4.ctg049861l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {
  "erm(F)": "#FF0000",
  "erm": "#FF0000",
  "Transposase": "#1E3A8A",
  "Transposase DDE domain-containing protein": "#3B82F6",
  "Transposase IS4 family protein": "#1D4ED8",
  "Transposase zinc-binding domain-containing protein": "#2563EB",
  "Mobile element protein": "#1E40AF",
  "Integrase": "#FFD700",
  "Integrase catalytic domain-containing protein": "#FBBF24",
  "Tyr recombinase domain-containing protein": "#FACC15",
  "DNA methylase": "#FF69B4",
  "Putative DNA methylase": "#FF85A1",
  "Methyltransferase": "#FF1493",
  "site-specific DNA-methyltransferase": "#FF6EB4",
  "Ribosomal RNA adenine dimethylase": "#FF4C8B",
  "DNA methylase N-4/N-6 domain-containing protein": "#FF7CA8",
  "Restriction endonuclease": "#ADD8E6",
  "GmrSD restriction endonucleases N-terminal": "#B0E0E6",
  "Type III restriction-modification enzyme EcoPI Mod": "#87CEFA",
  "MvaI/BcnI restriction endonuclease family protein": "#AFEEEE",
  "IS30 family transposase": "#3742FA",
  "PF03961 family protein": "#9467BD",
  "Putative DNA-binding protein": "#8C564B",
  "catB": "#FF7F0E",
  "brxL": "#1F77B4",
  "lpoB": "#2CA02C",
  "nucS": "#5380A6",
  "Band 7 domain-containing protein": "#7F7F7F",
  "Glycosyltransferase 2-like domain-containing": "#BCBD22",
  "ABC transporter domain-containing protein": "#17BECF",
  "ABC transporter permease": "#17BECF",
  "XRE family transcriptional regulator": "#C5B0D5",
  "HTH arsR-type domain-containing protein": "#C49C94",
  "HTH cro/C1-type domain-containing protein": "#E377C2",
  "HTH deoR-type domain-containing protein": "#FF9896",
  "DNA recombination-mediator protein A": "#AEC7E8",
  "DNA repair photolyase": "#98DF8A",
  "DEAD/DEAH box helicase family protein": "#C7C7C7",
  "Helicase ATP-binding domain-containing protein": "#E5D99E",
  "Helicase/UvrB N-terminal domain-containing protein": "#C49C94",
  "AAA+ ATPase domain-containing protein": "#DBDB8D",
  "Tetratricopeptide repeat protein": "#E7CB94",
  "DUF368 domain-containing protein": "#F7B6D2",
  "DUF262 domain-containing protein": "#C5B0D5",
  "DUF4440 domain-containing protein": "#9EDAE5",
  "DUF4252 domain-containing protein": "#FFBB78",
  "DUF418 domain-containing protein": "#98DF8A",
  "DUF47 domain-containing protein": "#C49C94",
  "Fido domain-containing protein": "#1B755F",
  "VOC domain-containing protein": "#8C564B",
  "GP-PDE domain-containing protein": "#17BECF",
  "phospholipase D-like domain-containing protein": "#9467BD",
  "nucleotidyl transferase AbiEii/AbiGii toxin family": "#FF7F50",
  "ParA family protein": "#e6194b",
  "mobC": "#3cb44b",
  "relaxase/mobilization nuclease domain-containing": "#ffe119",
  "tet(X)": "#46158A",
  "winged helix-turn-helix domain-containing protein": "#f58231",
  "helix-turn-helix transcriptional regulator": "#911eb4",
  "Putative small multidrug resistance protein": "#46f0f0",
  "Multidrug transporter": "#f032e6",
  "BstEII": "#87CEFA",
  "Antitoxin HigA-1": "#bcf60c",
  "Plasmid maintenance system killer family protein": "#fabebe",
  "blaOXA": "#008080",
  "aadS": "#e6beff",
  "xerD": "#9a6324",
  "sigW": "#fffac8",
  "trnV": "#800000",
  "3-demethylubiquinone-9 3-methyltransferase": "#aaffc3",
  "Oxidoreductase": "#808000",
  "Amidophosphoribosyltransferase": "#ffd8b1",
  "Xylose isomerase-like TIM barrel domain-containing": "#000075",
  "Putative RNA 2'-phosphotransferase": "#808080",
  "Integration host factor subunit beta": "#a9a9a9",
  "N-acetyltransferase domain-containing protein": "#ff69b4",
  "phosphate transporter": "#b0c4de",
  "Phosphate transporter": "#299190"
}

#default_color = "#7d7d7c"
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
#align_coords = Blast(gbff_list, seqtype="protein").run()
align_coords = AlignCoord.filter(align_coords, length_thr=100, identity_thr=30)

# Plot BLAST alignment links
if len(align_coords) > 0:
    min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
    color, inverted_color = "#97a5c0", "red"
    for ac in align_coords:
        gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident, curve=True)
    gv.set_colorbar([color, inverted_color], vmin=min_ident)

#gv.savefig("genbank_comparison_by_blast_prot.png")
#gv.savefig("genbank_comparison_by_blast_nucl.png")

gv.savefig("genbank_comparison_by_blast_nucl.svg")
