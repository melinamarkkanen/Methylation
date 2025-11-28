import os
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.align import Blast, AlignCoord

import textwrap

input_dir = "input_dir_mod"

desired_order = [
"M17808.1.gbff",
"s10.ctg005298l.gbff",
"CP054002.1.gbff",
"rev_s252.ctg000287l.gbff",
"rev_s3603.ctg005124l.gbff",
"s11680.ctg016811l.gbff",
"rev_s1.ctg011090l.gbff",
"rev_s3.ctg008624l.gbff",
"rev_s2030.ctg002865l.gbff",
"s1.ctg049768l.gbff",
"rev_s45902.ctg073735l.gbff",
"rev_s1.ctg052735l.gbff",
"rev_s1.ctg001235l.gbff",
"rev_s1.ctg021279l.gbff",
"s4.ctg000209l.gbff",
"s4.ctg034248l.gbff",
"s1.ctg023398l.gbff",
"s3.ctg005828l.gbff",
"s30909.ctg045958l.gbff",
"rev_s4.ctg016419l.gbff",
"rev_s854.ctg000934l.gbff",
"rev_s29774.ctg043727l.gbff",
"rev_s3.ctg021958l.gbff",
"rev_s2396.ctg003397l.gbff",
"rev_s498.ctg000705l.gbff",
"rev_s3.ctg007810l.gbff",
"rev_s25641.ctg036837l.gbff",
"s1.ctg010434l.gbff",
"s35793.ctg053218l.gbff",
"rev_s4.ctg005206l.gbff",
"s4.ctg049861l.gbff",
"s3.ctg002540l.gbff",
"s40323.ctg061682l.gbff",
"rev_s3.ctg034885l.gbff"
]

# Load and parse Genbank files in desired order
gbff_list = [Genbank(os.path.join(input_dir, fname)) for fname in desired_order]

gv = GenomeViz(track_align_type="left", fig_track_height=0.5, feature_track_ratio=0.7)
gv.set_scale_bar()
gv.set_scale_xticks(labelsize=0)

manual_colors = {
# ARGs
"erm(F)": "#F54927",
"erm": "#F54927",
"tet(X)": "#309651",
"Putative small multidrug resistance protein": "#686487",
"Multidrug transporter": "#304D18",
"catB": "#FF7800",
"blaOXA": "#92D957",
"aadS": "#6E9699",

# Mobilization
"mobC": "#F576CF",
"relaxase/mobilization nuclease domain-containing": "#9D6AAB",
"tnp": "#76A8F5",
"Plasmid maintenance system killer family protein": "#F2CBE7",
"Mobile element protein": "#DAC3DB",
"IS30 family transposase": "#4B76D1",
"xerD": "#E0CD80",
"Transposase DDE domain-containing protein": "#093696",
"Integrase catalytic domain-containing protein": "#ff0066",
"Transposase": "#ff3300",
"Integrase": "#ffff66",
"Transposase IS4 family protein": "#0249DE",
"Tyr recombinase domain-containing protein": "#ECF745",
"Transposase zinc-binding domain-containing protein": "#069BD6",

# Methylation
"DNA methylase N-4/N-6 domain-containing protein": "#F09560",
"Ribosomal RNA adenine dimethylase": "#DB5E14",
"site-specific DNA-methyltransferase": "#C93C04",
"DNA methylase": "#F7CFBE",
"Methyltransferase": "#D9803D",
"3-demethylubiquinone-9 3-methyltransferase": "#F56E49",
"Putative DNA methylase": "#F54949",

# RM
"GmrSD restriction endonucleases N-terminal": "#8ABCBF",
"Restriction endonuclease": "#5FBCC2",
"Type III restriction-modification enzyme EcoPI Mod": "#41A3BA",


"ParA family protein": "#ff0000",
"winged helix-turn-helix domain-containing protein": "#55ff00",
"nucleotidyl transferase AbiEii/AbiGii toxin family": "#00ff0a",
"helix-turn-helix transcriptional regulator": "#00ffb7",
"Band 7 domain-containing protein": "#5E8275",
"Putative DNA-binding protein": "#7A5126",
"Helicase/UvrB N-terminal domain-containing protein": "#b700ff",
"AAA+ ATPase domain-containing protein": "#ff00fa",
"ubiE": "#ff00a6",
"DNA repair photolyase": "#ff0052",
"Antitoxin HigA-1": "#d46a00",
"DEAD/DEAH box helicase family protein": "#BDF2E0",
"PF03961 family protein": "#00d46a",
"DUF368 domain-containing protein": "#0000d4",
"VOC domain-containing protein": "#6a00d4",
"DUF262 domain-containing protein": "#d4006a",
"Acid-resistance membrane protein": "#8c0000",
"doxX": "#8c3f00",
"DUF4440 domain-containing protein": "#8c8c00",
"Potassium channel domain-containing protein": "#3f8c00",
"lpoB": "#008c8c",
"Lipoprotein": "#00468c",
"Tetratricopeptide repeat protein": "#944C01",
"BstEII": "#46008c",
"TIR domain-containing protein": "#8c0046",
"Helicase ATP-binding domain-containing protein": "#525200",
"MvaI/BcnI restriction endonuclease family protein": "#005232",
"Cadmium, cobalt and zinc/H(+)-K(+) antiporter": "#005252",
"N-acetyltransferase domain-containing protein": "#520052",
"trnV": "#520032",
"ATPase components of ABC transporters with": "#00ff00",
"Glycosyltransferase 2-like domain-containing": "#00ffcc",
"Sodium:proton antiporter": "#0066ff",
"HTH arsR-type domain-containing protein": "#0000ff",
"brxL": "#6600ff",
"Amidophosphoribosyltransferase": "#cc00ff",
"DNA recombination-mediator protein A": "#ff00cc",
"ABC transporter permease": "#ff9933",
"ABC transporter domain-containing protein": "#ffff33",
"GH92 family glycosyl hydrolase": "#99ff33",
"HTH cro/C1-type domain-containing protein": "#33ffff",
"DUF47 domain-containing protein": "#3333ff",
"Phosphate transporter": "#9933ff",
"DNA polymerase III subunit gamma/tau": "#ff33ff",
"Elongation factor Ts": "#ff3399",
"histidine kinase": "#ff6666",
"Elongation factor Tu GTP binding domain protein": "#ffcc66",
"DUF418 domain-containing protein": "#66ff66",
"GP-PDE domain-containing protein": "#66ccff",
"phospholipase D-like domain-containing protein": "#6666ff",
"Nuclease": "#cc66ff",
"HTH deoR-type domain-containing protein": "#ff66cc",
"Integration host factor subunit beta": "#ff6699",
"Xylose isomerase-like TIM barrel domain-containing": "#99ffcc",
"nucS": "#99ccff",
"putative RNA 2'-phosphotransferase": "#9999ff",
"XRE family transcriptional regulator": "#cc99ff",
"Oxidoreductase": "#ff99ff",
"Abortive phage infection protein C-terminal": "#ff99cc",
"Fido domain-containing protein": "#ff9999",
"Cadmium-exporting ATPase": "#ffff99",
"Intein N-terminal splicing region": "#ccffcc",
"Large ribosomal subunit protein uL3": "#ccffff",
"Small ribosomal subunit protein uS10": "#99ccff",
"Elongation factor G": "#9999cc",
"Small ribosomal subunit protein uS7": "#cc99cc",
"30S ribosomal protein S12": "#ffccff",
"sigW": "#ffcccc",
"DUF4252 domain-containing protein": "#ffe6cc",
"oxaloacetate decarboxylase (Na(+) extruding)": "#ffffcc",
"Methylmalonyl-CoA decarboxylase subunit beta": "#ccffe6"
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
#align_coords = Blast(gbff_list, seqtype="protein").run()
align_coords = AlignCoord.filter(align_coords, length_thr=100, identity_thr=30)

# Plot BLAST alignment links
if len(align_coords) > 0:
    min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
    color, inverted_color = "#83B7C9", "#83B7C9"
    for ac in align_coords:
        gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident, curve=True)
    gv.set_colorbar([color, inverted_color], vmin=min_ident)


gv.savefig("genbank_comparison_by_blast_nucl.png")
gv.savefig("genbank_comparison_by_blast_nucl.svg")
