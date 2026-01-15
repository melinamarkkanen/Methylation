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
# ARGs
    "blaOXA": "#F54927",
    "blaOXA-464": "#F54927",
    "msrA": "#A327F5",

    #"Putative membrane-bound redox modulator Alx": "#ebdf21",
    #"Diacylglycerol kinase": "#f96f9c",
    "rhlE": "#21179A",
    "Aldose 1-epimerase": "#709C6E",
    #"MOSC domain-containing protein": "#46d9c0",
    #"N-succinyldiaminopimelate-aminotransferase /": "#172b86",		# N-succinyldiaminopimelate-aminotransferase/acetylornithine transaminase
    #"Pyridoxamine 5'-phosphate oxidase family protein": "#2702d5",
    #"histidine kinase": "#c5c06c",
    #"BaiN-like insert domain-containing protein": "#0376fd",
    #"Metallophosphoesterase": "#062b3f",
    #"N-acetylmuramoyl-L-alanine amidase": "#3b5c87",
    #"DNA repair protein Rad50": "#3f71fa",
    #"Potassium transporter": "#cbd097",
    #"DUF3010 domain-containing protein": "#4a4c2b",
    #"hipA": "#7be79a",
    #"HTH cro/C1-type domain-containing protein": "#2cc4f2",
    #"DUF4376 domain-containing protein": "#7cd093",
    #"DUF302 domain-containing protein": "#dcb11a",
    #"ycaO": "#74ed3b",
    #"Glutamine amidotransferase, class I": "#33d5f5",
    #"Dihydrolipoyl dehydrogenase": "#bb6c26",
    #"DUF155 domain-containing protein": "#8479cb",
    #"AAA+ ATPase domain-containing protein": "#416add",
    "Resolvase/invertase-type recombinase catalytic domain-containing protein": "#F9E45E",	# Resolvase/invertase-type recombinase catalytic domain-containing protein
    #"Nicotinamide mononucleotide transporter": "#2f79ad",
    #"Glutathione peroxidase": "#429c83",
    #"Mutator family transposase": "#89bbd4",
    #"30S ribosomal protein S6": "#5ccf11",
    #"Acetyltransferase component of pyruvate": "#4056ef",		# Acetyltransferase component of pyruvate dehydrogenase complex
    #"HD-GYP domain-containing protein": "#aa5e69",
    #"Nucleotidyl transferase AbiEii/AbiGii toxin family": "#e8c909",	# Nucleotidyl transferase AbiEii/AbiGii toxin family protein
    "putative DNA invertase": "#e1ad4c",
    #"Benzoate transporter": "#02cce9",
    "istB": "#5E98F9",
    #"htpG": "#774879",
    #"ciaB": "#207bcb",
    #"HNH domain-containing protein": "#b64834",
    #"Bro-N domain-containing protein": "#bdde66",
    #"lysE": "#2fbe0b",
    #"DEAD/DEAH box helicase": "#83c548",
    #"algI": "#76cbcb",
    #"DUF2721 domain-containing protein": "#87f39d",
    #"PAS domain-containing protein": "#f2992c",

# MGEs
    "IS21 family transposase": "#7AC1F5",
    "IS5 family transposase": "#0C7ACC",
    "Transposase IS4-like domain-containing protein": "#A2DEF2",
    "Transposase": "#6FFCF3",
    "Phage protein": "black",
    "Site-specific tyrosine recombinase, phage integrase family (INT_XerDC domain)": "#8A8A79",	# Site-specific tyrosine recombinase, phage integrase family (INT_XerDC domain)

# Methylation
    #"site-specific DNA-methyltransferase": "#ED98D3",			# site-specific DNA-methyltransferase (adenine-specific)
    #"Serine hydroxymethyltransferase": "#AB87A0",
    #"Methyl-accepting transducer domain-containing": "#FAA07D"		# Methyl-accepting transducer domain-containing protein
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

