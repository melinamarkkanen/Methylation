import os
from Bio import SeqIO

input_dir = "input_dir"
output_dir = "input_dir_mod"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Process each .gbff file in the input directory
for filename in os.listdir(input_dir):
    if not filename.endswith(".gbff"):
        continue

    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, filename)

    print(f"Processing: {filename}")
    
    records = []
    for record in SeqIO.parse(input_path, "genbank"):
        for feature in record.features:
            # Only modify features that lack a 'gene' qualifier and have a 'product'
            if "gene" not in feature.qualifiers and "product" in feature.qualifiers:
                product_value = feature.qualifiers["product"][0].strip()

                # Skip if product is "hypothetical protein" (case-insensitive)
                if product_value.lower() == "hypothetical protein":
                    continue

                # Assign product as gene
                feature.qualifiers["gene"] = [product_value]
                print(f" - Added gene='{product_value}' from product")

        records.append(record)

    # Write the modified GenBank records to the output file
    SeqIO.write(records, output_path, "genbank")

print("\n✅ All .gbff files processed and saved to 'output_dir'")