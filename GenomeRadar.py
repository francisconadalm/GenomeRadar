import argparse
import os
from Bio import SeqIO
import gzip

ALLOWED_TYPES = [
    'CDS', 'tRNA', 'operon', 'regulatory', 'mobile_element', 'ncRNA', 'mRNA', 'rRNA',
    'oriT', 'protein_bind', 'repeat_region', 'tmRNA'
]

def load_positions(positions_file):
    positions = {}
    with open(positions_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                filename = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    print(
                        f"Warning: Invalid position '{parts[1]}' for file '{filename}'. Skipping."
                    )
                    continue
                if filename not in positions:
                    positions[filename] = []
                positions[filename].append(pos)
    return positions

def find_elements_around_position(
    record,
    position,
    allowed_types,
    upstream_count=20,
    downstream_count=20
):
    features = [f for f in record.features if f.type in allowed_types]
    sorted_features = sorted(features, key=lambda f: f.location.start)

    target_feature = None
    for feature in sorted_features:
        if feature.location.start <= position < feature.location.end:
            target_feature = feature
            break

    if target_feature:
        downstream_features = [target_feature]
        upstream_features = []
        remaining_features = [f for f in sorted_features if f != target_feature]
    else:
        downstream_features = []
        upstream_features = []
        remaining_features = sorted_features

    for feature in remaining_features:
        if feature.location.end <= position:
            upstream_features.append(feature)
        elif feature.location.start > position:
            downstream_features.append(feature)

    upstream_features = sorted(
        upstream_features,
        key=lambda f: f.location.end,
        reverse=True
    )
    downstream_features = sorted(
        downstream_features,
        key=lambda f: f.location.start
    )

    upstream_features = upstream_features[:upstream_count]
    downstream_features = downstream_features[:downstream_count]

    return upstream_features, downstream_features, target_feature

def process_record(record, genbank_filename, position, outfile):
    organism = record.annotations.get('organism', 'Unknown')
    lineage = ";".join(record.annotations.get('taxonomy', ['Unknown']))

    upstream_features, downstream_features, target_feature = (
        find_elements_around_position(
            record,
            position,
            ALLOWED_TYPES,
            upstream_count=20,
            downstream_count=20
        )
    )

    # If a feature contains the target position, write it first as downstream
    if target_feature:
        element_start = target_feature.location.start
        element_end = target_feature.location.end
        element = target_feature.type
        product = target_feature.qualifiers.get('product', ['Unknown'])[0]
        protein_id = target_feature.qualifiers.get('protein_id', ['Unknown'])[0]
        translation = target_feature.qualifiers.get('translation', ['Unknown'])[0]
        strand = target_feature.location.strand
        strand_symbol = '+' if strand == 1 else '-'
        locus_tag = target_feature.qualifiers.get('locus_tag', ['Unknown'])[0]

        outfile.write(
            f"{record.id}\t{genbank_filename}\t{organism}\t{lineage}\t"
            f"{position}\t{locus_tag}\t.\t.\t{element_start}\t{element_end}\t"
            f"{element}\t{product}\t{protein_id}\t{translation}\t"
            f"0\t1\t{strand_symbol}\n"
        )

    # Write upstream features
    for i, feature in enumerate(upstream_features, start=1):
        element_start = feature.location.start
        element_end = feature.location.end
        element = feature.type
        product = feature.qualifiers.get('product', ['Unknown'])[0]
        protein_id = feature.qualifiers.get('protein_id', ['Unknown'])[0]
        translation = feature.qualifiers.get('translation', ['Unknown'])[0]
        strand = feature.location.strand
        strand_symbol = '+' if strand == 1 else '-'
        distance = position - feature.location.end
        order = -i
        locus_tag = feature.qualifiers.get('locus_tag', ['Unknown'])[0]

        outfile.write(
            f"{record.id}\t{genbank_filename}\t{organism}\t{lineage}\t"
            f"{position}\t{locus_tag}\t.\t.\t{element_start}\t{element_end}\t"
            f"{element}\t{product}\t{protein_id}\t{translation}\t"
            f"{distance}\t{order}\t{strand_symbol}\n"
        )

    # Write downstream features (excluding target feature if already written)
    start_order = 2 if target_feature else 1
    for i, feature in enumerate(downstream_features, start=start_order):
        element_start = feature.location.start
        element_end = feature.location.end
        element = feature.type
        product = feature.qualifiers.get('product', ['Unknown'])[0]
        protein_id = feature.qualifiers.get('protein_id', ['Unknown'])[0]
        translation = feature.qualifiers.get('translation', ['Unknown'])[0]
        strand = feature.location.strand
        strand_symbol = '+' if strand == 1 else '-'
        distance = feature.location.start - position
        order = i
        locus_tag = feature.qualifiers.get('locus_tag', ['Unknown'])[0]

        outfile.write(
            f"{record.id}\t{genbank_filename}\t{organism}\t{lineage}\t"
            f"{position}\t{locus_tag}\t.\t.\t{element_start}\t{element_end}\t"
            f"{element}\t{product}\t{protein_id}\t{translation}\t"
            f"{distance}\t{order}\t{strand_symbol}\n"
        )

def process_local_file(filepath, positions_map, outfile):
    genbank_filename = os.path.basename(filepath)
    positions = positions_map.get(genbank_filename)

    if not positions:
        print(f"No positions found for {genbank_filename}. Skipping.")
        return

    try:
        if filepath.endswith(".gz"):
            with gzip.open(filepath, 'rt', encoding='utf-8') as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for pos in positions:
                        process_record(record, genbank_filename, pos, outfile)
        else:
            with open(filepath, 'r', encoding='utf-8') as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for pos in positions:
                        process_record(record, genbank_filename, pos, outfile)
    except Exception as e:
        print(f"Error processing {genbank_filename}: {e}")

def main(input_folder, output_file, positions_file):
    positions_map = load_positions(positions_file)

    with open(output_file, 'w', encoding='utf-8') as outfile:
        header = (
            "Contig\tGenBank_file\tOrganism\tLineage\tTarget_Position\tLocus_Tag\t"
            "Target_Start\tTarget_End\tElement_start\tElement_end\tElement\t"
            "Product\tProtein_ID\tTranslation\tDistance\tOrder\tStrand\n"
        )
        outfile.write(header)

        for filename in os.listdir(input_folder):
            if filename.endswith((".gb", ".gbff", ".gb.gz", ".gbff.gz")):
                filepath = os.path.join(input_folder, filename)
                print(f"Processing file: {filename}")
                process_local_file(filepath, positions_map, outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extracts flanking genomic features around specific positions "
            "from local GenBank files."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Folder containing GenBank files."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file."
    )
    parser.add_argument(
        "-p", "--positions",
        required=True,
        help="TSV file with GenBank filename and genomic position (one per line)."
    )

    args = parser.parse_args()
    main(args.input, args.output, args.positions)
