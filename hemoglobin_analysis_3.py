from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation, SeqFeature
from io import StringIO
import gzip

# Set your email for Entrez
Entrez.email = "your@email.com."


# Defining a function to fetch and save sequences
def fetch_and_save(protein_id, file_format, output_filename):
    handle = Entrez.efetch(db="protein", id=protein_id, rettype=file_format, retmode="text")
    with open(output_filename, "w") as out_file:
        out_file.write(handle.read())
    handle.close()
    print(f"Saved {output_filename}")

# 1. Fetching and saving hemoglobin sequences from NCBI
print("1. Fetching hemoglobin sequences from NCBI")
hemoglobin_id = "NP_000509"  # Beta globin
fetch_and_save(hemoglobin_id, "gb", "hemoglobin.gb")
fetch_and_save(hemoglobin_id, "fasta", "hemoglobin.fasta")

# 2. Parsing and iterating through sequence files fasta and genbank
print("\n2. Parsing and iterating through sequence files")
print("GenBank file:")
for seq_record in SeqIO.parse("hemoglobin.gb", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

print("\nFASTA file:")
for seq_record in SeqIO.parse("hemoglobin.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


print("\nFASTA file:")
for seq_record in SeqIO.parse("hemoglobin.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

# 3. Using iteration to extract information 
print("\n3. Extracting information using list comprehension")
identifiers = [seq_record.id for seq_record in SeqIO.parse("hemoglobin.fasta", "fasta")]
print(f"Identifiers: {identifiers}")

# 4. Getting the first record 
print("\n4. Getting the first record")
first_record = next(SeqIO.parse("hemoglobin.gb", "genbank"))
print(f"First record ID: {first_record.id}")

# 5. Getting a list of records and last record
print("\n5. Getting a list of records")
records = list(SeqIO.parse("hemoglobin.fasta", "fasta"))
print(f"Number of records: {len(records)}")
print(f"Last record ID: {records[-1].id}")

# 6. Extracting annotations 
print("\n6. Extracting annotations")
gb_iterator = SeqIO.parse("hemoglobin.gb", "genbank")
for seq_record in gb_iterator:
    print(f"ID: {seq_record.id}")
    print(f"Sequence: {seq_record.seq[:10]}...")
    print(f"Description: {seq_record.description}")
    print("Annotations:")
    for key, value in seq_record.annotations.items():
        print(f"  {key}: {value}")
    print("Features:")
    for feature in seq_record.features:
        print(f"  {feature.type}: {feature.location}")

# 7. Modifying annotations (5.1.5)
print("\n7. Modifying annotations")
for seq_record in SeqIO.parse("hemoglobin.fasta", "fasta"):
    seq_record.annotations["molecule_type"] = "protein"
    seq_record.annotations["organism"] = "Homo sapiens"
    print(f"Modified record: {seq_record.id}")
    print(f"  Annotations: {seq_record.annotations}")

# 8. Converting file formats
print("\n8. Converting GenBank to FASTA")
count = SeqIO.convert("hemoglobin.gb", "genbank", "hemoglobin_converted.fasta", "fasta")
print(f"Converted {count} records")

# 9. Using SeqIO.to_dict()
print("\n9. Using SeqIO.to_dict()")
fasta_dict = SeqIO.to_dict(SeqIO.parse("hemoglobin.fasta", "fasta"))
print(f"Number of sequences in dictionary: {len(fasta_dict)}")

# 10. Using SeqIO.index()
print("\n10. Using SeqIO.index()")
gb_idx = SeqIO.index("hemoglobin.gb", "genbank")
print(f"Sequences in index: {len(gb_idx)}")
print(f"First sequence ID: {next(iter(gb_idx))}")

# 11. Writing sequences
print("\n11. Writing sequences")
modified_records = []
for record in SeqIO.parse("hemoglobin.fasta", "fasta"):
    record.seq = record.seq[10:30]  # Trim sequence
    record.id = "trimmed_" + record.id
    record.description = "Trimmed " + record.description
    modified_records.append(record)

SeqIO.write(modified_records, "hemoglobin_trimmed.fasta", "fasta")

# 12. Working with compressed files
print("\n12. Working with compressed files")
with gzip.open("hemoglobin.gb.gz", "wt") as handle:
    SeqIO.write(first_record, handle, "genbank")

with gzip.open("hemoglobin.gb.gz", "rt") as handle:
    compressed_record = SeqIO.read(handle, "genbank")

print(f"Compressed record ID: {compressed_record.id}")

# 13. Using low-level parsers
print("\n13. Using low-level parsers")
from Bio.SeqIO.FastaIO import SimpleFastaParser

total_len = 0
with open("hemoglobin.fasta") as handle:
    for title, seq in SimpleFastaParser(handle):
        total_len += len(seq)

print(f"Total sequence length: {total_len}")

# 14. Sequence manipulation and writing to StringIO
print("\n14. Sequence manipulation and writing to StringIO")
output_handle = StringIO()
for record in SeqIO.parse("hemoglobin.fasta", "fasta"):
    # Create a hypothetical DNA sequence from the protein sequence
    # This is for demonstration and isn't biologically accurate
    dna_seq = Seq("".join(['ATG' if aa == 'M' else 'GCT' for aa in record.seq]))
    rev_comp_record = SeqRecord(dna_seq.reverse_complement(), 
                                id="rc_"+record.id, 
                                description="reverse complement of hypothetical DNA")
    SeqIO.write(rev_comp_record, output_handle, "fasta")

print(output_handle.getvalue()[:200] + "...")  # Print first 200 characters

print("Script completed successfully!")
