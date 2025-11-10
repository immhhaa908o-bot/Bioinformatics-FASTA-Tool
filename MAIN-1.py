import os
from TOL import bio_seq, read_fasta

# إذا وُجد ملف FASTA باسم input.fasta نقرأه تلقائياً، وإلا نستخدم التسلسل المضمَّن.
fasta_file = "input.fasta"

# Try current working directory first, then the script's directory
fasta_path = fasta_file
if not os.path.exists(fasta_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    alt_path = os.path.join(script_dir, fasta_file)
    if os.path.exists(alt_path):
        fasta_path = alt_path

if os.path.exists(fasta_path):
    label, seq = read_fasta(fasta_path)
    test_seq = seq
    test_label = label
else:
    test_label = "No Label"
    test_seq = """"""

# AUTO: يحدد النوع تلقائياً (DNA/RNA) بناءً على T/U
test_dna = bio_seq(test_seq, "AUTO", test_label)

print(test_dna.get_seq_info())
print(test_dna.nucleotide_frequency())
print(test_dna.transcription())
print(test_dna.reverse_complement())
print(test_dna.gc_content())
print(test_dna.gc_content_subsec())
print(test_dna.translate_seq())
print(test_dna.codon_usage('L'))

for rf in test_dna.gen_reading_frames():
    print(rf)

print(test_dna.all_proteins_from_orfs())
