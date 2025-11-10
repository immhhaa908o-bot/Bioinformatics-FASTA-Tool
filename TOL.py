from STR2323 import *
from collections import Counter
import random
import re
from pathlib import Path


def _infer_seq_type(seq: str, default: str = "DNA") -> str:
    """Infer sequence type (DNA/RNA) from characters.
    - Contains 'U' and not 'T' -> RNA
    - Contains 'T' and not 'U' -> DNA
    - Contains neither -> default
    - Contains both -> error
    """
    s = set(seq)
    has_u = 'U' in s
    has_t = 'T' in s
    if has_u and not has_t:
        return "RNA"
    if has_t and not has_u:
        return "DNA"
    if not has_u and not has_t:
        return default
    raise ValueError("Sequence appears to contain both T and U; cannot infer type.")


class bio_seq:
    """DNA sequnece class. Defalt value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation.
        seq_type may be "DNA", "RNA", or "AUTO"/None to infer from sequence.
        """
        # Remove whitespace/newlines and normalize case to support multi-line inputs
        self.seq = re.sub(r"\s+", "", str(seq)).upper()
        self.label = label
        # Allow automatic type detection
        if seq_type is None or str(seq_type).upper() == "AUTO":
            self.seq_type = _infer_seq_type(self.seq, default="DNA")
        else:
            self.seq_type = str(seq_type).upper()
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    # DNA Toolkit functions:

    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def nucleotide_frequency(self):
        """Count nucleotides in a given sequence. Return a dictionary"""
        return dict(Counter(self.seq))

    def transcription(self):
        """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    def reverse_complement(self):
        """
        Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC Content in a DNA/RNA sequence"""
        # Handle empty sequences to avoid division by zero
        if not self.seq:
            return 0
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)

    def gc_content_subsec(self, k=20):
        """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res


def from_fasta(path, seq_type="AUTO") -> bio_seq:
    """Read first FASTA record and return a bio_seq.
    If seq_type is AUTO/None, infer DNA/RNA from sequence.
    """
    label, sequence = read_fasta(path)
    return bio_seq(sequence, seq_type, label)


# ------------------------------
# FASTA/Text IO utilities
# ------------------------------
def read_fasta(path):
    """
    Read the first record from a FASTA file and return (label, sequence).
    - Ignores blank lines and trims whitespace.
    - Concatenates multiline sequences.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    label = None
    seq_chunks = []
    with p.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if label is None:
                    label = line[1:].strip()
                else:
                    # Stop at next header (only first record)
                    break
            else:
                # keep only non-whitespace chars
                seq_chunks.append(re.sub(r"\s+", "", line))

    sequence = ''.join(seq_chunks).upper()
    if not sequence:
        raise ValueError(f"No sequence data found in FASTA: {path}")
    return (label or 'No Label'), sequence


def write_fasta(path, label, sequence, width=60):
    """Write a FASTA file with the given label and sequence, wrapping lines to width."""
    p = Path(path)
    seq = re.sub(r"\s+", "", str(sequence)).upper()
    with p.open('w', encoding='utf-8') as f:
        f.write(f">{label}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")


def read_text_file(path):
    """Read plain text file contents (UTF-8)."""
    p = Path(path)
    return p.read_text(encoding='utf-8')


def write_text_file(path, text):
    """Write plain text to file (UTF-8)."""
    p = Path(path)
    p.write_text(str(text), encoding='utf-8')


def read_fasta_all(path):
    """
    Read all records from a (multi-)FASTA file and return a list of (label, sequence).
    - Ignores blank lines and trims whitespace.
    - Concatenates multiline sequences per record.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records = []
    label = None
    seq_chunks = []
    with p.open('r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if label is not None:
                    sequence = ''.join(seq_chunks).upper()
                    if sequence:
                        records.append((label or 'No Label', sequence))
                    seq_chunks = []
                label = line[1:].strip()
            else:
                seq_chunks.append(re.sub(r"\s+", "", line))
    # flush last record
    if label is not None or seq_chunks:
        sequence = ''.join(seq_chunks).upper()
        if sequence:
            records.append((label or 'No Label', sequence))

    if not records:
        raise ValueError(f"No sequence data found in FASTA: {path}")
    return records
