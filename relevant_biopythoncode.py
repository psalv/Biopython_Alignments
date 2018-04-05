

# Setup: Add biopython to the interpreter packages
# Biopython docs: http://biopython.org/DIST/docs/tutorial/Tutorial.html


### SEQUENCE ALIGNMENT ###
#
# To show that channelrhodopsin and rhodopsin do not align very well because they are not very closely related
#

def sequence_alignment():

    from Bio import pairwise2
    from Bio import SeqIO
    from Bio.SubsMat.MatrixInfo import blosum62


    ### Nucleotide alignment
    seq1 = SeqIO.read("JN223416.1.faa", "fasta")
    seq2 = SeqIO.read("NC_000003.12.faa", "fasta")

    # Aligning seq1 with itself to show what an optimal alignment would look like:
    alignments = pairwise2.align.globalds(seq1.seq, seq1.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))

    # Performing a sequence alignment of two sequences:
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    print(pairwise2.format_alignment(*alignments[0]))

    # Altering the penalties (gap open penalty of 10 and a gap extension penalty of 0.5), to give better alignment:
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))


    ### Protein alignment
    seq1 = SeqIO.read("crho_prot.fasta", "fasta")
    seq2 = SeqIO.read("rho_prot.fasta", "fasta")

    # Aligning sequence with itself to show perfect alignment
    alignments = pairwise2.align.globalds(seq1.seq, seq1.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))

    # Aligning channelrhodopsin and rhodopsin
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))


    ### Nucleotide alignment with newly found sequences (from BLAST)
    seq1 = SeqIO.read("found_from_blast.faa", "fasta")
    seq2 = SeqIO.read("NC_000003.12.faa", "fasta")

    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))


### BLAST ###
#
# To find another sequence that is closely related to rhodopsin
#


def prot_blast():

    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import SeqIO

    seq2 = SeqIO.read("rho_prot.fasta", "fasta")

    # Perform a blast of sequence 1 and write it to xml, then read it again to be able to manipulate it
    result_handle = NCBIWWW.qblast("blastp", "nr", seq2.format("fasta"))
    with open("my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()

    E_VALUE_THRESH = 0.04
    result_handle = open("my_blast.xml")
    for blast_record in NCBIXML.parse(result_handle):
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    print('\n****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')


def nt_blast():

    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import SeqIO


    ### Nucleotide blast
    seq2 = SeqIO.read("NC_000003.12.faa", "fasta")

    # Perform a blast of sequence 1 and write it to xml, then read it again to be able to manipulate it
    result_handle = NCBIWWW.qblast("blastn", "nt", seq2.format("fasta"))
    with open("my_blast.xml", "w") as out_handle:
         out_handle.write(result_handle.read())

    result_handle.close()

    E_VALUE_THRESH = 0.04
    result_handle = open("my_blast.xml")
    for blast_record in NCBIXML.parse(result_handle):
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    print('\n****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('e value:', hsp.expect)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')


def transcription_translation():

    from Bio import pairwise2
    from Bio import SeqIO
    from Bio.SubsMat.MatrixInfo import blosum62



    seq1 = SeqIO.read("found_from_blast.faa", "fasta")
    seq2 = SeqIO.read("NC_000003.12.faa", "fasta")

    coding_dna1 = seq1.seq.transcribe()
    coding_dna2 = seq2.seq.transcribe()

    alignments = pairwise2.align.globalxx(coding_dna1, coding_dna2)
    print(pairwise2.format_alignment(*alignments[0]))

    protein_seq1 = coding_dna1.translate()
    protein_seq2 = coding_dna2.translate()

    alignments = pairwise2.align.globalxx(protein_seq1, protein_seq2)
    print(pairwise2.format_alignment(*alignments[0]))


transcription_translation()
