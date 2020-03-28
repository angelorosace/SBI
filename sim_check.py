from Bio import pairwise2

def sequence_similarity (chain1, chain2, sequence_dict):
    """
    Takes two protein chain sequences and returns true if their sequences have > 95% similarity and false if not
    """
    # Get the sequences of the chains passed
    if sequence_dict:
        seq1 = sequence_dict[chain1]
        seq2 = sequence_dict[chain2]
     if seq1 and seq2:
        # Align the aa sequences to see if there is similiarity
        alignment = pairwise2.align.globalxx(seq1, seq2)
        score = alignment[0][2]
        length = max(len(seq1), len(seq2))
        sim_score = score / length

        if sim_score > 0.95:
            return True
        else:
            return False
