from Bio import pairwise2

def get_novel_sequences(sequence_dict):
    """
    The function performs a pairwise comparison between protein chains in order to select the ones that have unique
    sequence. A chain sequence is considered to be unique if the identity score (percentage of identical aa
    in the alignment) between it self and any other sequence is below 95%.

    :param sequence_dict: dictionary containing the information for each protein chain and their sequences
    :return: list of all unique/novel chain IDs
    """
    novel_sequences = []
    similar_sequences = []
    # TODO - consider case in which we have DNA and RNA
    #sequence_keys = [key for key in sequence_dict.keys() if key != "DNA" or key != "RNA"]
    sequence_keys = [key for key in sequence_dict.keys()]
    for current, first_key in enumerate(sequence_keys):
        for second_key in sequence_keys[:current]:
            if sequence_similarity(first_key, second_key, sequence_dict):
                if second_key not in similar_sequences:
                    similar_sequences.append(second_key)
                if first_key not in novel_sequences:
                    novel_sequences.append(first_key)
                else:
                    novel_sequences.remove(first_key)
            else:
                if first_key not in novel_sequences and first_key not in similar_sequences:
                    novel_sequences.append(first_key)
                if second_key not in novel_sequences and second_key not in similar_sequences:
                    novel_sequences.append(second_key)
    return sorted(novel_sequences)

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
            return sim_score > 0.95

