from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)

    dp = [[0.0 for j in range(m + 1)] for i in range(n + 1)]
    traceback = [[0 for j in range(m + 1)] for i in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + scoring_function(seq1[i - 1], "-")
        traceback[i][0] = 1

    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + scoring_function("-", seq2[j - 1])
        traceback[0][j] = 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diagonal = (
                dp[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            up = (
                dp[i - 1][j]
                + scoring_function(seq1[i - 1], "-")
            )

            left = (
                dp[i][j - 1]
                + scoring_function("-", seq2[j - 1])
            )

            scores = [diagonal, up, left]

            best_direction = max(range(3), key=lambda k: scores[k])

            dp[i][j] = scores[best_direction]
            traceback[i][j] = best_direction

    aligned1 = []
    aligned2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        direction = traceback[i][j]

        if direction == 0:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == 1:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1

        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1 = "".join(reversed(aligned1))
    aligned2 = "".join(reversed(aligned2))

    return aligned1, aligned2, float(dp[n][m])

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    n = len(seq1)
    m = len(seq2)

    dp = [[0.0 for j in range(m + 1)] for i in range(n + 1)]
    traceback = [[0 for j in range(m + 1)] for i in range(n + 1)]

    maxx = 0.0
    max_i = 0
    max_j = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diagonal = (
                dp[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            up = (
                dp[i - 1][j]
                + scoring_function(seq1[i - 1], "-")
            )

            left = (
                dp[i][j - 1]
                + scoring_function("-", seq2[j - 1])
            )

            scores = [diagonal, up, left, 0.0]

            best_direction = max(range(4), key=lambda k: scores[k])

            dp[i][j] = scores[best_direction]
            traceback[i][j] = best_direction

            if dp[i][j] > maxx:
                maxx = dp[i][j]
                max_i = i
                max_j = j

    aligned1 = []
    aligned2 = []

    i = max_i
    j = max_j

    while i > 0 and j > 0 and dp[i][j] > 0.0:

        direction = traceback[i][j]

        if direction == 0:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == 1:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1

        elif direction == 2:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

        else:
            break

    aligned1 = "".join(reversed(aligned1))
    aligned2 = "".join(reversed(aligned2))

    return aligned1, aligned2, float(maxx)

def scoring_function(a, b):
    """Calculate alignment score for 2 sequences, using BLOSUM62.
    
        Parameters
        ----------
        a, b: str
            aligned sequences for scoring.
    
        Returns
        -------
        n/a: float
            alignment score.
    
        Examples
        --------
        >>> scoring_function("A", "A")
        4.0
    """
    if a == "-" or b == "-":
        return -5.0

    return float(blosum62[a, b])

def get_sequence_from_accession(accession_code):
    """Get a sequence from an accession code.

    Parameters
    ----------
    accession_code: str
        Accession code of the sequence to be retrieved.

    Returns
    -------
    n/a: str
        Sequence corresponding to the accession code.

    Examples
    --------
    >>> get_sequence_from_accession("NC_045512.2")
    'ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA...'

    """
    from Bio import Entrez, SeqIO

    Entrez.email = "z.brookes@student.unsw.edu.au"

    with Entrez.efetch(
        db="nucleotide",
        id=accession_code,
        rettype="fasta",
        retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "fasta")
        
    return str(record.seq)[21562:25384]
    
