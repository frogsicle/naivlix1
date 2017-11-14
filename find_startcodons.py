from genome2simplesttrain import atcg2numbers, fasta2seqs, cds_from_gff
import sys
import getopt
import random
import copy


def y_formater(numseq, cds_start_or_not):
    lineout = ','.join([str(x) for x in numseq + [cds_start_or_not]]) + '\n'
    return lineout


def start_range(cds, fasta, to_start=-98, to_end=98):
    if cds['strand'] == '+':
        start = cds['begin'] + to_start
        end = cds['begin'] + to_end
        out = fasta[cds['chr']][start:end]
    elif cds['strand'] == '-':
        start = cds['end'] - to_start
        end = cds['end'] - to_end
        out = fasta[cds['chr']][end:start]
        out = rc(out)
    else:
        raise ValueError('strand character not recognized (only + or - allowed)')
    return out


def reverse(x):
    x = list(x)
    x.reverse()
    x = ''.join(x)
    return x


def complement(x):
    repkey = {}
    pairings = (('a', 't'),
                ('g', 'c'),
                ('m', 'k'),
                ('r', 'y'),
                ('w', 'w'),
                ('s', 's'),
                ('v', 'b'),
                ('h', 'd'),
                ('n', 'n'))

    for a, b in pairings:
        A = a.upper()
        B = b.upper()
        repkey[a] = b
        repkey[b] = a
        repkey[A] = B
        repkey[B] = A

    out = []
    for char in x:
        out += repkey[char]
    return ''.join(out)


def rc(x):
    x = complement(x)
    x = reverse(x)
    return x


def notstart_range(cds_old, cds, fasta, to_start=-98, to_end=98):
    """get random subsequence between two genes"""
    ok_strand = ['+', '-']
    if cds_old['strand'] not in ok_strand or cds['strand'] not in ok_strand:
        raise ValueError('strand character not recognized (only + or - allowed)')
    # prevent large overlaps
    if cds_old['strand'] == '+':
        min_start = cds_old['begin'] + to_end
    else:
        min_start = cds_old['end'] + to_end
    if cds['strand'] == '+':
        max_start = cds['begin'] + to_start
    else:
        max_start = cds['end'] + to_start
    start = random.uniform(min_start, max_start)
    # take random point between old and new CDS, then treat it as start
    # for next non-gene range
    start = int(start)

    mockcds = copy.deepcopy(cds)
    mockcds['begin'] = start
    mockcds['end'] = start  # because start range might take either, depending on orientation of cds

    out = start_range(mockcds, fasta, to_start, to_end)
    return out


def seq2line(subseq, y, expected_length):
    if len(subseq) != expected_length:
        print('subsequence: {}'.format(subseq))
        print('y: {}'.format(y))
        raise ValueError('length of subsequence does not match expected: {} bp'.format(expected_length))
    numseq = atcg2numbers(subseq)
    lineout = y_formater(numseq, y)
    return lineout


def process_cds2line(cds, genome, to_start, to_end, size, oldcds=None):
    if oldcds is None:
        subseq = start_range(cds, genome, to_start, to_end)
        y = 1
    else:
        subseq = notstart_range(oldcds, cds, genome, to_start, to_end)
        y = 0
    try:
        lineout = seq2line(subseq, y, size)
    except ValueError:
        lineout = ''
    return lineout


def range_from_size(size=196):
    to_end = size // 2
    to_start = to_end * -1
    to_end += size % 2  # any remainder for odd numbers gets tacked on to end
    return to_start, to_end


def usage(x=''):
    usestr = """ find_startcodons.py -f genome.fa -g genes.gff3 -s piece_size
#########
returns genome slices (of piece_size) as vector of [0s and 1s] and a score for % coding

requires:
-f|--fasta=     genome file in fasta format
-g|--gff=       genes / features in gff3 format (including CDS), [skipping -g assumes no CDS]
-o|--out=       output file prefix (produces prefix.x.csv and prefix.y.csv)
optional:
-s|--size=      desired bp per genome 'piece' (default = 196)
-h|--help       prints this
    """
    if len(x) > 0:
        usestr = x + '\n#########\n' + usestr

    print(usestr)
    sys.exit(1)


def main():
    # getopt, if file > return sequences from file
    fasta = None
    gff = None
    fileout = None
    size = 196
    train_prob = 0.6
    val_prob = 0.2

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:g:s:o:h",
                                       ["fasta=", "gff=", "size=", "help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-g", "--gff"):
            gff = a
        elif o in ("-s", "--size"):
            size = int(a)
        elif o in ("-o", "--out"):
            fileout = a
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if fasta is None or fileout is None or gff is None:
        usage("required input missing")

    # get
    to_start, to_end = range_from_size(size)

    # output files
    train_out = open(fileout + '.train.csv', 'w')
    val_out = open(fileout + '.val.csv', 'w')
    test_out = open(fileout + '.test.csv', 'w')

    # import fasta sequences
    genome = fasta2seqs(fasta)

    cdsgen = cds_from_gff(gff)

    # handle first element seperately, so the rest can be done in start, not-start pairs
    oldcds = next(cdsgen)
    lineout = process_cds2line(oldcds, genome, to_start, to_end, size=size)
    train_out.write(lineout)

    for cds in cdsgen:
        ns_lineout = ''
        s_lineout = '' # so that nothing is written / happens when subseqs can't be found (e.g. between chromosomes)
        # the space in between the previous and current cds, not a start
        if oldcds['chr'] == cds['chr']:
            ns_lineout = process_cds2line(cds, genome, to_start, to_end, size=size, oldcds=oldcds)

        # the current is a start
        s_lineout = process_cds2line(cds, genome, to_start, to_end, size=size)

        linesout = ns_lineout + s_lineout

        # chose a file, and export to it
        a_number = random.uniform(0, 1)
        if a_number <= train_prob:
            train_out.write(linesout)
        elif a_number <= (train_prob + val_prob):
            val_out.write(linesout)
        else:
            test_out.write(linesout)

        oldcds = cds

    train_out.close()
    test_out.close()
    val_out.close()


if __name__ == "__main__":
    main()
