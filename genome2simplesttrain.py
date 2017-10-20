import intervaltree
import sys
import getopt
import random


def fasta2seqs(fastafile, headerchar=' '):
    seqs = {}
    running_seq = ''
    running_id = None
    with open(fastafile) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if len(running_seq) > 0:
                    seqs[running_id] = running_seq
                    running_seq = ''
                running_id = line.split(headerchar)[0][1:]
            else:
                running_seq += line
    # save last sequence too
    seqs[running_id] = running_seq
    return seqs


def gff2trees(gfffile, feature='CDS'):
    trees = {}
    i = 0
    for cds in cds_from_gff(gfffile, feature=feature):
        chr = cds['chr']
        if chr not in trees:
            trees[chr] = intervaltree.IntervalTree()
        # merge with any existing interval
        olaps = trees[chr][cds['begin']:cds['end']]
        olaps_begins = [x[0] for x in olaps]
        olaps_ends = [x[1] for x in olaps]
        newmin = min([cds['begin']] + olaps_begins)
        newmax = max([cds['end']] + olaps_ends)
        # clean up and add
        trees[chr].chop(newmin, newmax)
        cdsint = intervaltree.Interval(newmin, newmax, i)
        trees[chr].add(cdsint)
        i += 1
    return trees


def cds_from_gff(gfffile, feature='CDS'):
    with open(gfffile) as f:
        for line in f:
            if not line.startswith('#'):
                sline = line.split('\t')
                if sline[2] == feature:
                    out = {'begin': int(sline[3]), 'end': int(sline[4]) + 1, 'chr': sline[0]}
                    yield out


def tree2train(tree, size, seq=None):
    size = int(size)
    starts = range(0, tree.end() - size, size)  # -size so no handling of end-cases is required
    for i in starts:
        end = i + size
        minitree = tree.search(i, end)
        minitree = intervaltree.IntervalTree(minitree)
        cds_or_not = [int(minitree.overlaps(x)) for x in range(i, end)]
        percent_coding = sum(cds_or_not) / len(cds_or_not)

        subseq = seq[i:end]
        numseq = atcg2numbers(subseq)
        yield((numseq, percent_coding))


def atcg2numbers(seq):
    decode = {'a': [1, 0, 0, 0],
              't': [0, 1, 0, 0],
              'c': [0, 0, 1, 0],
              'g': [0, 0, 0, 1],
              'y': [0, 0.5, 0.5, 0],
              'r': [0.5, 0, 0, 0.5],
              'w': [0.5, 0.5, 0, 0],
              's': [0, 0, 0.5, 0.5],
              'k': [0, 0.5, 0, 0.5],
              'm': [0.5, 0, 0.5, 0],
              'd': [0.33, 0.33, 0, 0.33],
              'v': [0.33, 0, 0.33, 0.33],
              'h': [0.33, 0.33, 0.33, 0],
              'b': [0, 0.33, 0.33, 0.33],
              'n': [0.25, 0.25, 0.25, 0.25]}
    onfail = [0, 0, 0, 0]
    seq = seq.lower()
    out = []
    for letter in seq:
        try:
            out += decode[letter]
        except KeyError:
            print('non DNA character: {0}, adding {1} for no match'.format(letter, onfail))
            out += onfail
    return out


def usage(x=''):
    usestr = """ genome2simplesttrain.py -f genome.fa -g genes.gff3 -s piece_size
#########
returns genome slices (of piece_size) as vector of [0s and 1s] and a score for % coding

requires:
-f|--fasta=     genome file in fasta format
-g|--gff=       genes / features in gff3 format (including CDS)
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

    if fasta is None or gff is None or fileout is None:
        usage("required input missing")

    test_prob = 1. - (train_prob + val_prob)
    if test_prob <= 0:
        usage("training and validation set probabilities set too high to have a test set")

    train_out = open(fileout + '.train.csv', 'w')
    val_out = open(fileout + '.val.csv', 'w')
    test_out = open(fileout + '.test.csv', 'w')

    # import fasta sequences
    genome = fasta2seqs(fasta)
    # gff -> interval tree
    trees = gff2trees(gff)
    for treek in trees:
        tree = trees[treek]
        for numseq, percent_coding in tree2train(tree, size, genome[treek]):
            lineout = ','.join([str(x) for x in numseq + [round(percent_coding, 2)]]) + '\n'
            a_number = random.uniform(0, 1)
            if a_number <= train_prob:
                train_out.write(lineout)
            elif a_number <= (train_prob + val_prob):
                val_out.write(lineout)
            else:
                test_out.write(lineout)

    train_out.close()
    test_out.close()
    val_out.close()


if __name__ == "__main__":
    main()

