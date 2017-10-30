import sys
import getopt
from genome2simplesttrain import Yformatr, atcg2numbers
import random
import numpy as np

def usage(x=''):
    usestr = """ make_dummy_data.py -o example_prefix
#########
returns genome slices (of piece_size) as vector of [0s and 1s] and a to-predict type

requires:
-o|--out=           output file prefix (produces prefix.(test|train|val).csv)
optional:
-e|--encode_amb=    encoding for ambiguous characters ('even' for evenly spread or '0'; default 'even')
-t|--type=          type of score {'%': percent coding, 'b': boolean "is coding" by pos} (default: 'b')
-m|--motif=         target motif (default ATG; may consist of only ATCGYRWSKMDVHBN characters)
-f|--frac_targ=     fraction of target motif (default 0.5)
-s|--size=          desired bp per genome 'piece' (default = 196)
-n|--num_total=     number of total X,Y pairs to produce (default = 10000)
-h|--help           prints this
    """
    if len(x) > 0:
        usestr = x + '\n#########\n' + usestr

    print(usestr)
    sys.exit(1)


def not_motif(len_motif):
    overweight_known = 'ATCG' * 500
    overweight_n = 'N' * 50
    choices = 'YRWSKMDVHB' + overweight_known + overweight_n
    out = ''
    for i in range(len_motif):
        out += random.choice(choices)
    return out


def maybe_toggle(x, check_flip_prob, target_fraction_true):
    if random.uniform(0, 1) < check_flip_prob:
        x = random.uniform(0, 1) <= target_fraction_true
    return x


def sub_mock2train(motif, spread_evenly, frac_targ, size, check_flip_prob=0.01):
    lm = len(motif)
    motif_on = random.uniform(0, 1) <= frac_targ
    motif_or_not = []
    out = ''
    for i in range(0, size, lm):
        if motif_on:
            out += motif
        else:
            out += not_motif(lm)
        motif_or_not += [motif_on] * lm
        motif_on = maybe_toggle(motif_on, check_flip_prob, frac_targ)
    out = out[:size]
    motif_or_not = motif_or_not[:size]
    out = atcg2numbers(out, spread_evenly)
    motif_or_not = [int(x) for x in motif_or_not]
    percent_motif = sum(motif_or_not) / len(motif_or_not)
    return out, percent_motif, motif_or_not


def mock2train(n, motif, spread_evenly, frac_targ, size, check_flip_prob=0.01):
    for i in range(n):
        yield sub_mock2train(motif, spread_evenly, frac_targ, size, check_flip_prob)


def main():
    fileout = None
    size = 196
    train_prob = 0.6
    val_prob = 0.2
    frac_targ = 0.5
    spread_evenly = True
    allowed_types = Yformatr.allowed_types
    out_type = 'b'
    motif = 'ATG'
    iterations = 10000
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:m:s:n:o:e:t:h",
                                       ["frac_targ=", "motif=", "size=", "num_total=", "out=",  "type=", "encode_amb=",
                                        "help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--frac_targ"):
            frac_targ = float(a)
        elif o in ("-m", "--motif"):
            motif = a
        elif o in ("-s", "--size"):
            size = int(a)
        elif o in ("-o", "--out"):
            fileout = a
        elif o in ('-t', '--type'):
            if a in allowed_types:
                out_type = a
            else:
                usage("invalid entry for -t/--type. Valid options are:{0}".format(allowed_types))
        elif o in ("-e", "--encode_amb"):
            if a == 'even':
                spread_evenly = True
            elif a == '0':
                spread_evenly = False
            else:
                usage("-e must be either 'even' or '0'")
        elif o in ("-n", "--num_total"):
            iterations = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if fileout is None:
        usage("required output file missing")

    test_prob = 1. - (train_prob + val_prob)
    if test_prob <= 0:
        usage("training and validation set probabilities set too high to have a test set")

    yformatr = Yformatr()
    yformatr.out_type = out_type

    train_out = open(fileout + '.train.csv', 'w')
    val_out = open(fileout + '.val.csv', 'w')
    test_out = open(fileout + '.test.csv', 'w')

    for numseq, percent_coding, cds_or_not in mock2train(iterations, motif, spread_evenly, frac_targ, size):
        lineout = yformatr.format(numseq, percent_coding, cds_or_not)
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
