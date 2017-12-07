import HTSeq
import getopt
import sys


def usage(x=''):
    usestr = """ make_mock_gff.py -o out.gff -f input.fa.fai -b input.bam
#########
average coverage across genome in steps of size

requires:
-o|--out=           output gff file
-f|--fai=           input fa.fai file (produced by e.g. samtools faidx genome.fa
-b|--bam=           input bam file

optional:
-s|--size=          desired bp per genome 'piece' (default = 32)
-h|--help           prints this
    """
    if len(x) > 0:
        usestr = x + '\n#########\n' + usestr

    print(usestr)
    sys.exit(1)


def read_fai(fai):
    out = {}
    with open(fai) as f:
        for line in f:
            line = line.rstrip()
            molecule, length, start, chars, with_endline = line.split()
            out[molecule] = int(length)
    return out


def coverage_on_chr(read_group, i, inext):
    bp = 0
    size = inext - i
    for read in read_group:
        len_inside = min(inext, read.iv.end) - max(i, read.iv.start)
        bp += len_inside
    return bp / size


def group_reads(bam, fai, size=32):

    bamgen = gen_bam(bam)
    latest = next(bamgen)
    read_group = []
    more_remaining = True

    while more_remaining:
        chromo = latest.iv.chrom
        end_at = fai[chromo]
        for i in range(0, end_at, size):
            i_next = i + size
            while (latest.iv.start <= i_next) & (latest.iv.chrom == chromo):
                read_group.append(latest)
                try:
                    latest = next(bamgen)
                except StopIteration:
                    more_remaining = False
            read_group = [x for x in read_group if x.iv.end >= i]  # drop alignments that no longer overlap
            yield read_group, i, i_next, chromo


def gen_bam(bam):
    for aln in bam:
        yield aln


def main():
    size = 32
    fai = None
    bam = None
    fileout = None
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "o:b:f:s:h",
                                       ["out=", "bam=", "fai=", "size=", "help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for o, a in opts:
        if o in ("-o", "--out"):
            fileout = a
        elif o in ("-b", "--bam"):
            bam = a
        elif o in ("-f", "--fai"):
            fai = a
        elif o in ("-s", "--size"):
            size = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if fai is None or fileout is None or bam is None:
        usage("required input missing")

    seq_lengths = read_fai(fai)
    alignment_file = HTSeq.BAM_Reader(bam)

    for key in seq_lengths:
        print('{}: {}'.format(key, seq_lengths[key]))

    j = 0  # temporary for faster testing

    print('chr', 'count', 'coverage')
    for read_group, i, i_next, chrom in group_reads(alignment_file, seq_lengths, size=size):
        print('{}\t{}\t{}'.format(chrom, len(read_group), coverage_on_chr(read_group, i, i_next)))
        if j >= 100000:
            break
        j += 1


if __name__ == "__main__":
    main()