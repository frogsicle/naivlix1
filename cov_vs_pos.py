import HTSeq
import getopt
import sys
from genome2simplesttrain import atcg2numbers, fasta2seqs
import random


class Genome:
    def __init__(self, bam, fasta):
        self.bam = bam
        self.fasta = fasta
        self.chromosomes = {}

    def fill_chromosomes(self, size):
        seq_lengths = {}
        for key in self.fasta:
            seq_lengths[key] = len(self.fasta[key])

        group_gen = group_reads(self.bam, seq_lengths, size=size)
        read_group, i, i_next, old_chrom = next(group_gen)
        coverage = [coverage_on_chr(read_group, i, i_next)]
        for read_group, i, i_next, chrom in group_gen:
            if chrom != old_chrom:
                self.chromosomes[old_chrom] = Chromosome(coverage, self.fasta[old_chrom], bin_size=size)
                coverage = []
            coverage.append(coverage_on_chr(read_group, i, i_next))
            old_chrom = chrom
        self.chromosomes[old_chrom] = Chromosome(coverage, self.fasta[old_chrom], bin_size=size)

    def slice_xy_pairs(self, n_pieces):
        for chromo in self.chromosomes:
            for seq, cov in self.chromosomes[chromo].slice(n_pieces):
                yield seq, cov


class Chromosome:
    def __init__(self, coverage, sequence, bin_size):
        self.coverage = coverage
        self.sequence = sequence
        self.bin_size = bin_size
        # I doubt ends are actually expressed, but padding just in case
        padding = "N" * (len(self.sequence) % bin_size)
        self.sequence += padding

    def slice(self, n_coverage_bins):
        n_bp = n_coverage_bins * self.bin_size

        for i in range(0, len(self.sequence), n_bp):
            i_cov = i // self.bin_size
            seq = self.sequence[i:(i + n_bp)]
            cov = self.coverage[i_cov:(i_cov + n_coverage_bins)]
            if len(cov) < n_coverage_bins:  # AKA, add padding if we ran off end of sequence
                cov += [0] * (n_coverage_bins - len(cov))
                seq += "N" * (n_bp - len(seq))
            yield seq, cov


class Exporter:
    sets = ["train", "val", "test"]

    def __init__(self):
        self.files = {}
        self.formatter = onewarm_formatter
        self.train_prob = 0.6
        self.val_prob = 0.2

    def which_set(self):
        r = random.uniform(0, 1)
        if r <= self.train_prob:
            out = 'train'
        elif r <= (self.train_prob + self.val_prob):
            out = 'val'
        else:
            out = 'test'
        return out

    def export(self, x, y):
        pass

    def open(self, prefix):
        pass

    def close(self):
        for key in self.files:
            self.files[key].close()


class TextExporter(Exporter):
    def open(self, prefix):
        for key in Exporter.sets:
            self.files[key] = open('{}.{}.csv'.format(prefix, key), 'w')

    def export(self, x, y):
        key = self.which_set()
        x = self.formatter(x)
        lineout = to_text_line(x, y)
        self.files[key].write(lineout)


def to_text_line(x, y):
    if isinstance(y, list):
        y = ','.join([str(w) for w in y])
    elif isinstance(y, int) or isinstance(y, float):
        y = str(y)
    if isinstance(x, list):
        x = ','.join([str(w) for w in x])
    if isinstance(x, str) and isinstance(y, str):
        return '{},{}\n'.format(x, y)
    else:
        raise NotImplementedError


def onewarm_formatter(x):
    numseq = atcg2numbers(x)
    return numseq


def onehot_formatter(x):
    numseq = atcg2numbers(x, False)
    return numseq


def raw_formatter(x):
    return x


def usage(x=''):
    usestr = """ cov_vs_pos -o prefix_out -f genome.fa -b input.sorted.bam
#########
average coverage across genome in steps of size

requires:
-o|--out=           output prefix
-f|--fasta=         input genome fasta file
-b|--bam=           input bam file

optional:
-s|--size=          desired bp per genome 'piece' (default = 32)
-n|--pieces=        desired number of pieces per training unit (default = 128)
-h|--help           prints this
    """
    if len(x) > 0:
        usestr = x + '\n#########\n' + usestr

    print(usestr)
    sys.exit(1)


def read_fai(fai):
    """from .fai file to dict[sequence_ID] = length"""
    out = {}
    with open(fai) as f:
        for line in f:
            line = line.rstrip()
            molecule, length, start, chars, with_endline = line.split()
            out[molecule] = int(length)
    return out


def coverage_on_chr(read_group, i, inext):
    """calculate coverage between i and inext of alignments in read_group"""
    bp = 0
    size = inext - i
    for read in read_group:
        len_inside = min(inext, read.iv.end) - max(i, read.iv.start)
        bp += len_inside
    return bp / size


def group_reads(bam, fai, size=32):
    """step through genome and yield all reads overlapping each step"""
    bamgen = gen_bam(bam)
    latest = next(bamgen)
    while latest.iv is None: # should skip passed all unaligned if they are sorted to the start
        latest = next(bamgen)
    read_group = []
    more_remaining = True

    while more_remaining:
        chromo = latest.iv.chrom
        end_at = fai[chromo]
        for i in range(0, end_at, size):  # should pad to next length divisible by size, I think...
            i_next = i + size
            if latest.iv is not None:  # stop looking for more reads if the last one in the file was not aligned
                while (latest.iv.start <= i_next) & (latest.iv.chrom == chromo):
                    read_group.append(latest)
                    try:
                        latest = next(bamgen)
                        while latest.iv is None:  # ignore non-aligned reads?
                            latest = next(bamgen)
                    except StopIteration:
                        more_remaining = False
                        if latest.iv is None:
                            print('Info, last non-none alignment')
                            print('chr: {}, i: {}, rg [-1]: {}'.format(chromo, i, read_group[-1]))
                        break
            read_group = [x for x in read_group if x.iv.end >= i]  # drop alignments that no longer overlap
            yield read_group, i, i_next, chromo


def gen_bam(bam):
    for aln in bam:
        yield aln


def main():
    size = 32
    fastafile = None
    bam = None
    fileout = None
    n_pieces = 128

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "o:b:f:s:n:h",
                                       ["out=", "bam=", "fasta=", "size=", "pieces=", "help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for o, a in opts:
        if o in ("-o", "--out"):
            fileout = a
        elif o in ("-b", "--bam"):
            bam = a
        elif o in ("-f", "--fai"):
            fastafile = a
        elif o in ("-s", "--size"):
            size = int(a)
        elif o in ("-n", "--pieces"):
            n_pieces = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if fastafile is None or fileout is None or bam is None:
        usage("required input missing")

    exporter = TextExporter()
    exporter.open(fileout)

    fasta = fasta2seqs(fastafile)
    alignment_file = HTSeq.BAM_Reader(bam)

    genome = Genome(alignment_file, fasta)
    genome.fill_chromosomes(size)

    for x, y in genome.slice_xy_pairs(n_pieces):
        exporter.export(x, y)

    exporter.close()

if __name__ == "__main__":
    main()