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

    def __init__(self, formatter):
        self.files = {}
        self.formatter = formatter
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

    def export(self, x, y, export_x=True):
        pass

    def open(self, prefix):
        pass

    def close(self):
        for key in self.files:
            self.files[key].close()


class TextExporter(Exporter):
    """exports x,y as csv"""
    def open(self, prefix):
        for key in Exporter.sets:
            self.files[key] = open('{}.{}.csv'.format(prefix, key), 'w')

    def export(self, x, y, export_x=True):
        key = self.which_set()
        if export_x:
            x = self.formatter(x)
            lineout = to_text_line(x, y)
        else:
            lineout = to_text_line(y)
        self.files[key].write(lineout)


class SplitTextExporter(Exporter):
    """exports x as csv, and y as csv"""
    def open(self, prefix):
        for key in Exporter.sets:
            for dat in ('x', 'y'):
                self.files[dat + key] = open('{}.{}_{}.csv'.format(prefix, dat, key), 'w')

    def export(self, x, y, export_x=True):
        key = self.which_set()
        if export_x:
            x = self.formatter(x)
            lineout_x = to_text_line(x)
            self.files['x' + key].write(lineout_x)
        lineout_y = to_text_line(y)
        self.files['y' + key].write(lineout_y)


def to_text_line(*args):
    to_join = []
    for sub_group in args:
        if isinstance(sub_group, str):
            as_string = sub_group
        elif isinstance(sub_group, list):
            as_string = ','.join([str(x) for x in sub_group])
        elif isinstance(sub_group, float) or isinstance(sub_group, int):
            as_string = str(sub_group)
        elif sub_group is None:
            as_string = ''
        else:
            raise NotImplementedError('while formatting: {}, of type: {}'.format(sub_group, type(sub_group)))
        to_join.append(as_string)
    return ','.join(to_join) + '\n'


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
--train_only        returns all x,y pairs in one file (train)
--skip_x            export only y
[the above two parameters are there for cases where you have multiple y per x, and want to join them later]
--seperate_x        x gets seperate file from y
-d|--deterministic= seed for deterministic "random" numbers (for later combination or simply a constant test set)
[another option to facilitate post processing]

-s|--size=          desired bp per genome 'piece' (default = 32)
-n|--pieces=        desired number of pieces per training unit (default = 128)
-F|--formatter=     one of [warm, hot, raw] (default = warm) to format x as follows:
                        warm ~ a one hot for ATCG, except with ambiguity codes respected
                        hot = a one hot vector for ATCG
                        raw = ATCG as is
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

    out = bp / size
    if out < 0:
        for read in read_group:
            print('chr {}, start {}, end {}'.format(read.iv.chrom, read.iv.start, read.iv.end))
        raise ValueError("coverage can't be {}, with i {}, inext {}, "
                         "bp {}, size {}".format(out, i, inext, bp, size))
    return bp / size


def group_reads(bam, fai, size=32):
    """step through genome and yield all reads overlapping each step"""
    bamgen = gen_bam(bam)
    latest = next(bamgen)
    counter = 0
    while latest.iv is None:  # should skip passed all unaligned if they are sorted to the start
        latest = next(bamgen)
        counter += 1
    print('{} unaligned reads skipped at start'.format(counter))
    more_remaining = True

    while more_remaining:
        read_group = []  # start every chromosome empty
        chromo = latest.iv.chrom
        print('processing chromosome: {}'.format(chromo))
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
            read_group = [x for x in read_group if x.iv.end > i]  # drop alignments that no longer overlap
            yield read_group, i, i_next, chromo


def gen_bam(bam):
    for aln in bam:
        yield aln

def choose_formatter(user_in):
    if user_in == 'warm':
        return onewarm_formatter
    elif user_in == 'hot':
        return onehot_formatter
    elif user_in == 'raw':
        return raw_formatter
    else:
        raise ValueError('only implemented formatters are: "warm", "hot", and "raw"')


def main():
    size = 32
    fastafile = None
    bam = None
    fileout = None
    n_pieces = 128
    train_only = False
    export_x = True
    seperate_x = False
    formatter = onewarm_formatter

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "o:b:f:s:n:d:F:h",
                                       ["out=", "bam=", "fasta=", "size=", "pieces=", "deterministic=", "help",
                                        "train_only", "skip_x", "seperate_x", "formatter="])
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
        elif o == "--train_only":
            train_only = True
        elif o == "--skip_x":
            export_x = False
        elif o == '--seperate_x':
            seperate_x = True
        elif o in ('-d', '--deterministic'):
            random.seed(a)
        elif o in ('-F', '--formatter'):
            formatter = choose_formatter(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if fastafile is None or fileout is None or bam is None:
        usage("required input missing")

    if seperate_x:
        exporter = SplitTextExporter(formatter)
    else:
        exporter = TextExporter(formatter)

    if train_only:
        exporter.train_prob = 1.0
        exporter.val_prob = 0.0  # technically irrelevant, but mentally nicer

    exporter.open(fileout)

    fasta = fasta2seqs(fastafile)
    alignment_file = HTSeq.BAM_Reader(bam)

    genome = Genome(alignment_file, fasta)
    genome.fill_chromosomes(size)

    for x, y in genome.slice_xy_pairs(n_pieces):
        exporter.export(x, y, export_x)

    exporter.close()


if __name__ == "__main__":
    main()
