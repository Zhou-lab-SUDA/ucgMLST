import click, subprocess, re, collections, numpy as np, gzip

# read in BAM file
# for each site
#    get proportion of bases
#    compare actual ratio and uncertainties
#    find SNPs that are reliable

def readFasta(fname) :
    seq = collections.OrderedDict()
    with open(fname, 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                n = line[1:].split()[0]
                seq[n] = []
            else :
                seq[n].extend(line.strip().split())
    for n, s in seq.items() :
        seq[n] = ''.join(s)
    return seq


@click.command()
@click.option('-r', '--reference')
@click.option('-b', '--bam')
@click.option('-o', '--outfile', help='output in FASTQ format. default: <bam>.fastq.gz', default=None)
@click.option('-p', '--proportion', help='proportion of the majority type. default: 1.', default=1., type=float)
@click.option('-u', '--uncertain', help='proportion of the majority type in uncertain sites. default: 0.5', default=0.5, type=float)
@click.option('-d', '--min_depth', help='minumin read depth for a site. default: 3', default=3, type=int)
def main(reference, bam, outfile, proportion, uncertain, min_depth) :
    if outfile == None :
        outfile = f'{bam}.fastq.gz'
    elif not outfile.endswith('.gz') :
        outfile = f'{outfile}.gz'

    ref = readFasta(reference)

    consensus = collections.OrderedDict({n:list(s) for n, s in ref.items() })
    qual = collections.OrderedDict({n:['!' for _ in s] for n, s in ref.items() })

    # subprocess.run(f'samtools index {bam}'.split(), check=True)

    with subprocess.Popen(f'samtools mpileup -AB -q 0 -Q 0 -f {reference} {bam}'.split(), stdout=subprocess.PIPE, universal_newlines=True) as proc :
        for line in proc.stdout :
            p = line.strip().split('\t')
            if p[0] not in consensus :
                continue

            bases = re.sub('[\$\*]', '', re.sub(r'\^.', '', p[4])).upper()
            bases = re.sub('[\.,]', p[2], bases)
            bases = ''.join([ block[int(n_indel):] for n_indel, block in re.findall('[\+\-](\d+)([A-Z]+)', '+0' + bases)])

            if len(bases) == 0 :
                continue
            base_types = collections.Counter(list(bases))
            major_type, major_p = max(base_types.items(), key=lambda t:t[1])
            consensus[p[0]][int(p[1])-1] = major_type
            if len(bases) < min_depth :
                continue
            elif len(base_types) == 1 :
                qual[p[0]][int(p[1])-1] = 'I'
            else :
                exp1, exp0 = 0., 0.
                for btype, bprop in base_types.items() :
                    if btype == major_type :
                        exp1 += np.log(proportion)*bprop
                        exp0 += np.log(uncertain)*bprop
                    else :
                        exp1 += np.log(1-proportion)*bprop
                        exp0 += np.log(1-uncertain)*bprop
                qual[p[0]][int(p[1])-1] = 'I' if exp0 - exp1 < np.log(0.05) else '!'

    # for name, seq in consensus.items() :
    #     seqLen = max(seq.keys()) + 1
    #     s = ['N' for _ in np.arange(seqLen)]
    #     q = ['!' for _ in np.arange(seqLen)]
    #     for site, base in seq.items() :
    #         s[site] = base
    #         q[site] = qual[name][site]
    #     consensus[name] = ''.join(s)
    #     qual[name] = ''.join(q)

    with gzip.open(outfile, 'wt') as fout :
        for name, seq in consensus.items() :
            s = ''.join(seq)
            q = ''.join(qual[name])
            fout.write(f'@{name}\n{s}\n+\n{q}\n')




if __name__ == '__main__' :
    main()
