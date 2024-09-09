#!/usr/bin/env python
import os, gzip, click, re, json
import subprocess, tempfile
import numpy as np, pandas as pd
from multiprocessing import Pool


try :
    from configure import executables, readFasta, get_md5
except :
    from .configure import executables, readFasta, get_md5


complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-',
              'a':'t', 't':'a', 'g':'c', 'c':'g', '-':'-'}
def rc(seq) :
    return ''.join([ complement.get(s, 'N') for s in seq[::-1] ])





def detranseq(d, orf_coords) :
    s, e = orf_coords[d[0]]
    d[0] = d[0].rsplit('_', 1)[0]
    
    if s < e :
        d[3], d[4] = s + (d[3]-1)*3, s + (d[4]*3 - 1)
        d[7] = 1
    else :
        d[3], d[4] =-(s - (d[3]-1)*3), -(s - (d[4]*3 - 1))
        d[7] = -1
    return d


def parseHits(data, c) :
    data.sort(key=lambda x:[-x[8], x[0], x[3], x[4]])
    outputs = []
    for d in data :
        ingroup = False
        for oo in outputs :
            if d[0] == oo[0][0] and (d[3] > 0) == (oo[0][3] > 0) and \
                    min(abs(oo[0][3] - d[4]), abs(oo[0][4] - d[3])) < 10000 :
                ingroup = True
                minus = 0
                for o in oo :
                    if (d[3] - o[3]) * (d[5] - o[5]) < 0 or (d[4] - o[4]) * (d[6] - o[6]) < 0 :
                        ingroup = False
                        break
                    o1 = min(o[4], d[4]) - max(o[3], d[3]) + 1
                    o2 = min(o[6], d[6]) - max(o[5], d[5]) + 1
                    if o1 >= 0.9 * min(o[4] - o[3]+1, d[4] - d[3]+1) or o2 >= 0.9 * min(o[6] - o[5]+1, d[6] - d[5]+1) :
                        ingroup = False
                        break
                    if o2 > 0 :
                        minus += o2/(d[6] - d[5]+1)*d[8]
                if ingroup :
                    d[8] = d[8] - minus if d[8] > minus else 0.
                    oo.append(d)
                    break
        if not ingroup :
            outputs.append([d])
    o2 = []
    for dat in outputs :
        score = 0
        cov = np.zeros(dat[0][10])
        for d in dat :
            cov[d[5]-1:d[6]] = 1
            score += d[8]
        n_cov = np.sum(cov)
        s, e = min([d[5] for d in dat]), max([d[6] for d in dat])
        if score >= c[0] and (min(c[2], e) - max(c[1], s)+1) >= (c[2]-c[1]+1) * 0.8 and n_cov >= (c[2]-c[1]+1) * 0.6 :
            o2.append(sorted(dat, key=lambda d:[d[3], d[4]]))
    return o2


def check_uscgs(fname, dbname, domain,  all_hits=False) :
    toMerge = {}
    with open(os.path.join(dbname, 'ortho_map')) as fin :
        for line in fin :
            p0, p1 = line.strip().split()
            toMerge[p1] = p0
    hmms = {}
    with open(os.path.join(dbname, 'ortho_group')) as fin :
        for line in fin :
            p0, p1 = line.strip().split()
            if domain in p1 :
                fn = os.path.join(dbname, 'hmms', f'{p0}.hmm')
                hmms[p0] = [fn, 1] if domain == p1 else [fn, 0]
    
    genome_profile = {}
    with gzip.open(fname, 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                p = line[1:].strip().split(' ')
                g1 = re.split('__', p[0])[0]
                g1 = toMerge.get(g1, g1)
                if g1 not in hmms :
                    continue
                
                genome_profile[g1] = genome_profile.get(g1, []) + [[line]]
            elif g1 in hmms :
                genome_profile[g1][-1].append(line)
    ids = {}
    for g, profile in genome_profile.items() :
        profile = [ p for p in profile if p[0] != '' ]
        genome_profile[g] = profile
        ids[g] = len(profile)

    return np.sum([hmms[hmm][1] != 1 for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int), np.sum([hmms[hmm][1] for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int)


def get_uscgs(query, dbname, acc, dirname, domain, all_hits=False) :
    outfile = '{0}.USCGs.ffn'.format(acc)
    
    subprocess.Popen('{getorf} -table 4 -minsize 100 -sequence {0} -nomethionine -outseq {1}'.format(
        query, os.path.join(dirname,'{0}.aa'.format(acc)), **executables).split(), stderr=subprocess.PIPE).communicate()

    orf_coords = {}
    with open(os.path.join(dirname,'{0}.aa'.format(acc)), 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                n, s, e = re.findall('>(\S+) \[(\d+) - (\d+)\]', line)[0]
                orf_coords[n] = [int(s), int(e)]

    dataset = []
    cutoffs, toMerge = {}, {}

    hmms = {}
    with open(os.path.join(dbname, 'ortho_group')) as fin :
        for line in fin :
            p0, p1 = line.strip().split()
            if (domain in p1) or all_hits :
                fn = os.path.join(dbname, 'hmms', f'{p0}.hmm')
                hmms[p0] = [fn, 1] if domain == p1 else [fn, 0]

    with open(os.path.join(dbname, 'ortho_map')) as fin :
        for line in fin :
            p0, p1 = line.strip().split()
            if p0 in hmms :
                toMerge[p1] = p0
                hmms[p1] = [os.path.join(dbname, 'hmms', f'{p1}.hmm'), -1]

    with open(os.path.join(dbname, 'scores_cutoff')) as fin :
        for line in fin :
            p = line.strip().split()
            cutoffs[p[0]] = [float(p[1]), 0, 0]
    with open(os.path.join(dbname, 'lengths_cutoff')) as fin :
        for line in fin :
            p = line.strip().split()
            cutoffs[p[0]][1:] = [float(p[2]), float(p[3])]

    for hmm, _ in sorted(hmms.values()) :
        data = []
        key = os.path.basename(hmm)[:-4]
        subprocess.Popen('{hmmsearch} --notextw --noali -T {3} --cpu {4} --domT {3} --domtblout {2} {0} {1} '.format(
            hmm, os.path.join(dirname,'{0}.aa'.format(acc)), os.path.join(dirname,'{0}.hmm'.format(acc)),
            cutoffs[key][0]*0.4, 2, **executables
        ).split(), stdout=subprocess.PIPE).communicate()
        with open(os.path.join(dirname, '{0}.hmm'.format(acc))) as fin :
            for line in fin :
                if line.startswith('#') :
                    continue
                p = line.strip().split()
                d = [p[0], toMerge.get(p[3], p[3]), float(p[21]), int(p[17]), int(p[18]), \
                     int(p[15]), int(p[16]), float(p[12]), float(p[13]), int(p[2]), int(p[5])]
                data.append(detranseq(d, orf_coords))
        if len(data) :
            dataset.extend(parseHits(data, cutoffs[key]))

    for d in dataset :
        if d[0][3] < 0 :
            for dd in d :
                dd[3:5] = -dd[4], -dd[3]
            d.sort(key=lambda dd:dd[3])

    dataset.sort(key=lambda d:[d[0][0], d[0][3]] )

    for i, d1 in enumerate(dataset) :
        if d1[0][0] == '' :
            continue
        for d2 in dataset[i+1:] :
            if d2[0][0] == '' or d1[0][0] != d2[0][0] :
                break
            if d1[0][1] != d2[0][1] :
                continue
            overlap = min(d1[-1][4], d2[-1][4]) - max(d1[0][3], d2[0][3]) + 1
            if overlap < 0 :
                break
            if overlap >= 0.8 * (d1[-1][4] - d1[0][3]+1) or overlap >= 0.8 * (d2[-1][4] - d2[0][3]+1) :
                if (d1[-1][4] - d1[0][3]+1) >= (d2[-1][4] - d2[0][3]+1) :
                    d2[0][0] = ''
                else :
                    d1[0][0] = ''
                    break
    dataset = sorted([ d for d in dataset if d[0][0] != '' ], key=lambda x:[x[0][1], x[0][0], x[0][3]])

    sequences = readFasta(query)
    ids = {}
    
    with open(os.path.join(dirname, outfile), 'wt') as ffn_out :
        for data in dataset :
            s = sequences[data[0][0]][data[0][3]-1:data[-1][4]]
            if data[0][7] < 0 :
                s = rc(s)
            ids[data[0][1]] = ids.get(data[0][1], 0) + 1
            n = '{0}__{1}__{2}'.format(data[0][1], acc, ids[data[0][1]])
            coding = ','.join([ '{0}-{1}'.format(d[3] - data[0][3]+1, d[4] - data[0][3]+1) for d in data]) if data[0][7] > 0 \
                else ','.join([ '{0}-{1}'.format(data[-1][4]-d[4] +1, data[-1][4]-d[3]+1) for d in data[::-1]])
            ffn_out.write('>{0} {2} {3} {4} {5} {6}\n{1}\n'.format(
                n, s, data[0][0], data[0][3], data[-1][4], data[0][7], coding))
    return os.path.abspath(os.path.join(dirname, outfile)), \
            np.sum([hmms[hmm][1] != 1 for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int), \
            np.sum([hmms[hmm][1] for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int)
    

def prepare_uscg(in_fna, db, acc, tmpdir, domain, all_hits=False) :
    if in_fna.lower().endswith('.gz') :
        subprocess.Popen('{gzip} -cd {0} > {1}'.format(
            in_fna,
            os.path.join(tmpdir, f'{acc}.fna'),
            **executables), shell=True).communicate()
        in_fna = os.path.join(tmpdir, f'{acc}.fna')

    if os.stat(in_fna).st_size >= 3.5*1000000000 :
        return '', -1, -1

    uscg_ffn, n_shared, n_specific = get_uscgs(in_fna, db, acc, tmpdir, domain, all_hits)
    return uscg_ffn, n_shared, n_specific


def process_query(data) :
    acc, fna, db, outdir, domain = data
    output = os.path.join(outdir, f'{acc}.USCGs.ffn.gz')
    if os.path.isfile(output) and domain.lower() not in ('virus', 'viral') :
        n_uscg, n_specific = check_uscgs(output, db, domain)
    else :
        with tempfile.TemporaryDirectory(dir='.') as tmpdir :        
            if domain.lower() in ('bacteria', 'archaea', 'eukaryota') :
                fas, n_uscg, n_specific = prepare_uscg(fna, db, acc, tmpdir, domain)
            else :
                f2, n_uscg, n_specific = prepare_uscg(fna, db, acc, tmpdir, domain, all_hits=True)
                fin = gzip.open(fna, 'rt') if fna.lower().endswith('.gz') else open(fna, 'rt')
                with open(f2, 'wt') as fout:
                    gene_id = 0
                    for line in fin :
                        if line.startswith('>') :
                            gene_id += 1
                            p = line[1:].strip().split(' ')
                            p[0] = f'>{p[0]}__{acc}__1'
                            fout.write(' '.join(p)+'\n')
                        else :
                            fout.write(line)
                fin.close()
                fas = f2
            subprocess.Popen('{gzip} -c {0} > {1}'.format(
                fas, output, 
                **executables), shell=True).communicate()
    return acc, output, n_uscg, n_specific



@click.command()
@click.option('-d', '--db', help='folder for the database')
@click.option('-D', '--domain', help='one of bacteria [default], archaea, eukaryota, or virus. Use virus to skip USCG step (only do low complex filtering)', default='bacteria')
@click.option('-p', '--representative', help='ANI99, ANI98 [default], ANI95, or ANI90', default='ANI98')
@click.option('-r', '--reference', help='uscg references. default ucgMLST/db/uscgs', default='/titan/softwares/ucgMLST/db/uscgs/')
@click.option('-m', '--min_shared', help='minimum number of unique, corss-domain uscg genes for quality control. default: 45', default=45, type=int)
@click.option('-M', '--min_specific', help='minimum number of unique, domain specific uscg genes for quality control. default: 5', default=5, type=int)
def main(db, reference, domain, representative, min_shared, min_specific) :
    assert db, 'db needs to be specified.'
    pool = Pool(4)
    
    if not db.endswith('.db') :
        db = os.path.join(db, os.path.basename(db)+'.db')
    metadata = pd.read_feather(db)
    fasta_files = []
    fasta_files = metadata.loc[metadata['accession'] == metadata[representative], ['accession', 'genome_path']].values
    
    out = db[:-3]
    if not os.path.isdir(out) :
        os.makedirs(out)

    profiles = {}
    alleles = {}
    with open(f"{out}.USCGs.status", "wt") as fout, gzip.open(f"{out}.USCGs.alleles.gz", 'wt') as allele_out :
        for acc, output, n_shared, n_specific in pool.imap_unordered(process_query, [[acc, fas, reference, out, domain] for acc, fas in fasta_files ]) :
            fout.write(f"{acc}\t{output}\t{n_shared}\t{n_specific}\n")
            if domain not in ('virus', 'viral') :
                if (n_shared >= min_shared) and (n_specific >= min_specific) :
                    parse_uscg(acc, output, alleles, profiles, allele_out)
            elif (n_shared < min_shared) and (n_specific < min_specific) :
                parse_uscg(acc, output, alleles, profiles, allele_out)
    
    with gzip.open(f"{out}.USCGs.profile.gz", "wt") as fout :
        json.dump(profiles, fout)        
    
    cmds = ['{pigz} -cd {0}.USCGs.alleles.gz > {0}.USCGs.alleles'.format(out, **executables), 
           '{samtools} faidx {0}.USCGs.alleles'.format(out, **executables), 
           '{minimap2} -x sr -T20 -d {0}.USCGs.alleles.mmi {0}.USCGs.alleles'.format(out, **executables), 
           'rm {0}.USCGs.alleles'.format(out)]
    for cmd in cmds :
        subprocess.Popen(cmd, shell=True).communicate()
    
    

def parse_uscg(acc, fn, alleles, profiles, fout) :
    seqs = readFasta(fn)
    profiles[acc] = {}
    for n, s in seqs.items() :
        gene = n.split('_', 1)[0]
        allele = get_md5(s)
        key = f'{gene}_{allele}'
        if key not in alleles :
            alleles[key] = 1
            fout.write(f'>{key}\n{s}\n')
        profiles[acc][key] = 1
    profiles[acc] = sorted(profiles[acc].keys())



if __name__ == '__main__' :
    main()
