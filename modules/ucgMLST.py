#!/usr/bin/env python
import click, os, numpy as np, pandas as pd, tempfile, subprocess, json, _collections, re
from configure import executables, logging
from SRA_Funcs import get_taxonomy


def get_db(module, formal_genus, formal_species, num_threads) :
    tag = os.path.basename(module)
    aa_db = os.path.join(module, f'{tag}.USCGs.faa.gz')
    if not os.path.isfile(aa_db) :
        na_db = os.path.join(module, f'{tag}.USCGs.alleles.gz')
        with tempfile.TemporaryDirectory(dir='.', prefix='adb_') as outdir :
            na_tmp = os.path.join(outdir, 'alleles.ffn')
            aa_tmp = os.path.join(outdir, 'alleles.faa')
            subprocess.run(f'{executables["pigz"]} -cd {na_db} > {na_tmp}', shell=True)
            subprocess.run(f'{executables["getorf"]} -table 4 -minsize 90 -sequence {na_tmp} -nomethionine -noreverse -outseq {aa_tmp}'.split())
            subprocess.run(f"sed 's/_[0-9]\{{1,\}} \[/_/g; s/]//; s/ - /_/' {aa_tmp} |{executables['pigz']} > {aa_db}", shell=True)
        subprocess.run(f'{executables["diamond"]} makedb --threads {num_threads} --in {aa_db} --db {aa_db}'.split())
    na_sizes = pd.read_csv(os.path.join(module, f'{tag}.USCGs.alleles.fai'), sep='\t', header=None, usecols=[0, 1])
    
    md_file = os.path.join(module, f'{tag}.db')
    metadata = pd.read_feather(md_file).set_index('accession')
    if formal_genus :
        metadata = metadata.loc[[tax.find('g__')>=0 for tax in metadata['taxonomy']]]
    if formal_species :
        metadata = metadata.loc[[tax.find('s__')>=0 and tax.find('__unc') < 0 and tax.find('n__environmental') < 0 for tax in metadata['taxonomy']]]

    profile_file = os.path.join(module, f'{tag}.USCGs.profile.gz')
    profile = json.load(subprocess.Popen(f'{executables["pigz"]} -cd {profile_file}'.split(), universal_newlines=True, stdout=subprocess.PIPE).stdout)
    return aa_db, {s:n for s, n in na_sizes.values}, {genome:profile[genome] for genome in metadata.index if genome in profile}, metadata


def prepare_query(query, outdir) :
    rcmd = f'{executables["pigz"]} -cd' if query.lower().endswith('gz') else 'cat'
    qcmd = "|awk '(NR-1) % 4 < 2'|sed 's/^@/>/'" if query.lower().endswith('q.gz') or query.lower().endswith('q') else ""
    qry_file = os.path.join(outdir, 'qry_na')
    subprocess.run(f'{rcmd} {query} {qcmd} > {qry_file}', shell=True)
    qry_aa = os.path.join(outdir, 'qry_aa')
    subprocess.run(f'{executables["getorf"]} -table 4 -minsize 100 -sequence {qry_file} -nomethionine -outseq {qry_aa}_tmp'.split())
    subprocess.run(f"sed 's/_[0-9]\{{1,\}} \[/_/g; s/]//; s/ - /_/' {qry_aa}_tmp > {qry_aa}", shell=True)
    os.unlink(f'{qry_aa}_tmp')
    return qry_aa


@click.command()
@click.option('-q', '--query', help='assembly genomes', required=True)
@click.option('-d', '--dbname', help='name of the databases [default: /titan/databases/ncbi_20240609z/]', default='/titan/databases/ncbi_20240609z/')
@click.option('-m', '--modules', help='name of the modules [default: bacteria,archaea,eukaryota]', default='bacteria,archaea,eukaryota')
@click.option('-o', '--outdir', help='folder name storing the output', required=True)
@click.option('-t', '--num_threads', help='number of threads [Default: 16]', default=16, type=int)
@click.option('-I', '--min_iden', help='minimum identity of predicted universal core genes [Default: 40]', default=40, type=float)
@click.option('-M', '--min_gene', help='minimum ratio of universal core genes [Default: 0.05]', default=0.05, type=float)
@click.option('-g', '--formal_genus', help='only accept formal genus designations [Default: False]', default=False, is_flag = True)
@click.option('-s', '--formal_species', help='only accept formal species designations [Default: False]', default=False, is_flag = True)
def query_asm(query, dbname, modules, outdir, num_threads, formal_genus, formal_species, min_iden, min_gene) :
    logging.info('Reading database ...')
    np.random.seed(42)
    query = os.path.abspath(query)
    dbname = os.path.abspath(dbname)
    modules = [ os.path.join(dbname, module) for module in modules.split(',') ]
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    qry_aa = prepare_query(query, outdir)
    uscgs, profile, metadata = [], {}, []
    for module in modules :
        logging.info(f'Matching module {module} ...')
        aa_db, na_sizes, prof, md = get_db(module, num_threads, formal_genus, formal_species)
        profile.update(prof)
        metadata.append(md)
        uscgs.append(retrieve_uscg(qry_aa, aa_db, na_sizes, min_iden, num_threads))
    logging.info(f'Merging results ...')
    uscgs = mergeUSCGs(uscgs)
    logging.info(f'Obtained {len(uscgs)} UCGs.')
    matches = mapToGenomes(uscgs, profile, min_gene)
    logging.info(f'Writing outputs ...')
    write_matches(matches, pd.concat(metadata), outdir)


def write_matches(matches, metadata, outdir) :
    json_out = os.path.join(outdir, 'uscg.json')
    res = {'OTU':[]}
    for genome, match in sorted(matches.items(), key=lambda m:m[1], reverse=True) :
        species, taxonomy = get_taxonomy(metadata, genome)
        res['OTU'].append( [int(match[0]/match[1]*100000+0.5)/1000., int(match[2]), int(match[1]*1000+0.5)/1000., 
                            species[0] + ('' if len(species) == 1 else ' ({0})'.format(','.join(species[1:]))), genome, 
                            {m[1]:[m[0], int(abs(m[6])), int(abs(m[7])), int(m[2]*100+0.5)/100., int(m[3]*100+0.5)/100.] for m in match[3]}, taxonomy] )
    json.dump(res, open(json_out, 'wt'))
    os.unlink(os.path.join(outdir, 'qry_aa'))
    os.unlink(os.path.join(outdir, 'qry_na'))
    os.unlink(os.path.join(outdir, 'qry_aa.bsp.gz'))


def mergeUSCGs(res) :
    res = pd.DataFrame(np.vstack(res))
    res = res.sort_values([0, 3, 11], ascending=False)
    results = {}
    for d in res.values.tolist() :
        if d[0] not in results :
            results[d[0]] = [[d]]
        else :
            oo = results[d[0]]
            ovl = False
            for o in oo :
                
                s, e = max(d[6], o[0][6]), min(d[7], o[0][7])
                if ((e - s + 1) >= 0.5 * (d[7] - d[6] + 1)) or ((e - s + 1) >= 0.5 * (o[0][7] - o[0][6] + 1)) :
                    o.append(d)
                    ovl = True
                    break
            if not ovl :
                results[d[0]].append([d])
    results = [rr for r in results.values() for rr in r]
    return results


def mapToGenomes(uscgs, profile, min_genes) :
    uscg_map = {}
    for ui, uscg in enumerate(uscgs) :
        for gi, gene in enumerate(uscg) :
            if gene[1] not in uscg_map :
                uscg_map[gene[1]] = [[gene[3], gene[2], ui, gi]]
            else :
                uscg_map[gene[1]].append([gene[3], gene[2], ui, gi])
    
    for gene, info in uscg_map.items() :
        uscg_map[gene] = sorted(info, reverse=True)
    
    gene_counts = {}
    genomes = _collections.defaultdict(dict)
    for genome, genes in profile.items() :
        for gene in genes :
            if gene in uscg_map :
                genomes[genome][gene] = uscg_map[gene]
        if genome in genomes :
            gene_counts[genome] = len(profile[genome])

    results = {}
    inUse = np.ones(len(uscgs), dtype=bool)
    while len(genomes) :
        genome_scores = []
        for genome, genes in list(genomes.items()) :
            idents = np.array([np.max([g[1] * inUse[g[2]] for g in gene]) for gene in genes.values()])
            if np.sum(idents > 0) < gene_counts[genome] * min_genes :
                genomes.pop(genome)
            else :
                genome_scores.append([np.sum([np.max([g[0]*inUse[g[2]] for g in gene]) for gene in genes.values()])/gene_counts[genome], \
                    np.mean(idents[idents > 0]), np.sum(idents > 0), genome])
        if not genome_scores :
            break
        score, ident, _, genome = max(genome_scores)
        matches = np.array([m + [gene] for gene, mat in genomes[genome].items() for m in mat], dtype=object)
        q1, q3 = np.quantile(matches.T[1], 0.25), np.quantile(matches.T[1], 0.75)
        delta_q = max(q3 - q1, 2)
        idx = (matches.T[1] >= q1 - 3 * delta_q)
        matches = matches[idx]
        pos_idx = (inUse[matches[:, 2].astype(int)] > 0)
        if np.sum(pos_idx) >= gene_counts[genome] * min_genes and (np.sum(pos_idx) > 0.8 * matches.shape[0] or np.mean(matches[pos_idx, 1]) > 2 + np.mean(matches[~pos_idx, 1])) :
            genes = []
            for m in sorted(matches.tolist(), reverse=True) :
                if inUse[m[2]] :
                    genes.append(uscgs[m[2]][m[3]])
                    inUse[m[2]] = False
            logging.info(f'  Identified {genome} - Identity: {ident:.3f}, Align_fraction: {100*score/ident:.3f}, UCGs: {len(genes)}. ')
            results[genome] = [score, ident, np.sum(pos_idx), genes]
        genomes.pop(genome)
    return results


def parseHits(param) :
    data, min_iden = param
    outputs = []
    for d in data[np.argsort(-data.T[11])] :
        ingroup = False
        for oo in outputs :
            if d[0] == oo[0][0] and (d[6] > 0) == (oo[0][6] > 0) and \
                    min(abs(oo[0][6] - d[7]), abs(oo[0][7] - d[6])) < 10000 :
                ingroup = True
                for o in oo :
                    if (d[6] - o[6]) * (d[8] - o[8]) < 0 or (d[7] - o[7]) * (d[9] - o[9]) < 0 :
                        ingroup = False
                        break
                if ingroup :
                    oo.append(d)
                    break
        if not ingroup :
            outputs.append([d])
    o2 = []
    for dat in outputs :
        dat = sorted(dat, key=lambda d:d[8])
        d0 = dat[0]
        for d1 in dat[1:] :
            s, sr, (e, er) = d1[8], d0[8], sorted([d0[9], d1[9]])
            sq, eq = min(d0[6], d1[6]), max(d0[7], d1[7])
            if e >= s :
                ss0, ss1 = d0[11]/(d0[9]-d0[8]+1), d1[11]/(d1[9]-d1[8]+1)
                score = d0[11] + d1[11] - min(ss0, ss1)*(e - s + 1)
                ident = (d0[2]*(d0[9] - d0[8] + 1) + d1[2]*(d1[9] - d1[8] + 1) - min(d0[2], d1[2])*(e - s + 1))/(er - sr + 1)
            else :
                score = d0[11] + d1[11]
                ident = (d0[2]*(d0[9] - d0[8] + 1) + d1[2]*(d1[9] - d1[8] + 1))/(er - sr + 1)
            d0[2], d0[11] = ident, score
            d0[6:10] = (sq, eq, sr, er)

        if d0[2] >= min_iden and (d0[9] - d0[8]+1) >= 0.5 * d0[10] :
            o2.append(d0)
    return o2


def retrieve_uscg(qry_aa, aa_db, na_sizes, min_iden, num_threads) :
    blastp = f'{qry_aa}.bsp.gz'
    subprocess.run(f'{executables["diamond"]} blastp --db {aa_db} --threads {num_threads} --compress 1 --max-target-seqs 100 --query {qry_aa} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --out {blastp}'.split(), stderr = subprocess.PIPE)
    bsp = pd.read_csv(blastp, sep='\t', header=None)
    bsp = bsp.loc[bsp[2] >= min_iden - 5]
    q = bsp[0].str.rsplit('_', n=2, expand=True)
    r = bsp[1].str.rsplit('_', n=2, expand=True)
    r[1] = r[1].astype(int)
    bsp[8], bsp[9] = bsp[8]*3 - 3 + r[1], bsp[9]*3 - 1 + r[1]
    q.values[:, 1:] = q.values[:, 1:].astype(int)
    pos = q[1] < q[2]
    bsp.loc[ pos, 6], bsp.loc[ pos, 7] = bsp.loc[ pos, 6]*3 - 3 + q.loc[ pos, 1], bsp.loc[ pos, 7]*3 - 1 + q.loc[ pos, 1]
    bsp.loc[~pos, 6], bsp.loc[~pos, 7] = bsp.loc[~pos, 6]*3 - 3 - q.loc[~pos, 1], bsp.loc[~pos, 7]*3 - 1 - q.loc[~pos, 1]
    bsp[0], bsp[1] = q[0], r[0]
    bsp[10] = np.vectorize(na_sizes.get)(r[0])
    
    bsp = bsp.sort_values(by=[1, 0, 6])
    for d in np.arange(1,10) :
        for _ in np.arange(100) :
            idx = np.where(np.all(bsp.values[:-d, :2] == bsp.values[d:, :2], 1) & ((bsp.values[:-d, 6] > 0) == (bsp.values[d:, 6] > 0)))[0] + d
            idx = idx[(bsp.values[idx-d, 7] >= bsp.values[idx, 7]) & (bsp.values[idx-d, 8] <= bsp.values[idx, 8]) & (bsp.values[idx-d, 9] >= bsp.values[idx, 9]) & (bsp.values[idx-d, 2] >= bsp.values[idx, 2])]
            if not idx.size :
                break
            bsp = bsp.reset_index(drop=True)
            bsp.loc[idx, 0] = ''
            bsp = bsp.loc[bsp[0] != '']

    idx = np.any(bsp.values[:-1, :2] != bsp.values[1:, :2], 1) | ((bsp.values[:-1, 6] > 0) != (bsp.values[1:, 6] > 0))
    scope = np.concatenate([[0], np.where(idx)[0]+1, [bsp.shape[0]]])

    b = bsp.values[scope[:-1][(scope[1:] - scope[:-1]) == 1]]
    res = [b[(b.T[2] >= min_iden) & ((b.T[9] - b.T[8]+1) >= 0.5 * b.T[10])]]

    b = np.array([bsp.values[scope[:-1][(scope[1:] - scope[:-1]) == 2]], bsp.values[scope[:-1][(scope[1:] - scope[:-1]) == 2]+1]])
    overlap = np.min(b[:, :, 9], 0) - np.max(b[:, :, 8], 0) + 1
    overlap[overlap < 0] = 0
    min_score = np.min(b[:, :, 11]/(b[:, :, 9] - b[:, :, 8] + 1), 0)
    b[0, :, 11] += b[1, :, 11] - min_score * overlap
    b[0, :, 2] = (np.sum((b[:, :, 9] - b[:, :, 8] + 1)*b[:, :, 2], 0)-np.min(b[:, :, 2], 0)*overlap)/(np.max(b[:, :, 9], 0) - np.min(b[:, :, 8], 0) + 1)
    b[0, :, 6], b[0, :, 8] = np.min(b[:, :, 6], 0), np.min(b[:, :, 8], 0)
    b[0, :, 7], b[0, :, 9] = np.max(b[:, :, 7], 0), np.max(b[:, :, 9], 0)
    res.append(b[0][(b[0].T[2] >= min_iden) & ((b[0].T[9] - b[0].T[8]+1) >= 0.5 * b[0].T[10])])

    def subset(bsp, start, end) :
        for ss, se in zip(start, end) :
            yield bsp.values[ss:se], min_iden

    for r in map(parseHits, subset(bsp, scope[:-1][(scope[1:] - scope[:-1]) > 2], scope[1:][(scope[1:] - scope[:-1]) > 2])) :
        res.extend(r)

    res = pd.DataFrame(np.vstack(res))
    res[3] = res[2]*(res[9]-res[8]+1)/res[10]
    logging.info(f'Identified {len(res)} potential matches.')
    return res
            


if __name__ == '__main__' :
    query_asm()