#!/usr/bin/env python
import os, click, numpy as np, pandas as pd
import subprocess, re, gzip
from multiprocessing import Pool
import json

from configure import executables, logging
from SRA_Funcs import generate_outputs, write_seq


def uscg2frag(read_maps, uscgs, block_size) :
    sites = np.array([(read_maps.T[3] // block_size), ((read_maps.T[4]-1) // block_size)], dtype=int).T
    
    new = []
    while sites.shape[0] > 0 :
        read_maps.T[0] = sites.T[0]
        new.append(read_maps)
        sites.T[0] += 1
        read_maps = read_maps[sites.T[0] <= sites.T[1]].copy()
        sites = sites[sites.T[0] <= sites.T[1]]

        new[-1].T[4] = np.min([(new[-1].T[0]+1) * block_size, new[-1].T[4]], 0) - np.max([new[-1].T[0] * block_size, new[-1].T[3]], 0)
        new[-1] = new[-1][new[-1].T[4] >= 10]
        
    read_maps = np.vstack(new)
    tmp = np.zeros(max(uscgs.keys())+1, dtype=np.uint32)
    for k, (s,e) in uscgs.items() :
        tmp[k] = s
    read_maps.T[0] += tmp[read_maps.T[1]]
    read_maps.T[3] = 0
    return read_maps


def find_cov_outlier(covs, block_size) :
    q1, q3 = np.quantile(covs.T[0], 0.25), np.quantile(covs.T[0], 0.75)
    delta_q = max(q3 - q1, 1.5/block_size)
    idx = (covs.T[1]/covs.T[2] <= q3 + 3*delta_q)
    return idx



def get_matches(metadata, genome_info, uscgs, read_maps, allowed_distance=2, min_gene_match=3, block_size=500) :
    '''frag_id, gene_id, read_id, start/cur_diff, end/size, mutation, diff
        0          1       2       3                4       5        6   '''
    allowed_distance *= 10000.
    pos_genes = set(read_maps[read_maps.T[6] <= allowed_distance, 0])
    genomes = {genome:len([g for g, s in genes if g in pos_genes]) for genome, genes in genome_info.items()}
    genomes = {genome:genes for genome, genes in genome_info.items() if genomes[genome] >= min_gene_match or genomes[genome]*3 >= len(genes)}
    
    uscgs = {g:s for gg in genomes.values() for g, s in gg}
    read_maps = read_maps[pd.Series(read_maps.T[0]).isin(uscgs)]
    if read_maps.shape[0] <= 0 :
        return [], []
    
    d = np.array(sorted(uscgs.items()))
    fragments = []
    while d.shape[0] > 0 :
        fragments.append(d.copy())
        fragments[-1][fragments[-1][:, 1] > block_size, 1] = block_size
        d = d[d.T[1] > block_size]
        d[:, 1] -= block_size
    fragments = np.vstack(fragments)
    fragments = fragments[np.argsort(fragments.T[0], kind='mergesort')]
    fragments.T[1] += int(0.1 * block_size)
    
    change_indices = np.concatenate([[0], np.where(np.diff(fragments.T[0]) > 0)[0]+1, [len(fragments)]])
    uscg_frag = {fragments[start, 0]: [start, end] for start, end in zip(change_indices[:-1], change_indices[1:])}
    fragments = fragments.T[1]
    
    genome2 = {genome:np.array([[i, fragments[i], g] for g, s in genes for i in range(*uscg_frag[g])], dtype=int) for genome, genes in genomes.items()}
    max_frag = np.max([np.max(gn.T[0]) for gn in genome2.values()])+1
    read_maps = uscg2frag(read_maps, uscg_frag, block_size)

    match_results = []
    summed_reads = np.zeros([int(np.max(read_maps.T[2]) + 1), 4], dtype=np.int32)
    summed_reads[:, 1:].fill(9999999)
    summed_reads[:, 0].fill(-1)
    coverages = [ [-1, -1, -1, genome, []] for genome in genome2.keys() ]
    
    while len(coverages) > 0 :
        frag_cov = np.bincount(read_maps.T[0], weights=np.power(0.398107171, read_maps.T[6].astype(np.float64) / 100.), minlength=max_frag)
        
        max_i = -1
        for i, (depth, n_gene, n_frag, genome, g_cov) in enumerate(coverages) :
            if depth == -1 or max_i < 0 or depth >= coverages[max_i][0] :
                genes = genome2[genome]
            else :
                break
            
            covs = np.array([ [((frag_cov[f]+.5)/(s+.5)), frag_cov[f], s, g] \
                                  if f < frag_cov.size else [0., 0, s, g] for f, s, g in genes ])
            idx = find_cov_outlier(covs, block_size)
            if len(g_cov) > 0 :
                if np.sum((covs[idx, 1] < 0.02 * g_cov[idx]) | ((covs[idx, 1] < 0.2 * g_cov[idx]) & (covs[idx, 1] < 0.5))) >= 0.8 * np.sum(g_cov[idx] > 0) :
                    coverages[i] = [0, 0, 0, genome, g_cov]
                    continue
            cov = np.sum(covs[idx, 1]) / np.sum(covs[idx, 2])
            n_frag = covs[(covs[:, 1] >= 1.) & idx, 3].shape[0]
            n_gene = np.unique(covs[idx & (covs[:, 1] >= 1.), 3]).size
            coverages[i] = [cov, n_gene, n_frag, genome, covs.T[1] if len(g_cov) == 0 else g_cov ]
            if (max_i < 0 or cov > coverages[max_i][0]) and \
                (n_gene >= min_gene_match or n_gene*3 >= len(genomes[genome])) and \
                (n_frag >= min_gene_match or n_frag*3 >= len(genome2[genome])) :
                max_i = i
        if max_i < 0 :
            break
        (depth, n_gene, n_frag, match, g_cov) = coverages[max_i]
        logging.info(f'    Ref: {match} with {n_gene} USCGs. ')
        coverages = [c for c in sorted(coverages, reverse=True) if c[0] > 0]
        
        match_genes = [g for g, s in genomes[match]]
        matches = read_maps[pd.Series(read_maps.T[1]).isin(match_genes)]
        match_covs = dict(zip(*np.unique(matches[matches.T[6] <= allowed_distance, 0], return_counts=True)))
        
        covs = np.array([[(match_covs.get(g, 0) + 0.5)/(s + 0.5), match_covs.get(g, 0), s, og] for g, s, og in genome2[match]])
        idx = find_cov_outlier(covs, block_size)

        # report hits
        match_results.append([match, np.sum(covs[idx].T[1])/np.sum(covs[idx].T[2])])

        frag1, frag2 = set(genome2[match][idx, 0]), set(genome2[match][~idx, 0])
        reads_ignored = set(matches[pd.Series(matches.T[0]).isin(frag2), 2]) - set(matches[pd.Series(matches.T[0]).isin(frag1), 2])
        
        matched_reads = np.zeros([summed_reads.shape[0], 3], dtype=np.int32)
        matched_reads[:] = 9999999
        
        _, ridx = np.unique(matches.T[2], return_index=True)
        matched_reads[matches[ridx, 2], :] = matches[ridx, :][:, (0,1,5)]
        matched_reads[list(reads_ignored)] = 9999999
        # report read matches
        idx = (matched_reads.T[2] < summed_reads.T[3])
        summed_reads[idx, 0] = len(match_results) - 1
        summed_reads[idx, 1:] = matched_reads[idx]
        
        reads_todrop = set(matches[matches.T[6] <= allowed_distance, 2]) - reads_ignored
        frag_todrop = set(matches.T[0])
        read_maps = read_maps[~pd.Series(read_maps.T[0]).isin(frag_todrop) & ~pd.Series(read_maps.T[2]).isin(reads_todrop)]
        
        coverages = [ c for c in coverages if c[3] != match ]

    for m in np.unique(summed_reads[summed_reads.T[0] >= 0, 0]):
        p = summed_reads[summed_reads.T[0] == m, 3]
        q1, q3 = np.quantile(p, 0.25), np.quantile(p, 0.75)
        delta_q = max(q3-q1, 100)
        scope = q3 + 3*delta_q
        summed_reads[(summed_reads.T[0] == m) & (summed_reads.T[3] > scope), :] = -1

    i1 = 0
    res = []
    for i0 in np.unique(summed_reads[summed_reads.T[0] >= 0, 0]) :
        ref = match_results[i0][0]
        x = summed_reads[summed_reads.T[0] == i0]
        n_frag = np.unique(x.T[1]).size
        n_gene = np.unique(x.T[2]).size
        if x.shape[0] >= min_gene_match and \
            (n_frag >= min_gene_match or n_frag*3 >= len(genome2[ref])) and \
            (n_gene >= min_gene_match or n_gene*3 >= len(genomes[ref])) :
                res.append(match_results[i0])
                if i0 != i1 :
                    summed_reads[summed_reads.T[0] == i0, 0] = i1
                i1 += 1
        else :
            summed_reads[summed_reads.T[0] == i0, :] = [-1, 9999999, 9999999, 9999999]

    return np.array(res, dtype=object), summed_reads



def parse_paf(data) :
    outfile, tmpdir, uscg_info = data
    rmaps = []
    p = subprocess.Popen(f"{executables['pigz']} -cd {outfile}".split(), cwd=tmpdir, stdout=subprocess.PIPE,
                         universal_newlines=True)

    for line in p.stdout:
        p = line.strip().split('\t')
        if p[5] not in uscg_info:
            continue
        p[0] = int(p[0], 16)
        p[1:4]  = [ int(pp) for pp in p[1:4]  ]
        p[6:11] = [ int(pp) for pp in p[6:11] ]
        p[5] = uscg_info[p[5]]
        if p[4] == '+' :
            s, e = min(p[2], p[7]), min(p[1]-p[3], p[6] - p[8])
        else :
            s, e = min(p[1]-p[3], p[7]), min(p[2], p[6] - p[8])
        mut = (p[10]-p[9])*10 + s + e
        rmaps.append([p[5], 0, p[0], p[7], p[8], int(mut*1000/(p[10]+s+e)+0.5)])

    if len(rmaps):
        np.savez_compressed( os.path.join(tmpdir, f'{outfile}.npz'), reads=np.array(rmaps, dtype=np.uint32))
        return os.path.join(tmpdir, f'{outfile}.npz')
    return ''


def map_to_uscgs(paf_files, uscgs, tmpdir, pool) :
    rmaps = []
    for rmap_file in pool.imap_unordered(parse_paf, [ [sfile, tmpdir, uscgs] for sfile in paf_files ]) :
        if rmap_file :
            dat = np.load(rmap_file)
            rmaps.append(dat['reads'])
            try :
                os.unlink(rmap_file)
            except :
                pass
    if len(rmaps) :
        rmaps = np.vstack(rmaps)
        
        read_rename, read_dist = {}, []
        for r in rmaps :
            if r[2] not in read_rename :
                read_rename[r[2]] = len(read_rename)
                r[2] = read_rename[r[2]]
                read_dist.append(r[5])
            else :
                r[2] = read_rename[r[2]]
                if read_dist[r[2]] > r[5] :
                    read_dist[r[2]] = r[5]
        
        read_dist = np.array(read_dist)
        rmaps.T[1] = rmaps.T[0]
        rmaps = np.hstack([rmaps, (rmaps.T[5] - read_dist[rmaps.T[2]]).reshape([-1, 1])]).astype(np.int32)
    return rmaps, np.array([r[0] for r in sorted(read_rename.items(), key=lambda r:r[1])])



def map_reads(query, dbname, mode, tmpdir, max_dist, num_threads) :
    outputs = []
    total_reads = 0

    for qid, qry in enumerate(query):
        if qry.lower().endswith('q.gz') :
            qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fastq.gz'))
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            with open(qry_file, 'wb') as fout2 :
                fout = subprocess.Popen('{pigz} -c'.format(**executables).split(), stdin=subprocess.PIPE, stdout=fout2, universal_newlines=True)
                for id, line in enumerate(fin.stdout) :
                    if id%4 == 0 :
                        fout.stdin.write(f'@{total_reads:X}\n')
                        total_reads += 1
                    else :
                        fout.stdin.write(line)
                fout.communicate()
            fin.communicate()
        elif qry.lower().endswith('q') :
            qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fastq'))
            with open(qry, 'rt') as fin, open(qry_file, 'wt') as fout:
                for id, line in enumerate(fin) :
                    if id%4 == 0 :
                        fout.write(f'@{total_reads:X}\n')
                        total_reads += 1
                    else :
                        fout.write(line)
        elif qry.lower().endswith('.gz') :
            qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fasta.gz'))
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            with open(qry_file, 'wb') as fout2 :
                fout = subprocess.Popen('{pigz} -c'.format(**executables).split(), stdin=subprocess.PIPE, stdout=fout2, universal_newlines=True)
                for line in fin.stdout :
                    if line.startswith('>') :
                        fout.stdin.write(f'>{total_reads:X}\n')
                        total_reads += 1
                    else :
                        fout.stdin.write(line)
                fout.communicate()
            fin.communicate()
        else :
            qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fasta'))
            with open(qry, 'rt') as fin, open(qry_file, 'wt') as fout:
                for line in fin :
                    if line.startswith('>') :
                        fout.write(f'>{total_reads:X}\n')
                        total_reads += 1
                    else :
                        fout.write(line)

        for rid, db in enumerate(dbname) :
            outfile = f'{qid}.{rid}.paf.gz'
            uscg_db = os.path.join(db, os.path.basename(db)+'.USCGs.alleles.mmi')
            p_dist = 0.6 if mode == 'sr' else 0.6

            subprocess.Popen(
                '{minimap2} -t{3} -cx {5} -T20 --frag=yes -p{6} -N90000 -Y --end-bonus 12 -2 --secondary=yes {0} {1}|{EnFlt} {4}|{pigz} -c > {2}'.format(
                    uscg_db, qry_file, outfile, num_threads, max_dist, mode, p_dist, **executables,
                ), cwd=tmpdir, shell=True).communicate()
                
            outputs.append(outfile)
        os.unlink(qry_file)

    return outputs, total_reads


def query_sra(query, dbname, metadata, genome_info, uscg_info, output, mode, max_dist, allowed_distance, min_depth, min_consensus, pool, debug=[False, False, False, False]) :
    if not debug[0] :
        logging.info('Running read mapping...')
        paf_files, n_reads = map_reads(query, dbname, mode, output, max_dist, len(pool._pool))
        # np.savez_compressed(os.path.join(output, 'uscg.npz'), n_reads=np.array(n_reads))
        logging.info('Done')
    else :
        paf_files = [os.path.abspath(os.path.join(output, f'{id}.{jd}.paf.gz')) for id, _ in enumerate(query) for jd, _ in enumerate(dbname)]
        n_reads = int(np.load(os.path.join(output, 'uscg.npz'))['n_reads'])
 
    if not debug[1] : 
        logging.info('Extracting USCG information...')
        read_maps, r_ids = map_to_uscgs(paf_files, uscg_info, output, pool)
        logging.info('Done')
        if len(read_maps) == 0 :
            return {'profile':[], 'OTU':[]}
        # np.savez_compressed(os.path.join(output, 'uscg.npz'), reads=read_maps, r_ids=r_ids, n_reads=np.array(n_reads))
    else :
        data = np.load(os.path.join(output, 'uscg.npz'), allow_pickle=True)
        if not debug[2] :
            read_maps = data['reads']
        r_ids = data['r_ids']
    
    if not debug[2] :
        logging.info('Extracting best aligned references...')
        matches, reads = get_matches(metadata, genome_info, uscg_info, read_maps, allowed_distance=allowed_distance)
        logging.info('Done')
        if len(matches) == 0 :
            return {'profile':[], 'OTU':[]}
        # np.savez_compressed(os.path.join(output, 'out1.npz'), matches=matches, reads=reads)
    else :
        data = np.load(os.path.join(output, 'out1.npz'), allow_pickle=True)
        matches, reads = data['matches'], data['reads']
    if not debug[3] :
        logging.info('Preparing outputs...')
        outputs, bam = generate_outputs(paf_files, query, metadata, matches, reads, r_ids, output, genome_info, uscg_info, n_reads)
        # json.dump(dict(outputs=outputs, bam=bam), open(os.path.join(output, 'out2.json'), 'wt'))
        logging.info('Done')
    else :
        data = json.load(open(os.path.join(output, 'out2.json'), 'rt'))
        outputs, bam = data['outputs'], data['bam']
    write_seq(output, outputs, bam, min_depth, min_consensus)
    return outputs


def read_metadata(modules, formal_genus, formal_species) :
    metadata = []
    for db in modules :
        md_file = os.path.join(db, os.path.basename(db) + '.db')
        md = pd.read_feather(md_file)
        if formal_genus :
            md = md.loc[[tax.find('g__')>=0 for tax in md['taxonomy']]]
        if formal_species :
            md = md.loc[[tax.find('s__')>=0 and tax.find('__unc') < 0 and tax.find('n__environmental') < 0 for tax in md['taxonomy']]]
        metadata.append(md)
    return pd.concat(metadata).set_index('accession')


def read_uscg(modules, metadata) :
    genomes = {}
    uscgs = {}
    for db in modules :
        fai_file = os.path.join(db, os.path.basename(db) + '.USCGs.alleles.fai')
        gene_sizes = dict(pd.read_csv(fai_file, header=None, sep='\t', usecols=[0,1]).values)
        n_uscgs = len(uscgs)
        uscgs.update({g:i+n_uscgs for i, g in enumerate([g for g in gene_sizes.keys() if g not in uscgs])})
        
        cg_file = os.path.join(db, os.path.basename(db) + '.USCGs.profile.gz')
        cg = json.load(gzip.open(cg_file))
        g = { acc:[[uscgs[g], gene_sizes[g]] for g in cg[acc]] for acc in metadata.index if acc in cg if cg[acc][0] in uscgs }
        genomes.update(g)
    return genomes, uscgs


@click.command()
@click.option('-q', '--query', help='fastq file(s), specify --query multiple times for additional reads', required=True, multiple=True)
@click.option('-d', '--dbname', help='path of the databasesc', required=True)
@click.option('-m', '--modules', help='name of the modules [default: bacteria,archaea,viral]', default='bacteria,archaea,viral')
@click.option('-o', '--outdir', help='folder name storing the output', required=True)
@click.option('-g', '--formal_genus', help='only accept formal genus designations [Default: False]', default=False, is_flag = True)
@click.option('-s', '--formal_species', help='only accept formal species designations [Default: False]', default=False, is_flag = True)
@click.option('-t', '--num_threads', help='number of threads [Default: 16]', default=16, type=int)
@click.option('-M', '--mode', help='One of sr [default], map-ont, map-hifi, map-pb, asm20', default='sr')
@click.option('-D', '--max_dist', help='maximum distance of alignment [Default: 0.05 or 0.20 for asm20]', default=None, type=float)
@click.option('-x', '--allowed_difference', help='allowed SNP difference for top hits [Default: 0.02 for sr and map-hifi or 0.04 for others]', default=None, type=float)
@click.option('--min_depth', help='minimum read depth to call a base reliably. [Default: 3]', default=3, type=int)
@click.option('--min_consensus', help='minimum proportion of consensus to call a base reliably [Default: 0.8]', default=0.8, type=float)
def main(query, dbname, modules, outdir, mode, formal_genus, formal_species, max_dist, allowed_difference, num_threads, min_depth, min_consensus) :
    if allowed_difference == None :
        if mode in ('sr', 'map-hifi') :
            allowed_difference = 0.02
        else :
            allowed_difference = 0.04
    if max_dist == None :
        if mode in ('asm20', ) :
            max_dist = 0.2
        else :
            max_dist = 0.05
    max_dist += allowed_difference

    pool = Pool(num_threads)
    logging.info('Reading database...')
    np.random.seed(42)
    query = [os.path.abspath(qry) for qry in query]
    dbname = os.path.abspath(dbname)
    modules = [ os.path.join(dbname, module) for module in modules.split(',') ]
    
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    metadata = read_metadata(modules, formal_genus, formal_species)
    genomes, uscgs = read_uscg(modules, metadata)
    logging.info('Done')
    query_sra(query, modules, metadata, genomes, uscgs, outdir, mode, max_dist, allowed_difference, min_depth, min_consensus, pool)



if __name__ == '__main__' :
    main()
