#!/usr/bin/env python
import os, click, numpy as np, pandas as pd
import subprocess, re, gzip
from multiprocessing import Pool
import json

from configure import executables, logging
from SRA_Funcs import generate_outputs, write_seq, get_matches


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
@click.option('-d', '--dbname', help='name of the databases [default: /titan/databases/ncbi_20240609z/]', default='/titan/databases/ncbi_20240609z/')
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
