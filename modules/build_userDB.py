import click, os, sys, numpy as np, pandas as pd
import subprocess, shutil, tempfile, gzip, re, json
from multiprocessing import Pool

try:
    from configure import executables, logging, readFasta, get_md5
except Exception:
    from .configure import executables, logging, readFasta, get_md5

bindash = executables['bindash']

# ============================================================================
# Utility Functions
# ============================================================================

def makedirs(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

def move(src, tgt):
    subprocess.Popen(['mv', '-f', src, tgt]).communicate()

def copy(src, tgt):
    subprocess.Popen(['cp', '-f', src, tgt]).communicate()

# ============================================================================
# Bindash Comparison Functions
# ============================================================================

def compare_with_existing(output_dir, user_genomes, existing_db, module, threads, cutoff=0.01):
    """
    Compare user assemblies with existing database to identify novel genomes.
    Returns DataFrame with comparison results and ANI assignments.
    """
    tag = os.path.basename(output_dir)
    # Create list file for user genomes
    existing_metadata = pd.read_feather(os.path.join(existing_db, module, f'{module}.db'))
    user_metadata = {os.path.basename(genome):{'accession': os.path.basename(genome), 'genome_path': genome,
                        'score':-idx, 'species_taxid':0, 'taxonomy': f'ANI95_{os.path.basename(genome)}',
                        'organism_name': os.path.basename(genome), 'status': 'DOWNLOADED',
                        'ANI90': os.path.basename(genome), 'ANI95': os.path.basename(genome), 'ANI98': os.path.basename(genome), 'ANI99': os.path.basename(genome),
                        'USCG_shared': 0, 'USCG_specific': 0, 'best_hit':['', 0, 0]} for idx, genome in enumerate(user_genomes)}
    user_list = os.path.join(output_dir, f'{tag}.files')
    with open(user_list, 'wt') as fout:
        for genome in user_genomes:
            fout.write(f'{genome}\n')
   
    # Sketch user genomes
    user_sketch = os.path.join(output_dir, f'{tag}.sketch')
    cmd = f'{bindash} sketch --nthreads={threads} --listfname={user_list} --kmerlen=21 --sketchsize64=160 --outfname={user_sketch}'
    subprocess.Popen(cmd.split()).communicate()
   
    # Compare with existing database
    existing_sketch = os.path.join(existing_db, module, f'{module}.sketch')
    dist_file = os.path.join(output_dir, f'{tag}.dist')

    cmd = f'{bindash} dist --nthreads={threads} --mthres=0.15 --outfname={dist_file} {user_sketch} {user_sketch}'
    subprocess.Popen(cmd.split()).communicate()

    dist_data = pd.read_csv(dist_file, sep='\t', header=None)
    for row in dist_data.values:
        query = os.path.basename(row[0])
        sbj = os.path.basename(row[1])
        if user_metadata[query]['score'] < user_metadata[sbj]['score']:
            ani = 100 - row[2]*100
            if user_metadata[query]['best_hit'][0] == '' or user_metadata[query]['best_hit'][1] < ani :
                user_metadata[query]['best_hit'] = [sbj, ani, 0]
   
    cmd = f'{bindash} dist --nthreads={threads} --mthres=0.15 --outfname={dist_file} {user_sketch} {existing_sketch}'
    subprocess.Popen(cmd.split()).communicate()
   
    dist_data = pd.read_csv(dist_file, sep='\t', header=None)
    for row in dist_data.values:
        query = os.path.basename(row[0])
        ani = 100 - row[2]*100
        if user_metadata[query]['best_hit'][0] == '' or user_metadata[query]['best_hit'][1] < ani :
            user_metadata[query]['best_hit'] = [existing_metadata.loc[existing_metadata['genome_path'] == row[1], 'accession'].values[0], ani, 1]

    for genome in user_genomes:
        query = os.path.basename(genome)
        best_hit, ani, old = user_metadata[query]['best_hit']
        if best_hit == '' :
            continue
        match = dict(zip(existing_metadata.columns, existing_metadata.loc[existing_metadata['accession'] == best_hit].values[0])) if old == 1 else user_metadata[best_hit]

        if ani >= 90 :
            user_metadata[query]['ANI90'] = match['ANI90']
            if ani >= 95 :
                user_metadata[query]['ANI95'] = match['ANI95']
                user_metadata[query]['taxonomy'] = match['taxonomy']
                if ani >= 98 :
                    user_metadata[query]['ANI98'] = match['ANI98']
                    if ani >= 99 :
                        user_metadata[query]['ANI99'] = match['ANI99']
                        user_metadata[query]['status'] = 'REDUNDANT'
            else :
                user_metadata[query]['taxonomy'] = re.split(r';s__', match['taxonomy'])[0]

    user_metadata = pd.DataFrame.from_dict(user_metadata.values())
    user_metadata = user_metadata.drop(columns=['best_hit'])
    user_metadata['score'] = 0
    return user_metadata

# ============================================================================
# USCG Extraction Functions (from second script)
# ============================================================================

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-',
              'a':'t', 't':'a', 'g':'c', 'c':'g'}

def rc(seq):
    return ''.join([complement.get(s, 'N') for s in seq[::-1]])

def detranseq(d, orf_coords):
    s, e = orf_coords[d[0]]
    d[0] = d[0].rsplit('_', 1)[0]
   
    if s < e:
        d[3], d[4] = s + (d[3]-1)*3, s + (d[4]*3 - 1)
        d[7] = 1
    else:
        d[3], d[4] = -(s - (d[3]-1)*3), -(s - (d[4]*3 - 1))
        d[7] = -1
    return d

def parseHits(data, c):
    data.sort(key=lambda x:[-x[8], x[0], x[3], x[4]])
    outputs = []
    for d in data:
        ingroup = False
        for oo in outputs:
            if d[0] == oo[0][0] and (d[3] > 0) == (oo[0][3] > 0) and \
                    min(abs(oo[0][3] - d[4]), abs(oo[0][4] - d[3])) < 10000:
                ingroup = True
                minus = 0
                for o in oo:
                    if (d[3] - o[3]) * (d[5] - o[5]) < 0 or (d[4] - o[4]) * (d[6] - o[6]) < 0:
                        ingroup = False
                        break
                    o1 = min(o[4], d[4]) - max(o[3], d[3]) + 1
                    o2 = min(o[6], d[6]) - max(o[5], d[5]) + 1
                    if o1 >= 0.9 * min(o[4] - o[3]+1, d[4] - d[3]+1) or o2 >= 0.9 * min(o[6] - o[5]+1, d[6] - d[5]+1):
                        ingroup = False
                        break
                    if o2 > 0:
                        minus += o2/(d[6] - d[5]+1)*d[8]
                if ingroup:
                    d[8] = d[8] - minus if d[8] > minus else 0.
                    oo.append(d)
                    break
        if not ingroup:
            outputs.append([d])
   
    o2 = []
    for dat in outputs:
        score = 0
        cov = np.zeros(dat[0][10])
        for d in dat:
            cov[d[5]-1:d[6]] = 1
            score += d[8]
        n_cov = np.sum(cov)
        s, e = min([d[5] for d in dat]), max([d[6] for d in dat])
        if score >= c[0] and (min(c[2], e) - max(c[1], s)+1) >= (c[2]-c[1]+1) * 0.8 and n_cov >= (c[2]-c[1]+1) * 0.6:
            o2.append(sorted(dat, key=lambda d:[d[3], d[4]]))
    return o2

def get_uscgs(query, dbname, acc, dirname, domain, all_hits=False):
    outfile = f'{acc}.USCGs.ffn'
   
    subprocess.Popen('{getorf} -table 4 -minsize 100 -sequence {0} -nomethionine -outseq {1}'.format(
        query, os.path.join(dirname, f'{acc}.aa'), **executables).split(),
        stderr=subprocess.PIPE).communicate()

    orf_coords = {}
    with open(os.path.join(dirname, f'{acc}.aa'), 'rt') as fin:
        for line in fin:
            if line.startswith('>'):
                n, s, e = re.findall(r'>(\S+) \[(\d+) - (\d+)\]', line)[0]
                orf_coords[n] = [int(s), int(e)]

    dataset = []
    cutoffs, toMerge = {}, {}
    hmms = {}
   
    with open(os.path.join(dbname, 'ortho_group')) as fin:
        for line in fin:
            p0, p1 = line.strip().split()
            if (domain in p1) or all_hits:
                fn = os.path.join(dbname, 'hmms', f'{p0}.hmm')
                hmms[p0] = [fn, 1] if domain == p1 else [fn, 0]

    with open(os.path.join(dbname, 'ortho_map')) as fin:
        for line in fin:
            p0, p1 = line.strip().split()
            if p0 in hmms:
                toMerge[p1] = p0
                hmms[p1] = [os.path.join(dbname, 'hmms', f'{p1}.hmm'), -1]

    with open(os.path.join(dbname, 'scores_cutoff')) as fin:
        for line in fin:
            p = line.strip().split()
            cutoffs[p[0]] = [float(p[1]), 0, 0]
   
    with open(os.path.join(dbname, 'lengths_cutoff')) as fin:
        for line in fin:
            p = line.strip().split()
            cutoffs[p[0]][1:] = [float(p[2]), float(p[3])]

    for hmm, _ in sorted(hmms.values()):
        data = []
        key = os.path.basename(hmm)[:-4]
        subprocess.Popen('{hmmsearch} --notextw --noali -T {3} --cpu {4} --domT {3} --domtblout {2} {0} {1}'.format(
            hmm, os.path.join(dirname, f'{acc}.aa'), os.path.join(dirname, f'{acc}.hmm'),
            cutoffs[key][0]*0.4, 2, **executables
        ).split(), stdout=subprocess.PIPE).communicate()
       
        with open(os.path.join(dirname, f'{acc}.hmm')) as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                p = line.strip().split()
                d = [p[0], toMerge.get(p[3], p[3]), float(p[21]), int(p[17]), int(p[18]),
                     int(p[15]), int(p[16]), float(p[12]), float(p[13]), int(p[2]), int(p[5])]
                data.append(detranseq(d, orf_coords))
       
        if len(data):
            dataset.extend(parseHits(data, cutoffs[key]))

    for d in dataset:
        if d[0][3] < 0:
            for dd in d:
                dd[3:5] = -dd[4], -dd[3]
            d.sort(key=lambda dd:dd[3])

    dataset.sort(key=lambda d:[d[0][0], d[0][3]])

    for i, d1 in enumerate(dataset):
        if d1[0][0] == '':
            continue
        for d2 in dataset[i+1:]:
            if d2[0][0] == '' or d1[0][0] != d2[0][0]:
                break
            if d1[0][1] != d2[0][1]:
                continue
            overlap = min(d1[-1][4], d2[-1][4]) - max(d1[0][3], d2[0][3]) + 1
            if overlap < 0:
                break
            if overlap >= 0.8 * (d1[-1][4] - d1[0][3]+1) or overlap >= 0.8 * (d2[-1][4] - d2[0][3]+1):
                if (d1[-1][4] - d1[0][3]+1) >= (d2[-1][4] - d2[0][3]+1):
                    d2[0][0] = ''
                else:
                    d1[0][0] = ''
                    break
   
    dataset = sorted([d for d in dataset if d[0][0] != ''], key=lambda x:[x[0][1], x[0][0], x[0][3]])

    sequences = readFasta(query)
    ids = {}
   
    with open(os.path.join(dirname, outfile), 'wt') as ffn_out:
        for data in dataset:
            s = sequences[data[0][0]][data[0][3]-1:data[-1][4]]
            if data[0][7] < 0:
                s = rc(s)
            ids[data[0][1]] = ids.get(data[0][1], 0) + 1
            n = f"{data[0][1]}__{acc}__{ids[data[0][1]]}"
            coding = ','.join([f'{d[3] - data[0][3]+1}-{d[4] - data[0][3]+1}' for d in data]) if data[0][7] > 0 \
                else ','.join([f'{data[-1][4]-d[4]+1}-{data[-1][4]-d[3]+1}' for d in data[::-1]])
            ffn_out.write(f'>{n} {data[0][0]} {data[0][3]} {data[-1][4]} {data[0][7]} {coding}\n{s}\n')
   
    return os.path.abspath(os.path.join(dirname, outfile)), \
           np.sum([hmms[hmm][1] != 1 for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int), \
           np.sum([hmms[hmm][1] for hmm, cnt in ids.items() if cnt == 1 or all_hits], dtype=int)

def prepare_uscg(in_fna, db, acc, tmpdir, domain, all_hits=False):
    if in_fna.lower().endswith('.gz'):
        subprocess.Popen('{gzip} -cd {0} > {1}'.format(
            in_fna, os.path.join(tmpdir, f'{acc}.fna'), **executables),
            shell=True).communicate()
        in_fna = os.path.join(tmpdir, f'{acc}.fna')

    if os.stat(in_fna).st_size >= 0.6*1000000000:
        return '', -1, -1

    uscg_ffn, n_shared, n_specific = get_uscgs(in_fna, db, acc, tmpdir, domain, all_hits)
    return uscg_ffn, n_shared, n_specific

def parse_uscg(acc, fn, alleles, profiles, fout):
    seqs = readFasta(fn)
    profiles[acc] = {}
    for n, s in seqs.items():
        gene = n.split('_', 1)[0]
        allele = get_md5(s)
        key = f'{gene}_{allele}'
        if key not in alleles:
            alleles[key] = 1
            fout.write(f'>{key}\n{s}\n')
        profiles[acc][key] = 1
    profiles[acc] = sorted(profiles[acc].keys())


def check_uscgs(fname, dbname, domain, all_hits=False) :
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


def process_user_genome(data):
    """Process a single user genome: extract USCGs and create profile"""
    acc, fna, uscg_db, outdir, domain = data
    output = os.path.join(outdir, f'{acc}.USCGs.ffn.gz')
    if os.path.isfile(output) and domain.lower() not in ('virus', 'viral'):
        n_uscg, n_specific = check_uscgs(output, uscg_db, domain)
    else :
        with tempfile.TemporaryDirectory(dir='.') as tmpdir:
            if domain.lower() in ('bacteria', 'archaea', 'eukaryota'):
                fas, n_uscg, n_specific = prepare_uscg(fna, uscg_db, acc, tmpdir, domain)
            else:
                f2, n_uscg, n_specific = prepare_uscg(fna, uscg_db, acc, tmpdir, domain, all_hits=True)
                fin = gzip.open(fna, 'rt') if fna.lower().endswith('.gz') else open(fna, 'rt')
                with open(f2, 'wt') as fout:
                    gene_id = 0
                    for line in fin:
                        if line.startswith('>'):
                            gene_id += 1
                            p = line[1:].strip().split(' ')
                            p[0] = f'>{p[0]}__{acc}__1'
                            fout.write(' '.join(p)+'\n')
                        else:
                            fout.write(line)
                fin.close()
                fas = f2
           
            subprocess.Popen('{gzip} -c {0} > {1}'.format(fas, output, **executables), shell=True).communicate()
   
    return acc, output, n_uscg, n_specific

# ============================================================================
# Main Processing Function
# ============================================================================

@click.command()
@click.option('-i', '--input_list', help='File containing paths to user assemblies (one per line)', required=True)
@click.option('-e', '--existing_db', help='Path to existing database folder. default: /titan/databases/ncbi_20251109/', default='/titan/databases/ncbi_20251109/')
@click.option('-m', '--module', help='Module name from existing database', default='bacteria')
@click.option('-o', '--output_dir', help='Output directory for user database', required=True)
@click.option('-r', '--uscg_reference', help='USCG reference database path', default='/titan/softwares/ucgMLST/db/uscgs/')
@click.option('-d', '--domain', help='Domain: bacteria [default], archaea, eukaryota, or virus', default='bacteria')
@click.option('-c', '--cutoff', help='ANI distance cutoff for novelty (default: 0.01 = 99% ANI)', default=0.01, type=float)
@click.option('-t', '--threads', help='Number of threads [default: 20]', default=20, type=int)
@click.option('-p', '--representative', help='ANI99 [default], ANI98, ANI95, or ANI90', default='ANI99,ANI98')
@click.option('--min_uscg', help='Minimum number of USCGs for QC [default: 60]', default=60, type=int)
@click.option('--n_proc', help='Number of parallel processes [default: 10]', default=10, type=int)
def process_user_assemblies(input_list, existing_db, module, output_dir, uscg_reference,
                           domain, cutoff, threads, representative, min_uscg, n_proc):
    """
    Process user assemblies: compare with existing database, extract USCGs,
    and build a temporary database in the same structure.
    """
    pool = Pool(n_proc)
   
    # Read user assembly paths
    representatives = sorted(representative.split(','), reverse=True)
    with open(input_list, 'rt') as fin:
        user_genomes = [os.path.abspath(line.strip()) for line in fin if line.strip()]
   
    logging.info(f'Loaded {len(user_genomes)} user assemblies')
   
    # Create output directory
    makedirs(output_dir)
    output_dir = os.path.abspath(output_dir)
    tag = os.path.basename(output_dir)
   
    # Step 1: Compare with existing database
    logging.info('Comparing user assemblies with existing database...')
   
    user_metadata = compare_with_existing(output_dir, user_genomes, existing_db, module, threads, cutoff)
    # user_metadata.to_feather(os.path.join(output_dir, f'{tag}.db'))
    # user_metadata.to_csv(os.path.join(output_dir, f'{tag}.csv'))
    # user_metadata = pd.read_feather(os.path.join(output_dir, f'user.db'))

    # Step 2: Extract USCGs for all user genomes
    logging.info('Extracting USCGs from user assemblies...')
    uscg_dir = os.path.join(output_dir, 'uscgs')
    makedirs(uscg_dir)
   
    process_data = [[accession, path, uscg_reference, uscg_dir, domain] for accession, path in user_metadata.loc[user_metadata['accession'].isin(user_metadata[representatives].values.ravel()), ['accession', 'genome_path']].values]
   
    with open(os.path.join(output_dir, 'uscg_extraction.log'), 'wt') as fout:
        for idx, (acc, output, n_shared, n_specific) in enumerate(pool.imap_unordered(process_user_genome, process_data)):
           
            if idx % 10 == 0:
                logging.info(f'Processed {idx}/{len(process_data)} genomes for USCG extraction')
           
            fout.write(f'{acc}\t{output}\t{n_shared}\t{n_specific}\n')
            user_metadata.loc[user_metadata['accession'] == acc, ['USCG_shared', 'USCG_specific']] = n_shared, n_specific
   
    pool.close()
    pool.join()
   
    # Step 3: Build USCG profiles
    logging.info('Building USCG profiles...')
    for representative in representatives :
        allele_file = os.path.join(output_dir, f'{tag}.{representative}.USCGs.alleles.gz')
        with gzip.open(allele_file, 'wt') as allele_out:
            profiles = {}
            alleles = {}
           
            for idx, row in user_metadata.loc[user_metadata['accession'] == user_metadata[representative]].iterrows() :
                acc = row['accession']
                n_shared = row['USCG_shared']
                n_specific = row['USCG_specific']
               
                if domain not in ('virus', 'viral'):
                    if n_shared + n_specific >= min_uscg:
                        output = os.path.join(uscg_dir, f'{acc}.USCGs.ffn.gz')
                        parse_uscg(acc, output, alleles, profiles, allele_out)
                    else:
                        logging.info(f'{acc} skipped: insufficient USCGs ({n_shared + n_specific} < {min_uscg})')
                else:
                    output = os.path.join(uscg_dir, f'{acc}.USCGs.ffn.gz')
                    parse_uscg(acc, output, alleles, profiles, allele_out)
        # Step 4: Index alleles for mapping
        logging.info('Indexing alleles...')
        allele_uncompressed = allele_file.replace('.gz', '')
       
        cmds = [
            f"{{pigz}} -cd {allele_file} > {allele_uncompressed}".format(**executables),
            f"{{samtools}} faidx {allele_uncompressed}".format(**executables),
            f"{{minimap2}} -x sr -T20 -d {allele_uncompressed}.mmi {allele_uncompressed}".format(**executables),
            f"rm {allele_uncompressed}"
        ]
        for cmd in cmds:
            subprocess.Popen(cmd, shell=True).communicate()
   
    # Save profiles
    with gzip.open(os.path.join(output_dir, f'{tag}.USCGs.profile.gz'), 'wt') as fout:
        json.dump(profiles, fout)
   
   
    # Step 5: Save metadata
    user_metadata.to_feather(os.path.join(output_dir, f'{tag}.db'))
    user_metadata.to_csv(os.path.join(output_dir, f'{tag}.csv'), index=False)
   
    logging.info('User database created successfully!')
    logging.info(f'Output directory: {output_dir}')


if __name__ == '__main__':
    process_user_assemblies()
