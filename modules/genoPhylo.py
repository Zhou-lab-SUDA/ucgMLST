#!/usr/bin/env python
import os, click, gzip, glob, json, numpy as np, ete3, pandas as pd, _collections
import subprocess, tempfile, re
from multiprocessing import Pool
from configure import executables, logging, rc, readFasta


def each_minimap(data) :
    tmpdir, target, min_iden, min_presence, ftag = data
    if target.lower().endswith('q.gz') :
        fn = os.path.basename(target).rsplit('.', 2)[0]
        with gzip.open(target, 'rt') as fin, open(os.path.join(tmpdir, f'{fn}.fas'), 'wt') as fout :
            if target.lower().endswith('q.gz') :
                for lid, line in enumerate(fin) :
                    if lid % 4 == 0 :
                        fout.write('>' + line[1:])
                    elif lid % 4 == 1 :
                        fout.write(line)
            elif target.lower().endswith('.gz') :
                for line in fin :
                    fout.write(line)
        target = f'{fn}.fas'
    elif target.lower().endswith('q') :
        with open(target, 'rt') as fin, open(os.path.join(tmpdir, f'{fn}.fas'), 'wt') as fout :
            for lid, line in enumerate(fin) :
                if lid % 4 == 0 :
                    fout.write('>' + line[1:])
                elif lid % 4 == 1 :
                    fout.write(line)
    else :
        fn = os.path.basename(target).rsplit('.', 1)[0]

    t_seq = readFasta(os.path.join(tmpdir, target), upper=False)

    map_cmd = f"{executables['minimap2']} -cx asm20 -t1 --frag=yes -P -m30 -s40 -N5000 -n2 --end-bonus 12 --secondary=yes query {target}"
    p = subprocess.Popen(map_cmd.split(), cwd=tmpdir, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    res = []

    for line in p.stdout :
        if line.startswith('[') :
            continue
        p = line.strip().split('\t')
        p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
        if p[10] <= min_iden :
            continue

        p[1:4] = [int(p[1]), int(p[2]), int(p[3])]
        p[6:9] = [int(p[6]), int(p[7]), int(p[8])]

        p[11], p[12] = p[-1][5:], int(p[14][5:])
        res.append(p[:13])
    res.sort(key=lambda x:(x[5], -x[12]))
    
    res_seqs = {}
    for p in res :
        if p[5] in res_seqs :
            s = res_seqs[p[5]]
            if np.sum(s[p[7]:p[8]] != '-') >= 0.2 * (p[8] - p[7]) :
                continue
        else :
            s = np.array(['-'] * p[6])
        i0, i1, d = (p[2], p[7], 1) if p[4] == '+' else (p[3], p[7], -1)

        for n, t in re.findall('(\d+)([MDI])', p[11]) :
            n = int(n)
            if t == 'M' :
                ss = np.array(list(t_seq[p[0]][i0:(i0+n)] if p[4] == '+' else rc(t_seq[p[0]][(i0-n):i0])))
                s[i1:(i1+n)][s[i1:(i1+n)] == '-'] = ss[s[i1:(i1+n)] == '-']
                i0, i1 = i0 + n*d, i1 + n
            elif t == 'D' :
                s[i1:(i1+n)][s[i1:(i1+n)] == '-'] = '.'
                i1 += n
            elif t == 'I' :
                i0 += n*d
        res_seqs[p[5]] = s
    res = {}
    for g, s in res_seqs.items() :
        s = re.sub('[nN\.]', '-', ''.join(s))
        if len(s.replace('-', '')) >= min_presence * len(s) :
            gene = g.rsplit('_', 1)[0]
            res[gene] = s
    return fn, res, (len(t_seq), np.sum([len(s) for s in t_seq.values()])), ftag



def minimap_align(tmpdir, metadata, profiles, uscgs, genomes, min_identity, min_presence, min_presence_ref, no_risky, no_db, pool) :
    logging.info('Building phylogeny.')
    ref_seqs = {g.rsplit('_', 1)[0]:s for g, s in profiles[metadata.index[0]].items()}
    genes = sorted(ref_seqs.keys())
    with open(os.path.join(tmpdir, 'query'), 'wt') as fout :
        for n, s in ref_seqs.items() :
            fout.write(f'>{n}\n{s}\n')

    if no_db :
        profiles = {}
    for acc, db_seqs in profiles.items() :
        with open(os.path.join(tmpdir, f"{acc}.fas"), 'wt') as fout :
            for n, s in db_seqs.items() :
                fout.write(f'>{n}\n{s}\n')

    for tag, items in uscgs.items() :
        with open(os.path.join(tmpdir, f"{tag}.fas"), 'wt') as fout:
            for n, s in items.items() :
                fout.write(f'>{n}\n{s}\n')

    queries = dict([[os.path.abspath(os.path.join(tmpdir, f"{acc}.fas")), 'ref'] for acc in profiles.keys()] + \
                     [[os.path.abspath(os.path.join(tmpdir, f"{acc}.fas")), 'qry'] for acc in uscgs.keys()] + \
                     [[os.path.abspath(g), 'ref'] for g in genomes])


    qryseqs = _collections.OrderedDict([[metadata.index[0], {}]])
    concatenated_seqs = _collections.OrderedDict([[metadata.index[0], '']])
    for qry, seqs, qry_stats, ftag in pool.imap_unordered(each_minimap, [ [tmpdir, fn, min_identity, {'ref':min_presence_ref, 'qry':min_presence}[ftag], ftag] \
                                                                        for fn, ftag in queries.items() ]) :
        min_p = {'ref':min_presence_ref, 'qry':min_presence}[ftag]
        if len(seqs) >= min_p * len(genes) and len(seqs) >= min_presence_ref * qry_stats[0] :
            lq_geneNum = np.sum([len(seqs.get(g, '').replace('-', ''))/len(ref_seqs[g]) > min_p for g in genes])
            hq_geneNum = np.sum([len(re.sub('[^ACGT]', '', seqs.get(g, '')))/len(ref_seqs[g]) > min_p for g in genes])

            qryseqs[qry] = seqs
            if hq_geneNum >= min_p * len(genes) and hq_geneNum >= min_presence_ref * qry_stats[0] :
                s = ''.join([ seqs.get(g, '-' * len(ref_seqs[g])) for g in genes ])
                s2 = re.sub('[^ACGT]', '-', s)
                concatenated_seqs[qry] = s2
            elif (not no_risky) and ftag == 'qry' and lq_geneNum >= min_p * len(genes) and lq_geneNum >= min_presence_ref * qry_stats[0] :
                s = ''.join([ seqs.get(g, '-' * len(ref_seqs[g])) for g in genes ])
                concatenated_seqs[f'{qry}|low_qual'] = s

            if len(qryseqs) % 100 == 0 :
                logging.info(f'Extracted {len(qryseqs)} USCGs from both the db and the samples.')
        
    logging.info(f'Identified {len(concatenated_seqs)} samples with good sequences.')
    aln_file = os.path.join(tmpdir, 'USCG_align.fas')
    with open(aln_file, 'w') as fout :
        for n, s in concatenated_seqs.items() :
            if len(s) :
                fout.write('>{0}\n{1}\n'.format(n, s))

    subprocess.Popen('{iqtree} -nt {0} -fast -redo -s USCG_align.fas -m GTR -pre USCG_align --runs 8 --polytomy'.format(
        len(pool._pool), **executables).split(), stdout=subprocess.PIPE, cwd=tmpdir).communicate()
    tre = ete3.Tree(os.path.join(tmpdir, 'USCG_align.treefile'), format=0)
    tre.set_outgroup(tre.get_midpoint_outgroup())
    for i, n in enumerate(tre.traverse('postorder')) :
        if n.name == '' :
            n.name = f'Node{i}'

    return aln_file, tre


def readJson(dat) :
    uscgs, ref = dat
    items = _collections.defaultdict(dict)
    for uscg in uscgs :
        if os.path.basename(uscg) == 'profile.json' or os.path.basename(uscg).startswith('resolved') :
            tag = os.path.basename(os.path.dirname(uscg))
        else :
            tag = os.path.basename(uscg).rsplit('.', 1)[0]
        data = json.load(open(uscg, 'rt'))['OTU']
        for d in data :
            if len(d) < 7 or not d[6] :
                continue
            test = '\t'.join([d[3], d[5], d[4]])
            if len(re.findall(ref.replace(' ', '_'), test)) > 0 :
                for n, s in d[6].items() :
                    key = f'{tag}|{d[4]}'
                    items[key][n] = s
    return items


def extract_reference(reference, dbname, module, genus) :
    if module == 'auto' :
        modules = [os.path.join(dname, os.path.basename(dname)) for dname in glob.glob(os.path.join(dbname, '*'))]
    else :
        modules = [os.path.join(dbname, dname, dname) for dname in module.split(',')]

    for module in modules :
        metadata = pd.read_feather(f'{module}.db').set_index('accession')
        metadata = metadata.loc[metadata.status == 'DOWNLOADED']
        r = metadata[metadata.index.str.contains(reference)]
        if len(r) <= 0 :
            r = metadata[metadata['taxonomy'].str.contains(reference)]
        if len(r) <= 0 :
            continue

        reference = r.index[0]
        logging.info(f'Use {reference} as the reference.')

        records = {reference:0}
        profiles = json.load(gzip.open(f'{module}.USCGs.profile.gz'))
        ref = metadata.loc[reference]
    
        for field in ('ANI99', 'ANI98', 'ANI95') :
            for rec in metadata.loc[metadata[field] == ref[field]].index :
                if rec in profiles and rec not in records :
                    records[rec] = len(records)
        spe_info = re.findall('s__[^;]+', ref['taxonomy'])
        if len(spe_info) :
            for rec in metadata.loc[metadata['taxonomy'].str.contains(spe_info[0])].index :
                if rec in profiles and rec not in records :
                    records[rec] = len(records)
        if genus :
            for rec in metadata.loc[metadata['ANI90'] == ref['ANI90']].index :
                if rec in profiles and rec not in records :
                    records[rec] = len(records)
            gen_info = re.findall('g__[^;]+', ref['taxonomy'])
            if len(gen_info) :
                for rec in metadata.loc[metadata['taxonomy'].str.contains(gen_info[0])].index :
                    if rec in profiles and rec not in records :
                        records[rec] = len(records)
        index = np.array([r[0] for r in sorted(records.items(), key=lambda r:r[1])])
        profile = {idx:profiles[idx] for idx in index}
        alleles = {allele:[] for alleles in profile.values() for allele in alleles}
    
        p = subprocess.Popen(f'pigz -cd {module}.USCGs.alleles.gz'.split(), stdout=subprocess.PIPE, universal_newlines=True)
        for line in p.stdout :
            if line.startswith('>') :
                n = line[1:].strip()
                if n not in alleles :
                    n = ''
            elif n != '' :
                alleles[n].append(line.strip())
        profile = {idx:{ g:''.join(alleles[g]) for g in prof } for idx, prof in profile.items() }

        return metadata.loc[index], profile
    return None, None
    


def get_reference(ref, uscgs) :
    references = {}
    for key, uscgs in uscgs.items() :
        if len(re.findall(ref.replace(' ', '_'), '\t'.join(key))) > 0 :
            references[key] = [sum([len(re.sub('[nN]', '', s)) for s in uscgs.values()]), len(uscgs.keys())]
    return max(references.items(), key=lambda r:r[1])[0]



@click.command()
@click.argument('uscg_files', nargs=-1)
@click.option('-d', '--dbname', help='absolute path of the database. [required]', default="/titan/databases/ncbi_20240609z")
@click.option('-m', '--module', help='modules in the database. [default: auto detect]', default='auto')
@click.option('-r', '--reference', help='reference genome, can be one of accession_code, tax_code, or species_name. [required]', required=True)
@click.option('--min_identity', help='minumum identity of a sequence comparing to ref_acc. default: 0.93 [0. - 1.]', default=None, type=float)
@click.option('-p', '--min_presence', help='minumum coverages of genes for a query strain to be evaluated. [default: 0.1]', default=0.1, type=float)
@click.option('-P', '--min_presence_ref', help='minumum coverages of genes for a refernce genome to be included. [default: 0.75]', default=None, type=float)
@click.option('-G', '--genus', help='Include genomes from the same genus, will also set min_presence_ref=0.45 and min_identity=0.8 by default.', default=False, is_flag=True)
@click.option('--no_risky', help='flag to prohibit using low quality bases. default: False', default=False, is_flag=True)
@click.option('-N', '--no_db', help='do not include genomes in the database', default=False, is_flag=True)
@click.option('-g', '--genome_list', help='additional genomes to be included in the analysis', default='')
@click.option('--resolve_db', help='add genomes in a resolve_db into the analysis', default='')
@click.option('-u', '--uscg_list', help='list of uscgs as a file [default: None]', default='')
@click.option('-o', '--outdir', help='folder storing outputs.', required=True)
@click.option('-n', '--n_proc', help='number of processes [default: 8]', default=8, type=int)
def main(uscg_files, dbname, module, reference, min_identity, min_presence, min_presence_ref, genus, no_db, genome_list, uscg_list, outdir, n_proc, no_risky, resolve_db) :
    genoPhylo(uscg_files, dbname, module, reference, min_identity, min_presence, min_presence_ref, genus, no_db, genome_list, uscg_list, outdir, n_proc, no_risky, resolve_db)

def genoPhylo(uscg_files, dbname, module, reference, min_identity, min_presence, min_presence_ref, genus, no_db, genome_list, uscg_list, outdir, n_proc, no_risky, resolve_db) :
    pool = Pool(n_proc)
    
    if min_identity == None :
        min_identity = 0.8 if genus else 0.93
    if min_presence_ref == None :
        min_presence_ref = 0.45 if genus else 0.75
    
    if min_identity <= 1 :
        min_identity *= 100.

    if uscg_list :
        with open(uscg_list, 'rt') as fin :
            uscg_files = list(uscg_files) + [line.strip().split()[0] for line in fin]
    uscg_files = [ os.path.join(fn, 'profile.json') if os.path.isfile(os.path.join(fn, 'profile.json')) else fn for fn in uscg_files ]

    uscgs = {}
    for uscg in pool.imap_unordered(readJson, [ [uscg_files[i::n_proc], reference] for i in np.arange(n_proc) ]) :
        uscgs.update(uscg)

    if len(uscgs) :
        uscg_cnt = {}
        for tag, uscg in uscgs.items() :
            acc = tag.rsplit('|', 1)[-1]
            uscg_cnt[acc] = uscg_cnt.get(acc, 0) + len(uscg)
        reference, _ = max(uscg_cnt.items(), key=lambda n:n[1])
        logging.info(f'Kept {len(uscgs)} sets of USCGs after initial filtering.')
    metadata, profiles = extract_reference(reference, dbname, module, genus)
    logging.info(f'Extract {len(profiles)} sets of USCGs from the main database.')

    genome_lists = [genome_list, os.path.join(resolve_db, 'genome.list')]
    genomes = {}
    for genome_list in genome_lists :
        if os.path.isfile(genome_list) :
            with open(genome_list, 'rt') as fin :
                genomes.update({os.path.abspath(line.strip().split()[0]):1 for line in fin if os.path.isfile(line.strip().split()[0])})
            logging.info(f'Read in {len(genomes)} additional genomes.')

    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    if os.path.isfile(os.path.join(outdir, 'tree_info.dump')) :
        os.unlink(os.path.join(outdir, 'tree_info.dump'))

    with tempfile.TemporaryDirectory(prefix='qry_', dir='.') as tmpdir :
        aln_file, tre = minimap_align(tmpdir, metadata, profiles, uscgs, sorted(genomes.keys()), min_identity, min_presence, min_presence_ref, no_risky, no_db, pool)
        with open(os.path.join(outdir, 'uscg.concat.fas'), 'wt') as aln_out, open(aln_file, 'rt') as fin :
            aln_out.write(fin.read())
    with open(os.path.join(outdir, 'uscg.nwk'), 'wt') as nwk_out :
        nwk_out.write(tre.write(format=1)+'\n')
    if len(genomes) :
        with open(os.path.join(outdir, 'genome.list'), 'wt') as fout :
            for genome in sorted(genomes.keys()) :
                fout.write(f'{genome}\n')
    logging.info(f'Tree "uscg.nwk" generated under folder {outdir}.')



if __name__ == '__main__' :
    main()
