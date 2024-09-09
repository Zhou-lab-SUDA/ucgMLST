import click, os, sys, numpy as np, pandas as pd
import subprocess, shutil
try :
    from configure import executables, logging
except Exception :
    from .configure import executables, logging

bindash = executables['bindash']

def rmtree(dirnames) :
    if not isinstance(dirnames, list) :
        dirnames = [dirnames]
    for i in range(0, len(dirnames), 1000) :
        subprocess.Popen(['rm', '-rf'] + dirnames[i:i+1000]).communicate()

def move(src, tgt) :
    subprocess.Popen(['mv', '-f', src, tgt]).communicate()

def copy(src, tgt) :
    subprocess.Popen(['cp', '-f', src, tgt]).communicate()

def makedirs(dirname) :
    if not os.path.isdir(dirname) :
        os.makedirs(dirname)


def get_sketches(metadata_tab, dbname, module, cutoff, threads, batch_num=4000) :
    metadata_tab['status'] = 'NEW'
    metadata_tab['ANI80'] = metadata_tab['accession']
    metadata_tab['ANI90'] = metadata_tab['accession']
    metadata_tab['ANI95'] = metadata_tab['accession']
    metadata_tab['ANI98'] = metadata_tab['accession']
    metadata_tab['ANI99'] = metadata_tab['accession']
    
    sketch_dir = os.path.join(dbname, 'sketches')
    makedirs(sketch_dir)
    
    ids = np.arange(0, metadata_tab.shape[0], batch_num)
    data_ids = {fn:id for fn, id in zip(metadata_tab['genome_path'], metadata_tab.index)}
    
    for idx, i in enumerate(ids) :
        results = run_bindash(metadata_tab[i:i+batch_num], ids[:idx+1], data_ids, dbname, module, threads)
        logging.info(f'Run bindash on genomes {i}')
        for acc, tgt, iden in results :
            if iden < 1- cutoff :
                metadata_tab.loc[acc, 'status'] = 'REDUNDANT'
            else :
                metadata_tab.loc[acc, 'status'] = 'DOWNLOADED'
            if iden <= 0.2 :
                metadata_tab.loc[acc, 'ANI80'] = metadata_tab.loc[tgt, 'ANI80']
                if iden <= 0.1:
                    metadata_tab.loc[acc, 'ANI90'] = metadata_tab.loc[tgt, 'ANI90']
                    if iden <= 0.05 :
                        metadata_tab.loc[acc, 'ANI95'] = metadata_tab.loc[tgt, 'ANI95']
                        if iden <= 0.02 :
                            metadata_tab.loc[acc, 'ANI98'] = metadata_tab.loc[tgt, 'ANI98']
                            if iden <= 0.01:
                                metadata_tab.loc[acc, 'ANI99'] = metadata_tab.loc[tgt, 'ANI99']

    metadata_tab.to_feather(os.path.join(dbname, f'{module}.db'))

    for i in ids :
        os.unlink(os.path.join(dbname, 'sketches', f'{module}.{i}.list'))
        os.unlink(os.path.join(dbname, 'sketches', f'{module}.{i}.sketch'))
        os.unlink(os.path.join(dbname, 'sketches', f'{module}.{i}.sketch.txt'))
        os.unlink(os.path.join(dbname, 'sketches', f'{module}.{i}.sketch.dat'))
    os.unlink(os.path.join(dbname, 'sketches', f'{module}.dist'))

    with open(os.path.join(dbname, f'{module}.files'), 'wt') as fout :
        for fname in metadata_tab.loc[metadata_tab['status'] == 'DOWNLOADED', 'genome_path'] :
            fout.write(fname + '\n')

    logging.info(f'Run bindash to generate final sketches')
    cmd = f'{bindash} sketch --nthreads={threads} --listfname={module}.files --kmerlen=21 --sketchsize64=160 --outfname={module}.sketch'
    subprocess.Popen(cmd.split(), cwd=dbname).communicate()

    return metadata_tab


def run_bindash(data, ids, data_ids, dbname, module, threads) :
    ref_dir = os.path.join(dbname, 'sketches')
    id = ids[-1]
    with open(os.path.join(ref_dir, f'{module}.{id}.list'), 'wt') as fout :
        for fname in data['genome_path'] :
            fout.write(fname + '\n')

    cmd = f'{bindash} sketch --nthreads={threads} --listfname={module}.{id}.list --kmerlen=21 --sketchsize64=100 --outfname={module}.{id}.sketch'
    subprocess.Popen(cmd.split(), cwd=ref_dir).communicate()

    links = {}
    for jd in ids :
        cmd = f'{bindash} dist --nthreads={threads} --mthres=0.2 --outfname={module}.dist {module}.{id}.sketch {module}.{jd}.sketch'
        subprocess.Popen(cmd.split(), cwd=ref_dir).communicate()

        try :
            for q, r, d in pd.read_csv(os.path.join(ref_dir, f'{module}.dist'), sep='\t', usecols=[0, 1, 2]).values :
                qi, ri = data_ids.get(q, 999999999), data_ids.get(r, 999999999)
                if qi > ri :
                    if qi not in links or links[qi][1] > d :
                        links[qi] = [ri, d]
        except :
            pass
    
    results = [ [i] + links.get(i, [-1, 1]) for i in data.index]
    return results



def read_queries(queries, genome) :
    genome_files = {acc:gfile for acc, gfile in pd.read_csv(genome, sep=',', header=None, dtype=str).values}
    q = []
    for qry in queries :
        query = pd.read_csv(qry, sep='\t', dtype=str, na_filter=None, skiprows=1)
        query.columns = ['accession'] + query.columns[1:].tolist()
        query['genome_path'] = [genome_files.get(acc, '') for acc in query['accession']]
        query = query.loc[query['genome_path'] != '']
        
        if 'excluded_from_refseq' in query.columns :
            problems = ('partial', 'contaminated', 'genome length too', 'sequence duplications', 'chimeric', 'hybrid', 'partial', 'mixed culture', 'completeness check')
            idx = [ not any([p in clause for p in problems]) for clause in query['excluded_from_refseq'] ]
            query = query.loc[idx]
        
        query['score'] = 0
        if 'relation_to_type_material' in query.columns :
            query.loc[query['relation_to_type_material'] != 'na', 'score'] += 16
        if 'refseq_category' in query.columns :
            query.loc[query['refseq_category'] == 'reference genome', 'score'] += 8
            query.loc[query['refseq_category'] == 'representative genome', 'score'] += 4
        if 'assembly_level' in query.columns :
            query.loc[query['assembly_level'] == 'Complete Genome', 'score'] += 2
            query.loc[query['assembly_level'] == 'Chromosome', 'score'] += 1
        
        if 'scaffold_count' in query.columns :
            query = query.sort_values(by=['score', 'scaffold_count'], ascending=False)
        else :
            query = query.sort_values(by=['score'], ascending=False)
        
        if 'species_taxid' not in query.columns :
            query['species_taxid'] = ''
        if 'taxonomy' not in query.columns :
            query['taxonomy'] = ''
        if 'organism_name' not in query.columns :
            query['organism_name'] = ''
        q.append(query[['accession', 'genome_path', 'score', 'species_taxid', 'taxonomy', 'organism_name']])
    return pd.concat(q).sort_values(by=['score'], ascending=False).reset_index(drop=True) #.set_index('accession')


def read_taxonomy(dbname) :
    def ite_taxa(desc, n, taxonomy) :
        for c, s in desc.get(n, []) :
            taxonomy[c] = taxonomy[n] + [[c, s]]
            taxonomy = ite_taxa(desc, c, taxonomy)
        return taxonomy
    tags = dict( [line.strip().split('\t') for line in '''class	c
        family	f
        genus	g
        infraclass	n
        infraorder	n
        kingdom	k
        order	o
        phylum	p
        species	s
        subclass	sc
        subfamily	sf
        subgenus	sg
        subkingdom	sk
        suborder	so
        subphylum	sp
        subspecies	ss
        strain	t
        superclass	uc
        superfamily	uf
        superkingdom	d
        superorder	uo
        superphylum	up'''.split('\n')])
    tax_db = os.path.join(dbname, 'taxonomy')
    name_tax = os.path.join(tax_db, 'names.dmp')
    node_tax = os.path.join(tax_db, 'nodes.dmp')
    nodes = pd.read_csv(node_tax, sep='\t', header=None, dtype=str, usecols=[0, 2, 4]).values[1:]

    desc = {}
    for n, p, s in nodes :
        if p not in desc :
            desc[p] = [[n, tags.get(s, 'n')]]
        else :
            desc[p].append([n, tags.get(s, 'n')])
    taxonomy = ite_taxa(desc, '1', {'1':[]})

    names = pd.read_csv(name_tax, sep='\t', header=None, dtype=str, usecols=[0, 2, 6]).values
    nn = names[names.T[2] == 'scientific name', :2].tolist() + names[names.T[2] == 'synonym', :2].tolist()
    names = {}
    for t, n in nn :
        if t not in names :
            names[t] = n

    output = {}
    for k, t in taxonomy.items() :
        if k in names :
            # x = names[k]
            if k in output : # and output[k].find('uk__Eukaryota') < 0 :
                continue
            output[k] = ';'.join([ '{0}__{1}'.format(s, names[n]).replace(' ', '_') for n, s in t if n in names ])
    return output



def prepare_taxa(dbname, genbank, gtdb, metadata_tab) :
    gtdb_taxa = {}
    ncbi_taxa = {}
    if gtdb :
        gtdb_taxa = {acc:tax for acc, tax in pd.read_csv(gtdb, sep='\t', header=None).values}
    if genbank :
        taxa_dir = os.path.join(dbname, 'taxonomy')
        makedirs(taxa_dir)
        subprocess.Popen('tar -vxzf {0}'.format(os.path.abspath(genbank)).split(), cwd=taxa_dir, stdout=subprocess.PIPE).communicate()
        ncbi_taxa = read_taxonomy(dbname)
        shutil.rmtree(taxa_dir)

    def try_gtdb(dat) :
        if dat['taxonomy'] :
            return dat['taxonomy']
        acc = dat['accession']
        if acc in gtdb_taxa :
            return gtdb_taxa[acc]
        elif 'RS_GCF' + acc[3:] in gtdb_taxa :
            return gtdb_taxa['RS_GCF' + acc[3:]]
        elif 'GB_GCA' + acc[3:] in gtdb_taxa :
            return gtdb_taxa['GB_GCA' + acc[3:]]
        return ncbi_taxa.get(dat['species_taxid'], dat['organism_name'])
    
    metadata_tab['taxonomy'] = [try_gtdb(dat).replace(' ', '_') for acc, dat in metadata_tab.iterrows()]
    return metadata_tab



@click.command()
@click.option('-d', '--dbname',  help='name of the database to build', required=True)
@click.option('-m', '--module',  help='name of the module to build', required=True)
@click.option('-c', '--cutoff',  help='cutoff for adding records into database [default: 0.99]', default=0.99, type=float)
@click.option('-t', '--threads', help='number of threads to use [Default:80]', type=int, default=80)
@click.option('-T', '--gtdb',    help='taxa in format of gtdb (https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz)', default=None)
@click.option('-G', '--genbank', help='https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz', default=None)
@click.option('-g', '--genome',  help='genome list in CSV format. accession,genome_file', required=True)
@click.argument('queries', nargs=-1)
def build(dbname, module, cutoff, threads, genome, queries, genbank, gtdb) :
    '''only one field is required in queries: accession, which needs to be matched to file specified in --genome.
    field "organism_name" is preferred for taxonomic designations'''
    makedirs(dbname)
    dbname = os.path.abspath(dbname)
    module_dir = os.path.join(dbname, module)
    makedirs(module_dir)
    with open(os.path.join(module_dir, 'build_command'), 'wt') as fout :
        fout.write(' '.join(sys.argv)+'\n')

    metadata_tab = read_queries(queries, genome)
    logging.info('Load {0} queries'.format(metadata_tab.shape[0]))
    metadata_tab = prepare_taxa(module_dir, genbank, gtdb, metadata_tab)
    logging.info('Added taxonomy information')
    metadata_tab.to_feather(os.path.join(module_dir, f'{module}.db'))
    get_sketches(metadata_tab, module_dir, module, cutoff, threads)
    shutil.rmtree(os.path.join(module_dir, 'sketch'))


if __name__ == '__main__' :
    build()
