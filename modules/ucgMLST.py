import click, os, gzip, shutil
from build_ucgDB import get_uscgs
from genoQuery import genoQuery



@click.command()
@click.option('-q', '--query', help='fastq file(s), specify --query multiple times for additional reads', required=True)
@click.option('-r', '--reference', help='uscg references. default ucgMLST/db/uscgs', default='/titan/softwares/ucgMLST/db/uscgs/')
@click.option('-d', '--dbname', help='name of the databases [default: /titan/databases/ncbi_20251109/]', default='/titan/databases/ncbi_20251109/')
@click.option('-m', '--modules', help='name of the modules [default: bacteria,archaea,eukaryota]', default='bacteria,archaea,eukaryota')
@click.option('-R', '--representative', help='name of the representatives [default: ANI99]', default='ANI99')
@click.option('-o', '--outdir', help='folder name storing the output', required=True)
@click.option('-t', '--num_threads', help='number of threads [Default: 16]', default=16, type=int)
@click.option('-I', '--min_iden', help='minimum identity of predicted universal core genes [Default: 40]', default=40, type=float)
@click.option('-M', '--min_gene', help='minimum ratio of universal core genes [Default: 0.05]', default=0.05, type=float)
@click.option('-g', '--formal_genus', help='only accept formal genus designations [Default: False]', default=False, is_flag = True)
@click.option('-s', '--formal_species', help='only accept formal species designations [Default: False]', default=False, is_flag = True)
def query_asm(query, reference, dbname, modules, representative, outdir, num_threads, formal_genus, formal_species, min_iden, min_gene) :
    # Create output directory if it doesn't exist
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
   
    # Step 1: Prepare query file - check if gzipped, uncompress if needed
    query_fasta = os.path.join(outdir, 'query.fasta')
    if query.lower().endswith('.gz'):
        print(f"Decompressing {query}...")
        with gzip.open(query, 'rb') as f_in, open(query_fasta, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        print(f"Copying {query} to output directory...")
        shutil.copy(query, query_fasta)
       
    # get_uscgs
    uscg_ffn, n_shared, n_specific = get_uscgs(query_fasta, reference, outdir, outdir, '', all_hits=True)
    # search closest references
    max_counts = 100000000
    genoQuery([uscg_ffn], dbname, modules, None, representative, outdir, 'asm20', formal_genus, formal_species, 0.5, max_counts, 0.02, num_threads, 1, 3, 0.5, 3)
    print()

if __name__ == '__main__' :
    query_asm()
