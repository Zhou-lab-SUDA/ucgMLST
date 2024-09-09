import click
import genoPhylo



@click.command()
@click.option('-d', '--dbname', help='absolute path of the database. [required]', default="/titan/databases/ncbi_20240609z")
@click.option('-m', '--module', help='modules in the database. [default: auto detect]', default='auto')
@click.option('-r', '--reference', help='reference genome, can be one of accession_code, tax_code, or species_name. [required]', required=True)
@click.option('--min_identity', help='minumum identity of a sequence comparing to ref_acc. default: 0.93 [0. - 1.]', default=None, type=float)
@click.option('-P', '--min_presence_ref', help='minumum coverages of genes for a refernce genome to be included. [default: 0.75]', default=None, type=float)
@click.option('-G', '--genus', help='Include genomes from the same genus, will also set min_presence_ref=0.45 and min_identity=0.83 by default.', default=False, is_flag=True)
@click.option('-N', '--no_main', help='do not include genomes in the main database', default=False, is_flag=True)
@click.option('-g', '--genome_list', help='additional genomes to be included in the analysis', default='')
@click.option('-o', '--outdir', help='folder storing the prepared database', required=True)
@click.option('-n', '--n_proc', help='number of processes [default: 8]', default=8, type=int)
def main(dbname, module, reference, min_identity, min_presence_ref, genus, no_main, genome_list, outdir, n_proc) :
    genoPhylo.genoPhylo([], dbname, module, reference, min_identity, 0.2, min_presence_ref, genus, no_main, genome_list, None, outdir, n_proc, True, '')

if __name__ == '__main__' :
    main()
