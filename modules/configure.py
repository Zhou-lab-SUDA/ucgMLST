import logging, os, _collections, shutil, hashlib, gzip

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

external_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'externals')
db_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'db')


def check_exe(exe, external_dir='/') :
    p = os.path.join(external_dir, exe)
    return p if os.path.isfile(p) else shutil.which(exe)


executables = dict(
    pigz      = check_exe('pigz') if check_exe('pigz') else check_exe('gzip'), 
    gzip      = check_exe('gzip'), 
    bindash   = check_exe('bindash', external_dir), 
    getorf    = check_exe('getorf', external_dir), #os.path.join(external_dir, 'getorf'),
    diamond   = check_exe('diamond', external_dir), 
    hmmsearch = check_exe('hmmsearch', external_dir), #os.path.join(external_dir, 'hmmsearch'),
    minimap2  = check_exe('minimap2', external_dir), #os.path.join(external_dir, 'minimap2'),
    samtools  = check_exe('samtools', external_dir), #os.path.join(external_dir, 'samtools'),
    EnFlt     = check_exe('_EnFlt.py', external_dir), #os.path.join(external_dir, '_EnFlt.py'),
    iqtree    = check_exe('iqtree', external_dir), #os.path.join(external_dir, 'iqtree'), 
)


comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n'}
def rc(seq) :
    return ''.join([comp.get(s, 'N') for s in list(seq)[::-1]])


def get_md5(value, integers=False):
    r = hashlib.md5(str(value).encode()).hexdigest()[-16:]
    if integers :
        return (int(r[:8], 16), int(r[8:], 16))
    return r


def readFasta(fname, upper=True) :
    seqs = _collections.OrderedDict()
    if fname.lower().endswith('.gz') :
        with gzip.open(fname, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    seqs[name] = []
                else :
                    seqs[name].extend(line.strip().split())
    else :
        with open(fname, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    seqs[name] = []
                else :
                    seqs[name].extend(line.strip().split())

    for n, s in seqs.items() :
        seqs[n] = ''.join(s).upper() if upper else ''.join(s)
    return seqs


if __name__ == '__main__' :
    main()
