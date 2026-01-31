import os, click, numpy as np, pandas as pd
import subprocess, re, gzip
from multiprocessing import Pool
import json

from configure import executables, logging

def write_seq(output, outputs, bam, min_depth=3, min_consensus=0.65):
    """Generate consensus sequences for matched genes from BAM alignment"""
    sequences = {}
    prev = 0
    p = subprocess.Popen(f"{executables['samtools']} mpileup -AB {bam}".split(), 
                        universal_newlines=True, stdout=subprocess.PIPE)
    
    for line in p.stdout:
        parts = line.strip().split('\t')
        parts[1] = int(parts[1])
        
        # Initialize new gene sequence
        if parts[0] not in sequences:
            sequences[parts[0]] = []
            prev = parts[1]
            dist = 0
        else:
            dist = parts[1] - prev - 1
            prev = parts[1]
        
        # Parse pileup format
        s = re.sub(r'\^.', '', parts[4]).upper()
        s = re.sub(r'\$', '', s)
        s = [re.findall(r'^(\d*)(.+)$', ss)[0] for ss in re.split(r'[+-]', s)]
        s = list(''.join([s2[int(s1):] if s1 else s2 for s1, s2 in s]))
        
        # Get consensus base
        base, cdp = sorted(zip(*np.unique(s, return_counts=True)), key=lambda x: -x[1])[0]
        depth = len(s)
        
        if base in ('*', 'N'):
            base = ''
        elif depth < min_depth or float(cdp) < min_consensus * float(depth):
            base = base.lower()
        
        sequences[parts[0]].append('n' * dist + base)
    
    # Organize sequences by gene
    seq_out = {}
    for name, seq in sequences.items():
        # Parse gene name from reference (format: gene__geneid)
        gene_id = re.split(r'__', name)[1] if '__' in name else name
        if gene_id not in seq_out:
            seq_out[gene_id] = {}
        seq_out[gene_id][name] = ''.join(seq).replace('*', '').replace('+', '').replace('-', '')
    
    # Add sequences to gene matches
    for gene_match in outputs.get('gene', []):
        gene_id = gene_match[3]  # gene identifier
        gene_match.append(seq_out.get(gene_id, {}))
    
    # Add sequences to function profile
    for func_profile in outputs.get('function', []):
        gene_ids = func_profile[3]  # list of gene identifiers
        s = {}
        for gene_id in gene_ids:
            s.update(seq_out.get(gene_id, {}))
        func_profile.append(s)
    
    json.dump(outputs, open(f'{output}/profile.json', 'wt'), indent=2)
    return outputs


def generate_outputs(paf_files, read_files, genes, matches, reads, r_ids, output, n_reads, left_cuts=None):
    """Generate gene alignment outputs with functional annotations"""
    if left_cuts is None:
        left_cuts = [0] * len(read_files)
    
    gene_dict = { mat:genes[mat] for mat, coverage in matches }
    
    # Map reads to their matched genes
    read_matches = {}
    for r in reads :
        if r[0] >= 0 :
            if r[4] not in read_matches :
                read_matches[r[4]] = [matches[r[0]][0]]
            else :
                read_matches[r[4]].append(matches[r[0]][0])
    
    # Prepare results for each gene match
    results = [[] for _ in matches]
    
    for i, (gene_name, coverage) in enumerate(matches):
        match_reads = reads[reads.T[0] == i]
        n_diffs = np.sum(match_reads.T[3]) / 100.
        gene_id, gene_size, gn2, annotation = gene_dict[gene_name]
        
        # Calculate metrics
        # RPKM-like normalization: (reads * 1000 / gene_length) * 1000000 / total_reads
        abundance = float(match_reads.shape[0]) * 1000. / gene_size * 1000000. / n_reads if gene_size > 0 else 0
        read_count = match_reads.shape[0]
        identity = int(100000 - 1000. * n_diffs / match_reads.shape[0] + 0.5) / 1000. if match_reads.shape[0] > 0 else 0
        
        results[i] = [
            abundance,           # normalized abundance
            read_count,          # number of reads
            identity,            # average identity
            gn2,           # gene identifier
            annotation,          # functional annotation
            coverage             # coverage depth
        ]

    # Filter and format results
    results = [r for r in results if r[0] > 0]
    
    outputs = {
        'gene': sorted(results, reverse=True)
    }
    
    json.dump(outputs, open(f'{output}/profile.json', 'wt'), indent=2)
    
    # Generate detailed read alignment information
    read_details = {}
    read_id = -1
    
    for fn_id, fn in enumerate(read_files):
        if fn.lower().endswith('gz'):
            p = subprocess.Popen(f"{executables['pigz']} -cd {fn}".split(), 
                               cwd=output, stdout=subprocess.PIPE, universal_newlines=True)
        else:
            p = subprocess.Popen(f'cat {fn}'.split(), 
                               cwd=output, stdout=subprocess.PIPE, universal_newlines=True)
        
        if fn.lower().endswith('q.gz') or fn.lower().endswith('q'):
            # FASTQ format
            for i, line in enumerate(p.stdout):
                if i % 4 == 0:
                    read_id += 1
                    rn = line[1:].strip().split()[0]
                    if read_id in read_matches:
                        read_details[read_id] = [f'{rn}']
                elif read_id in read_matches:
                    if i % 4 in (1, 3):
                        read_details[read_id].append(line.strip()[left_cuts[fn_id]:])
        else:
            # FASTA format
            x = 0
            for line in p.stdout:
                if line.startswith('>'):
                    read_id += 1
                    rn = line[1:].strip().split()[0]
                    if read_id in read_matches:
                        read_details[read_id] = [f'{rn}', [], '']
                        x = 0
                elif read_id in read_matches:
                    if x == 0:
                        read_details[read_id][1].append(line.strip()[left_cuts[fn_id]:])
                        x = 1
                    else:
                        read_details[read_id][1].append(line.strip())
            
            # Join sequences for FASTA
            for read_id, read_info in read_details.items():
                read_info[1] = ''.join(read_info[1])
                read_info[2] = 'I' * len(read_info[1])  # Dummy quality
        
        p.communicate()
    
    # Generate SAM file
    with gzip.open(f'{output}/primary.sam.gz', 'wt') as pout:
        pout.write('@HD\tVN:1.6\tSO:unsorted\tGO:query\n')
        
        # Write gene references to SAM header
        for gene_name, (gene_id, size, gn2, annotation) in sorted(gene_dict.items(), key=lambda x: x[1][0]):
            pout.write(f'@SQ\tSN:{gn2}\tLN:{size}\n')
        
        # Process PAF files and write alignments
        for fn_id, fname in enumerate(paf_files):
            p = subprocess.Popen(f"{executables['pigz']} -cd {fname}".split(), 
                               cwd=output, stdout=subprocess.PIPE, universal_newlines=True)
            
            for i, line in enumerate(p.stdout):
                parts = line.strip().split('\t')
                r_id = int(parts[0], 16)
                
                refs = read_matches.get(r_id, [])
                if len(refs) <= 0 or parts[5] not in refs :
                    continue
                
                cigar = ['', parts[-1][5:], '']
                rn, rs, rq = read_details.get(r_id, [str(r_id), '', ''])
                
                if parts[4] == '+':
                    flag = '0'
                    if parts[2] != '0':
                        cigar[0] = f'{parts[2]}S'
                    if parts[1] != parts[3]:
                        cigar[2] = '{0}S'.format(int(parts[1]) - int(parts[3]))
                else:
                    flag = '16'
                    rs, rq = reverse_complement(rs), rq[::-1]
                    if parts[2] != '0':
                        cigar[2] = f'{parts[2]}S'
                    if parts[1] != parts[3]:
                        cigar[0] = '{0}S'.format(int(parts[1]) - int(parts[3]))
                
                res = [
                    rn, flag, f'{parts[5]}', str(int(parts[7]) + 1), 
                    parts[11], ''.join(cigar), '*', '0', '0', rs, rq
                ] + parts[12:-1]
                pout.write('\t'.join(res) + '\n')
            
            try:
                os.unlink(os.path.join(output, fname))
            except:
                pass
    
    # Convert SAM to sorted BAM
    subprocess.Popen(
        f"{executables['pigz']} -cd {output}/primary.sam.gz | "
        f"{executables['samtools']} sort -m 4G -@ 8 -O bam -l 0 -T {output}/tmp - > {output}/primary.bam",
        shell=True
    ).communicate()
    
    os.unlink(f'{output}/primary.sam.gz')
    
    return outputs, f'{output}/primary.bam'

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


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


def find_cov_outlier(covs, block_size, delta_fold=3) :
    q1, q3 = np.quantile(covs.T[0], 0.25), np.quantile(covs.T[0], 0.75)
    delta_q = max(q3 - q1, 1.5/block_size)
    idx = (covs.T[1]/covs.T[2] <= q3 + delta_fold*delta_q)
    return idx


def get_matches(gene_info, read_maps, allowed_distance=0.005, min_frag_match=2, delta_fold=3, block_size=300) :
    '''gfrag_id, gene_id, read_id/rfrag_id, start/read_id, end/size, mutation, diff
        0          1       2       3                4       5        6   '''
    allowed_distance *= 10000.
    read_fragments, read_frag_idx = np.unique(read_maps[:, 1:3], axis=0, return_inverse=True)
    fragment2read = {frag_id:read for frag_id, (x, read) in enumerate(read_fragments) }
    read_maps[:, 2] = read_frag_idx
    read_maps[:, 1] = read_maps[:, 0]
    pos_genes = set(read_maps[:, 0])
    genes = { idx:s for g, (idx, s, name, metadata) in gene_info.items() if idx in pos_genes }
    
    read_maps = read_maps[pd.Series(read_maps.T[0]).isin(genes)]
    if read_maps.shape[0] <= 0 :
        return [], []
    
    d = np.array(sorted(genes.items()))
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
    
    gene2 = {g:np.array([[i, fragments[i], g] for i in range(*uscg_frag[g])], dtype=int) for g, s in genes.items()}
    max_frag = np.max([np.max(gn.T[0]) for gn in gene2.values()])+1
    read_maps = uscg2frag(read_maps, uscg_frag, block_size)
    read_maps.T[3] = np.vectorize(fragment2read.get)(read_maps.T[2])
    match_results = []
    summed_reads = np.zeros([int(np.max(read_maps.T[2]) + 1), 5], dtype=np.int32)
    summed_reads[:, 1:].fill(9999999)
    summed_reads[:, 0].fill(-1)
    coverages = [ [-1, -1, -1, gene, []] for gene in gene2.keys() ]
    
    while len(coverages) > 0 :
        frag_cov = np.bincount(read_maps.T[0], weights=np.power(0.33333333, read_maps.T[6].astype(np.float64) / 100.), minlength=max_frag)
        
        max_i = -1
        for i, (depth, n_gene, n_frag, gene, g_cov) in enumerate(coverages) :
            if depth == -1 or max_i < 0 or depth >= coverages[max_i][0] :
                fragments = gene2[gene]
            else :
                break
            
            covs = np.array([ [((frag_cov[f]+.5)/(s+.5)), frag_cov[f], s, g] \
                                  if f < frag_cov.size else [0., 0, s, g] for f, s, g in fragments ])
            idx = find_cov_outlier(covs, block_size, delta_fold)
            if len(g_cov) > 0 :
                if np.sum((covs[idx, 1] < 0.1 * g_cov[idx]) | ((covs[idx, 1] < 0.2 * g_cov[idx]) & (covs[idx, 1] < 1))) >= 0.75 * np.sum(g_cov[idx] > 0) :
                    coverages[i] = [0, 0, 0, gene, g_cov]
                    continue
            cov = np.sum(covs[idx, 1]) / np.sum(covs[idx, 2])
            n_frag = covs[(covs[:, 1] >= 1.) & idx, 3].shape[0]
            n_gene = np.unique(covs[idx & (covs[:, 1] >= 1.), 3]).size
            coverages[i] = [cov, n_gene, n_frag, gene, covs.T[1] if len(g_cov) == 0 else g_cov ]
            if (max_i < 0 or cov > coverages[max_i][0]) and \
                (n_frag >= min_frag_match or n_frag*3 >= len(gene2[gene])) :
                max_i = i
        if max_i < 0 :
            break
        (depth, n_gene, n_frag, match, g_cov) = coverages[max_i]
        logging.info(f'    Ref: {match} with {n_gene} USCGs. ')
        coverages = [c for c in sorted(coverages, reverse=True) if c[0] > 0]
        
        matches = read_maps[pd.Series(read_maps.T[1]) == match]
        match_covs = dict(zip(*np.unique(matches[:, 0], return_counts=True)))
        
        covs = np.array([[(match_covs.get(g, 0) + 0.5)/(s + 0.5), match_covs.get(g, 0), s, og] for g, s, og in gene2[match]])
        idx = find_cov_outlier(covs, block_size, delta_fold)

        # report hits
        match_results.append([match, np.sum(covs[idx].T[1])/np.sum(covs[idx].T[2])])

        frag1, frag2 = set(gene2[match][idx, 0]), set(gene2[match][~idx, 0])
        reads_ignored = set(matches[pd.Series(matches.T[0]).isin(frag2), 2]) - set(matches[pd.Series(matches.T[0]).isin(frag1), 2])
        
        matched_reads = np.zeros([summed_reads.shape[0], 4], dtype=np.int32)
        matched_reads[:] = 9999999
        
        _, ridx = np.unique(matches.T[2], return_index=True)
        matched_reads[matches[ridx, 2], :] = matches[ridx, :][:, (0,1,5,3)]
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
        if x.shape[0] >= 1 and \
            (n_frag >= min_frag_match or n_frag*3 >= len(gene2[ref])) :
                res.append(match_results[i0])
                if i0 != i1 :
                    summed_reads[summed_reads.T[0] == i0, 0] = i1
                i1 += 1
        else :
            summed_reads[summed_reads.T[0] == i0, :] = [-1, 9999999, 9999999, 9999999]

    gene_map = {info[0]:g for g, info in gene_info.items()}
    res = [[gene_map[r[0]], r[1]] for r in res]
    return np.array(res, dtype=object), summed_reads



def parse_paf(data) :
    outfile, tmpdir, gene_info, allowed_distance = data
    rmaps = []
    p = subprocess.Popen(
        f"{executables['pigz']} -cd {outfile}".split(), cwd=tmpdir, 
        stdout=subprocess.PIPE, universal_newlines=True, bufsize=1024*1024  # 1MB buffer for faster I/O
    )

    for line in p.stdout:
        p = line.split('\t', 11)
        if p[5] not in gene_info:
            continue
        p[0] = int(p[0], 16)
        p[1:4]  = [ int(pp) for pp in p[1:4]  ]
        p[6:11] = [ int(pp) for pp in p[6:11] ]
        p[5] = gene_info[p[5]][0]
        if p[4] == '+' :
            s, e = min(p[2], p[7]), min(p[1]-p[3], p[6] - p[8])
        else :
            s, e = min(p[1]-p[3], p[7]), min(p[2], p[6] - p[8])
        mut = (p[10]-p[9])*10 + s + e
        rmaps.append([p[5], int(p[2]/100), p[0], p[7], p[8], int(mut*1000/(p[10]+s+e)+0.5)])

    if len(rmaps):
        rmaps = np.array(rmaps, dtype=np.uint32)
        rmaps = rmaps[rmaps.T[5] <= allowed_distance * 10000]
        np.savez_compressed( os.path.join(tmpdir, f'{outfile}.npz'), reads=rmaps)
        return os.path.join(tmpdir, f'{outfile}.npz')
    return ''


def map_to_genes(paf_files, genes, tmpdir, allowed_distance, pool) :
    rmaps = []
    for rmap_file in pool.imap_unordered(parse_paf, [ [sfile, tmpdir, genes, allowed_distance] for sfile in paf_files ]) :
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
        # rmaps.T[1] = rmaps.T[0]
        rmaps = np.hstack([rmaps, (rmaps.T[5] - read_dist[rmaps.T[2]]).reshape([-1, 1])]).astype(np.int32)
    return rmaps, np.array([r[0] for r in sorted(read_rename.items(), key=lambda r:r[1])])


def read_filter(qry, tmpdir, total_reads) :
    left_cut = 0
    if qry.lower().endswith('q') or qry.lower().endswith('q.gz') :
        if qry.lower().endswith('q.gz') :
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True).stdout
        elif qry.lower().endswith('q') :
            fin = open(qry, 'rt')
        reads = []
        for id, line in enumerate(fin) :
            if id > 2000000 :
                break
            if id % 40 == 1 :
                reads.append(list(line.strip())[:12])
        fin.close()
        reads = np.array(reads).T
        for ix in range(max(reads.shape[0], 12)) :
            site = reads[ix]
            _, cnt = np.unique(site, return_counts=True)
            max_cnt = np.max(cnt)/site.size
            if max_cnt >= 0.8 :
                left_cut = ix + 1

        if left_cut > 0 : 
            logging.info(f'Trimmed {left_cut} bases from the beginning of reads in {os.path.basename(qry)} based on base composition bias.')
            
        if qry.lower().endswith('q.gz') :
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True).stdout
        elif qry.lower().endswith('q') :
            fin = open(qry, 'rt')
        qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fastq.gz'))
        with open(qry_file, 'wb') as fout2 :
            fout = subprocess.Popen('{pigz} -c'.format(**executables).split(), stdin=subprocess.PIPE, stdout=fout2, universal_newlines=True)
            for id, line in enumerate(fin) :
                if id % 4 == 0 :
                    fout.stdin.write(f'@{total_reads:X}\n')
                    total_reads += 1
                elif id % 4 == 2 :
                    fout.stdin.write(line)
                else :
                    fout.stdin.write(line[left_cut:])
            fout.communicate()
        fin.close()
    else :
        if qry.lower().endswith('.gz') :
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True).stdout
        else :
            fin = open(qry, 'rt')
        reads = []
        n_read = 0
        for line in fin :
            if line.startswith('>') :
                n_read += 1
                if n_read > 500000 :
                    break
                elif n_read % 10 == 1 :
                    reads.append([])
            elif n_read % 10 == 1 :
                if len(reads[-1]) < 12 :
                    reads[-1].extend(list(line.strip()))
                    reads[-1] = reads[-1][:12]
        fin.close()
        reads = np.array(reads).T
        for ix in range(max(reads.shape[0], 12)) :
            site = reads[ix]
            _, cnt = np.unique(site, return_counts=True)
            max_cnt = np.max(cnt)/site.size
            if max_cnt >= 0.8 :
                left_cut = ix + 1

        if left_cut > 0 : 
            logging.info(f'Trimmed {left_cut} bases from the beginning of reads in {os.path.basename(qry)} based on base composition bias.')
            
        if qry.lower().endswith('.gz') :
            fin = subprocess.Popen('{pigz} -cd {0}'.format(qry, **executables).split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True).stdout
        else :
            fin = open(qry, 'rt')
        qry_file = os.path.abspath(os.path.join(tmpdir, 'r.fasta.gz'))
        with open(qry_file, 'wb') as fout2 :
            x = 0
            fout = subprocess.Popen('{pigz} -c'.format(**executables).split(), stdin=subprocess.PIPE, stdout=fout2, universal_newlines=True)
            for line in fin :
                if line.startswith('>') :
                    fout.stdin.write(f'>{total_reads:X}\n')
                    total_reads += 1
                    x = 0
                elif x == 0 :
                    fout.stdin.write(line[left_cut:])
                    x = 1
                else :
                    fout.stdin.write(line)
            fout.communicate()
        fin.close()
    return qry_file, left_cut, total_reads



def map_reads(query, dbname, mode, tmpdir, max_dist, num_threads) :
    outputs = []
    total_reads = 0
    left_cuts = [0 for qry in query]
    for qid, qry in enumerate(query):
        qry_file, left_cuts[qid], total_reads = read_filter(qry, tmpdir, total_reads)
        outfile = f'{qid}.0.paf.gz'
        p_dist = 0.6 if mode == 'sr' else 0.6

        subprocess.Popen(
            '{minimap2} -t{3} -cx {5} -T20 --frag=yes -p{6} -N90000 -Y --end-bonus 12 -2 --secondary=yes {0} {1}|{EnFlt} {4}|{pigz} -c > {2}'.format(
                dbname, qry_file, outfile, num_threads, max_dist, mode, p_dist, **executables,
            ), cwd=tmpdir, shell=True).communicate()
            
        outputs.append(outfile)
        os.unlink(qry_file)

    return outputs, total_reads, np.array(left_cuts, dtype=np.int32)


def query_sra(query, dbname, genes, output, mode, max_dist, min_depth, min_frag_match, min_consensus, cover_fold, pool, debug=[False, False, False, False]) :
    if not debug[0] :
        logging.info('Running read mapping...')
        paf_files, n_reads, left_cuts = map_reads(query, dbname, mode, output, max_dist, len(pool._pool))
        # np.savez_compressed(os.path.join(output, 'uscg.npz'), n_reads=np.array(n_reads), left_cuts=left_cuts)
        logging.info('Done')
    else :
        paf_files = [os.path.abspath(os.path.join(output, f'{id}.{jd}.paf.gz')) for id, _ in enumerate(query) for jd, _ in enumerate(dbname)]
        # n_reads = int(np.load(os.path.join(output, 'uscg.npz'))['n_reads'])
        # left_cuts = np.load(os.path.join(output, 'uscg.npz'))['left_cuts']
 
    if not debug[1] : 
        logging.info('Extracting USCG information...')
        read_maps, r_ids = map_to_genes(paf_files, genes, output, max_dist, pool)
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
        matches, reads = get_matches(genes, read_maps, allowed_distance=0.005, min_frag_match=min_frag_match, delta_fold=cover_fold)
        logging.info('Done')
        if len(matches) == 0 :
            return {'profile':[], 'OTU':[]}
        # np.savez_compressed(os.path.join(output, 'out1.npz'), matches=matches, reads=reads)
    else :
        data = np.load(os.path.join(output, 'out1.npz'), allow_pickle=True)
        matches, reads = data['matches'], data['reads']
        
    if not debug[3] :
        logging.info('Preparing outputs...')
        outputs, bam = generate_outputs(paf_files, query, genes, matches, reads, r_ids, output, n_reads, left_cuts)
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


def read_gene(dbname) :
    genes = {}
    with gzip.open(dbname, 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                parts = line[1:].strip().split()
                if parts[0].startswith('VFG') :
                    n = parts[0]
                    annotations = ' '.join(parts[1:]) if len(parts) > 1 else ''
                else :
                    p = parts[0].split('|')
                    n = p[1]
                    annotations = f'({p[5]}) {p[7]}'
                genes[parts[0]] = [len(genes), 0, n, annotations]
            else :
                genes[parts[0]][1] += len(line.strip())
    return genes
@click.command()
@click.option('-q', '--query', help='fastq file(s), specify --query multiple times for additional reads', required=True, multiple=True)
@click.option('-d', '--dbname', help='name of the databases [default: AMRfinder + VFDB]', default=None)
@click.option('-o', '--outdir', help='folder name storing the output', required=True)
@click.option('-t', '--num_threads', help='number of threads [Default: 16]', default=16, type=int)
@click.option('-M', '--mode', help='One of sr [default], map-ont, map-hifi, map-pb, asm20', default='sr')
@click.option('-D', '--max_dist', help='maximum distance of alignment [Default: 0.10 for map-ont and 0.05 for others]', default=None, type=float)
@click.option('-f', '--coverage_fold_diff', help='allowed coverage fold differene (relative to std) for identifying nonspecific matches. [default: 3]', default=3, type=float)
@click.option('--min_frag_match', help='minimum fragments to call the presence of a gene. [Default: 2]', default=2, type=int)
@click.option('--min_depth', help='minimum read depth to call a base reliably. [Default: 3]', default=None, type=int)
@click.option('--min_consensus', help='minimum proportion of consensus to call a base reliably [Default: 0.8]', default=0.8, type=float)
def main(query, dbname, outdir, mode, max_dist, num_threads, min_depth, min_frag_match, min_consensus, coverage_fold_diff) :
    genoEffectors(query, dbname, outdir, mode, max_dist, num_threads, min_depth, min_frag_match, min_consensus, coverage_fold_diff)

def genoEffectors(query, dbname, outdir, mode, max_dist, num_threads, min_depth, min_frag_match, min_consensus, coverage_fold_diff) :
    if dbname == None :
        dbname = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/db/effectors/effectors.centroid.fas.gz'
    if max_dist == None :
        max_dist = 0.1 if mode in ('map-ont', 'map-pb') else 0.05

    pool = Pool(num_threads)
    np.random.seed(42)
    query = [os.path.abspath(qry) for qry in query]
    dbname = os.path.abspath(dbname)
    genes = read_gene(dbname)
    
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    logging.info('Done')
    query_sra(query, dbname, genes, outdir, mode, max_dist, min_depth, min_frag_match, min_consensus, coverage_fold_diff, pool)



if __name__ == '__main__' :
    main()
