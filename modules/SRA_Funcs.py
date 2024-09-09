import os, pandas as pd, subprocess, re, numpy as np, gzip, _collections

from configure import executables, logging, rc
import json


def parse_species(taxa) :
    species = _collections.defaultdict(float)
    s2t = {}
    for tax in taxa :
        for t, n in tax :
            s = re.split('s__', t)[-1]
            w = (4. if 'g__' in t else 2.) if 's__' in t else 1.
            species[s] += n*w
            if s not in s2t :
                s2t[s] = {}
            if t not in s2t[s] :
                s2t[s][t] = 0
            s2t[s][t] += n*w

    species = sorted(species.items(), key=lambda s:[-s[1], s[0]])
    taxonomy, _ = max(s2t[species[0][0]].items(), key=lambda x:x[1])
    spec = [s[0] for s in species if s[1] >= 0.5 * species[0][1]]
    spec[1:] = sorted(spec[1:])
    return spec, taxonomy


def get_taxonomy(metadata, acc) :
    top_taxa = [(metadata.loc[acc].taxonomy, 1)]
    ani99_taxa = [[t, np.power(n, 1.   )] for t, n in zip(*np.unique(metadata.loc[metadata.ANI99 == metadata.loc[acc].ANI99].taxonomy, return_counts=True))]
    ani98_taxa = [[t, np.power(n, 2./3.)] for t, n in zip(*np.unique(metadata.loc[metadata.ANI98 == metadata.loc[acc].ANI98].taxonomy, return_counts=True))]
    ani95_taxa = [[t, np.power(n, 1./3.)] for t, n in zip(*np.unique(metadata.loc[metadata.ANI95 == metadata.loc[acc].ANI95].taxonomy, return_counts=True))]
    species, taxonomy = parse_species([top_taxa, ani99_taxa, ani98_taxa, ani95_taxa])
    species[0] += ' [TAX_{0}]'.format(metadata.loc[acc].ANI95[4:])
    return species, taxonomy


def generate_outputs(paf_files, read_files, metadata, matches, reads, r_ids, tmpdir, genome_info, uscg_info, n_reads) :
    uscg_list = [g for g, _ in sorted(uscg_info.items(), key=lambda g:g[1])]
    uscgs = {}
    for mat, _ in matches :
        for g, s in genome_info[mat] :
            if uscg_list[g] not in uscgs :
                uscgs[uscg_list[g]] = [mat, s]

    read_matches = { r_ids[i]:matches[r[0]][0] for i, r in enumerate(reads) if r[0] >= 0 }
    results = [[] for match in matches]
    species_res = {}
    for i, match in enumerate(matches) :
        match_reads = reads[reads.T[0] == i]
        n_diffs = np.sum(match_reads.T[3])/100.
        
        species, taxonomy = get_taxonomy(metadata, match[0])

        results[i] = [float(match_reads.shape[0])*1000./sum([ s for g, s in genome_info[match[0]]])*1000000./n_reads,
                      match_reads.shape[0],
                      int(100000 - 1000.*n_diffs/match_reads.shape[0]+0.5)/1000.,
                      species[0] + ('' if len(species) == 1 else ' ({0})'.format(','.join(species[1:]))),
                      match[0], taxonomy]
        if species[0] not in species_res :
            species_res[species[0]] = results[i][:3] + [species[0], [match[0]], taxonomy]
        else :
            species_res[species[0]][0] += results[i][0]
            species_res[species[0]][1] += results[i][1]
            species_res[species[0]][2] = (species_res[species[0]][2]*species_res[species[0]][1] + results[i][2]*results[i][1])/(species_res[species[0]][1]+results[i][1])
            species_res[species[0]][4].append(match[0])
    results = [ r for r in results if r[0] > 0 ]
    species_res = { s:r for s, r in species_res.items() if r[0] > 0 }

    outputs = {'profile':sorted(species_res.values(), reverse=True), 'OTU':sorted(results, reverse=True)}
    json.dump(outputs, open(f'{tmpdir}/profile.json', 'wt'))
    
    read_details = {}
    read_id = -1
    for fn in read_files :
        if fn.lower().endswith('gz') :
            p = subprocess.Popen(f"{executables['pigz']} -cd {fn}".split(), cwd=tmpdir, stdout=subprocess.PIPE, universal_newlines=True)
        else :
            p = subprocess.Popen(f'cat {fn}'.split(), cwd=tmpdir, stdout=subprocess.PIPE, universal_newlines=True)

        if fn.lower().endswith('q.gz') or fn.lower().endswith('q') :
            for i, line in enumerate(p.stdout) :
                if i % 4 == 0 :
                    read_id += 1
                    rn = line[1:].strip().split()[0]
                    if read_id in read_matches :
                        ref = read_matches[read_id]
                        read_details[read_id] = [f'{rn}']
                elif read_id in read_matches :
                    if i % 4 in (1, 3) :
                        read_details[read_id].append(line.strip())
        else :
            for line in p.stdout :
                if line.startswith('>') :
                    read_id += 1
                    rn = line[1:].strip().split()[0]
                    if read_id in read_matches :
                        ref = matches[read_matches[read_id]][0]
                        read_details[read_id] = [f'{rn}', [], '']
                elif read_id in read_matches :
                    read_details[read_id][1].append(line.strip())
        p.communicate()

    with gzip.open(f'{tmpdir}/primary.sam.gz', 'wt') as pout :
        pout.write('@HD\tVN:1.6\tSO:unsorted\tGO:query\n')
        for g, (r, s) in sorted(uscgs.items(), key=lambda g:[g[1], g[0]]) :
            pout.write(f'@SQ\tSN:{g}__{r}\tLN:{s}\n')
        for fn_id, fname in enumerate(paf_files) :
            p = subprocess.Popen(f"{executables['pigz']} -cd {fname}".split(), cwd=tmpdir, stdout=subprocess.PIPE,
                                 universal_newlines=True)
            for i, line in enumerate(p.stdout) :
                p = line.strip().split('\t')
                r_id = int(p[0], 16)
                if (ref := read_matches.get(r_id, '')) != uscgs.get(p[5], ['-'])[0] :
                    continue
                cigar = ['', p[-1][5:], '']
                rn, rs, rq = read_details.get(r_id, [str(r_id), '', ''])
                if p[4] == '+' :
                    flag = '0'
                    if p[2] != '0' :
                        cigar[0] = f'{p[2]}S'
                    if p[1] != p[3] :
                        cigar[2] = '{0}S'.format(int(p[1]) - int(p[3]))
                else :
                    flag = '16'
                    rs, rq = rc(rs), rq[::-1]
                    if p[2] != '0' :
                        cigar[2] = f'{p[2]}S'
                    if p[1] != p[3] :
                        cigar[0] = '{0}S'.format(int(p[1]) - int(p[3]))
                        
                res = [rn, flag, f'{p[5]}__{ref}', str(int(p[7])+1), p[11], ''.join(cigar), '*', '0', '0', rs, rq] + p[12:-1]
                pout.write('\t'.join(res)+'\n')

            try :
                os.unlink(os.path.join(tmpdir, fname))
            except :
                pass
    subprocess.Popen(f"{executables['pigz']} -cd {tmpdir}/primary.sam.gz | {executables['samtools']} sort -m 4G -@ 8 -O bam -l 0 -T {tmpdir}/tmp - > {tmpdir}/primary.bam", 
                     shell=True).communicate()

    os.unlink(f'{tmpdir}/primary.sam.gz')
    return outputs, f'{tmpdir}/primary.bam'



def write_seq(output, outputs, bam, min_depth=3, min_consensus=0.8) :
    sequences = {}
    prev = 0
    p = subprocess.Popen(f"{executables['samtools']} mpileup -AB {bam}".split(), universal_newlines=True, stdout=subprocess.PIPE)
    for line in p.stdout :
        p = line.strip().split('\t')
        p[1] = int(p[1])
        if p[0] not in sequences :
            sequences[p[0]] = []
            prev = p[1]
            dist = 0
        else :
            dist = p[1] - prev - 1
            prev = p[1]
        s = re.sub(r'\^.', '', p[4]).upper()
        s = re.sub(r'\$', '', s)
        s = [re.findall('^(\d*)(.+)$', ss)[0] for ss in re.split(r'[+-]', s)]
        s = list(''.join([ s2[int(s1):] if s1 else s2 for s1, s2 in s ]))
        base, cdp = sorted(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])[0]
        depth = len(s)
        if base in ('*', 'N') :
            base = ''
        elif cdp < min_depth or float(cdp) < min_consensus * float(depth) :
            base = base.lower()
        sequences[p[0]].append('n'*dist + base)
    seq_out = {}
    for n, s in sequences.items() :
        t = re.split(r'__', n)[1]
        if t not in seq_out :
            seq_out[t] = {}
        seq_out[t][n] = ''.join(s).replace('*', '').replace('+', '').replace('-', '')
    for otu in outputs.get('OTU', []) :
        ref = otu[4]
        otu.append(seq_out.get(ref, ''))
    for profile in outputs.get('profile', []) :
        refs = profile[4]
        s = {}
        for ref in refs :
            s.update(seq_out.get(ref, ''))
        profile.append(s)
    json.dump(outputs, open(f'{output}/profile.json', 'wt'))
    return outputs
