#!/usr/bin/env python
import os, ete3, numpy as np, subprocess, click, tempfile, pickle, json, re, _collections
import configure
from sklearn.linear_model import Ridge

executables = configure.executables
rc = configure.rc

def get_sites(aln_fas, nwk) :
    names = list(aln_fas.keys())
    fas = np.array([list(aln_fas[n]) for n in names]).T
    sites, vseqs = [], []
    
    for i, s in enumerate(fas) :
        stype = np.unique(s)
        if stype.size > 2 or (stype.size == 2 and stype[0] != '-') :
            sites.append(i+1)
            vseqs.append(s)
    vseqs = [''.join(s) for s in np.array(vseqs).T]
    
    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir :
        with open(os.path.join(tmpdir, 'aln'), 'wt') as fout :
            for n, s in zip(names, vseqs) :
                fout.write('>{0}\n{1}\n'.format(n, s))
        
        subprocess.Popen('{0} --ancestral -s {1} -te {2} -m GTR+G4 -redo --prefix {3}'.format(configure.executables['iqtree'], os.path.join(tmpdir, 'aln'), nwk, os.path.join(tmpdir, 'aln')).split(),
                         stdout=subprocess.PIPE).communicate()
        tre = ete3.Tree(os.path.join(tmpdir, 'aln.treefile'), format=1)
        nodes = {n.name:[] for n in tre.traverse() if not n.is_leaf() }
        
        with open(os.path.join(tmpdir, 'aln.state'), 'rt') as fin :
            for line in fin :
                if line.startswith('#') or line.startswith('Node\t') :
                    continue
                p = line.strip().split()
                nodes[p[0]].append([ float(v) for v in p[3:] ])
        base_encoding = {'A':[1., 0., 0., 0.], 'C':[0., 1., 0., 0.], 'G':[0., 0., 1., 0.], 'T':[0., 0., 0., 1.]}
        for n, s in zip(names, vseqs) :
            nodes[n] = [base_encoding.get(b, [0., 0., 0., 0.]) for b in s ]
        for n, s in nodes.items() :
            nodes[n] = np.array(s)

        for n in tre.iter_descendants('postorder') :
            if n.is_leaf() :
                mask = (np.sum(nodes[n.name], 1) == 0)
                nodes[n.name][mask] = nodes[n.up.name][mask]
        # mutations = np.zeros(len(sites))
        # for n in tre.iter_descendants('postorder') :
        #     mutations += np.sum(np.abs(nodes[n.name] - nodes[n.up.name]), 1)
        # mut_q1, mut_q3 = np.quantile(mutations, 0.25), np.quantile(mutations, 0.75)
        # idx = (mutations <= mut_q3 + 3.*max(mut_q3-mut_q1, 3))
        # for n in tre.traverse('postorder') :
        #     nodes[n.name] = nodes[n.name][idx]
    return nodes, [[s, []] for s in sites], tre

def map_qry(aln_fas, otu, sites) :
    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir :
        with open(os.path.join(tmpdir, 'ref'), 'wt') as fout :
            n, s = list(aln_fas.items())[0]
            s = s.upper()
            fout.write('>ref\n{1}\n'.format(n, s))
            ref_seq = {'ref':s}
        qry_seq = {}
        with open(os.path.join(tmpdir, 'qry'), 'wt') as fout :
            for n, s in otu[6].items() :
                s = s.upper()
                qry_seq[n] = list(s)
                fout.write('>{0}\n{1}\n'.format(n, s))

        map_cmd = f"{executables['minimap2']} -k13 -w5 -c -t1 --frag=yes --rmq -A1 -B4 -O8,16 -E2,1 -r20k,40k -g10k -P -N5000 -f1000,5000 -n2 -m20 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes ref qry"
        p = subprocess.Popen(map_cmd.split(), cwd=tmpdir, stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.PIPE, )
        for line in p.stdout :
            if line.startswith('[') :
                continue
            p = line.strip().split('\t')
            p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
            p[1:4] = [int(p[1]), int(p[2])+1, int(p[3])]
            p[6:9] = [int(p[6]), int(p[7])+1, int(p[8])]

            if p[4] == '+' :
                qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1 
            else :
                qi, ri, cigar, d = p[3], p[7], p[-1][5:], -1
            xi = 0
            while xi < len(sites) and sites[xi][0] < ri :
                xi += 1
            for s, t in re.findall('(\d+)([MDI])', cigar) :
                s = int(s)
                if t != 'I' :
                    rj = ri + s
                if t != 'D' :
                    qj = qi + s*d
                while xi < len(sites) and sites[xi][0] >= ri and sites[xi][0] < rj :
                    rd = sites[xi][0] - ri
                    rx = ri + (rd if t != 'I' else 0) - 1
                    qx = qi + (rd*d if t != 'D' else 0) - 1
                    rseq = ref_seq[p[5]][rx]
                    qseq = qry_seq[p[0]][qx] if d > 0 else configure.rc(qry_seq[p[0]][qx])
                    if qseq != '-' :
                        rr = sites[xi][1]
                        if len(rr) == 0 or rr[5] < p[10] :
                            sites[xi][1] = [p[0], qx+d, rseq, qseq, p[9], p[10], p[4]]
                    xi += 1
                ri, qi = rj, qj
            qry_seq[p[0]][p[2]-1:p[3]] = ['-']* (p[3] - p[2] + 1)
    return [s[:2] for s in sites]

def parse_bam(bam, ref_acc, sites) :
    base_comp = {}
    if bam :
        qry_sites = {}
        for site, var in sites :
            if len(var) :
                qry_sites[(var[0], var[1])] = ('ref', site)

        p = subprocess.Popen('{0} mpileup -AB -q 0 -Q 0 {1}'.format(executables['samtools'], bam).split(),
                            universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        contigs = {}

        for line in p.stdout :
            p = line.strip().split('\t')
            if re.split('__', p[0])[1] != ref_acc :
                continue
            if p[0] not in contigs :
                contigs[p[0]] = int(p[1]) - 1

            s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
            s = list(''.join([b[int(n):] for n, b in re.findall('[+-](\d+)(.+)', '+0' + s)]))
            base, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
            if base in ('*', 'N') :
                contigs[p[0]] -= 1
                continue
            
            site = int(p[1]) - contigs[p[0]]
            key = (p[0], site)
            if key not in qry_sites :
                continue
            
            bases = dict(zip(*np.unique(s, return_counts=True)))
            base_comp[key] = [int(bases.get(b, 0.)) for b in ('A', 'C', 'G', 'T')]

    for i, (site, v) in enumerate(sites) :
        if len(v) :
            x = base_comp.get((v[0], v[1]), [0, 0, 0, 0])
            sites[i].append(x if v[6] == '+' else list(x[::-1]))
        else :
            sites[i].append([0, 0, 0, 0])
    return sites


def estimate(nodes, tre, max_nGenotype, min_score, min_freq) :
    sample = nodes['__query__'].astype(float)
    cov = np.sum(sample, 1)
    if np.sum(cov) == 0 :
        return []
    q1, q3 = np.quantile(cov, 0.25), np.quantile(cov, 0.75)
    max_cov = q3 + 3. * np.max([q3 - q1, 2.])
    min_cov = q1 - 3. * np.max([q3 - q1, 2.])

    sample_idx = (cov > max(0, min_cov)) & (cov <= max_cov)
    sample = sample.T/(cov+1e-8)
    mean_cov = np.mean(cov[sample_idx])
    cov *= sample_idx/mean_cov

    polymorphic = ((np.sum(sample, 0) - np.max(sample, 0)) > min_freq) * cov
    if np.sum(polymorphic) < 1 :
        polymorphic = cov
        max_nGenotype = 1
    sample -= min_freq

    explained = np.zeros(sample.shape, dtype=bool)
    used = _collections.OrderedDict()
    for nGenotype in np.arange(max_nGenotype) :
        node_map = []
        for node in tre.iter_descendants('postorder') :
            if node.name in used :
                continue
            g0, g1 = np.argmax(nodes[node.up.name], 1), np.argmax(nodes[node.name], 1)
            s2 = (sample * (1-explained))
            ss = s2[g0, np.arange(sample.shape[1])] - s2[g1, np.arange(sample.shape[1])]
            s_coord = (np.sum(-ss[ss<0])+0.1)/(np.sum(np.abs(ss))+0.1)
            genotypes = np.argmax([s2[g0, np.arange(sample.shape[1])], s2[g1, np.arange(sample.shape[1])]], 0)
            genotypes = np.array([g0, g1])[genotypes, np.arange(g0.size)]

            bases = np.zeros(sample.shape, dtype=bool)
            bases[genotypes, np.arange(sample.shape[1])] = True
            additional_sites = np.sum((np.sum(explained | bases, 0) - np.sum(explained, 0))*cov)
            shared_poly = np.sum(((np.sum(explained | bases, 0) - np.sum(explained, 0)) * (np.sum(s2, 0) > 0) * polymorphic)) if nGenotype > 0 else additional_sites
            
            poly_sites = np.sum((np.sum(s2, 0) > 0)*polymorphic)
            poly_base = np.max([np.sum(np.sum(s2 * bases, 0)*polymorphic), 0])
            poly_res = np.sum((sample*polymorphic)[(explained | bases)])
            node_map.append([np.sqrt(poly_base/(poly_sites+.5))*shared_poly/(additional_sites+.5), poly_res, poly_base/(poly_sites+.5), additional_sites, shared_poly/(additional_sites+.5), s_coord, node, bases])

        poly_score, poly_res, poly_base, additional_sites, shared_poly, s_coord, node, bases = max(node_map, key=lambda n:n[:2])
        configure.logging.info(f'   Trying: {node.name}: {poly_score:.3f} \t- explained polymorphsm: {poly_base:.3f}\t- inaccurate polymorphism: {1-shared_poly:.3f}')
        if nGenotype > 0 and poly_score < min_score :
            break
        if node.name not in used :
            used[node.name] = []
        used[node.name].append([s_coord, node.name, bases])
        explained |= bases
    
    data = [[node, d] for node, dat in used.items() for d in dat]
    if len(data) > 1 :
        sample += min_freq
        explained = np.zeros(sample.shape, dtype=bool)
        for node, dat in data :
            explained |= dat[2]

        polymorphic = np.where(np.sum(explained, 0) > 1)[0]
        x, y, w = [], [], []
        for site in polymorphic :
            for b, presence in enumerate(zip(*[dat[2][:, site] for node, dat in data])) :
                if any(presence) :
                    x.append(np.array(presence, dtype=float))
                    y.append(sample[b, site])
                    w.append(cov[site])
        x, y, w = np.array(x), np.array(y), np.array(w)
        param = Ridge(positive=True, fit_intercept=False, random_state=0, alpha=0.1).fit(x, y, sample_weight=w)
        delta = np.abs(param.predict(np.array(x)) - y)
        idx = delta < (np.quantile(delta, 0.75) + 3* max(np.quantile(delta, 0.75)-np.quantile(delta, 0.25),0.15))
        param = Ridge(positive=True, fit_intercept=False, random_state=0, alpha=0.1).fit(x[idx], y[idx], sample_weight=w[idx])
        coef = param.coef_/np.sum(param.coef_)
        coef[coef < min_freq] = 0
        coef = coef/np.sum(coef)
    else :
        coef = [1.]
    bestModel = [[p]+d for (n, d), p in zip(data, coef) if p > 0]
    for m in bestModel :
        configure.logging.info(f'{m[2]}:{m[1]:.3f} - p={m[0]:.3f}')
    return bestModel



def seperate_strains(bam, ref_acc, sites, best_model) :
    snps = np.array([np.argmax(m[3], 0) for m in best_model ])
    pGenotype = np.array([m[0] for m in best_model])
    
    nGenotype = len(best_model)
    
    qry_sites = {}
    for idx, (site, var, bases) in enumerate(sites) :
        if len(var) :
            qry_sites[(var[0], var[1])] = snps[:, idx] if var[6] == '+' else 3 - snps[:, idx]

    p = subprocess.Popen('{0} mpileup -AB -q 0 -Q 0 {1}'.format(executables['samtools'], bam).split(),
                        universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    contigs = {}
    seqs = [{} for _ in np.arange(nGenotype)]
    for line in p.stdout :
        p = line.strip().split('\t')
        if re.split('__', p[0])[1] != ref_acc :
            continue

        if p[0] not in contigs :
            contigs[p[0]] = int(p[1]) - 1
            for s in seqs :
                s[p[0]] = []

        s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
        s = list(''.join([b[int(n):] for n, b in re.findall('[+-](\d+)(.+)', '+0' + s)]))
        cons, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
        if cons in ('*', 'N') :
            contigs[p[0]] -= 1
            continue

        bases = dict(zip(*np.unique(s, return_counts=True)))
        bases = np.array([int(bases.get(b, 0.)) for b in ('A', 'C', 'G', 'T')])
        n_bases = np.sum(bases)
        bi = set(np.where(bases>0)[0])

        site = int(p[1]) - contigs[p[0]]
        key = (p[0], site)
        if key in qry_sites :
            bi |= set(qry_sites[key].tolist())
        bi = np.array(list(bi))
        
        max_comb = None
        for i in np.arange(bi.size**nGenotype) :
            pp_idx = np.array([ bi[int(i/(bi.size**pi))%bi.size] for pi in np.arange(nGenotype) ])
            base2 = np.zeros(4)
            for pp, pt in zip(pGenotype, pp_idx) :
                base2[pt] += pp
            base2 *= n_bases
            diff = np.sum(np.abs(bases-base2))*0.5
            same = n_bases - diff
            prob = np.log(0.02)*diff + np.log(0.94)*same
            if key in qry_sites :
                for gt, pt in zip(qry_sites[key], pp_idx) :
                    if gt == pt :
                        prob += np.log(9997.)

            if not max_comb or prob > max_comb[0] :
                max_comb = [prob, pp_idx, same/n_bases]
        for i, (b, c) in enumerate(zip(max_comb[1], pGenotype * n_bases)) :
            if max_comb[2] < 0.8 :
                c = 0
            elif bases[b] < c :
                c = bases[b]
            seqs[i][p[0]].append('ACGT'[b] if c>=3 else 'acgt'[b])
            
    return [ {n:''.join(s) for n, s in seq.items()} for seq in seqs ]



@click.command()
@click.option('-d', '--resolve_db', help='resolve_db generated by build_resolveDB')
@click.option('-q', '--query', help='query results generated by genoQuery')
@click.option('-r', '--ref', help='reference. could be accession, tax_id, or taxonomy')
@click.option('-o', '--outdir', help='folder storing the outputs. default: same as query', default=None)
@click.option('-n', '--num_genotype', help='maximum number of genotypes. default:20', default=20)
@click.option('-s', '--poly_score', help='minimum level of poly_score. default 0.14', default=0.14)
@click.option('-f', '--min_freq', help='maximum level of frequency of a genotype. default 0.02', default=0.02)
def explore(resolve_db, query, ref, outdir, num_genotype, poly_score, min_freq): 
    resolve_db = os.path.abspath(resolve_db)
    
    if query.endswith('profile.json') or query.endswith('primary.bam') :
        query = os.path.dirname(query)
    
    if not outdir :
        outdir = query
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    prefix = os.path.join(outdir, f'resolved.{ref}')
    if os.path.islink(f'{prefix}.db') :
        os.unlink(f'{prefix}.db')
    subprocess.run(f'ln -s {resolve_db} {prefix}.db'.split())
    
    nwk, aln = os.path.join(resolve_db, 'uscg.nwk'), os.path.join(resolve_db, 'uscg.concat.fas')
    uscg, bam = os.path.join(query, 'profile.json'), os.path.join(query, 'primary.bam')
    
    aln_fas = configure.readFasta(aln)
    configure.logging.info('Reading database.')
    if not os.path.isfile((dump_file := os.path.join(resolve_db, 'tree_info.dump'))) :
        nodes, sites, tre = get_sites(aln_fas, nwk)
        pickle.dump([nodes, sites, tre], open(dump_file, 'wb'))
    else :
        nodes, sites, tre = pickle.load(open(dump_file, 'rb'))
    
    uscgs = json.load(open(uscg))
    otus = [o for o in uscgs['OTU'] if '\t'.join(o[3:6]).find(ref) >= 0]

    res = {'OTU':[], 'profile':[]}

    strains = {}
    for otu in otus :
        ref_acc = otu[4]
        configure.logging.info('Retrieving BAM alignment.')
        sites = map_qry(aln_fas, otu, sites)
        sites = parse_bam(bam, ref_acc, sites)
        nodes['__query__'] = np.array([s[2] for s in sites])

        best_model = estimate(nodes, tre, num_genotype, poly_score, min_freq)# explained, inaccurate)
        configure.logging.info(f'Writing {len(best_model)} OTUs.')
        if len(best_model) == 0 :
            continue
        elif len(best_model) <= 1 :
            m = best_model[0]
            otu[2] = m[1]
            otu[4] = m[2]
            res['OTU'].append(otu)
        else :
            seqs = seperate_strains(bam, ref_acc, sites, best_model)
            otus = []
            for m, seq in zip(best_model, seqs) :
                seq = { re.split('__',n)[0] + f'__{m[2]}':s for n, s in seq.items() }
                otus.append([otu[0]*m[0], int(otu[1]*m[0]+0.5), m[1], otu[3], m[2], otu[5], seq])
            res['OTU'].extend(otus)

        for i, otu in enumerate(res['OTU']) :
            if otu[4] not in strains :
                strains[otu[4]] = []
            strains[otu[4]].append([otu[2], int(otu[0]*1000+0.5)/1000., f'OTU{i}'])

    json.dump(res, open(prefix + '.json', 'wt'))

    for node in tre.iter_descendants('preorder') :
        if node.name in strains :
            n_dist = node.dist
            for loc, depth, name in sorted(strains[node.name]) :
                parent = node.up
                n_dist1 = n_dist * (1-loc)
                new0 = ete3.TreeNode(dist=node.dist - n_dist1)
                new0.up = parent
                new1 = ete3.TreeNode(dist=0., name=f'{name}|{node.name}_{depth}')

                node.dist = n_dist1
                node.up = new0
                new1.up = new0

                parent.remove_child(node)
                parent.add_child(new0)
                new0.add_child(new1)
                new0.add_child(node)
    tre.write(format=1, outfile=f'{prefix}.nwk')
    configure.logging.info('Done.')

      
if __name__ == '__main__' :
    explore()
