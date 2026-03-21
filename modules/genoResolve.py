import os, ete3, numpy as np, subprocess, click, tempfile, pickle, json, re
import configure
from scipy.special import gammaln, xlogy, logsumexp
import warnings
warnings.filterwarnings('ignore')
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

executables = configure.executables
rc = configure.rc


def get_log_comb(a, b) :
    """Log combinatorial term for Binomial distribution"""
    n = a + b
    log_comb = np.zeros_like(a, dtype=float)
    mask = n > 0
    log_comb[mask] = (gammaln(n[mask] + 1) -
                    gammaln(a[mask] + 1) -
                    gammaln(b[mask] + 1))
    return log_comb


def log_binom(a, b, p, log_comb=None):
    """Complete log Binomial(a | a+b, p) with combinatorial term"""
    p = np.clip(p, 1e-12, 1 - 1e-12)
   
    res = xlogy(a, p) + xlogy(b, 1.0 - p)
    res += log_comb if log_comb is not None else get_log_comb(a, b)
    return res


def get_sites(aln_fas, nwk):
    """Extract variable sites and perform ancestral state reconstruction"""
    names = list(aln_fas.keys())
    fas = np.array([list(aln_fas[n]) for n in names]).T
    sites, vseqs = [], []

    for i, s in enumerate(fas):
        stype = np.unique(s)
        if stype.size > 2 or (stype.size == 2 and stype[0] != '-'):
            sites.append(i+1)
            vseqs.append(s)
    vseqs = [''.join(s) for s in np.array(vseqs).T]

    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir:
        with open(os.path.join(tmpdir, 'aln'), 'wt') as fout:
            for n, s in zip(names, vseqs):
                fout.write(f'>{n}\n{s}\n')

        subprocess.Popen('{0} --ancestral -s {1} -te {2} -m GTR+G4 -redo --prefix {3}'.format(
            configure.executables['iqtree'], os.path.join(tmpdir, 'aln'), nwk, os.path.join(tmpdir, 'aln')).split(),
            stdout=subprocess.PIPE).communicate()
        tre = ete3.Tree(os.path.join(tmpdir, 'aln.treefile'), format=1)
        nodes = {n.name:[] for n in tre.traverse() if not n.is_leaf()}

        with open(os.path.join(tmpdir, 'aln.state'), 'rt') as fin:
            for line in fin:
                if line.startswith('#') or line.startswith('Node\t'):
                    continue
                p = line.strip().split()
                base = np.argmax([float(v) for v in p[3:]])
                x = [0., 0., 0., 0.]
                x[base] = 1.
                nodes[p[0]].append(x)
        base_encoding = {'A':[1., 0., 0., 0.], 'C':[0., 1., 0., 0.],
                         'G':[0., 0., 1., 0.], 'T':[0., 0., 0., 1.]}
        for n, s in zip(names, vseqs):
            nodes[n] = [base_encoding.get(b, [0., 0., 0., 0.]) for b in s]
        for n, s in nodes.items():
            nodes[n] = np.array(s)

        for n in tre.iter_descendants('postorder'):
            if n.is_leaf():
                mask = (np.sum(nodes[n.name], 1) == 0)
                nodes[n.name][mask] = nodes[n.up.name][mask]
    return nodes, [[s, []] for s in sites], tre

def map_qry(aln_fas, profiles):
    """Map query sequences to reference alignment using minimap2"""
    aligns = {}
   
    with tempfile.TemporaryDirectory(prefix='se_', dir='.') as tmpdir:
        with open(os.path.join(tmpdir, 'ref'), 'wt') as fout:
            n, s = list(aln_fas.items())[0]
            s = s.upper()
            fout.write(f'>ref\n{s}\n')
            ref_seq = {'ref':s}
        qry_seq = {}
        with open(os.path.join(tmpdir, 'qry'), 'wt') as fout:
            for otu in profiles['OTU']:
                taxon = otu[3].rsplit('[', 1)[-1].split(']', 1)[0]
                for n, s in otu[6].items():
                    s = s.upper()
                    qry_seq[n] = list(s)
                    fout.write(f'>{taxon}|{otu[4]}|{n}\n{s}\n')

        alns = []
        map_cmd = f"{executables['minimap2']} -k13 -w5 -c -t1 --frag=yes --rmq -A1 -B4 -O8,16 -E2,1 -r20k,40k -g10k -P -N5000 -f1000,5000 -n2 -m20 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes ref qry"
        p = subprocess.Popen(map_cmd.split(), cwd=tmpdir, stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.PIPE)
        for line in p.stdout:
            if line.startswith('['):
                continue
            p = line.strip().split('\t')
            p[9:11] = int(p[9]), 100. * float(p[9])/float(p[10])
            if p[10] >= 90 :
                p[1:4] = [int(p[1]), int(p[2]) + 1, int(p[3])]
                p[6:9] = [int(p[6]), int(p[7]) + 1, int(p[8])]
                alns.append(p)
       
        x0 = 0
        for p in sorted(alns, key=lambda x:x[7]):
            taxon, ref, contig = p[0].split('|', 2)
            if (taxon, ref) not in aligns:
                aligns[(taxon, ref)] = {}

            if p[4] == '-':
                continue
            qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1
           
            for s, t in re.findall(r'(\d+)([MDI])', cigar):
                s = int(s)
                if t != 'I' and t != 'D' :
                    for x in range(0, s) :
                        rx = ri + x
                        qx = qi + x
                        rseq = ref_seq[p[5]][rx-1]
                        qseq = qry_seq[contig][qx-1]
                        if qseq != '-':
                            rr = aligns[(taxon, ref)].get(rx, [])
                            if len(rr) == 0 or rr[5] < p[10]:
                                aligns[(taxon, ref)][rx] = [contig, qx, rseq, qseq, p[9], p[10], p[4]]
                if t != 'I' :
                    ri += s
                if t != 'D' :
                    qi += s
            qry_seq[contig][p[2]-1:p[3]] = ['-'] * (p[3] - p[2] + 1)

    taxon = max([[len(aln), key[0], key[1]] for key, aln in aligns.items()])[1]
    sites = [[site, a] for key, aln in aligns.items() if key[0] == taxon for site, a in aln.items()]
    return taxon, sites


def parse_bam(bam, sites):
    """Parse BAM file to extract base compositions at all sites"""
    base_comp = {}
    contigs = {}
    if bam:
        qry_sites = {}
        for site, var in sites:
            if len(var):
                contigs[var[0]] = -1
                qry_sites[(var[0], var[1])] = ('ref', site)

        p = subprocess.Popen(f'{executables["samtools"]} mpileup -AB -q 0 -Q 0 {bam}'.split(),
                            universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        for line in p.stdout:
            p = line.strip().split('\t')
            if p[0] not in contigs:
                continue
            elif contigs[p[0]] < 0 :
                contigs[p[0]] = int(p[1]) - 1

            s = re.sub(r'\^.', '', p[4]).upper().replace('$', '')
            s = list(''.join([b[int(n):] for n, b in re.findall(r'[+-](\d+)(.+)', '+0' + s)]))
            base, cdp = min(zip(*np.unique(s, return_counts=True)), key=lambda x:-x[1])
            if base in ('*', 'N'):
                contigs[p[0]] -= 1
                continue

            site = int(p[1]) - contigs[p[0]]
            key = (p[0], site)
            if key not in qry_sites:
                continue

            bases = dict(zip(*np.unique(s, return_counts=True)))
            base_comp[key] = [int(bases.get(b, 0.)) for b in ('A', 'C', 'G', 'T')]

    results = {s[0]:np.zeros(4, dtype=int) for s in sites}
    for site, v in sites:
        if len(v):
            x = base_comp.get((v[0], v[1]), [0, 0, 0, 0])
            if v[6] == '-' :
                x = list(x[::-1])
            results[site] += x
    return sorted(results.items())


def fit_model_em(samples, states, lineages, pi_flip=0.05, max_iter=100, tol=1e-6):
    K = len(lineages)
    logit_prior = np.log(pi_flip) - np.log(1.0 - pi_flip)
   
    site_reads = np.sum(samples, 1)
    geno_states = np.stack([states[l] for l in lineages], axis=0).astype(np.float32)
    remain_reads = site_reads - (samples*geno_states.any(0)).sum(1) # sites not covered by any known lineage

    p_vec = np.clip(np.arange(K+1, 0, -1)/((K+2)*(K+1)*0.5), tol, 1 - tol)
    r_vec = np.zeros(K+1)

    ll = -np.inf
    for ite in range(max_iter) :
        genotypes = geno_states * p_vec[:-1, None, None]
        per_base_sum = np.clip(genotypes.sum(0), 1e-12, None)
        genotypes /= per_base_sum
       
        geno_reads = np.einsum("ij,kij->ki", samples, genotypes, optimize=True)
        total_reads = np.vstack([geno_reads, [remain_reads]])

        other_reads = site_reads-total_reads
        log_comb = get_log_comb(total_reads, other_reads)
        logL_normal = log_binom(total_reads, other_reads, p_vec[:, None], log_comb)
        logL_flipped = log_binom(total_reads, other_reads, 1.0 - p_vec[:, None], log_comb)

        logLR = logL_flipped - logL_normal
        # Posterior probability of flip given data
        rate = np.clip(-(logit_prior + logLR), -300, 300)
        p_flip = 1 / (1 + np.exp(rate))

        # log-sum-exp for mixture likelihood
        r_vec = np.mean(p_flip, 1)
        log_mix = np.logaddexp(np.log(1.0 - r_vec)[:, None] + logL_normal, np.log(r_vec)[:, None] + logL_flipped)

        A_corr = (1 - p_flip) * total_reads + p_flip * other_reads
        p_corr = np.sum(A_corr, 1)/np.sum(site_reads)
        p_vec = np.clip(p_corr/p_corr.sum(), tol, 1 - tol)

        # ---------- log-likelihood ----------
        new_ll = np.sum(log_mix)
        if abs(new_ll - ll) < tol :
            break
        ll = new_ll

    x = (geno_states * (1-p_flip[:-1])[:, :, None])
    base_content = x + (1-geno_states)*(1-np.sum(x, 2)[:, :, None])/3
    results = []
    for k, lineage in enumerate(lineages):
        results.append([p_vec[k], r_vec[k], lineage])
    # print(results, new_ll)
    return results, new_ll, base_content, p_corr[-1]


def estimate(genotypes, states, max_nGenotype, min_rate=0.02, beam_width=3, pi_flip=0.05, min_bic_improvement=10):

    cov = np.sum(genotypes, axis=1)
    if np.sum(cov) == 0:
        return []
   
    median_cov = np.median(cov)
    site_weights = np.clip(cov / (median_cov + 1e-10), 0, 3) * median_cov

    mad = np.median(np.abs(cov - median_cov))
    robust_z = 0.6745 * (cov - median_cov) / (mad + 1e-10)
    mask = (np.abs(robust_z) < 5) & (site_weights > 0)

    genotypes = (genotypes[mask] / cov[mask][:, None] * site_weights[mask][:, None]+0.5).astype(int)
    states = {k: v[mask] for k, v in states.items()}
   
    n_sites = genotypes.shape[0]
    logger.info(f"Retained {np.sum(mask)}/{len(mask)} sites after filtering")

    beam = [{
        "lineages": [],
        "results": [],
        "bic": np.inf,
        "loglik": -np.inf
    }]

    # === 3. ITERATIVE GENOTYPE SELECTION ===
    for iteration in range(max_nGenotype):
        logger.info(f"Iteration {iteration + 1}/{max_nGenotype}")
        new_beam = []
       
        for model in beam:
            accepted_genotypes = model["lineages"]

            for idx, (lineage, state) in enumerate(states.items()):
                if idx % 100 == 0:
                    logger.info(f"  Evaluating lineage {idx + 1}/{len(states)}: {lineage}")
                if lineage in accepted_genotypes:
                    continue
                # if lineage.startswith('Node') : continue
                test_lineages = accepted_genotypes + [lineage]

                results, loglik, base_content, miss_rate = fit_model_em(genotypes, states, test_lineages, pi_flip=pi_flip)

                K = len(results)
                bic = -2 * loglik + (2 * K - 1) * np.log(n_sites)

                new_beam.append({
                    "lineages": test_lineages,
                    "results": results,
                    "loglik": loglik,
                    "bic": bic,
                    "base_content": base_content,
                    "miss_rate": miss_rate
                })
        # ---------- beam pruning ----------
        new_beam.sort(key=lambda x: x["bic"])
        best = new_beam[0]

        improvement = beam[0]["bic"] - best["bic"]
        logger.info(f"Best BIC improvement: {improvement:.2f}")

        if improvement < min_bic_improvement:
            break

        beam = new_beam[:beam_width]
    bases = beam[0]["base_content"]
    base0 = np.zeros([bases.shape[0], mask.size, bases.shape[2]], dtype=np.float32)
    base0[:, mask, :] = bases
    beam[0]["base_content"] = base0
    return beam[0]


@click.command()
@click.option('-d', '--resolve_db', help='resolve_db generated by build_resolveDB')
@click.option('-q', '--query', help='query results generated by genoQuery')
@click.option('-o', '--outdir', help='folder storing the outputs. default: same as query', default=None)
@click.option('-n', '--num_genotype', help='maximum number of genotypes per hits. default:10', default=10)
@click.option('-f', '--min_freq', help='minimum frequency of a genotype. default 0.02', default=0.02)
@click.option('--beam_width', help='beam width for beam search. default 1', default=1)
@click.option('--p_flip', help='probability of flipping reads. default 0.05', default=0.05)
def explore(resolve_db, query, outdir, num_genotype, min_freq, beam_width, p_flip):
    """
    Merged strain resolution combining greedy phylogenetic approach with NMF/Lasso refinement.
    """
    resolve_db = os.path.abspath(resolve_db)

    if query.endswith('profile.json') or query.endswith('primary.bam'):
        query = os.path.dirname(query)

    if not outdir:
        outdir = query
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    prefix = os.path.join(outdir, f'resolved')
    if os.path.islink(f'{prefix}.db'):
        os.unlink(f'{prefix}.db')
    subprocess.run(f'ln -s {resolve_db} {prefix}.db'.split())

    nwk, aln = os.path.join(resolve_db, 'uscg.nwk'), os.path.join(resolve_db, 'uscg.concat.fas')
    uscg, bam = os.path.join(query, 'profile.json'), os.path.join(query, 'primary.bam')

    aln_fas = configure.readFasta(aln)
    configure.logging.info('Reading database.')
    if not os.path.isfile((dump_file := os.path.join(resolve_db, 'tree_info.dump'))):
        nodes, mut_sites, tre = get_sites(aln_fas, nwk)
        pickle.dump([nodes, mut_sites, tre], open(dump_file, 'wb'))
    else:
        nodes, mut_sites, tre = pickle.load(open(dump_file, 'rb'))

    uscgs = json.load(open(uscg))

    res = {'OTU':[], 'profile':[]}

    strains = {}
    taxon, all_sites = map_qry(aln_fas, uscgs)
    if taxon is None:
        logger.warning("No valid taxon found")
        json.dump(res, open(prefix + '.json', 'wt'))
        return

    all_sites = parse_bam(bam, all_sites)
   
    sites_map = dict(all_sites)
   
    genotypes = np.array([sites_map.get(s[0], [0, 0, 0, 0]) for s in mut_sites])
   
    # best_model = json.load(open(f'{prefix}_best_model.json', 'rt'))
    # best_model['base_content'] = np.array(best_model['base_content'], dtype=np.float32)
    best_model = estimate(genotypes, nodes, num_genotype, min_freq, beam_width, p_flip)
    # best_model['base_content'] = best_model['base_content'].tolist()
    # json.dump(best_model, open(f'{prefix}_best_model.json', 'wt'))

    strains = {n: [p, r, query] for p, r, n in best_model['results']}
    # Update tree with strain placements
    logger.info('Updating phylogenetic tree...')
    for node in tre.iter_descendants('postorder'):
        if node.name in strains:
            p, r, query = strains[node.name]
            new_node = ete3.TreeNode(name=node.name, dist=1e-8)
            node.add_child(new_node)
            dist = len(mut_sites)/len(all_sites) * r
            new_node = ete3.TreeNode(name=f'{query}|{node.name}|{p:.2f}', dist=dist)
            node.add_child(new_node)
            node.name = ''
    tre.write(format=1, outfile=f'{prefix}.nwk')
    logger.info(f'Tree written to {prefix}.nwk')

    if len(best_model) == 0:
        logger.warning("No strains identified")
        json.dump(res, open(prefix + '.json', 'wt'))
        return
   
    logger.info(f'Identified {len(best_model)} strains')
    for i, (prop, flip_rate, lineage) in enumerate(best_model['results']):
        logger.info(f'  Strain {i+1}: {lineage} (proportion: {prop:.3f}, flip_rate: {flip_rate:.3f})')
   
    logger.info('Reconstructing genotype sequences...')
    seqs = reconstruct_genotype_sequences(all_sites, mut_sites, best_model)
   
    taxon_profile = [profile for profile in uscgs['profile'] if profile[3].find(taxon) > 0][0]
    taxon_profile[4] = [f'node__{n}' for n in best_model['lineages']]
    taxon_profile[6] = {f'concatenated__{n}':s for n, s in zip(best_model['lineages'], seqs) }
   
    res = {"profile": [taxon_profile], "OTU": []}
   
    for (p, r, n), s in zip(best_model['results'], seqs) :
        otu = [taxon_profile[0] * p, int(taxon_profile[1] * p + 0.5),
               100-100*r, taxon_profile[3], f'node__{n}',
               taxon_profile[5], {f'concatenated__{n}': s}]
        res['OTU'].append(otu)

    # Save results
    logger.info(f'Writing {len(res["OTU"])} OTUs to {prefix}.json')
    json.dump(res, open(prefix + '.json', 'wt'))
   
    logger.info('Done.')



def reconstruct_genotype_sequences(all_bases, mut_sites, best_model, min_posterior=0.75):
    """
    Vectorized reconstruction of genotype-specific consensus sequences.

    Exhaustively evaluates all genotype-base combinations at each site
    using multinomial likelihood + mutation uncertainty + phylogenetic prior.
    """

    BASES = np.array(['A', 'C', 'G', 'T'])
    EPS = 1e-6
    mut_rate = max(best_model.get('miss_rate', 0.01), 0.01)

    mut_site_idx = {site: idx for idx, (site, _) in enumerate(mut_sites)}

    if not best_model or not best_model.get('results'):
        return []

    nG = len(best_model['results'])
    base_content = best_model['base_content']  # (nG, n_mut_sites, 4)
    proportions = np.array([m[0] for m in best_model['results']])
    proportions = proportions / proportions.sum()

    logger.info(f"Vectorized reconstruction for {nG} genotypes")

    L = max(site for site, _ in all_bases)
    seqs = [np.full(L, '-', dtype='U1') for _ in range(nG)]

    # ---------- precompute all base combinations ----------
    # Shape: (n_config, nG)
    grids = np.meshgrid(*([np.arange(4)] * nG), indexing='ij')
    all_cfg = np.stack(grids, axis=-1).reshape(-1, nG)
    n_cfg = all_cfg.shape[0]

    # ---------- mutation noise template ----------
    eye4 = np.eye(4)

    for site, obs in all_bases:

        total = sum(obs)

        if total == 0:
            continue

        obs = np.array(obs, dtype=float)
        obs_log = np.log(obs + EPS)

        # ---------- expected proportions ----------
        # exp shape: (n_cfg, 4)
        exp = np.zeros((n_cfg, 4))

        # vectorized genotype contribution
        for g in range(nG):
            bases = all_cfg[:, g]  # (n_cfg,)
            exp += proportions[g] * (
                (1 - mut_rate) * eye4[bases] + mut_rate / 3.0 * (1 - eye4[bases])
            )

        exp = np.clip(exp, EPS, None)
        exp /= exp.sum(axis=1, keepdims=True)

        # ---------- multinomial log-likelihood ----------
        # obs dot log(exp)
        ll = obs @ np.log(exp.T)

        # ---------- phylogenetic prior ----------
        if site in mut_site_idx:
            idx_site = mut_site_idx[site]
            prior = np.zeros(n_cfg)

            for g in range(nG):
                prior += np.log(np.clip(base_content[g, idx_site, all_cfg[:, g]], EPS, 1))
            ll += prior

        # ---------- choose best ----------
        ll -= logsumexp(ll)
        lll = np.exp(ll)
       
        onehot_all = eye4[all_cfg]
        possibility = np.einsum('c,cgk->gk', lll, onehot_all)
        best_base = np.argmax(possibility, axis=1)
        posterior = possibility[np.arange(nG), best_base] / (possibility.sum(axis=1) + EPS)
       
        for g in range(nG):
            if posterior[g] >= 0.3 :
                if posterior[g] >= min_posterior:
                    seqs[g][site-1] = BASES[best_base[g]]
                else:
                    seqs[g][site-1] = BASES[best_base[g]].lower()
       
    result = [''.join(s) for s in seqs]

    logger.info(f"Reconstructed {len(result)} genotype sequences")
    for i, seq in enumerate(result):
        logger.info(f"  Genotype {i}: {len(seq)} bp")

    return result


if __name__ == '__main__':
    explore()

