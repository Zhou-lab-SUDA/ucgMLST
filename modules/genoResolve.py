import os, ete3, numpy as np, subprocess, click, tempfile, pickle, json, re
import configure
from scipy.special import gammaln, xlogy, logsumexp
from typing import List, Tuple, Dict
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


def log_binom(a, b, p):
    """Complete log Binomial(a | a+b, p) with combinatorial term"""
    n = a + b
    # Avoid log(0!) issues
    log_comb = np.zeros_like(a, dtype=float)
    # mask = n > 0
    # log_comb[mask] = (gammaln(n[mask] + 1) - 
    #                   gammaln(a[mask] + 1) - 
    #                   gammaln(b[mask] + 1))
    return log_comb + xlogy(a, p) + xlogy(b, 1.0 - p)

# def log_binom(a, b, p):
#     """Stable log Binomial(a | a+b, p) without combinatorial term"""
#     return xlogy(a, p) + xlogy(b, 1.0 - p)


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

def map_qry(aln_fas, profiles, sites):
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
                aligns[(taxon, ref)] = [[s[0], []] for s in sites]            

            if p[4] == '+':
                qi, ri, cigar, d = p[2], p[7], p[-1][5:], 1
            else:
                qi, ri, cigar, d = p[3], p[7], p[-1][5:], -1
            
            while x0 < len(sites) and sites[x0][0] < ri:
                x0 += 1
            xi = x0
            for s, t in re.findall(r'(\d+)([MDI])', cigar):
                s = int(s)
                if t != 'I':
                    rj = ri + s
                if t != 'D':
                    qj = qi + s * d
                while xi < len(sites) and sites[xi][0] >= ri and sites[xi][0] < rj:
                    rd = sites[xi][0] - ri
                    rx = ri + (rd if t != 'I' else 0) - 1
                    qx = qi + (rd * d if t != 'D' else 0) - 1
                    rseq = ref_seq[p[5]][rx]
                    qseq = qry_seq[contig][qx] if d > 0 else configure.rc(qry_seq[contig][qx])
                    if qseq != '-':
                        rr = aligns[(taxon, ref)][xi][1]
                        if len(rr) == 0 or rr[5] < p[10]:
                            aligns[(taxon, ref)][xi][1] = [contig, qx + d, rseq, qseq, p[9], p[10], p[4]]
                    xi += 1
                ri, qi = rj, qj
            qry_seq[contig][p[2]-1:p[3]] = ['-'] * (p[3] - p[2] + 1)

    taxon = max([[sum([len(a[1]) > 0 for a in aln]), key[0], key[1]] for key, aln in aligns.items()])[1]
    sites = [a[:2] for key, aln in aligns.items() if key[0] == taxon for a in aln]
    return taxon, sites


def parse_bam(bam, taxon, sites):
    """Parse BAM file to extract base compositions at variant sites"""
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


def identify_differentiating_sites(new, existing, states, genotypes):
    if not existing:
        return np.sum(genotypes, axis=1) > 0

    new_state = states[new[0]].argmax(1)
    existing_state = np.vstack([states[g].argmax(1) for g in existing])

    return np.all(existing_state != new_state, axis=0) & (np.sum(genotypes, axis=1) > 0)


def analyze_read_support_distribution(genotypes, new_genotypes, existing_genotypes, states, diff_sites, pi_flip = 0.05, max_iter = 100, tol = 1e-6):
    """
    Analyze the distribution of read support at differentiating sites to detect recombination.
    """
    if np.sum(diff_sites) == 0:
        return -np.inf, 0.0, 0.0

    new_states = np.array([states[new_genotype] for new_genotype in new_genotypes])
    # existing_states = (1. - new_states) if len(existing_genotypes) == 0 else np.array([states[g] for g in existing_genotypes])

    genotypes = genotypes[diff_sites]
    new_states = new_states[:, diff_sites].max(0)
    existing_states = (1. - new_states) if len(existing_genotypes) == 0 else np.array([states[g] for g in existing_genotypes])[:, diff_sites].max(0)
    
    new_state_reads = (genotypes * new_states).sum(axis = 1)
    other_state_reads = (genotypes * existing_states).sum(axis = 1)
    total_reads = new_state_reads + other_state_reads
    
    if total_reads.sum() < 1 :
        return -np.inf, 0.0, 0.0
    
    pA = np.sum(new_state_reads)/np.sum(total_reads)
    pA = 0.9 if pA > 0.5 else 0.1
        
    logit_prior = np.log(pi_flip) - np.log(1.0 - pi_flip)
    
    for i in range(max_iter) :
        logL_normal = log_binom(new_state_reads, other_state_reads, pA)
        logL_flipped = log_binom(new_state_reads, other_state_reads, 1.0 - pA)

        logLR = logL_flipped - logL_normal
        
        # Posterior probability of flip given data
        rate = np.clip(-(logit_prior + logLR), -300, 300)
        p_flip = 1 / (1 + np.exp(rate))

        A_corr = (1 - p_flip) * new_state_reads + p_flip * other_state_reads
        pA_corr = np.sum(A_corr)/np.sum(total_reads)
        if abs(pA_corr - pA) < tol:
            break
        pA = np.clip(pA_corr, tol, 1 - tol)
    
    # log-sum-exp for mixture likelihood
    mean_p_flip = np.mean(p_flip)
    log_mix = np.logaddexp(np.log(1.0 - mean_p_flip) + logL_normal, np.log(mean_p_flip) + logL_flipped)
    log_likelihood = np.sum(log_mix)
    
    return log_likelihood, pA, np.mean(p_flip)


def fit_model_em(samples, states, lineages, pi_flip=0.05, tol=1e-6):
    results = refine_proportions_em(samples, states, lineages, pi_flip, tol=tol)

    ll = -np.inf
    for ite in range(100) :
        genotypes = [ states[lineage] for lineage in lineages ] * np.array([r[1] for r in results])[:, None, None]
        per_base_sum = genotypes.sum(0)
        per_base_sum[per_base_sum == 0] = 1e-10
        genotypes = genotypes/per_base_sum
        
        for lineage, genotype, res in zip(lineages, genotypes, results) :
            geno_reads = np.sum(samples * genotype, 1)
            other_reads = np.sum(samples, 1) - geno_reads

            logit_prior = np.log(pi_flip) - np.log(1.0 - pi_flip)

            # for i in range(100) :
            logL_normal = log_binom(geno_reads, other_reads, res[1])
            logL_flipped = log_binom(geno_reads, other_reads, 1.0 - res[1])

            logLR = logL_flipped - logL_normal
                
            # Posterior probability of flip given data
            rate = np.clip(-(logit_prior + logLR), -300, 300)
            p_flip = 1 / (1 + np.exp(rate))

            A_corr = (1 - p_flip) * geno_reads + p_flip * other_reads
            p_corr = np.sum(A_corr)/np.sum(samples)
            res[1] = np.clip(p_corr, tol, 1 - tol)
            
            # log-sum-exp for mixture likelihood
            mean_p_flip = np.mean(p_flip)
            log_mix = np.logaddexp(np.log(1.0 - mean_p_flip) + logL_normal, np.log(mean_p_flip) + logL_flipped)
            res[0] = np.sum(log_mix)
            res[2] = np.mean(p_flip)
        new_ll = np.sum([r[0] for r in results])
        if abs(new_ll - ll) < 1e-6 :
            break
        ll = new_ll
    print(results)
    return results, new_ll


def estimate(genotypes, states, tre, max_nGenotype, min_rate, beam_width=3, pi_flip=0.05, min_bic_improvement=10):

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

            for lineage, state in states.items():
                if lineage in accepted_genotypes:
                    continue
                
                test_lineages = accepted_genotypes + [lineage]

                results, loglik = fit_model_em(genotypes, states, test_lineages, pi_flip=pi_flip)

                K = len(results)
                bic = -2 * loglik + (2 * K - 1) * np.log(n_sites)

                new_beam.append({
                    "lineages": test_lineages,
                    "results": results,
                    "loglik": loglik,
                    "bic": bic,
                })
        # ---------- beam pruning ----------
        new_beam.sort(key=lambda x: x["bic"])
        best = new_beam[0]

        improvement = beam[0]["bic"] - best["bic"]
        logger.info(f"Best BIC improvement: {improvement:.2f}")

        if improvement < min_bic_improvement:
            break

        beam = new_beam[:beam_width]
    return beam[0]["results"]


def refine_proportions_em(genotypes, states, lineages, pi_flip: float = 0.05, max_iterations: int = 100, tol: float = 1e-6):
    """
    EM refinement of genotype proportions with
    genotype-specific informative site sets.
    """

    # lineages = [r[3] for r in results]
    K = len(lineages)

    results = []
    for lineage in lineages:
        others = [x for x in lineages if x != lineage]
        diff_sites = identify_differentiating_sites([lineage], others, states, genotypes)

        # Fit pA once per genotype (as you already do)
        ll, pA, pR = analyze_read_support_distribution(genotypes, [lineage], others, states, diff_sites, pi_flip=pi_flip, max_iter=max_iterations, tol=tol)
        results.append([ll, pA, pR, lineage])

    p_sum = np.sum([r[1] for r in results])
    for r in results:
        r[1] = np.clip(r[1] / p_sum, tol, 1 - tol) if p_sum > 0 else 1.0 / K
    return results 


def reconstruct_genotype_sequences(bam, ref_acc, sites, best_model, min_posterior=0.9):
    """
    Reconstruct genotype-specific consensus sequences using
    posterior base probabilities (no read assignment).
    """

    nG = len(best_model)
    proportions = np.array([m[0] for m in best_model])
    genotype_states = [np.argmax(m[3], axis=0) for m in best_model]

    seqs = [{} for _ in range(nG)]
    contig_pos = {}

    p = subprocess.Popen(
        f"{executables['samtools']} mpileup -AB -q 0 -Q 0 {bam}".split(),
        universal_newlines=True,
        stdout=subprocess.PIPE
    )

    for line in p.stdout:
        p0 = line.strip().split('\t')
        if re.split('__', p0[0])[1] != ref_acc:
            continue

        contig = p0[0]
        pos = int(p0[1]) - 1

        if contig not in contig_pos:
            contig_pos[contig] = pos
            for s in seqs:
                s[contig] = []

        # parse bases
        s = re.sub(r'\^.', '', p0[4]).upper().replace('$', '')
        s = list(''.join([b[int(n):] for n, b in re.findall(r'[+-](\d+)(.+)', '+0' + s)]))

        bases, counts = np.unique(s, return_counts=True)
        base_counts = dict(zip(bases, counts))
        total = sum(base_counts.get(b, 0) for b in "ACGT")

        if total == 0:
            for s in seqs:
                s[contig].append('N')
            continue

        obs = np.array([base_counts.get(b, 0) for b in "ACGT"])

        # genotype-wise posterior
        for gi in range(nG):
            prior = np.zeros(4) + 1e-3
            prior[genotype_states[gi][pos]] = 1.0

            post = obs * prior
            post = post / post.sum()

            b = np.argmax(post)
            if post[b] >= min_posterior:
                seqs[gi][contig].append("ACGT"[b])
            elif post[b] >= 0.6:
                seqs[gi][contig].append("acgt"[b])
            else:
                seqs[gi][contig].append("N")

    return [{k: ''.join(v) for k, v in s.items()} for s in seqs]


@click.command()
@click.option('-d', '--resolve_db', help='resolve_db generated by build_resolveDB')
@click.option('-q', '--query', help='query results generated by genoQuery')
@click.option('-r', '--ref', help='reference. could be accession, tax_id, or taxonomy')
@click.option('-o', '--outdir', help='folder storing the outputs. default: same as query', default=None)
@click.option('-n', '--num_genotype', help='maximum number of genotypes per hits. default:10', default=10)
@click.option('-f', '--min_freq', help='minimum frequency of a genotype. default 0.02', default=0.02)
@click.option('-beam_width', help='beam width for beam search. default 3', default=3)
@click.option('-p_flip', help='probability of flipping reads. default 0.05', default=0.05)
def explore(resolve_db, query, ref, outdir, num_genotype, min_freq, beam_width, p_flip):
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
        nodes, sites, tre = get_sites(aln_fas, nwk)
        pickle.dump([nodes, sites, tre], open(dump_file, 'wb'))
    else:
        nodes, sites, tre = pickle.load(open(dump_file, 'rb'))

    uscgs = json.load(open(uscg))

    res = {'OTU':[], 'profile':[]}

    strains = {}
    taxon, sites = map_qry(aln_fas, uscgs, sites)
    if taxon is None:
        logger.warning("No valid taxon found")
        json.dump(res, open(prefix + '.json', 'wt'))
        return

    sites = parse_bam(bam, taxon, sites)
    
    genotypes = np.array([s[1] for s in sites])
    best_model = estimate(genotypes, nodes, tre, num_genotype, min_freq, beam_width, p_flip)
    
    logger.info(f'Writing {len(best_model)} OTUs.')
        
    if len(best_model) == 0:
        json.dump(res, open(prefix + '.json', 'wt'))
        return
    elif len(best_model) <= 1:
        m = best_model[0]
        # Find matching OTU and update
        for otu in uscgs['OTU']:
            otu_copy = list(otu)
            otu_copy[2] = m[1]
            # Store genotype identifier
            genotype_id = id(m[3])
            for name, state in nodes.items():
                if np.array_equal(state, m[3]):
                    genotype_id = name
                    break
            otu_copy[4] = genotype_id
            res['OTU'].append(otu_copy)
    else:
        seqs = reconstruct_genotype_sequences(bam, taxon, sites, best_model)
        for otu in uscgs['OTU']:
            otus = []
            for m, seq in zip(best_model, seqs):
                seq = {re.split('__',n)[0] + f'__{m[2]}': s for n, s in seq.items()}
                genotype_id = id(m[3])
                for name, state in nodes.items():
                    if np.array_equal(state, m[3]):
                        genotype_id = name
                        break
                otus.append([
                    otu[0]*m[1], int(otu[1]*m[1]+0.5),
                    m[2], otu[3], genotype_id, otu[5],
                    seq])
            res['OTU'].extend(otus)
    
    for i, otu in enumerate(res['OTU']):
        if otu[4] not in strains:
            strains[otu[4]] = []
        strains[otu[4]].append([otu[2], int(otu[0]*1000+0.5)/1000., f'OTU{i}'])

    json.dump(res, open(prefix + '.json', 'wt'))

    for node in tre.iter_descendants('preorder'):
        if node.name in strains:
            n_dist = node.dist
            for loc, depth, name in sorted(strains[node.name]):
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
    logger.info('Done.')


if __name__ == '__main__':
    explore()

