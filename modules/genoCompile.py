#!/usr/bin/env python
import sys, os, numpy as np, click, tempfile, subprocess, re, gzip
try :
    import ujson as json
except :
    import json

@click.command()
@click.option('-o', '--output', help='prefix for the outputs. Default: genoComplie', default='genoComplie')
@click.option('-m', '--min_rpkm', help='minimum level of RPKM to report. Default: 0.01', default=0.01, type=float)
@click.option('-M', '--min_ani', help='minimum level of ANI to report. Default: 0.95', default=0.95, type=float)
@click.option('-p', '--no_profile', help='do NOT show comparison of taxa profiles.', default=False, is_flag=True)
@click.option('-u', '--otu', help='show comparison of OTUs.', default=False, is_flag=True)
@click.argument('infiles', nargs=-1)
def main(infiles, output, no_profile, otu, min_rpkm, min_ani) :
    data = []
    for infile in infiles :
        if os.path.isfile(os.path.join(infile, 'profile.json')) :
            d = json.load(open(os.path.join(infile, 'profile.json'), 'rt'))
        else :
            d = json.load(open(infile, 'rt'))
        data.append(d)
    
    if not no_profile :
        profiles = {}
        for fname, d in zip(infiles, data) :
            for p in d['profile'] :
                if p[0] < min_rpkm or p[2] < min_ani:
                    continue
                if p[3] not in profiles :
                    profiles[p[3]] = {'taxonomy':p[5]}
                                        
                profiles[p[3]][fname] = [p[0], p[2]]
        profiles = sorted(profiles.items(), key=lambda p:-max([ p[1].get(fn, [0,0])[0] for fn in infiles ]))
        with open(f'{output}.profile', 'wt') as fout :
            fout.write('#Species(RPKM)\t{0}\t#Taxonomy\n'.format('\t'.join(infiles)))
            for s, p in profiles :
                fout.write('{0}\t{1}\t{2}\n'.format(s, '\t'.join( [ '{0:.4f}'.format(p.get(fn, [0,0])[0]) for fn in infiles ] ), p.get('taxonomy', '')))
    if otu :
        otu = {}
        for fname, d in zip(infiles, data) :
            for p in d['OTU'] :
                if p[0] < min_rpkm :
                    continue
                r = p[4][0]
                if r not in profiles :
                    otu[r] = {'taxonomy':p[5]}
                otu[r][fname] = [p[0], p[2]]
        otu = sorted(otu.items(), key=lambda p:-max([ p[1].get(fn, [0,0])[0] for fn in infiles ]))
        with open(f'{output}.otus', 'wt') as fout :
            fout.write('#Reference(RPKM)\t{0}\t#Taxonomy\n'.format('\t'.join(infiles)))
            for s, p in otu :
                fout.write('{0}\t{1}\t{2}\n'.format(s, '\t'.join( [ '{0:.4f}'.format(p.get(fn, [0,0])[0]) for fn in infiles ] ), p.get('taxonomy', '')))



if __name__ == '__main__' :
    main()