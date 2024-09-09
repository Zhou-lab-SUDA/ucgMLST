#!/usr/bin/env python

import sys, re

refLen = {}

if __name__ == '__main__' :
    cutoff = float(sys.argv[1])
    cin = sys.stdin if len(sys.argv) < 3 else open(sys.argv[2])
    
    for line in cin :
        part = line.strip().split('\t')
        part[1:4] = [ int(p) for p in part[1:4] ]
        part[6:11] = [ int(p) for p in part[6:11] ]
        if part[9] < part[10]*(1-cutoff) :
            continue

        if part[4] == '+' :
            start_gap = min((part[2]), (part[7]))
            end_gap = min((part[1])-(part[3]), (part[6])-(part[8]))
        else :
            start_gap = min((part[1])-(part[3]), (part[7]))
            end_gap = min((part[2]), (part[6])-(part[8]))

        gap = [start_gap, end_gap]

        if min(gap) >= 15 or min(gap) >= part[10]*0.15 :
            continue
        if part[10] < 75 and sum(gap) >= part[10] :
            continue
        sys.stdout.write(line)
        sys.stdout.flush()
