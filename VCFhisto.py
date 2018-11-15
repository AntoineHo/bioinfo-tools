#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def main() :
    """Gathers arguments and runs functions"""
    parser = argparse.ArgumentParser(description='Creates a bargraph of positions from .vcf file.')
    parser.add_argument('input',nargs=1,type=str,
                        help='A valid VCF file path.')
    parser.add_argument('target',nargs=1,type=str,
                        help='Header of the sequence to plot.')
    parser.add_argument('--output','-o',required=False,nargs=1,type=str,default=["Filename"],
                        help='Prefix of output files.')
    parser.add_argument('--BinSize','-b',required=False,nargs=1,type=int,default=[1000],
                        help='Size of the bin used to count substitutions. Default: 1000\nWARNING: while creating VCF you should set INDELs to false, to keep only substitutions.')
    parser.add_argument('--MinScore','-m',required=False,nargs=1,type=int,default=[25],
                        help='Minimum phred score to consider. Default: 25.0')
    args = parser.parse_args()

    infile = args.input[0]
    target = args.target[0]
    binsize = int(args.BinSize[0])
    minscore = float(args.MinScore[0])

    prefix = args.output[0]
    if prefix == "Filename" :
        prefix = infile

    bins, rbins, rbinsize, tlen = read(infile, target, binsize, minscore)
    x = open("relativebins.txt", "a+")
    x.write(target + "\t" + str(rbinsize) + "\t" + ",".join(str(b) for b in rbins) + "\t" + str(tlen) + "\n")
    x.close()
    doplot(bins, binsize, rbins, rbinsize, target, prefix)

def read(input, target, binsize, minscore) :
    bins = []
    relativebins = []
    x = []
    print("## Reading metadata of VCF file...")
    f = open(input, 'r')
    cur = f.readline()
    while re.search(target, cur) == None :
        cur = f.readline()
        continue
    listOfInfo = cur[int(re.search("<", cur).span()[1]):int(re.search(">", cur).span()[0])].split(",")
    f.close()
    dictInfo = {}
    for el in listOfInfo :
        equalpos = int(re.search("=", el).span()[0])
        key = el[:equalpos]
        val = el[equalpos+1:]
        dictInfo[key] = val
    targetLength = int(dictInfo["length"])
    print("Target name:\t\t" + str(target))
    print("Length of target:\t" + str(targetLength))
    print("## Computing bins ...")
    bins = [0 for x in range((targetLength//binsize))] + [0]
    print("Created " + str(len(bins)) + " bins of " + str(binsize) + "bp")
    print("## Computing relative bins ...")
    relbinsize = 0.01*targetLength
    relativebins = [0 for x in range(100)]
    print("Created 100 bins (1%) of " + str(relbinsize) + "bp")
    f = open(input, 'r')
    for line in f :
        if line[0] == "#" :
            continue
        else :
            s = line.split('\t')
            if s[0] != target or float(s[5]) < minscore :
                continue
            else :
                cpos = int(s[1])
                #print(cpos)
                for x in range(len(bins)) :
                    if cpos > x*binsize :
                        continue
                    else :
                        bins[x] += 1
                        break
                for x in range(len(relativebins)) :
                    if cpos > x*relbinsize :
                        continue
                    else :
                        relativebins[x] += 1
                        break

    f.close()
    return bins, relativebins, relbinsize, targetLength

def doplot(numSUB, binsize, relbins, relbinsize, target, prefix) :
    print("## Plotting bar graph of substitutions in bins...")
    headers = [str(x*binsize) for x in range(len(numSUB))]

    fig, ax = plt.subplots(figsize=(30,10))
    #sns.set_color_codes("pastel")
    ind = np.arange(len(numSUB))
    ax.bar(ind, numSUB, -0.9, align="edge", edgecolor="black", linewidth=0.1)
    ax.set_xlim(-1,len(headers))
    ax.set_xlabel('Position along chromosome (bp)', fontsize=14)
    ax.set_ylabel('#Substitutions found', fontsize=14)
    ax.set_title('Heterozigosity along {}'.format(target), fontsize=20)
    xlabs = []
    for n in range(len(headers)) :
        if n % 5 == 0 :
            xlabs.append(headers[n])
        else :
            xlabs.append("")
    ax.set_xticks(range(len(headers)))
    ax.set_xticklabels(xlabs, rotation=90)
    plt.tight_layout()
    plt.savefig(prefix + "." + target + "." + str(binsize) + "bp" + ".png", format="png", dpi=600)

    print("## Plotting bar graph of substitutions in relative bins...")
    headers = [str(x*relbinsize) for x in range(len(relbins))]

    fig, ax = plt.subplots(figsize=(30,10))
    #sns.set_color_codes("pastel")
    ind = np.arange(len(relbins))
    ax.bar(ind, relbins, -0.9, align="edge", edgecolor="black", linewidth=0.1)
    ax.set_xlim(-1,len(headers))
    ax.set_xlabel('Position along chromosome (bp)', fontsize=14)
    ax.set_ylabel('#Substitutions found', fontsize=14)
    ax.set_title('Heterozigosity along {}'.format(target), fontsize=20)
    xlabs = []
    for n in range(len(headers)) :
        if n % 5 == 0 :
            xlabs.append(headers[n])
        else :
            xlabs.append("")
    ax.set_xticks(range(len(headers)))
    ax.set_xticklabels(xlabs, rotation=90)
    plt.tight_layout()
    plt.savefig(prefix + "." + target + ".relative" + ".png", format="png", dpi=600)

if __name__ == "__main__" :
    main()
