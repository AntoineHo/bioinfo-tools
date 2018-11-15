#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main() :
    dic = {}
    f = open("relativebins.txt", "r")
    l = []
    for line in f :
        s = line.split("\t")
        l.append([s[0], float(s[1]), s[2], int(s[3])])
    f.close()

    l.sort(key=lambda x: x[3])

    for el in l :
        for n in el[2].split(",") :
            if el[0] not in dic :
                dic[el[0]] = [int(n)]
            else :
                dic[el[0]].append(int(n))

    names = {}
    for el in l :
        ln = str(el[3])
        nln = ""
        lnln = len(ln)
        bfr = lnln//3

        for p in range(len(ln)) :
            if p == bfr or (p-2)%3 == 0 :
                nln+="."+ln[p]
            else :
                nln+=ln[p]
        nln = nln + "bp"

        curname = el[0] + ": " + nln
        names[el[0]] = curname

    df0 = pd.DataFrame(dic)
    df = df0.rename(columns=names)
    df.index += 1
    f, ax = plt.subplots(figsize=(20, 10))
    title = "Relative SNV frequency along {} contigs".format(len(df.columns))
    df.plot(kind='bar', stacked=True, ax=ax, title=title, position=1.5)
    ax.set_xlabel("Relative bin (1%) of chromosome length")
    ax.set_ylabel("Cumulated # SNVs & INDELs")
    plt.savefig("relative.heterozigosity.allvariants" + ".png", format="png", dpi=600)

if __name__ == "__main__" :
    main()
