#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cairo
#import pickle

class Turtle :
    def __init__(self, vector, width, height, finalname) :
        self.a = open("t.txt", "w")
        self.vector = vector
        self.width = width
        self.height = height
        self.finalname = finalname
        self.vmin = min(vector)
        self.vmax = max(vector)
        self.maxsteps = len(vector)
        self.stps = 0
        self.x = 0
        self.y = 0
        # 0 is UP ; 1 is DOWN ; 2 is RIGHT ; 3 is LEFT
        self.o = 1

        # INITS CAIRO
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.ctx = cairo.Context(self.surface)
        self.ctx.scale(self.width, self.height)
        self.ctx.move_to(0, 0)
        self.ctx.set_line_width(0.02)
        self.ctx.set_source_rgb(1, 1, 1)
        self.ctx.rectangle(0, 0, 1, 1)
        self.ctx.fill()
        self.ctx.set_source_rgb(0, 0, 0)
        self.ctx.rectangle(0.001, 0.001, 0.998, 0.998)
        self.ctx.fill()
        self.ctx.set_source_rgb(1, 0, 0)
        self.ctx.rectangle(0, 0, 10/self.width, 10/self.height)  # Rectangle(x0, y0, x1, y1)
        self.ctx.fill()
    def turn_left(self) :
        if self.o == 0 :    # IF UP (0) THEN NEW IS LEFT (3)
            self.o = 3
        elif self.o == 1 :
            self.o = 2
        elif self.o == 2 :
            self.o = 0
        elif self.o == 3 :  # IF LEFT (3) THEN NEW IS DOWN (1)
            self.o = 1
    def turn_right(self) :
        if self.o == 0 :    # IF UP (0) THEN NEW IS RIGHT (2)
            self.o = 2
        elif self.o == 1 :
            self.o = 3
        elif self.o == 2 :
            self.o = 1
        elif self.o == 3 :
            self.o = 0
    def forward(self) :
        self.step()
    def step(self) :
        xy0 = (self.x, self.y)
        if self.o == 0 :    # IF UP (0) THEN Y + 1
            self.y += 1
        elif self.o == 1 :
            self.y -= 1
        elif self.o == 2 :
            self.x += 1
        elif self.o == 3 :
            self.x -= 1
        xy = (self.x, self.y)
        self.stps += 1
        if self.stps < self.maxsteps :
            self.draw(xy0, xy)
            self.a.write(str(xy) + "\n")
        else :
            self.a.close()
    def stop(self) :
        self.ctx.rectangle(1-10/self.width, 1-10/self.height, 10/self.width, 10/self.height)  # Rectangle(x0, y0, x1, y1)
        self.ctx.fill()
        self.ctx.stroke()
        self.surface.write_to_png(self.finalname)  # Output to PNG
    def draw(self, xy0, xy) :
        #print(self.stps)
        self.ctx.save()
        nxy0 = self.convert_coordinates(xy0)
        self.ctx.move_to(nxy0[0]+20/self.width, nxy0[1]+20/self.height)
        # COLOR DEPENDS OF THE VALUE IN VECTOR
        color = self.rgb(self.vmin, self.vmax, self.vector[self.stps])
        nxy = self.convert_coordinates(xy)
        self.ctx.line_to(nxy[0]+20/self.width, nxy[1]+20/self.height)
        rgb = [float((1/255)*cur) for cur in color]
        self.ctx.set_source_rgb(rgb[0], rgb[1], rgb[2])  # Solid color
        self.ctx.set_line_width(float(2/self.width))
        self.ctx.stroke()
        self.ctx.restore()
    def convert_coordinates(self, xy) :
        px = float(xy[0]*5/self.width)
        py = -float(xy[1]*5/self.height)
        return (px, py)
    def rgb(self, minimum, maximum, value):
        minimum, maximum = float(minimum), float(maximum)
        ratio = 2 * (value-minimum) / (maximum - minimum)
        b = int(max(0, 255*(1 - ratio)))
        r = int(max(0, 255*(ratio - 1)))
        g = 255 - b - r
        return r, g, b

def main() :
    """Gathers arguments and runs functions"""
    parser = argparse.ArgumentParser(description='Creates a heatmap for targets from .coverage file.')
    parser.add_argument('input',nargs=1,type=str,
                        help='A valid samtools depth coverage file path.')
    parser.add_argument('--TargetsNumber','-n',required=False,nargs=1,type=int,default=[1],
                        help='Number of targets to plot. Default: 1 (longest). To plot all targets set to 0.')
    parser.add_argument('--Specific','-s',required=False,nargs=1,type=int,default=[-1],
                        help='Specific target number to plot. Default: None (plot multiple contigs by default).')
    parser.add_argument('--Target','-t',required=False,nargs=1,type=str,default=["None"],
                        help='Specific target name to plot.')
    parser.add_argument('--output','-o',required=False,nargs=1,type=str,default=["covisu"],
                        help='Prefix of output files.')
    parser.add_argument('--plot','-p',required=False,nargs=1,type=str,default=["point"],choices=["line","heatmap","scatter","hilbert"],
                        help='Plot type of output.')
    parser.add_argument('--TargetList','-l',required=False,nargs=1,type=str,default=["None"],
                        help='Coma separated list of targets to plot.')
    args = parser.parse_args()

    infile = args.input[0]
    prefix = args.output[0]
    n_max = int(args.TargetsNumber[0])
    specific = int(args.Specific[0])
    typeplot = args.plot[0]
    target_list = list(args.TargetList[0].split(','))
    target_name = args.Target[0]

    # CHECK SPECIFIC -> always on top of other args
    print("## Reading file: {} ##".format(infile))
    if specific != -1 : # IF SPECIFIC ID IS SET
        matrix, headers = read_id(infile, specific) # READS A SPECIFIC ID
        print("## Plotting {}th target coverage...\ttype: {} ##".format(specific, typeplot))
        plot(matrix, headers, typeplot, prefix)
    else :
        if target_name != "None" : # IF SPECIFIC NAME IS SET
            matrix, headers = read_name(infile, target_name)
            print("Plotting {} coverage...\ttype: {} ##".format(target_name, typeplot))
            plot(matrix, headers, typeplot, prefix)
        else :
            if target_list[0] != "None" : # IF TARGET LIST IS SET
                matrix, headers = read_list(infile, target_list)
                """
                f1 = open("mx.pkl", "wb")
                f2 = open("hd.pkl", "wb")
                pickle.dump(matrix, f1)
                pickle.dump(headers, f2)
                f1.close()
                f2.close()
                f1 = open("mx.pkl", "rb")
                f2 = open("hd.pkl", "rb")
                matrix = pickle.load(f1)
                headers = pickle.load(f2)
                """
                print("## Plotting targets from list coverage...\ttype: {} ##".format(typeplot))
                plot(matrix, headers, typeplot, prefix)
            else :
                matrix, headers = read_0_to_id(infile, n_max) # DEFAULT n_max IS 1 IF NOTHING WAS SET
                if n_max not in [1,2,3] :
                    print("## Plotting first to {}th target coverage...\ttype: {} ##".format(n_max, typeplot))
                elif n_max == 1 :
                    print("## Plotting first target coverage...\ttype: {} ##".format(typeplot))
                elif n_max == 2 :
                    print("## Plotting first and second targets coverage...\ttype: {} ##".format(typeplot))
                elif n_max == 3 :
                    print("## Plotting first, second and third targets coverage...\ttype: {} ##".format(typeplot))
                plot(matrix, headers, typeplot, prefix)

def read_id(infile, p) :
    f = open(infile, "r")
    current_header = ""
    matrix = [[]]
    target_number = 0
    for line in f :
        if line[0] == "#" :
            continue
        s = line.split('\t')
        if s[0] != current_header :
            current_header = s[0]
            target_number += 1

        if target_number < p :
            continue
        elif target_number == p :
            matrix[0].append(int(s[2]))
        elif target_number > p :
            break

    return matrix, [current_header]

def read_name(infile, target_name) :
    f = open(infile, "r")
    headers = [target_name]
    matrix = [[]]
    for line in f :
        if line[0] == "#" :
            continue
        s = line.split('\t')
        if s[0] != target_name :
            if len(matrix[0]) > 1 :
                break
            else :
                continue
        else :
            matrix[0].append(int(s[2]))

    return matrix, headers

def read_0_to_id(infile, n_max) :
    f = open(infile, "r")
    headers = []
    matrix = []
    p = 0
    for line in f :
        if line[0] == "#" :
            continue
        s = line.split('\t')
        if s[0] not in headers and p < n_max :
            print("Next target: {}".format(s[0]))
            headers.append(s[0])
            matrix.append([])
            p += 1
        elif s[0] not in headers and p == n_max :
            break
        else :
            matrix[p-1].append(int(s[2]))

    return matrix, headers

def read_list(infile, target_list) :
    f = open(infile, "r")
    headers = target_list
    matrix = [[] for target in target_list]
    for line in f :
        if line[0] == "#" :
            continue
        s = line.split('\t')
        if s[0] not in target_list :
            continue
        else :
            idx = headers.index(s[0])
            matrix[idx].append(int(s[2]))

    return matrix, headers

def plot(matrix, headers, typeplot, prefix) :

    # SORT HEADERS & COV MATRIX BY LENGTH
    zipped = zip(headers, matrix)
    zipped = sorted(zipped, key= lambda x: len(x[1]), reverse=False)
    headers = []
    matrix = []
    for i,t in enumerate(zipped) :
        headers.append(t[0])
        matrix.append(t[1])

    plt.rcParams.update({'figure.autolayout': True})
    if typeplot == "scatter" :
        p = 0
        for target in matrix :
            x = [i for i in range(0, len(target))]
            plt.figure(figsize=(40,10))
            plt.rcParams['axes.axisbelow'] = True
            plt.grid(color='gray', linestyle='--')
            plt.ylim(-0.05*max(target),max(target)+0.05*max(target))
            plt.xlim(-0.05*len(target),len(target)+0.05*len(target))
            ax = plt.scatter(x=x, y=target, marker=".")
            plt.ylabel("Coverage")
            plt.xlabel("Position on target (bp)")
            plt.title(headers[p])
            plt.savefig(prefix + "." + str(headers[p]) + ".png", format="png", dpi=600)
            p += 1

    elif typeplot == "line" :
        p = 0
        for target in matrix :
            plt.figure(figsize=(40,10)) # prevents plotting on top of previous figure
            x = [i for i in range(0, len(target))]
            ax = sns.lineplot(x=x, y=target, palette="deep")
            ax.set_ylabel("Coverage")
            ax.set_xlabel("Position on target (bp)")
            ax.set_title(headers[p])
            plt.savefig(prefix + "." + str(headers[p]) + ".png", format="png", dpi=600)
            p += 1

    elif typeplot == "heatmap" :
        # DEFINES CHUNKS
        lengths = [np.mean([len(x) for x in matrix])]
        mean_length = np.mean(lengths)
        units = "bp"
        if mean_length > 1000000 :
            units = ["1Mbp", "100Kbp", "10Kbp"]
            chunksizes = [1000000, 100000, 10000]
        elif mean_length > 1000 :
            units = ["1Kbp", "100bp", "10bp"]
            chunksizes = [1000, 100, 10]
        else :
            units = ["100bp", "10bp", "1bp"]
            chunksizes = [100, 10, 1]

        # REMOVE TOO SMALL TARGETS
        pos = 0
        for target in matrix :
            if len(target) < 1000 :
                print("Cannot draw the heatmap of target as it is smaller than 1 Kb")
                del headers[pos]
                matrix.remove(target)
                continue
            pos += 1

        # CREATE A BUNCH PER CHUNKSIZE
        cp = 0
        for chksize in chunksizes :
            total_matrix = []
            for target in matrix :
                chunk_list = list(chunks(target, chksize))
                #print(len(chunk_list[-2]), len(chunk_list[-1]))
                mean_list = [int(np.mean(chunk)) for chunk in chunk_list]
                if len(chunk_list[-1]) < 0.5*chksize and len(chunk_list) != 1 :
                    mean_list[-1] = np.nan
                total_matrix.append(mean_list)
                #print(mean_list[-1])


            # Gets bunches in total_matrix
            headbuns, meanbuns = get_bunches(total_matrix, headers)
            longestbun = len(max(headbuns, key=len)) # Max number of lists in bunch -> vertical size of total plot
            #print(len(headbuns))

            # ADDS NaN IN MISSING POSITIONS OF EACH TARGET
            for bun in meanbuns :
                #print(bun)
                max_length = longest(bun) # Find longest list in current bunch and fill other lists with nan to this length
                for ml in bun :
                    while len(ml) < max_length :                        # Sets the same number of elements for every mean list
                        ml.append(np.nan)
                while len(bun) < longestbun :                           # Sets the same number of lists for every bunch
                    bun.append([np.nan for x in range(0, max_length)])

            maxmean = 0 # Finds vmax of all subplots
            for bun in meanbuns :
                for meanlist in bun :
                    if max(meanlist) > maxmean :
                        maxmean = max(meanlist)
                        #print(maxmean)

            # Computes nb of columns and lines according to nb of plots
            nplots = len(headbuns)
            nblin = 0
            nbcol = 0
            for x in range(1, nplots+1) :
                if x**2 >= nplots :
                    nblin = x
                    break
            for x in range(1, nblin+1) :
                if x*nblin >= nplots :
                    nbcol = x
                    break

            fig, axes = plt.subplots(nrows = nblin, ncols = nbcol, figsize=(20*nblin, 10*nbcol), sharex=False, sharey=False)

            if nbcol == 1 and nblin == 1 :
                curheaders = headbuns[0]
                curmeanlis = meanbuns[0]
                sns.heatmap(curmeanlis, cmap="RdBu_r", ax = axes, vmin=0, vmax=maxmean)
                axes.set_yticklabels(curheaders, rotation=0)

            else :
                sub = 0
                for ax in axes.flat :
                    if sub < len(headbuns) :
                        curheaders = headbuns[sub]
                        curmeanlis = meanbuns[sub]
                        sns.heatmap(curmeanlis, cmap="RdBu_r", ax = ax, vmin=0, vmax=maxmean)
                        ax.set_yticklabels(curheaders, rotation=0)
                    #ax.set_xlabel("Position on target (1:{})".format(units[cp]))
                    sub += 1

            plt.ylabel("Coverage")
            plt.xlabel("Position on target (1:{})".format(units[cp]))
            fig.savefig(prefix + ".{}".format(units[cp]) + ".png")
            cp += 1

    elif typeplot == "hilbert" :
        i = 0
        for vector in matrix :
            currentplotname = prefix + "." + str(headers[i]) + ".png"
            print("Plotting {}\t{} segments".format(headers[i], len(vector)))
            depth = 0
            while 4**depth < len(vector) :
                depth += 1
            if depth > 13 :
                raise Exception("Maximal recursive depth reached! Cannot output hilbert plot")
            print("Recursive depth={} #WARNING: File may take several Mbs!".format(depth))
            pxs = [50,60,80,120,200,360,680,1320,2600,5160,10280,20520]
            WIDTH = pxs[depth-1]
            HEIGHT = pxs[depth-1]
            turtle = Turtle(vector, WIDTH, HEIGHT, currentplotname)
            hilbert("a", depth, turtle)
            turtle.stop()
            i += 1

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def longest(l):
    if(not isinstance(l, list)): return(0)
    return(max([len(l),] + [len(subl) for subl in l if isinstance(subl, list)] + [longest(subl) for subl in l]))

def get_bunches(matrix, headers):
    bunches = []
    zipped = zip(headers, matrix)
    current = None
    maxlen = 0
    for n, v in enumerate(zipped) :

        if current == None :
            bunches.append([])
            current = v[0] # Current headers
            maxlen = len(v[1]) # list len of means
            #print("new max length = {}".format(maxlen))
            bunches[-1].append(v)
        else :
            curlen = len(v[1]) # list len of means
            #print("current length = {}".format(curlen))
            if curlen < int(1.25*maxlen) : # Current bunch
                current = v[0]
                #print("{} < {} smaller -> no new bunch".format(curlen, int(1.25*maxlen)))
                bunches[-1].append(v)
            else : # New bunch
                current = None
                #print("{} > {} bigger -> NEW bunch".format(curlen, int(1.25*maxlen)))

    headbunches = []
    meanbunches = []
    for bunch in bunches :
        headbunches.append([])
        meanbunches.append([])
        for v in bunch :
            headbunches[-1].append(v[0])
            meanbunches[-1].append(v[1])

    return headbunches, meanbunches

def hilbert(rule, depth, turtle) :
    if depth > 0 :
        a = lambda: hilbert("a", depth - 1, turtle)
        b = lambda: hilbert("b", depth - 1, turtle)

        left = lambda: turtle.turn_left()
        right = lambda: turtle.turn_right()
        forward = lambda: turtle.forward()

        if rule == "a":
            left()
            b()
            forward()
            right()
            a()
            forward()
            a()
            right()
            forward()
            b()
            left()
        if rule == "b":
            right()
            a()
            forward()
            left()
            b()
            forward()
            b()
            left()
            forward()
            a()
            right()


if __name__ == "__main__" :
    # main function
    main()
