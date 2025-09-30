import os
import re
import statistics
from collections import defaultdict
import argparse
import textwrap
import sys


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


parser = argparse.ArgumentParser(
    prog="triage.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber;
    Middle Author Bionformatics, LLC.
    Please send comments and inquiries to ark@midauthorbio.com

    *******************************************************
    '''))


parser.add_argument('-a', type=str, help="antiSMASH parsed directory")
parser.add_argument('-f', type=str, help="FeGenie output directory")
parser.add_argument('-b', type=str, help="BLAST output directory")
parser.add_argument('-m', type=str, help="MiBIG directory")
parser.add_argument('-o', type=str, help="output directory")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

antismashDir = args.a
feGenieDir = args.f
blastDir = args.b
mibigDir = args.m
outDir = args.o

productDict = defaultdict(lambda: 0)
products = open("naturalproducts.csv")
for i in products:
    productDict[i.rstrip().lower()] = 1

knownDict = defaultdict(lambda: 0)
known = open("knownclusterblast.csv")
for i in known:
    knownDict[i.rstrip().lower()] = 1

tfbsDict = defaultdict(lambda: 0)
tfbs = open("tfbs.csv")
for i in tfbs:
    tfbsDict[i.rstrip().lower()] = 1

print("Parsing MiBIG files for known natural products...")
mibigDict = defaultdict(lambda: defaultdict(list))
mibigDict2 = defaultdict(lambda: defaultdict(list))
mibigDir = os.listdir(f"{antismashDir}/mibig")
for i in mibigDir:
    if re.search(r'csv', i):
        mibig = open(f"{antismashDir}/mibig/{i}")
        for j in mibig:
            ls = j.rstrip().split(",")
            product = "-"
            for k in productDict.keys():
                if re.search(k, ls[4].lower()):
                    product = "+"
                    break
            if product == "+":
                orf = ls[0].split("_")[1]
                if orf not in mibigDict[i.split(".")[0]][ls[0].split("_")[0]]:
                    mibigDict[i.split(".")[0]][ls[0].split("_")[0]].append(orf)
                    mibigDict2[i.split(".")[0]][ls[0].split("_")[0]].append(ls[4])


print("Parsing BLAST results for contig matches...")
blastDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '-')))
blastDir = os.listdir(f"{blastDir}")
for i in blastDir:
    if re.search("blast", i):
        blast = open(f"/Users/agarber4/Desktop/MABanalysis/UNMC/blast/{i}")
        for j in blast:
            ls = j.rstrip().split("\t")
            sample = i.split(".")[0]
            contig = ls[0].split("_")[0] + "_" + ls[0].split("_")[1]
            contig2 = ls[1].split("_")[0]
            blastDict[sample][contig] = contig2
            blastDict[sample][contig2] = contig

print("Parsing antiSMASH results...")
antiDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '0')))
antiDir = os.listdir(f"{antismashDir}")
for i in antiDir:
    if re.search("csv", i):
        file = open(f"{antismashDir}/{i}")
        for j in file:
            ls = j.rstrip().split(",")
            if ls[0] != "File":
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["category"] = ls[3]
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["product"] = ls[2]
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["tfbs"] = ls[12]

                tfbsList = []
                tfbsList2 = []
                for k in ls[12].split(";"):
                    tfbs = k.split("(")[0]
                    if tfbsDict[tfbs.lower()] == 1:
                        tfbsList.append(k.split("(")[2].split(")")[0].split("_")[1])
                        tfbsList2.append(k.split("(")[0])
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["tfbsFiltNums"] = tfbsList
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["tfbsFilt"] = tfbsList2

                if ls[7] != "":
                    antiDict[i.split(".")[0]][ls[0].split(".")[0]]["knownclusterblastSource"] = ls[7]
                else:
                    antiDict[i.split(".")[0]][ls[0].split(".")[0]]["knownclusterblastSource"] = "-"
                antiDict[i.split(".")[0]][ls[0].split(".")[0]]["knownclusterblastType"] = ls[8]

print("Parsing FeGenie results...")
feDict = defaultdict(lambda: defaultdict(list))
feLocusDict = defaultdict(lambda: defaultdict(list))
fegenie = open(f"{antismashDir}/FeGenie-geneSummary.csv")
for i in fegenie:
    if not re.match(r'#', i):
        ls = i.rstrip().split(",")
        if re.findall(r'iron_aquisition', ls[0]):
            contig = ls[2].split("_")[0] + "_" + ls[2].split("_")[1]
            feDict[ls[1].split("_")[0]][contig].append(ls[3])
            feLocusDict[ls[1].split("_")[0]][contig].append(lastItem(ls[2].split("_")))


print("Comparing antiSMASH and FeGenie results...")
synth = re.compile(r"synth")
sidDict = defaultdict(lambda: defaultdict(lambda: '-'))
sidSynthDict = defaultdict(lambda: defaultdict(list))
sidLocusDict = defaultdict(lambda: defaultdict(list))
for i in feDict.keys():
    for j in feDict[i]:
        present = [item for item in feDict[i][j] if synth.search(item)]
        if len(present) > 0:
            sidDict[i.split("_")[0]][j] = "+"
            sidSynthDict[i.split("_")[0]][j] = feLocusDict[i][j]
            sidLocusDict[i.split("_")[0]][j] = feDict[i][j]

out1 = open(f"{outDir}/antiSMASH-FeGenie.csv", "w")
out1.write("sample,contig,antiSMASH_category,antiSMASH_product,antiSMASH_knownclusterblastSource,antiSMASH_knownclusterblastType,TFBS\n")

out2 = open(f"{outDir}/antiSMASH-only.csv", "w")
out2.write("sample,contig,antiSMASH_category,antiSMASH_product,antiSMASH_knownclusterblastSource,antiSMASH_knownclusterblastType,TFBS\n")

out3 = open(f"{outDir}/antiSMASH-FeGenie-rescued-detailed.csv", "w")
out3.write("sample,contig,locus_tag,feature\n")

out4 = open(f"{outDir}/antiSMASH-FeGenie-rescued.csv", "w")
out4.write("sample,contig,locus_tags,features\n")

antiDict2 = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '0')))
for i in antiDict.keys():
    for j in antiDict[i]:
        fegenie = sidDict[i][j]
        try:
            product = "-"
            known = "-"

            for k in knownDict.keys():
                if re.search(k, antiDict[i][j]["knownclusterblastSource"].lower()):
                    known = "+"
                    antiDict2[i][j] = antiDict[i][j]
                    if fegenie == "+":
                        out1.write(f"{i},{j},{antiDict[i][j]['category']},{antiDict[i][j]['product']},{antiDict[i][j]['knownclusterblastSource']},{antiDict[i][j]['knownclusterblastType']},{antiDict[i][j]['tfbs']}\n")
                        break
                    else:
                        tfbs = "-"
                        for k in antiDict[i][j]['tfbs'].split(";"):
                            if tfbsDict[k.split("(")[0].lower()] == 1:
                                tfbs = "+"
                                break
                        if tfbs == "+":
                            out2.write(f"{i},{j},{antiDict[i][j]['category']},{antiDict[i][j]['product']},{antiDict[i][j]['knownclusterblastSource']},{antiDict[i][j]['knownclusterblastType']},{antiDict[i][j]['tfbs']}\n")
                        break
            if known == "+":
                continue

            for k in productDict.keys():
                if re.search(k, antiDict[i][j]["product"].lower()):
                    product = "+"
                    antiDict2[i][j] = antiDict[i][j]
                    if fegenie == "+":
                        out1.write(f"{i},{j},{antiDict[i][j]['category']},{antiDict[i][j]['product']},{antiDict[i][j]['knownclusterblastSource']},{antiDict[i][j]['knownclusterblastType']},{antiDict[i][j]['tfbs']}\n")
                        break
                    else:
                        tfbs = "-"
                        for k in antiDict[i][j]['tfbs'].split(";"):
                            if tfbsDict[k.split("(")[0].lower()] == 1:
                                tfbs = "+"
                                break
                        if tfbs == "+":
                            out2.write(f"{i},{j},{antiDict[i][j]['category']},{antiDict[i][j]['product']},{antiDict[i][j]['knownclusterblastSource']},{antiDict[i][j]['knownclusterblastType']},{antiDict[i][j]['tfbs']}\n")
                        break
            if product == "+":
                continue

        except ZeroDivisionError:
            pass

out1.close()
out2.close()

rescueDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: '-')))
for i in sidDict.keys():
    sample = i
    for j in sidDict[i]:
        contig = j
        if j not in antiDict2[i]:
            if len(sidSynthDict[i][j]) > 1 and antiDict[i][j]['tfbs'] != '0':
                tfbs = 0
                for k in antiDict[i][j]['tfbs'].split(";"):
                    if tfbsDict[k.split("(")[0].lower()] == 1:
                        tfbs = 1
                if tfbs == 1:
                    for k in range(0, len(sidSynthDict[i][j])):
                        rescueDict[sample][contig][int(sidSynthDict[i][j][k])] = sidLocusDict[i][j][k]
                    for k in range(0, len(antiDict[i][j]['tfbsFiltNums'])):
                        rescueDict[sample][contig][int(antiDict[i][j]['tfbsFiltNums'][k])] = antiDict[i][j]['tfbsFilt'][k]
                    for k in range(0, len(mibigDict[i][blastDict[i][j]])):
                        rescueDict[sample][contig][int(mibigDict[i][blastDict[i][j]][k])] = mibigDict2[i][blastDict[i][j]][k]

for i in rescueDict.keys():
    for j in rescueDict[i]:
        clusters = (cluster(sorted(rescueDict[i][j].keys()), 10))
        longestCluster = 0
        for k in clusters:
            if len(k) > longestCluster:
                longestCluster = len(k)
        for k in clusters:
            if len(k) == longestCluster:
                features = []
                for l in k:
                    out3.write(f"{i},{j},{l},{rescueDict[i][j][l]}\n")
                    features.append(rescueDict[i][j][l])
                out4.write(f"{i},{j},{';'.join([str(x) for x in k])},{';'.join(features)}\n")
    out3.write("\n")
out3.close()
out4.close()
