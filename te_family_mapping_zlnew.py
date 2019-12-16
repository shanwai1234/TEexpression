import sys

fh = open(sys.argv[1], 'r')  # read sam file output from htseq

mdict = {}

for line in fh:
    if 'not_aligned' in line:
        continue  # removing reads that not aligned to reference genome
    if 'too_low_aQual' in line:
        continue  # removing reads that too low qual to reference genome
#    if 'no_feature' in line:
#        continue  # removing reads that no feature to reference genome
    if 'alignment_not_unique' in line:
        print ("Please editting NH:i in the input sam file! It is essential for calculating TE family based expression!")
        break  # removing reads that no feature to reference genome
    new = line.strip().split('\t')
    x = new[0].split('.')
    read = x[0] + '.' + x[1]
    if read not in mdict:
        mdict[read] = []
    if 'ambiguous' in new[-1]:
        ann = new[-1].split('ambiguous')[-1].replace('[', '').replace(']', '')
    else:
        if 'gene' not in new[-1]:
            ann = new[-1].split(':')[-1]
        else:
            y = new[-1].split(':')
            ann = y[-2] + ':' + y[-1]
    count = new[-2].split(':')[-1]
    rh_ann = ann + '|' + count
    mdict[read].append(rh_ann)  # storing all labels for each pair of read, i.e: gene:Zm00001g23594|20 or RLX00001Zm0001|3

fh.close()

# with open('temp-sam.pickle', 'wb') as out:
#    pickle.dump(mdict, out)

# with open('temp-sam.pickle', 'rb') as handle:
#    mdict = pickle.load(handle)

# notice that all reads mapping to te are restricted to family level but not individual level

sum_read = 0  # total read counter for future control
uniq_nd = 0  # counter for reads that are unique-mapping but are not defined as gene or te
multi_nd = 0  # reads with no feature define and NH>1
uniq_gene = 0  # read A in Sarah G3 paper
uniq_gene_te = 0  # read B in Sarah G3 paper
multi_gene_te = 0  # read C in Sarah G3 paper
multi_te = 0  # read D in Sarah G3 paper
uniq_te = 0  # read E in Sarah G3 paper


def single(mdict):

    global uniq_nd, multi_nd, uniq_gene, uniq_gene_te, multi_gene_te, multi_te, uniq_te

    te = 0
    fte = 0
    gene = 0
    gname = ''
    tname = ''
    fname = ''

    name = list(mdict.keys())
    item = [mdict[x] for x in name][0]

    try:
        if len(item) == 1:
            label = item[0]
            if '+' in label:
                num = int(label.split('|')[1])
                if num == 1:
                    uniq_gene_te += 1
                else:
                    multi_gene_te += 1
            else:
                if "no_feature" in label:
                    num = int(label.split('|')[1])
                    if num == 1:
                        uniq_nd += 1
                    else:
                        multi_nd += 1
                elif 'gene' in label:
                    num = int(label.split('|')[1])
                    gname = label.split(':')[1].split('|')[0]
                    if num == 1:
                        uniq_gene += 1
                        gene += 1
                else:
                    num = int(label.split('|')[1])
                    name = label.split('|')[0]
                    tname = name
                    fname = name[:8]
                    if num == 1:
                        uniq_te += 1
                        te += 1  # unique hits could only be added to single TE
                        fte += 1
                    elif num > 1:
                        multi_te += 1
                        fte += 1  # multiple hits could be added for one TE family

    except ValueError:
        print ("This (paired) read(s) contains multiple labels!\n")

    return gname, gene, tname, fname, te, fte


def double(mdict):
    '''
    priority for double strand, gene/te > ambiguous > no_feature
    '''
    global uniq_nd, multi_nd, uniq_gene, uniq_gene_te, multi_gene_te, multi_te, uniq_te

    te = 0
    fte = 0
    gene = 0
    gname = ''
    tname = ''
    fname = ''

    name = list(mdict.keys())
    item = [mdict[x] for x in name][0]

    ft = 0
    plus = 0

    try:
        if len(item) == 2:
            for x in item:
                m1 = x.count('no_feature')
                m2 = x.count('+')
                ft += m1
                plus += m2
            if ft == 2:  # both of strands are no_feature
                label1 = item[0]
                label2 = item[1]
                num1 = int(label1.split('|')[1])
                num2 = int(label2.split('|')[1])
                num = max(num1, num2)
                if num == 1:
                    uniq_nd += 1
                else:
                    multi_nd += 1
            elif ft == 1 and plus == 0:  # only the single strand is no_feature; ignoring this strand
                nlist = [i for i, e in enumerate(item) if 'no_feature' not in e]
                mindex = nlist[0]
                label = item[mindex]
                if 'gene' in label:
                    num = int(label.split('|')[1])
                    gname = label.split(':')[1].split('|')[0]
                    if num == 1:
                        uniq_gene += 1
                        gene += 1
                else:
                    num = int(label.split('|')[1])
                    name = label.split('|')[0]
                    tname = name
                    fname = name[:8]
                    if num == 1:
                        uniq_te += 1
                        te += 1  # unique hits could only be added to single TE
                        fte += 1
                    elif num > 1:
                        multi_te += 1
                        fte += 1  # multiple hits could be added for one TE family
            elif ft == 1 and plus >= 1:  # [no_feature, gene:ZmXXX+DHH0001]
                nlist = [i for i, e in enumerate(item) if 'no_feature' not in e]
                mindex = nlist[0]
                label = item[mindex]
                num = int(label.split('|')[1])
                if num == 1:
                    uniq_gene_te += 1
                else:
                    multi_gene_te += 1
            elif ft == 0 and plus == 0:  # [gene:ZmXXX|1, DHHXXX|1]
                if item[0] == item[1]:
                    label = item[0]
                    if 'gene' in item[0]:
                        gname = label.split(':')[1].split('|')[0]
                        num = int(label.split('|')[1])
                        if num == 1:
                            uniq_gene += 1
                            gene += 1
                    else:
                        name = label.split('|')[0]
                        num = int(label.split('|')[1])
                        tname = name
                        fname = name[:8]
                        if num == 1:
                            uniq_te += 1
                            te += 1  # unique hits could only be added to single TE
                            fte += 1
                        elif num > 1:
                            multi_te += 1
                            fte += 1  # multiple hits could be added for one TE family
                else:
                    num1 = int(item[0].split('|')[1])
                    num2 = int(item[1].split('|')[1])
                    num = max(num1, num2)
                    if num == 1:
                        uniq_gene_te += 1
                    else:
                        multi_gene_te += 1
            elif ft == 0 and plus >= 1:
                nlist = [i for i, e in enumerate(item) if '+' not in e]
                if len(nlist) != 0:
                    mindex = nlist[0]  # find '+' not in the item
                    label = item[mindex]
                    if 'gene' in label:
                        gname = label.split(':')[1].split('|')[0]
                        num = int(label.split('|')[1])
                        if num == 1:
                            uniq_gene += 1
                            gene += 1
                    else:
                        name = label.split('|')[0]
                        num = int(label.split('|')[1])
                        tname = name
                        fname = name[:8]
                        if num == 1:
                            uniq_te += 1
                            te += 1  # unique hits could only be added to single TE
                            fte += 1
                        elif num > 1:
                            multi_te += 1
                            fte += 1  # multiple hits could be added for one TE family
                else:  # both of them are ambiguous
                    num1 = int(item[0].split('|')[1])
                    num2 = int(item[1].split('|')[1])
                    num = max(num1, num2)
                    if num == 1:
                        uniq_gene_te += 1
                    else:
                        multi_gene_te += 1

    except ValueError:
        print ("This (paired) read(s) does not contain multiple labels!\n")

    return gname, gene, tname, fname, te, fte


sample = sys.argv[2]

out1 = open('{0}_library_count_summary.txt'.format(sample), 'w')
out2 = open('{0}_gene_count.txt'.format(sample), 'w')
out3 = open('{0}_te_count.txt'.format(sample), 'w')
out4 = open('{0}_te_family_count.txt'.format(sample), 'w')

gdict = {}
tdict = {}
fdict = {}

for x in sorted(mdict):

    sum_read += 1  # total read add 1

    ndict = {}
    if x not in ndict:
        ndict[x] = mdict[x]
    print (mdict[x])
    if len(mdict[x]) == 1:
        gene_name, gene_count, te_name, fam_name, te_count, fte_count = single(ndict)  # mdict['SRRXXX'] = ['gene:Zm000|20']
        if gene_count == 1 and te_count == 0 and fte_count == 0:
            if gene_name not in gdict:
                gdict[gene_name] = 0
            gdict[gene_name] += 1
        elif gene_count == 0 and te_count == 1 and fte_count == 1:
            if te_name not in tdict:
                tdict[te_name] = 0
            if fam_name not in fdict:
                fdict[fam_name] = 0
            tdict[te_name] += 1
            fdict[fam_name] += 1
        elif gene_count == 0 and te_count == 0 and fte_count == 1:
            if fam_name not in fdict:
                fdict[fam_name] = 0
            fdict[fam_name] += 1

    elif len(mdict[x]) == 2:
        gene_name, gene_count, te_name, fam_name, te_count, fte_count = double(ndict)  # mdict['SRRXXX'] = ['gene:Zm000|20']
        if gene_count == 1 and te_count == 0 and fte_count == 0:
            if gene_name not in gdict:
                gdict[gene_name] = 0
            gdict[gene_name] += 1
        elif gene_count == 0 and te_count == 1 and fte_count == 1:
            if te_name not in tdict:
                tdict[te_name] = 0
            if fam_name not in fdict:
                fdict[fam_name] = 0
            tdict[te_name] += 1
            fdict[fam_name] += 1
        elif gene_count == 0 and te_count == 0 and fte_count == 1:
            if fam_name not in fdict:
                fdict[fam_name] = 0
            fdict[fam_name] += 1

out1.write(str(sum_read) + '\t' + str(uniq_nd) + '\t' + str(multi_nd) + '\t' + str(uniq_gene) + '\t' + str(uniq_te) + '\t' + str(uniq_gene_te) + '\t' + str(multi_te) + '\t' + str(multi_gene_te) + '\n')
for i in sorted(gdict):
    out2.write(i + '\t' + str(gdict[i]) + '\n')
for i in sorted(tdict):
    out3.write(i + '\t' + str(tdict[i]) + '\n')
for i in sorted(fdict):
    out4.write(i + '\t' + str(fdict[i]) + '\n')
out1.close()
out2.close()
out3.close()
out4.close()
