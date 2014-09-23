
__author__ = 'mjsul_000'

import pysam
import sys
import subprocess


class variation:
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.qual = qual

def groupflow(samfile, vcffile, outfile, minreadsvar=2, minreadsflow=2, minsnpqual=0):
    try:
        import pysam
    except:
        return 0
    snps = []
    vcf = open(vcffile)
    sys.stdout.write('Reading vcf file...')
    sys.stdout.flush()
    # get variants from VCF file
    for line in vcf:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filt, info, form, unknown = line.split()
            aninstance = variation(chrom, pos, ref, alt, qual)
            for i in info.split(';'):
                if i.startswith('DP='):
                    aninstance.depth = int(i.split('=')[1])
                if i.startswith('AO='):
                    if ',' in i:
                        aninstance.altcount = int(i.split('=')[1].split(',')[0])
                    else:
                        aninstance.altcount = int(i.split('=')[1])
                if i.startswith('RO='):
                    aninstance.refcount = int(i.split('=')[1])
                if i.startswith('AB='):
                    if ',' in i:
                        aninstance.altrat = float(i.split('=')[1].split(',')[0])
                    else:
                        aninstance.altrat = float(i.split('=')[1])
            if aninstance.qual >= minsnpqual:
                snps.append(aninstance)
    sys.stdout.write(str(len(snps)) + ' snps found.\n')
    sam = pysam.Samfile(samfile, 'rb')
    variantmax = 0
    lastfirreads = set()
    lastsecreads = set()
    firblock = None
    secblock = None
    refblock = None
    printcount = 0
    blocks = []
    breaks = []
    novar, threevar, recombinants, forked, channelled, properflow, singflow = 0, 0, 0, 0, 0, 0, 0
    recombpos = []
    thecovlist = []
    sys.stdout.write('Finding flows...')
    sys.stdout.flush()
    lastthree = False
    depths = []
    for pileupcolumn in sam.pileup():
        depths.append(pileupcolumn.n)
    mediancov = depths[len(depths)/2]
    global halfmed, onehalfmed
    halfmed = mediancov / 2
    onehalfmed = mediancov * 3 / 2
    if halfmed < 5:
        halfmed = 5
    for snp in snps:
        covcount = 0
        for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1): # Find all pileup columns assosciated with position
            if pileupcolumn.pos == snp.pos - 1: # find the specific pileup column at the position of interest
                varlength = len(snp.ref) # find the amount of pilup columns assosciated with the variation called by freebayes
                variants = snp.alt.split(',') + [snp.ref] # list all potential variants at this site
                variantreads = [set() for i in range(len(variants))] # sets for putting reads assosciated with each variant in
                for pileupread in pileupcolumn.pileups: # iterate over all reads aligning to that base
                    covcount += 1
                    count = 0
                    rvar = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                           pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)] # read sequence at variant position
                    for i in variants:
                        if rvar == i:
                            variantreads[count].add(pileupread.alignment.qname)
                        count += 1
                flowa = []
                flowb = []
                count = 0
                variations = 0
                max1c = 0
                max1 = None
                max2c = 0
                max2 = None
                maxflow = 0
                for i in variantreads:
                    if len(i) >= minreadsvar:
                        if len(i) > max1c:
                            max2c = max1c
                            max2 = max1
                            max1c = len(i)
                            max1 = count
                        elif len(i) > max2c:
                            max2c = len(i)
                            max2 = count
                        variations += 1
                        firflow = len(lastfirreads.intersection(i))
                        if firflow >= minreadsflow:
                            flowa.append('fir')
                            flowb.append(count)
                        secflow = len(lastsecreads.intersection(i))
                        if firflow > maxflow:
                            maxflow = firflow
                        if secflow >= maxflow:
                            maxflow = secflow
                        if secflow >= minreadsflow:
                            flowa.append('sec')
                            flowb.append(count)
                    count += 1
                breakit = False
                if variations >= variantmax:
                    variantmax = variations
                if variations == 0: # if there is anything other than two variants at this site skip to next site
                    novar += 1
                    break
                elif variations == 1:
                    break
                elif variations == 2:
                    if lastthree:
                        breakit = True
                        lastthree = False
                elif variations >= 3:
                    threevar += 1
                    lastthree = True
                    break
                if len(flowa) == 0 or (len(flowa) == 1 and maxflow <= 5):
                    if not firblock is None:
                        blocks.append((firblock, secblock, refblock))
                        breaks.append('gap')
                    base = variants[max1]
                    firblock = [(snp.pos, base, len(variantreads[max1]))]
                    base = variants[max2]
                    secblock = [(snp.pos, base, len(variantreads[max2]))]
                    refblock = [snp.ref]
                    lastfirreads = variantreads[max1]
                    lastsecreads = variantreads[max2]
                elif len(flowa) == 1 and not breakit:
                    singflow += 1
                    if flowa[0] == 'fir':
                        base = variants[flowb[0]]
                        firblock.append((snp.pos, base, len(variantreads[flowb[0]])))
                        if flowb[0] == max1:
                            base = variants[max2]
                            secblock.append((snp.pos, base, len(variantreads[max2])))
                            refblock.append(snp.ref)
                            lastfirreads = variantreads[max1]
                            lastsecreads = variantreads[max2]
                        else:
                            base = variants[max1]
                            secblock.append((snp.pos, base, len(variantreads[max1])))
                            refblock.append(snp.ref)
                            lastfirreads = variantreads[max2]
                            lastsecreads = variantreads[max1]
                    else:
                        base = variants[flowb[0]]
                        secblock.append((snp.pos, base, len(variantreads[flowb[0]])))
                        refblock.append(snp.ref)
                        if flowb[0] == max1:
                            base = variants[max2]
                            firblock.append((snp.pos, base, len(variantreads[max2])))
                            lastfirreads = variantreads[max2]
                            lastsecreads = variantreads[max1]
                        else:
                            base = variants[max1]
                            firblock.append((snp.pos, base, len(variantreads[max1])))
                            lastfirreads = variantreads[max1]
                            lastsecreads = variantreads[max2]
                elif len(flowa) == 2 and not len(set(flowa)) == 1 and not len(set(flowb)) == 1 and not breakit:
                    properflow += 1
                    if flowa[0] == 'fir':
                        base = variants[flowb[0]]
                        firblock.append((snp.pos, base, len(variantreads[flowb[0]])))
                        base = variants[flowb[1]]
                        secblock.append((snp.pos, base, len(variantreads[flowb[1]])))
                        refblock.append(snp.ref)
                        lastfirreads = variantreads[flowb[0]]
                        lastsecreads = variantreads[flowb[1]]
                    else:
                        base = variants[flowb[0]]
                        secblock.append((snp.pos, base, len(variantreads[flowb[0]])))
                        refblock.append(snp.ref)
                        base = variants[flowb[1]]
                        firblock.append((snp.pos, base, len(variantreads[flowb[1]])))
                        lastfirreads = variantreads[flowb[1]]
                        lastsecreads = variantreads[flowb[0]]
                else: # more than 3 flows means recombination is happening at significant levels at this site
                    recombpos.append(snp.pos)
                    recombinants += 1
                    blocks.append((firblock, secblock, refblock))
                    breaks.append('recomb')
                    base = variants[max1]
                    firblock = [(snp.pos, base, len(variantreads[max1]))]
                    base = variants[max2]
                    secblock = [(snp.pos, base, len(variantreads[max2]))]
                    refblock = [snp.ref]
                    lastfirreads = variantreads[max1]
                    lastsecreads = variantreads[max2]
        thecovlist.append(covcount)
    sys.stdout.write(str(singflow + properflow) + ' flows found.\n')
    sys.stdout.write(str(recombinants) + ' recombinants found.\n')
    sys.stdout.write(str(threevar) + ' sites with three variants found.\n')
    sys.stdout.write(str(variantmax) + ' maximum variants at a single site.\n')
    thecovlist.sort()
    sys.stdout.write(str(mediancov) + ' median coverage, ' + str(halfmed) + '-' + str(onehalfmed)\
                     + ' coverage cutoffs being used.\nFlows outside of these coverages will not be used.\n')
    covcut = thecovlist[len(thecovlist) / 2] /2
    dominantreads = set()
    secondaryreads = set()
    out = open(outfile + '.txt', 'w')
    global recombblocklist
    recombblocklist = []
    mids, lowcovblock, goodblocks = 0, 0, 0
    ratios = []
    coveragelist = []
    coverageblocklist = []
    sys.stdout.write('Assigning blocks to strains...')
    sys.stdout.flush()
    for i in range(len(blocks)):
        firstreads, secreads = set(), set()
        blockcov = []
        blockcova = []
        blockcovb = []
        for j in range(len(blocks[i][0])):
            firvar = blocks[i][0][j][1]
            secvar = blocks[i][1][j][1]
            blockcov.append(blocks[i][0][j][2] + blocks[i][1][j][2])
            blockcova.append(blocks[i][0][j][2])
            blockcovb.append(blocks[i][1][j][2])
            coveragelist.append(blocks[i][0][j][2] + blocks[i][1][j][2])
            for pileupcolumn in sam.pileup(snp.chrom, blocks[i][0][j][0], blocks[i][0][j][0] + 1):
                if pileupcolumn.pos == blocks[i][0][j][0] - 1:
                    varlength = len(blocks[i][2][j])
                    for pileupread in pileupcolumn.pileups:
                        rvar = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                               pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
                        if rvar == firvar:
                            firstreads.add(pileupread.alignment.qname)
                        elif rvar == secvar:
                            secreads.add(pileupread.alignment.qname)
            outlist = blocks[i][0][j] + blocks[i][1][j][1:] + (blocks[i][2][j],)
            out.write('\t'.join(map(str, outlist)) + '\n')
        recombblock = False
        coverageblocklist.append(max(blockcov))
        if breaks[i] == 'recomb' and i != 0 and breaks[i-1] == 'recomb':
            if blocks[i][0][-1][0] - blocks[i][0][0][0] <= 100: # if the block is less than 50 bp long
                recombblock = True
                recombblocklist.append((blocks[i][0][0][0], blocks[i][0][-1][0]))
        if not recombblock and max(blockcov) >= covcut:
            highcovs = set()
            for j in range(len(blockcov)):
                if onehalfmed >= blockcov[j] and (blockcova[j] >= halfmed or blockcovb[j] >= halfmed):
                    if blockcova[j] > blockcovb[j]:
                        highcovs.add('a')
                        ratios.append(blockcova[j] * 1.0 / blockcov[j])
                    else:
                        highcovs.add('b')
                        ratios.append(blockcovb[j] * 1.0 / blockcov[j])
            if len(highcovs) > 1:
                print 'whatthe?'
                print blocks[i]
                sys.exit()
            else:
                if 'a' in highcovs:
                    dominantreads.update(firstreads)
                    secondaryreads.update(secreads)
                    goodblocks += 1
                elif 'b' in highcovs:
                    dominantreads.update(secreads)
                    secondaryreads.update(firstreads)
                    goodblocks += 1
                else:
                    mids += 1
        else:
            if not recombblock:
                lowcovblock += 1
        out.write(breaks[i] + '\n')
    sys.stdout.write(str(mids) + ' blocks removed due to indistinguishable coverage between strains.\n')
    sys.stdout.write(str(len(recombblocklist)) + ' blocks under active recombination removed.\n')
    sys.stdout.write(str(lowcovblock) + ' low coverage blocks removed.\n')
    sys.stdout.write(str(goodblocks) + ' blocks used to assign reads.\n')
    for i in range(len(firblock)):
        outlist = firblock[i] + secblock[i][1:] + (refblock[i],)
        out.write('\t'.join(map(str, outlist)) + '\n')
    out.write('end\n')
    out.close()
    totalreads = 0
    domsam = pysam.Samfile(outfile + '.dom.bam', 'wb', template=sam)
    secsam = pysam.Samfile(outfile + '.sec.bam', 'wb', template=sam)
    ratios.sort()
    try:
        print ratios[0], ratios[len(ratios)/4], ratios[len(ratios)/2], ratios[len(ratios) * 3 / 4], ratios[-1]
    except:
        pass
    sys.stdout.write('Writing dominant and secondary strain bam files....\n')
    sys.stdout.flush()
    for read in sam.fetch():
        if read.qname in dominantreads:
            domsam.write(read)
        elif read.qname in secondaryreads:
            secsam.write(read)
        totalreads += 1
    domsam.close()
    secsam.close()
    sys.stdout.write(str(len(dominantreads)) + ' assigned to the dominant strain.\n')
    sys.stdout.write(str(len(secondaryreads)) + ' assigned to the secondary strain.\n')
    sys.stdout.write('Out of a possible ' + str(totalreads) + ' reads.\n')
    sys.stdout.write('Sorting and indexing BAM files.\n')
    pysam.sort(outfile + '.dom.bam', outfile + '.dom.sorted')
    pysam.index(outfile + '.dom.sorted.bam')
    pysam.sort(outfile + '.sec.bam', outfile + '.sec.sorted')
    pysam.index(outfile + '.sec.sorted.bam')

def callStrains(reference, outfile, sam):
    global halfmed, recombblocklist, onehalfmed
    subprocess.Popen('samtools mpileup -uf ' + reference + ' ' + outfile + '.dom.sorted.bam | bcftools view -cg -  > ' + outfile + '.dom.pu', shell=True).wait()
    subprocess.Popen('samtools mpileup -uf ' + reference + ' ' + outfile + '.sec.sorted.bam | bcftools view -cg -  > ' + outfile + '.sec.pu', shell=True).wait()
    samfile = pysam.Samfile(sam)
    domseq = ['n' for i in range(samfile.lengths[0])]
    for pileupcolumn in samfile.pileup():
        varcounts = {}
        coverage = 0
        for pileupread in pileupcolumn.pileups:
            var = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                               pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + 1)]
            coverage += 1
            if var in varcounts:
                varcounts[var] += 1
            else:
                varcounts[var] = 1
        if coverage >= halfmed:
            for i in varcounts:
                if varcounts[i] * 1.0 / coverage >= 0.98:
                    domseq[pileupcolumn.pos] = i
    print len(domseq)
    domfile = open(outfile + '.dom.pu')
    lastpos = 0
    count1, count2, count3, count4 = 0, 0, 0, 0
    for line in domfile:
        if not line.startswith('#'):
            chrom, pos, eyed, ref, alt, qual, zefilter, info = line.split()[:8]
            pos = int(pos)
            if alt == '.':
                alt = ref
            for i in info.split(';'):
                if i.startswith('FQ='):
                    quality = - float(i.split('=')[1])
            if pos == lastpos:
                domseq[pos-1] = 'n'
                count1 += 1
            else:
                if quality >= 30:
                    domseq[pos-1] = alt
            lastpos = pos
    print len(domseq)
    domfile.close()
    secseq = ['n' for i in range(len(domseq))]
    secfile = open(outfile + '.sec.pu')
    lastpos = 0
    for line in secfile:
        if not line.startswith('#'):
            chrom, pos, eyed, ref, alt, qual, zefilter, info = line.split()[:8]
            pos = int(pos)
            if alt == '.':
                alt = ref
            for i in info.split(';'):
                if i.startswith('FQ='):
                    quality = - float(i.split('=')[1])
            if pos == lastpos:
                secseq[pos-1] = 'n'
                count2 += 1
            else:
                if quality >= 30:
                    secseq[pos-1] = alt
                    count3 += 1
                else:
                    count4 += 1
            lastpos = pos
    print len(domseq)
    print count1, count2, count3, count4
    for i in recombblocklist:
        for j in range(i[0] - 1, i[1]):
            domseq[j] = 'n'
            secseq[j] = 'n'
    domseq = ''.join(domseq)
    secseq = ''.join(secseq)
    domseqs = []
    secseqs = []
    domslen = []
    secslen = []
    for i in domseq.split('nnnnnnnnnn'):
        seq = i.strip('n')
        if len(seq) >= 100:
            domseqs.append(seq)
            domslen.append(len(seq))
    for i in secseq.split('nnnnnnnnnn'):
        seq = i.strip('n')
        if len(seq) >= 100:
            secseqs.append(seq)
            secslen.append(len(seq))
    sys.stdout.write('\tMax\tnumber\ttotal\n')
    sys.stdout.write('DOM: ' + str(max(domslen)) + ' ' + str(len(domslen)) + ' ' + str(sum(domslen)) + '\n')
    try:
        sys.stdout.write('SEC: ' + str(max(secslen)) + ' ' + str(len(secslen)) + ' ' + str(sum(secslen)) + '\n')
    except:
        sys.stdout.write('Not enough coverage to create sequence for secondary strain\n ')
    out = open(outfile + '.dom.fa', 'w')
    count = 1
    for i in domseqs:
        out.write('>dom_' + str(count) + '\n')
        for j in range(0, len(i), 60):
            out.write(i[j:j+60] + '\n')
        count += 1
    out.close()
    out = open(outfile + '.sec.fa', 'w')
    count = 1
    for i in secseqs:
        out.write('>sec_' + str(count) + '\n')
        for j in range(0, len(i), 60):
            out.write(i[j:j+60] + '\n')
        count += 1
    out.close()





            







mincov = 5
minratio = 0.9


if len(sys.argv) < 4 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    sys.stdout.write('''
groupFlows.py
Written by Mitchell Sullivan (mjsull@gmail.com)
LICENCE: GPLv3
USAGE: python groupFlows.py <args> bam_file vcf_file ref_file outfile_prefix
Bam file must be ordered and indexed, vcf file should be made using freebayes
Arguments
-m minimum coverage minimum coverage for base to be used in outfile
-r minimum ratio for base to be used in outfile

OUTPUT:
<prefix>.dom.bam               Alignement of reads assosciated with the dominant strain
<prefix>.dom.sorted.bam        sorted alignment
<prefix>.dom.sorted.bam.bai    index of sorted alignment
<prefix>.dom.fa                fasta of aligned regions of reference

<prefix>.sec.bam               Identical files but for the secondary strain
<prefix>.sec.sorted.bam
<prefix>.sec.sorted.bam.bai
<prefix>.sec.fa

<prefix>.txt                   list of flows seperated into blocks
''')
    sys.exit()
returned = groupflow(sys.argv[-4], sys.argv[-3], sys.argv[-1])
if returned == 0:
     sys.stderr.write('Pysam not found, please install.')
     sys.exit()
callStrains(sys.argv[-2], sys.argv[-1], sys.argv[-4])