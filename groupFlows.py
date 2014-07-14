
__author__ = 'mjsul_000'

import pysam
import sys

class variation:
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.qual = qual

def groupflow(samfile, vcffile, outfile, minreadsvar=2, minfractvar=0.01, minreadsflow=2, minfractflow=0.01):
    try:
        import pysam
    except:
        return 0
    snps = []
    vcf = open(vcffile)
    sys.stdout.write('Reading vcf file...')
    sys.stdout.flush()
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
            snps.append(aninstance)
    sys.stdout.write(str(len(snps)) + ' snps found.\n')
    sam = pysam.Samfile(samfile, 'rb')
    variantmax = 0
    lastfirreads = set()
    lastsecreads = set()
    firblock = None
    secblock = None
    refblock = None
    blocks = []
    breaks = []
    novar, threevar, recombinants, forked, channelled, properflow, singflow = 0, 0, 0, 0, 0, 0, 0
    recombpos = []
    thecovlist = []
    sys.stdout.write('Finding flows...')
    sys.stdout.flush()
    for snp in snps:
        covcount = 0
        for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1):
            if pileupcolumn.pos == snp.pos - 1:
                varlength = len(snp.ref)
                variants = snp.alt.split(',') + [snp.ref]
                variantreads = [set() for i in range(len(variants))]
                for pileupread in pileupcolumn.pileups:
                    covcount += 1
                    count = 0
                    rvar = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                           pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
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
                        if secflow >= minreadsflow:
                            flowa.append('sec')
                            flowb.append(count)
                    count += 1
                if variations >= variantmax:
                    variantmax = variations
                if variations == 0:
                    novar += 1
                    break
                elif variations == 1:
                    break
                elif variations == 2:
                    pass
                elif variations >= 3:
                    threevar += 1
                    break
                if len(flowa) == 0:
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
                if len(flowa) == 1:
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
                if len(flowa) == 2:
                    if len(set(flowa)) == 1:
                        forked += 1
                        break
                    elif len(set(flowb)) == 1:
                        channelled += 1
                        break
                    else:
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
                            lastsecreads = variantreads[flowb[1]]
                if len(flowa) > 2:
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
    sys.stdout.write(str(thecovlist[len(thecovlist) / 2]) + ' median coverage, ' + str(thecovlist[len(thecovlist) / 2] / 2)\
                     + ' coverage cutoff being used.\nFlows under this coverage will not be used.\n')
    covcut = thecovlist[len(thecovlist) / 2] /2
    dominantreads = set()
    secondaryreads = set()
    out = open(outfile + '.txt', 'w')
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
        for j in range(len(blocks[i][0])):
            firvar = blocks[i][0][j][1]
            secvar = blocks[i][1][j][1]
            blockcov.append(blocks[i][0][j][2] + blocks[i][1][j][2])
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
            if blocks[i][0][-1][0] - blocks[i][0][0][0] <= 50: # if the block is less than 50 bp long
                recombblock = True
                recombblocklist.append(len(blocks[i][0]))
        if not recombblock and max(blockcov) >= covcut:
            ratios.append(len(secreads) * 1.0 / (len(secreads) + len(firstreads)))
            if len(firstreads) * 1.0 / (len(secreads) + len(firstreads)) >= 0.7:
                dominantreads.update(firstreads)
                secondaryreads.update(secreads)
                goodblocks += 1
            elif len(secreads) * 1.0 / (len(secreads) + len(firstreads)) >= 0.7:
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
    domreads = 0
    secreads = 0
    domsam = pysam.Samfile(sys.argv[3] + '.dom.bam', 'wb', template=sam)
    secsam = pysam.Samfile(sys.argv[3] + '.sec.bam', 'wb', template=sam)
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
    pysam.sort(sys.argv[3] + '.dom.bam', sys.argv[3] + '.dom.sorted')
    pysam.index(sys.argv[3] + '.dom.sorted.bam')
    pysam.sort(sys.argv[3] + '.sec.bam', sys.argv[3] + '.sec.sorted')
    pysam.index(sys.argv[3] + '.sec.sorted.bam')

def callStrains(bamfile, outname, mincov=5, mincovorig=70):
    sys.stdout.write('Creating fasta files..\n')
    origbam = pysam.Samfile(bamfile, 'rb')
    dombam = pysam.Samfile(outname + '.dom.sorted.bam', 'rb')
    secbam = pysam.Samfile(outname + '.sec.sorted.bam', 'rb')
    reflen = int(origbam.lengths[0])
    refname = origbam.references[0]
    domlist = ['n' for i in range(reflen)]
    seclist = ['n' for i in range(reflen)]
    varbasedom, varbasesec = 0, 0
    for pileupcolumn in origbam.pileup(refname):
        varlist = []
        for pileupread in pileupcolumn.pileups:
            varlist.append(pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                               pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + 1)])
        coverage = len(varlist)
        if len(varlist) >= mincovorig:
            dombase = max(set(varlist), key=varlist.count)
            dombasecount = varlist.count(dombase)
            ratio = dombasecount * 1.0 / coverage
            if ratio >= 0.95 or dombasecount >= coverage - 1:
                domlist[pileupcolumn.pos] = dombase
    for pileupcolumn in dombam.pileup(refname):
        varlist = []
        for pileupread in pileupcolumn.pileups:
            varlist.append(pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                               pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + 1)])
        coverage = len(varlist)
        if len(varlist) >= mincov:
            dombase = max(set(varlist), key=varlist.count)
            dombasecount = varlist.count(dombase)
            ratio = dombasecount * 1.0 / coverage
            if ratio >= 0.95 or dombasecount >= coverage - 1:
                domlist[pileupcolumn.pos] = dombase
            else:
                domlist[pileupcolumn.pos] = 'n'
                varbasedom += 1
    print set(domlist)
    print domlist.count('n')
    for pileupcolumn in secbam.pileup(refname):
        varlist = []
        for pileupread in pileupcolumn.pileups:
            varlist.append(pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                               pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + 1)])
        coverage = len(varlist)
        if len(varlist) >= mincov:
            dombase = max(set(varlist), key=varlist.count)
            dombasecount = varlist.count(dombase)
            ratio = dombasecount * 1.0 / coverage
            if ratio >= 0.95 or dombasecount >= coverage - 1:
                seclist[pileupcolumn.pos] = dombase
    domlen, domcount, seclen, seccount = 0, 0, 0, 0
    domseq = ''.join(domlist)
    domseq = domseq.split('n')
    out = open(outname + '.dom.fa', 'w')
    for i in domseq:
        if len(i) >= 100:
            domlen += len(i)
            domcount += 1
            out.write('>dom_' + str(domcount) + '\n')
            for j in range(0, len(i), 60):
                out.write(i[j:j+60] + '\n')
    out.close()
    secseq = ''.join(seclist)
    secseq = secseq.split('n')
    out = open(outname + '.sec.fa', 'w')
    for i in secseq:
        if len(i) >= 100:
            seclen += len(i)
            seccount += 1
            out.write('>dom_' + str(seccount) + '\n')
            for j in range(0, len(i), 60):
                out.write(i[j:j+60] + '\n')
    sys.stdout.write('Dominant strain has ' + str(domcount) + ' sequences totaling ' + str(domlen) + ' base pairs.\n')
    sys.stdout.write(str(varbasedom) + ' sites were removed due to variability at the site.\n')
    sys.stdout.write('Secondary strain has ' + str(seccount) + ' sequences totaling ' + str(seclen) + ' base pairs.\n')
    sys.stdout.write(str(varbasesec) + ' sites were removed due to variability at the site.\n')
    sys.stdout.write('groupFlows.py finished.')




mincov = 5
minratio = 0.95


if len(sys.argv) < 4 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    sys.stdout.write('''
groupFlows.py
Written by Mitchell Sullivan (mjsull@gmail.com)
LICENCE: GPLv3
USAGE: python groupFlows.py <args> bam_file vcf_file outfile_prefix
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
returned = groupflow(sys.argv[-3], sys.argv[-2], sys.argv[-1])
if returned == 0:
     sys.stderr.write('Pysam not found, please install.')
     sys.exit()
callStrains(sys.argv[-3], sys.argv[-1])