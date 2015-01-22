__author__ = 'mjsull'
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox
import random


try:
    import pysam
except:
    pass
import sys

class variation:
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.qual = qual


class App:
    def __init__(self, master):
        self.menubar = Menu(master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Create flow file", command=self.create_flow)
        self.filemenu.add_command(label="Load flow file", command=self.load_flow)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.quit)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.toolmenu = Menu(self.menubar, tearoff=0)
        self.toolmenu.add_command(label="Goto base", command=self.goto_base)
        self.toolmenu.add_command(label="Create image", command=self.create_image)
        self.menubar.add_cascade(label="Tools", menu=self.toolmenu)
        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="About", command=self.about)
        self.helpmenu.add_command(label="Help", command=self.help)
        self.helpmenu.add_command(label="Citing Haploflow", command=self.cite)
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)
        master.config(menu=self.menubar)
        self.currxscroll = 1000
        self.curryscroll = 5000
        self.fontsize = 10 # When zooming in/out font does scale, we create a custom font and then change the font size when we zoom
        self.customFont = tkFont.Font(family="Courier", size=self.fontsize)
        self.mainframe = Frame(master)
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure(0, weight=1)
        self.mainframe.grid(row=0, column=0, sticky=NSEW)
        xscrollbar = Scrollbar(self.mainframe, orient=HORIZONTAL)
        xscrollbar.grid(row=1, column=0, sticky=E+W)
        yscrollbar = Scrollbar(self.mainframe)
        yscrollbar.grid(row=0, column=1, sticky=N+S)
        self.canvas = Canvas(self.mainframe, bd=0, bg='#FFFAF0', scrollregion=(0, 0, self.currxscroll, self.curryscroll),
                        xscrollcommand=xscrollbar.set,
                        yscrollcommand=yscrollbar.set)
        self.mainframe.grid_rowconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(0, weight=1)
        self.canvas.grid(row=0, column=0, sticky=N+S+E+W)
        xscrollbar.config(command=self.xscroll)
        yscrollbar.config(command=self.canvas.yview)
        self.reflength = None
        self.ypossnp = 60
        self.yposref = 20
        self.ypos1 = 160
        self.ypos2 = 225
        self.ypos3 = 290
        self.ypos4 = 355
        self.min_flow = 2
        self.block_width = 80
        self.xmod = 20
        self.ymod = 2
        self.gap_size = 8
        self.block_height = 60
        self.flowfile = None
        self.maxvar = 2
        if len(sys.argv) > 1:
            self.flowfile = sys.argv[1]
            self.load_flow()

    def update_frame(self, x2=None):
        x1 = self.canvas.canvasx(0)
        if x2 is None:
            x2 = self.canvas.canvasx(self.canvas.winfo_width())
        #print self.canvas.winfo_width(), x1, x2
        self.canvas.delete(ALL)
        indexa = int(x1 / self.xmod/4)
        indexb = int(x2 / self.xmod/4)
        if indexb >= len(self.flowlist):
            indexb = len(self.flowlist) - 1
        positions = []
        suppositions = []
        for i in range(min([indexa, 0]), indexb + 1):
            for j in self.flowlist[i]:
                if j[0] == 0: # if single forward
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0))
                elif j[0] == 1: # if single reverse
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0))
                elif j[0] == 2: # if pair F F
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][0]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][0])) for val in pair]
                    coords2 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][1])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0))
                    self.canvas.create_line(coords2, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0))
                    self.canvas.create_line((coords1[-2], coords1[-1], coords1[-2] + 1 * self.xmod, coords1[-1], coords2[0] - 1 * self.xmod, coords2[1], coords2[0], coords2[1]), smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2))
                elif j[0] == 3: # if pair F R
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][0]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][0])) for val in pair]
                    coords2 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][1])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0))
                    self.canvas.create_line(coords2, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0))
                    self.canvas.create_line((coords1[-2], coords1[-1], coords1[-2] + 1 * self.xmod, coords1[-1], coords2[0] - 1 * self.xmod, coords2[1], coords2[0], coords2[1]), smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2))
                elif j[0] == 4: # if pair R F
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][0]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][0])) for val in pair]
                    coords2 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][1])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0))
                    self.canvas.create_line(coords2, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0))
                    self.canvas.create_line((coords1[-2], coords1[-1], coords1[-2] + 1 * self.xmod, coords1[-1], coords2[0] - 1 * self.xmod, coords2[1], coords2[0], coords2[1]), smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2))
                elif j[0] == 5: # if pair R R
                    coords1 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][0]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][0])) for val in pair]
                    coords2 = [val for pair in zip(map(lambda x: int(x * self.xmod), j[1][1]),map(lambda x: int(x * self.ymod + self.ypos1), j[2][1])) for val in pair]
                    self.canvas.create_line(coords1, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0))
                    self.canvas.create_line(coords2, smooth=True, width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0))
                    self.canvas.create_line((coords1[-2], coords1[-1], coords1[-2] + 1 * self.xmod, coords1[-1], coords2[0] - 1 * self.xmod, coords2[1], coords2[0], coords2[1]), smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2))
            xpos = (i + 1) * self.xmod * 4
            if i >= indexa:
                suppositions.append(xpos)
                positions.append(self.poslist[i][0])
                thetext = str(self.poslist[i][0]) + '\n' + '\n'.join(self.poslist[i][1:])
                self.canvas.create_text(xpos + 4, self.ymod * self.maxcount * (self.maxvar + 2) + self.ypos1 + 5, anchor=NW, text=thetext, font=self.customFont, tags='top')
        self.canvas.create_rectangle(x1+10, self.ypossnp + 20, x2 - 10, self.ypossnp, tags='top', fill='#E1974C')
        #positions = positions[1:]
        #suppositions = suppositions[1:]
        #print positions
        if len(positions) > 1:
            for i in range(len(positions)):
                xpos = x1+15 + int((positions[i] - positions[0]) * 1.0 / (positions[-1] - positions[0]) * (x2 - x1 - 30))
                self.canvas.create_line(suppositions[i], self.ymod * self.maxcount * (self.maxvar + 2) + self.ypos1, suppositions[i], self.ypos1 - 5, xpos, self.ypossnp + 20, xpos, self.ypossnp, tags='top', width=2)
            self.canvas.create_rectangle(x1 + 10, self.yposref + 20, x2 - 10, self.yposref, tags='top', fill='#7293CB')
            starto = x1 + 10 + int(positions[0] * 1.0 / self.reflength * (x2 - x1 - 20))
            endo = x1 + 10 + int(positions[-1] * 1.0 / self.reflength * (x2 - x1 - 20))
            self.canvas.create_rectangle(starto, self.yposref + 20, endo, self.yposref, tags='top', fill='#E1974C')
            self.canvas.create_text(x1 + 10, self.ypossnp - 2, anchor=SW, text='SNP block length: ' + str(positions[-1] - positions[0]), font=self.customFont, tags='top')

    def create_flow(self):
        pass

    def quit(self):
        root.quit()

    def goto_base(self):
        base = tkSimpleDialog.askint('Goto base', 'Base number')
        pass

    def create_image(self):
        try:
            saveas = tkFileDialog.asksaveasfilename(parent=root)
        except IOError:
            tkMessageBox.showerror('File not valid', 'Please choose another file.')
        self.canvas.postscript(file=saveas, colormode='color')

    def about(self):
        pass

    def help(self):
        pass

    def cite(self):
        pass

    def load_flow(self):
        if self.flowfile is None:
            filename = tkFileDialog.askopenfilename()
            if filename == '' or filename == ():
                return
            self.flowfile = filename
        self.poslist = []
        self.flowlist = []
        self.verticalgap = 40
        flowfile = open(self.flowfile)
        lf = None
        initiallist = []
        for line in flowfile:
            if line.startswith('V'):
                self.poslist.append([int(line.split()[1])] + line.split()[2].split(','))
            elif line.startswith('F'):
                pos = int(line.split()[1].split(',')[0])
                count = int(line.split()[1].split(',')[-1])
                flow = line.split()[1].split(',')[1:-1]
                if lf != pos:
                    while self.poslist[len(initiallist)][0] != pos:
                        initiallist.append(None)
                    initiallist.append([[count, flow]])
                    lf = pos
                else:
                    initiallist[-1].append([count, flow])
            elif line.startswith('I'):
                maxvar, self.medcount, self.maxcount = map(int, line.split()[1].split(','))
                if self.maxvar is None:
                    self.maxvar = maxvar
        stacks = [[i for i in range(0, int(self.maxcount) * (self.maxvar + 1) + 1, int(self.maxcount))] for i in range(len(initiallist))]
        for i in range(len(initiallist)):
            self.flowlist.append([])
            if not initiallist[i] is None:
                for q in initiallist[i]:
                    thedirs = ''
                    x1s = []
                    x2s = []
                    y1s = []
                    y2s = []
                    count = 0
                    start = False
                    start2 = False
                    for j in q[1]:
                        if j == '+s':
                            thedirs += '+'
                            start2 = True
                        elif j == '+':
                            start = True
                            if len(thedirs) == 0:
                                x1s.append((i + count) * 4 + 3)
                            else:
                                x2s.append((i + count) * 4 + 3)
                            thedirs += '+'
                        elif j == 'e':
                            if len(thedirs) == 1:
                                x1s = x1s[:-1]
                                y1s = y1s[:-1]
                            else:
                                x2s = x2s[:-1]
                                y2s = y2s[:-1]
                        elif j == '-s':
                            thedirs += '-'
                            start2 = True
                        elif j == '-':
                            start = True
                            if len(thedirs) == 0:
                                x1s.append((i + count) * 4 + 3)
                            else:
                                x2s.append((i + count) * 4 + 3)
                            thedirs += '-'
                        elif j == 'r' and len(thedirs) == 1:
                            y1s.append(stacks[i + count][0] + q[0]/2)
                            y1s.append(stacks[i + count][0] + q[0]/2)
                            if start:
                                y1s.append(stacks[i + count][0] + q[0]/2)
                                start = False
                            elif not start2:
                                y1s.append(stacks[i + count][0] + q[0]/2)
                                x1s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x1s.append((i+count + 1) * 4)
                            x1s.append((i+count + 1) * 4 + 2)
                            stacks[i + count][0] += q[0]
                            count += 1
                        elif j == 'r' and len(thedirs) == 2:
                            y2s.append(stacks[i + count][0] + q[0]/2)
                            y2s.append(stacks[i + count][0] + q[0]/2)
                            if start:
                                y2s.append(stacks[i + count][0] + 0 + q[0]/2)
                                start = False
                            elif not start2:
                                y2s.append(stacks[i + count][0] + q[0]/2)
                                x2s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x2s.append((i+count+1) * 4)
                            x2s.append((i+count+1) * 4 + 2)
                            stacks[i + count][0] += q[0]
                            count += 1
                        elif j == '_':
                            count += 1
                        elif (j == 'x' or int(j) >= self.maxvar) and len(thedirs) == 1:
                            y1s.append(stacks[i + count][-1] + q[0]/2)
                            y1s.append(stacks[i + count][-1] + q[0]/2)
                            if start:
                                y1s.append(stacks[i + count][-1] + q[0]/2)
                                start = False
                            elif not start2:
                                y1s.append(stacks[i + count][-1] + q[0]/2)
                                x1s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x1s.append((i+count+1) * 4)
                            x1s.append((i+count+1) * 4 + 2)
                            stacks[i + count][-1] += q[0]
                            count += 1
                        elif (j == 'x' or int(j) >= self.maxvar) and len(thedirs) == 2:
                            y2s.append(stacks[i + count][-1] + q[0]/2)
                            y2s.append(stacks[i + count][-1] + q[0]/2)
                            if start:
                                y2s.append(stacks[i + count][-1] + q[0]/2)
                                start = False
                            elif not start2:
                                y2s.append(stacks[i + count][-1] + q[0]/2)
                                x2s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x2s.append((i+count+1) * 4)
                            x2s.append((i+count+1) * 4 + 2)
                            stacks[i + count][-1] += q[0]
                            count += 1
                        elif len(thedirs) == 1:
                            y1s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                            y1s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                            if start:
                                y1s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                                start = False
                            elif not start2:
                                y1s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                                x1s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x1s.append((i+count+1) * 4)
                            x1s.append((i+count+1) * 4 + 2)
                            stacks[i + count][int(j) + 1] += q[0]
                            count += 1
                        else:
                            y2s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                            y2s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                            if start:
                                y2s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                                start = False
                            elif not start2:
                                y2s.append(stacks[i + count][int(j) + 1] + q[0]/2)
                                x2s.append((i+count + 1) * 4 - 2)
                            else:
                                start2 = False
                            x2s.append((i+count+1) * 4)
                            x2s.append((i+count+1) * 4 + 2)
                            stacks[i + count][int(j) + 1] += q[0]
                            count += 1
                    rc = lambda: random.randint(100,255)
                    color = '#%02X%02X%02X' % (rc(), rc(), rc())
                    if x1s[-1] % 4 == 2:
                        x1s[-1] = x1s[-1] - 1
                    if len(thedirs) == 2 and x2s[-1] % 4 == 2:
                        x2s[-1] = x2s[-1] - 1
                    if thedirs == '+':
                        self.flowlist[-1].append([0, x1s, y1s, q[0], color])
                    elif thedirs == '-':
                        self.flowlist[-1].append([1, x1s, y1s, q[0], color])
                    elif thedirs == '++':
                        self.flowlist[-1].append([2, (x1s, x2s), (y1s, y2s), q[0], color])
                    elif thedirs == '+-':
                        self.flowlist[-1].append([3, (x1s, x2s), (y1s, y2s), q[0], color])
                    elif thedirs == '-+':
                        self.flowlist[-1].append([4, (x1s, x2s), (y1s, y2s), q[0], color])
                    elif thedirs == '--':
                        self.flowlist[-1].append([5, (x1s, x2s), (y1s, y2s), q[0], color])
        self.currxscroll = len(self.poslist) * self.xmod * 4 + 100
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        if self.reflength is None:
            self.reflength = self.poslist[-1][0]
        self.update_frame(600)



    def xscroll(self, stuff, stuff2, stuff3=None):
        if not stuff3 is None:
            self.canvas.xview_scroll(stuff2, 'units')
        else:
            self.canvas.xview_moveto(stuff2)
        self.update_frame()



def getflow(samfile, vcffile, outfile, maxdist=1000):
    try:
        import pysam
    except:
        return 0
    snps = []
    vcf = open(vcffile)
    out = open(outfile, 'w')
    depthcount = []
    maxvar = 0
    for line in vcf:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filt, info, form, unknown = line.split()
            aninstance = variation(chrom, pos, ref, alt, qual)
            for i in info.split(';'):
                if i.startswith('DP='):
                    aninstance.depth = int(i.split('=')[1])
                    depthcount.append(int(i.split('=')[1]))
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
            out.write('V ' + str(pos) + ' ' + ref + ',' + alt + '\n')
            snps.append(aninstance)
    sam = pysam.Samfile(samfile, 'rb')
    reads = {}
    for snp in snps:
        newflows = {}
        removereads = []
        for i in reads:
            if reads[i][0] < snp.pos - maxdist:
                rstring = str(reads[i][0]) + ',' + ','.join(reads[i][2:])
                if rstring in newflows:
                    newflows[rstring] += 1
                else:
                    newflows[rstring] = 1
                removereads.append(i)
        for i in removereads:
            del reads[i]
        newflows2 = []
        for i in newflows:
            newflows2.append(i.strip(',_') + ',' + str(newflows[i]))
        newflows2.sort(key=lambda x: int(x.split(',')[-1]), reverse=True)
        newflows2.sort(key=lambda x: int(x.split(',')[0]))
        for i in newflows2:
            out.write('F ' + i + '\n')
        for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1):
            if pileupcolumn.pos == snp.pos - 1:
                vardict = {}
                count = 0
                for i in snp.alt.split(','):
                    vardict[i] = str(count)
                    count += 1
                    if len(vardict) > maxvar:
                        maxvar = len(vardict)
                varlength = len(snp.ref)
                gottenreads = set() # prevent dovetailed paired-end reads recording two variants at a single position
                for pileupread in pileupcolumn.pileups:
                    readname = pileupread.alignment.query_name
                    if not readname in gottenreads: # ignore the start of the second pair in dovetailed reads, there might be a better solution to this but would take a long time to implement and wouldn't provide much additional clarity.
                        rvar = pileupread.alignment.seq[pileupread.query_position:pileupread.query_position +
                               pileupread.alignment.get_overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
                        if readname in reads:
                            if reads[readname][1] != pileupread.alignment.is_read1:
                                if pileupread.alignment.is_reverse and pileupread.query_position == 0:
                                    reads[readname].append('-s')
                                elif pileupread.alignment.is_reverse:
                                    reads[readname].append('-')
                                elif pileupread.query_position == 0:
                                    reads[readname].append('+s')
                                else:
                                    reads[readname].append('+')
                                reads[readname][1] = pileupread.alignment.is_read1
                        else:
                            reads[readname] = [snp.pos, pileupread.alignment.is_read1]
                            if pileupread.alignment.is_reverse and pileupread.query_position == 0:
                                reads[readname].append('-s')
                            elif pileupread.alignment.is_reverse:
                                reads[readname].append('-')
                            elif pileupread.query_position == 0:
                                reads[readname].append('+s')
                            else:
                                reads[readname].append('+')
                        if rvar == snp.ref:
                            reads[readname].append('r')
                        elif rvar in vardict:
                            reads[readname].append(vardict[rvar])
                        else:
                            reads[readname].append('x')
                        gottenreads.add(readname)
                        if pileupread.query_position == pileupread.alignment.query_alignment_end:
                            reads[readname].append('e')
                for i in reads:
                    if not i in gottenreads:
                        reads[i].append('_')
    newflows = {}
    for i in reads:
        rstring = str(reads[i][0]) + ',' + ','.join(reads[i][2:])
        if rstring in newflows:
            newflows[rstring] += 1
        else:
            newflows[rstring] = 1
        removereads.append(i)
    newflows2 = []
    for i in newflows:
        newflows2.append(i.strip(',_') + ',' + str(newflows[i]))
    newflows2.sort(key=lambda x: int(x.split(',')[-1]), reverse=True)
    newflows2.sort(key=lambda x: int(x.split(',')[0]))
    for i in newflows2:
        out.write('F ' + i + '\n')
    depthcount.sort()
    out.write('I ' + str(maxvar) + ',' + str(depthcount[3*len(depthcount)/4]) + ',' + str(depthcount[-1]) + '\n')
    out.close()






if len(sys.argv) > 1 and sys.argv[1] == '-cl':
    getflow(sys.argv[2], sys.argv[3], sys.argv[4])
else:
    root = Tk()
    root.title('HaploFlow')
    root.option_add('*Font', 'Courier 10')
    root.option_add("*Background", "#E0E0FF")
    root.option_add("*Foreground", "#2B3856")
    root.option_add("*Listbox.Background", '#FFFFFF')
    root.option_add("*Scrollbar.Background", "#C0C0FF")
    root.option_add("*Entry.Background", "#FFFFFF")
    root.geometry("600x400")
    app = App(root)
    root.mainloop()