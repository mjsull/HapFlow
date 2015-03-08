__author__ = 'mjsull'
from Tkinter import *
import threading
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox
import random
import time
import os
import Queue
import platform
import webbrowser

try:
    import pysam
except:
    pass
import sys

class clqueue:
    def put(self, theval):
        sys.stdout.write(theval + '\n')

class variation:
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.qual = qual


class App:
    def __init__(self, master):
        if master is None:
            try:
                maxdist = int(sys.argv[5])
            except:
                maxdist = 1000
            self.getflow(sys.argv[2], sys.argv[3], sys.argv[4], maxdist)
            return
        self.menubar = Menu(master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Create flow file", command=self.create_flow)
        self.filemenu.add_command(label="Load flow file", command=self.get_flow_name)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.quit)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.toolmenu = Menu(self.menubar, tearoff=0)
        self.toolmenu.add_command(label="Goto base", command=self.goto_base)
        self.toolmenu.add_command(label="Create image", command=self.create_image)
        self.menubar.add_cascade(label="Tools", menu=self.toolmenu)
        self.viewmenu = Menu(self.menubar, tearoff=0)
        self.viewmenu.add_command(label="Hide gapped", command=self.hide_gapped)
        self.viewmenu.add_command(label="Show gapped", command=self.show_gapped)
        self.viewmenu.add_command(label="Stretch X", command=self.stretch_x)
        self.viewmenu.add_command(label="Shrink X", command=self.shrink_x)
        self.viewmenu.add_command(label="Stretch Y", command=self.stretch_y)
        self.viewmenu.add_command(label="Shrink Y", command=self.shrink_y)
        self.menubar.add_cascade(label="View", menu=self.viewmenu)
        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="About", command=self.about)
        self.helpmenu.add_command(label="Help", command=self.help)
        self.helpmenu.add_command(label="Support", command=self.support)
        self.helpmenu.add_command(label="Citing", command=self.cite)
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
        self.bamfile = StringVar(value='')
        self.vcffile = StringVar(value='')
        self.flowfile = StringVar(value='')
        self.maxdist = IntVar(value=1000)
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
        self.maxvar = 2
        self.currflows = set()
        self.gapped_state = NORMAL
        if len(sys.argv) == 2:
            self.flowfile.set(sys.argv[1])
            self.load_flow()
        self.queue = Queue.Queue()
        self.rcmenu = Menu(root, tearoff=0)
        self.rcmenu.add_command(label="Details", command=self.details)
        self.rcmenu.add_command(label="Write flow readnames", command=self.write_flow_names)
        self.rcmenu.add_command(label="Write group readnames", command=self.write_group_names)
        if platform.system() in ['Windows', 'Linux']:
            self.canvas.tag_bind('map', '<Button-3>', self.rightClick)
        else:
            self.canvas.tag_bind('map', '<Button-2>', self.rightClick)
        self.canvas.tag_bind('map', '<Button-1>', self.select_flow)
        self.canvas.bind('<Button-1>', self.remove_rc)
        self.selected = [None, None, None]


    def rightClick(self, event):
        self.rctag = self.canvas.gettags(CURRENT)
        self.rcmenu.unpost()
        self.rcmenu.post(event.x_root, event.y_root)
        self.rcpos = (event.x_root, event.y_root)

    def select_flow(self, event):
        pos, num, amap, zecurrent = self.canvas.gettags(CURRENT)
        thecol = self.canvas.itemcget(CURRENT, 'fill')
        if not self.selected[0] is None:
            for i in self.canvas.find_withtag(self.selected[0]):
                if self.canvas.gettags(i)[1] == self.selected[1]:
                   self.canvas.itemconfig(i, fill=self.selected[2])
        if pos == self.selected[0] and num == self.selected[1]:
            self.selected = [None, None, None]
        else:
            for i in self.canvas.find_withtag(pos):
                if self.canvas.gettags(i)[1] == num:
                   self.canvas.itemconfig(i, fill='#000000')
            self.selected = [pos, num, thecol]
        self.rcmenu.unpost()

    def remove_rc(self, event):
        self.rcmenu.unpost()

    def find_flow(self, pos, num):
        flowfile = open(self.flowfile.get())
        count = 0
        for line in flowfile:
            if line.startswith('F'):
                if line[2:].startswith(pos + ','):
                    count += 1
                    if count == num:
                        flowfile.close()
                        return line.split()[1]

    def get_names(self, flows, positions):
        try:
            import pysam
        except ImportError:
            tkMessageBox.showerror('Pysam not installed.', 'This functionality is only available if pysam is isntalled.')
            return
        if self.bamfile.get() == '':
            filename = tkFileDialog.askopenfilename(title='Please select alignment file (BAM)')
            self.bamfile.set(filename)
        sam = pysam.Samfile(self.bamfile.get(), 'rb')
        minsnp = float('inf')
        maxsnp = 0
        for i in range(len(flows)):
            minisnp = int(positions[i])
            maxisnp = int(positions[i]) + len(filter(lambda x: not x in ['+s', '+', '-s', '-', 'e'], flows[i].split(',')[1:-2])) -1
            if minisnp < minsnp:
                minsnp = minisnp
            if maxisnp > maxsnp:
                maxsnp = maxisnp
        snplist = []
        getit = True
        minsnp = max([0, minsnp-1])
        i = minsnp
        while getit:
            position = self.poslist[i][0]
            if position >= self.poslist[maxsnp][0] + self.maxdist.get():
                getit = False
            alt = ''
            for j in self.poslist[i][1:]:
                if j.startswith('*'):
                    ref = j[1:]
                else:
                    alt += ',' + j
            alt = alt[1:]
            aninstance = variation(self.chrom, position, ref, alt, 0)
            snplist.append(aninstance)
            i += 1
        reads = {}
        for snp in snplist: # for each variant in the vcf file
            for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1):
                if pileupcolumn.pos == snp.pos - 1:
                    vardict = {}
                    varorder = []
                    varcount = {}
                    currerr = 0
                    varlength = len(snp.ref)
                    for i in [snp.ref] + snp.alt.split(','):
                        varcount[i] = 0
                    for pileupread in pileupcolumn.pileups:
                        rvar = pileupread.alignment.seq[pileupread.query_position:pileupread.query_position +
                               pileupread.alignment.get_overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
                        if rvar in varcount:
                            varcount[rvar] += 1
                        else:
                            currerr += 1
                    for i in varcount:
                        varorder.append((varcount[i], i))
                    varorder.sort(reverse=True)
                    count = 0
                    for i in varorder:
                        vardict[i[1]] = str(count)
                        count += 1
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
                            if rvar in vardict:
                                reads[readname].append(vardict[rvar])
                            else:
                                reads[readname].append('x')
                            gottenreads.add(readname)
                            if pileupread.query_position == pileupread.alignment.query_alignment_end:
                                reads[readname].append('e')
                    for i in reads:
                        if not i in gottenreads:
                            reads[i].append('_')
        outset = set()
        for i in reads:
            for j in flows:
                if int(j.split(',')[0]) == reads[i][0] and ','.join(reads[i][2:]).strip(',_') == ','.join(j.split(',')[1:-2]):
                    outset.add(i)
        outfile = tkFileDialog.asksaveasfilename(title='Write to BAM')
        if outfile[-4:] != '.bam':
            outfile += '.bam'
        newsam = pysam.Samfile(outfile, 'wb', template=sam)
        for read in sam.fetch():
            if read.qname in outset:
                newsam.write(read)
        newsam.close()
        sam.close()


    def details(self):
        pos, num, amap, securrent = self.rctag
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        try:
            self.detail_window.destroy()
        except:
            pass
        self.detail_window = Toplevel()
        splitline = flow.split(',')
        self.hl1 = Label(self.detail_frame, text='Position:', anchor=E)
        self.hl1.grid(column=0, row=1)
        self.he1 = Entry(self.detail_frame, textvariable=StringVar(value=splitline[0]), state='readonly')
        self.he1.grid(column=1, row=1)
        self.hl2 = Label(self.detail_frame, text='Flow:', anchor=E)
        self.hl2.grid(column=0, row=2)
        self.he2 = Entry(self.detail_frame, textvariable=StringVar(value=', '.join(splitline[1:-2])), state='readonly')
        self.he2.grid(column=1, row=2)
        self.hl3 = Label(self.detail_frame, text='Count:', anchor=E)
        self.hl3.grid(column=0, row=3)
        self.he3 = Entry(self.detail_frame, textvariable=StringVar(value=splitline[-2]), state='readonly')
        self.he3.grid(column=1, row=3)
        self.hl4 = Label(self.detail_frame, text='Group:', anchor=E)
        self.hl4.grid(column=0, row=4)
        self.he4 = Entry(self.detail_frame, textvariable=StringVar(value=splitline[-1]), state='readonly')
        self.he4.grid(column=1, row=4)
        self.detail_frame.grid(padx=5, pady=5)


    def write_flow_names(self):
        pos, num, amap, securrent = self.rctag
        posnum = int(pos[1:])
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        self.get_names([flow], [posnum])

    def get_group(self, group):
        flowfile = open(self.flowfile.get())
        count = 0
        flows, posnums = [], []
        for line in flowfile:
            if line.startswith('F'):
                if line[-len(group)-1:] == ',' + group:
                    flows.append(line.split()[1])
                    posnums.append(count)
                count += 1
        return flows, posnums


    def write_group_names(self):
        pos, num, amap, securrent = self.rctag
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        group = flow.split(',')[-1]
        flows, posnums = self.get_group(group)
        self.get_names(flows, posnums)

    def hsl_to_rgb(self, h, s, l):
        c = (1 - abs(2*l - 1)) * s
        x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
        m = l - c/2
        if h < 60:
            r, g, b = c + m, x + m, 0 + m
        elif h < 120:
            r, g, b = x + m, c+ m, 0 + m
        elif h < 180:
            r, g, b = 0 + m, c + m, x + m
        elif h < 240:
            r, g, b, = 0 + m, x + m, c + m
        elif h < 300:
            r, g, b, = x + m, 0 + m, c + m
        else:
            r, g, b, = c + m, 0 + m, x + m
        r = int(r * 255)
        g = int(g * 255)
        b = int(b * 255)

        return '#%02x%02x%02x' % (r, g, b)


    def update_frame(self, x2=None):
        x1 = self.canvas.canvasx(0)
        if x2 is None:
            x2 = self.canvas.canvasx(self.canvas.winfo_width())
        self.canvas.delete('top')
        indexa = int(x1 / self.xmod/4)
        indexb = int(x2 / self.xmod/4)
        if indexb >= len(self.flowlist):
            indexb = len(self.flowlist) - 1
        positions = []
        suppositions = []
        toremove = set()
        for i in self.currflows:
            if i < indexa -3 or i > indexb + 3:
                self.canvas.delete('p' + str(i))
                toremove.add(i)
        thetime = time.time()
        for i in toremove:
            self.currflows.remove(i)
        thetime = time.time()
        for i in range(max([indexa, 0]), indexb + 1):
            if not i in self.currflows:
                count = 0
                self.currflows.add(i)
                for j in self.flowlist[i]:
                    count += 1
                    if j[0] == 0: # if single forward
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1] - 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - 1) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k] + 1) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 1: # if single reverse
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1] - 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - 1) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k] + 1) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 2: # if pair F F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 3: # if pair F R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 4: # if pair R F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 5: # if pair R R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
            xpos = (i + 1) * self.xmod * 4
            if i >= indexa:
                suppositions.append(xpos)
                positions.append(self.poslist[i][0])
                thetext = str(self.poslist[i][0]) + '\n' + '\n'.join(self.poslist[i][1:])
                self.canvas.create_text(xpos + 4, self.ymod * self.flowend + self.ypos1 + 5, anchor=NW, text=thetext, font=self.customFont, tags='top')
        self.canvas.create_rectangle(x1+10, self.ypossnp + 20, x2 - 10, self.ypossnp, tags='top', fill='#E1974C')
        self.canvas.tag_raise('gapped')
        if len(positions) > 1:
            for i in range(len(positions)):
                xpos = x1+15 + int((positions[i] - positions[0]) * 1.0 / (positions[-1] - positions[0]) * (x2 - x1 - 30))
                self.canvas.create_line(suppositions[i], self.ypos1 + self.ymod * self.flowend, suppositions[i], self.ypos1 - 5, xpos, self.ypossnp + 20, xpos, self.ypossnp, tags='top', width=2)
            self.canvas.create_rectangle(x1 + 10, self.yposref + 20, x2 - 10, self.yposref, tags='top', fill='#7293CB')
            starto = x1 + 10 + int(positions[0] * 1.0 / self.reflength * (x2 - x1 - 20))
            endo = x1 + 10 + int(positions[-1] * 1.0 / self.reflength * (x2 - x1 - 20))
            self.canvas.create_rectangle(starto, self.yposref + 20, endo, self.yposref, tags='top', fill='#E1974C')
            self.canvas.create_text(x1 + 10, self.ypossnp - 2, anchor=SW, text='SNP block start..stop: ' + str(positions[-1]) + '..' + str(positions[0]), font=self.customFont, tags='top')

    def create_flow(self):
        self.create_flow_top = Toplevel()
        self.create_flow_top.grab_set()
        self.create_flow_top.wm_attributes("-topmost", 1)
        self.create_flow_top.geometry('+20+30')
        self.create_flow_top.title('Create Flow')
        self.create_flow_frame = Frame(self.create_flow_top)
        self.bamfilelabel = Label(self.create_flow_frame, text='BAM file:')
        self.bamfilelabel.grid(row=0, column=0, sticky=E)
        self.bamfileentry = Entry(self.create_flow_frame, textvariable=self.bamfile, justify=RIGHT, width=30)
        self.bamfileentry.grid(row=0, column=1)
        self.bamfileentrybutton = Button(self.create_flow_frame, text='...', command=self.loadbam)
        self.bamfileentrybutton.grid(row=0, column=2)
        self.vcffilelabel = Label(self.create_flow_frame, text='VCF file:')
        self.vcffilelabel.grid(row=1, column=0, sticky=E)
        self.vcffileentry = Entry(self.create_flow_frame, textvariable=self.vcffile, justify=RIGHT, width=30)
        self.vcffileentry.grid(row=1, column=1)
        self.vcffileentrybutton = Button(self.create_flow_frame, text='...', command=self.loadvcf)
        self.vcffileentrybutton.grid(row=1, column=2)
        self.flowfilelabel = Label(self.create_flow_frame, text='Output file:')
        self.flowfilelabel.grid(row=2, column=0, sticky=E)
        self.flowfileentry = Entry(self.create_flow_frame, textvariable=self.flowfile, justify=RIGHT, width=30)
        self.flowfileentry.grid(row=2, column=1)
        self.flowfileentrybutton = Button(self.create_flow_frame, text='...', command=self.loadflow)
        self.flowfileentrybutton.grid(row=2, column=2)
        self.maxdistlabel = Label(self.create_flow_frame, text='Max. distance (bp):', width=30, anchor=E)
        self.maxdistlabel.grid(row=3, column=0)
        self.maxdistentry = Entry(self.create_flow_frame, textvariable=self.maxdist)
        self.maxdistentry.grid(row=3, column=1)
        self.okflow = Button(self.create_flow_frame, text='Ok', command=self.ok_flow)
        self.okflow.grid(row=4, column=2, sticky=E)
        self.create_flow_frame.grid(padx=10, pady=10)

    def loadbam(self):
        filename = tkFileDialog.askopenfilename(parent=self.create_flow_top)
        if filename == '':
            return
        self.bamfile.set(filename)


    def loadvcf(self):
        filename = tkFileDialog.askopenfilename(parent=self.create_flow_top)
        if filename == '':
            return
        self.vcffile.set(filename)

    def loadflow(self):
        filename = tkFileDialog.asksaveasfilename(parent=self.create_flow_top)
        if filename == '':
            return
        else:
            self.flowfile.set(filename)

    def ok_flow(self):
        if not os.path.exists(self.bamfile.get()) or not os.path.exists(self.vcffile.get()):
            tkMessageBox.showerror('File missing', 'Please include a contig and read file.', parent=self.create_flow_top)
            return
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        try:
            import pysam
        except:
            tkMessageBox.showerror('Pysam not found',
                                   'Creating a flow file requires pysam, please install.')
            return
        self.create_flow_top.destroy()
        self.run_flow_top = Toplevel()
        self.run_flow_top.grab_set()
        self.run_flow_top.wm_attributes("-topmost", 1)
        self.run_flow_top.geometry('+20+30')
        self.run_flow_top.title('Running flow creation tool')
        self.run_flow_frame = Frame(self.run_flow_top)
        self.consoletext = StringVar(value='Creating Flow File.')
        self.consolelabel = Label(self.run_flow_frame, bg='#FFFF99', relief=SUNKEN, textvariable=self.consoletext, width=35, height=10)
        self.consolelabel.grid(row=0, column=0)
        self.consolebutton = Button(self.run_flow_frame, text='Ok', command=self.ok_console, state=DISABLED)
        self.consolebutton.grid(row=1, column=0, sticky=E)
        self.run_flow_frame.grid(padx=10, pady=10)
        self.thethread = threading.Thread(target=self.getflow)
        self.thethread.start()
        self.update_flow()

    def ok_console(self):
        self.run_flow_top.destroy()

    def update_flow(self):
        self.dot_console()
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text == 'Flow file successfully created.':
                    self.consolebutton.config(state=NORMAL)
                else:
                    root.after(1000, self.update_flow)
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            self.update_console('Flow file creation failed,\n please check console output.')
            self.consolebutton.config(state=NORMAL)
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_flow)
            return

    # push text here to update the console
    def update_console(self, text):
        if self.consoletext.get() == '':
            self.consoletext.set(text)
        else:
            self.consoletext.set(self.consoletext.get().strip('.') + '.\n' + text) # strip excess dots once process has stopped running

    # the . .. ... . .. ... effect in the console
    def dot_console(self):
        text = self.consoletext.get()
        if text[-3:] == '...':
            self.consoletext.set(text[:-2])
        else:
            self.consoletext.set(text + '.')


    def hide_gapped(self):
        self.canvas.itemconfig('gapped', state=HIDDEN)
        self.gapped_state = HIDDEN

    def show_gapped(self):
        self.canvas.itemconfig('gapped', state=NORMAL)
        self.gapped_state = NORMAL

    def stretch_x(self):
        self.xmod = max([int(self.xmod * 1.1), self.xmod + 1])
        self.update_frame()

    def shrink_x(self):
        self.xmod = min([int(self.xmod * 0.9), self.xmod - 1])
        self.xmod = max([self.xmod, 1])
        self.update_frame()

    def stretch_y(self):
        self.ymod = max([int(self.ymod * 1.1), self.ymod + 1])
        self.update_frame()

    def shrink_y(self):
        self.ymod = min([int(self.ymod * 0.9), self.ymod - 1])
        self.ymod = max([self.ymod, 1])
        self.update_frame()

    def quit(self):
        root.quit()

    def goto_base(self):
        base = tkSimpleDialog.askinteger('Goto base', 'Base number')
        for i in range(len(self.poslist)):
            if base < self.poslist[i][0]:
                fraction = (i-1) * 1.0 / len(self.poslist)
                break
        self.canvas.xview_moveto(fraction)
        self.update_frame()

    def create_image(self):
        try:
            saveas = tkFileDialog.asksaveasfilename(parent=root)
        except IOError:
            tkMessageBox.showerror('File not valid', 'Please choose another file.')
        self.canvas.postscript(file=saveas, colormode='color')

    def help(self):
        webbrowser.open_new('https://github.com/mjsull/HaploFlow/blob/master/README.md')

    def about(self):
        try:
            self.helppanel.destroy()
        except:
            pass
        self.helppanel = Toplevel()
        self.helppanel.wm_attributes("-topmost", 1)
        self.helppanel.geometry('+20+30')
        self.helppanel.title('About')
        self.frame7 = Frame(self.helppanel)
        self.about1label = Label(self.frame7, text='HaploFlow', font='TkDefaultFont 13 bold')
        self.about1label.grid(row=0, column=0)
        self.about2label = Label(self.frame7, text='HaploFlow is a Python application for visualising\n\
haplotypes present in sequencing data.\n\n\
Version 0.1\n')
        self.about2label.grid(row=1, column=0)
        self.frame7.grid(padx=10, pady=10)

    def support(self):
        try:
            self.helppanel.destroy()
        except:
            pass
        self.helppanel = Toplevel()
        self.helppanel.wm_attributes("-topmost", 1)
        self.helppanel.geometry('+20+30')
        self.helppanel.title('Support')
        self.frame9 = Frame(self.helppanel)
        self.about1label1 = Label(self.frame9, text='HaploFlow', font='TkDefaultFont 13 bold')
        self.about1label1.grid(row=0, column=0)
        self.supportlabel2 = Label(self.frame9, text='written by Mitchell Sullivan - mjsull@gmail.com\n\
Please do not hesitate to email with issues or bug reports.')
        self.supportlabel2.grid(row=1, column=0)
        self.frame9.grid(padx=10, pady=10)

    def cite(self):
        try:
            self.helppanel.destroy()
        except:
            pass
        self.helppanel = Toplevel()
        self.helppanel.wm_attributes("-topmost", 1)
        self.helppanel.geometry('+20+30')
        self.helppanel.title('Citing')
        self.frame9 = Frame(self.helppanel)
        self.about1label1 = Label(self.frame9, text='HaploFlow', font='TkDefaultFont 13 bold')
        self.about1label1.grid(row=0, column=0)
        self.supportlabel2 = Label(self.frame9, text='No citation information yet.')
        self.supportlabel2.grid(row=1, column=0)
        self.frame9.grid(padx=10, pady=10)

    def get_flow_name(self):
        filename = tkFileDialog.askopenfilename()
        if filename == '' or filename == ():
            return
        self.flowfile.set(filename)
        self.load_flow()

    def load_flow(self):
        self.poslist = []
        self.flowlist = []
        self.verticalgap = 40
        flowfile = open(self.flowfile.get())
        lf = None
        initiallist = []
        for line in flowfile:
            if line.startswith('V'):
                self.poslist.append([int(line.split()[1].split(',')[0])] + line.split()[1].split(',')[1:])
        flowfile.close()
        flowfile = open(self.flowfile.get())
        for line in flowfile:
            if line.startswith('F'):
                pos = int(line.split()[1].split(',')[0])
                count = int(line.split()[1].split(',')[-2])
                group = int(line.split()[1].split(',')[-1])
                flow = line.split()[1].split(',')[1:-2]
                if lf != pos:
                    while self.poslist[len(initiallist)][0] != pos:
                        initiallist.append(None)
                    initiallist.append([[count, flow, group]])
                    lf = pos
                else:
                    initiallist[-1].append([count, flow, group])
            elif line.startswith('I'):
                maxvar, self.medcount, self.maxcount = map(int, line.split()[1].split(','))
                if self.maxvar is None:
                    self.maxvar = maxvar
            elif line.startswith('G'):
                count = 0
                cumalative = 0
                errcount = 0
                stacker = [0 for i in range(self.maxvar + 1)]
                for i in map(int, line[1:].rstrip().split(',')):
                    count += 1
                    if count > self.maxvar:
                        errcount += i
                    else:
                        cumalative += i
                        stacker[count] = cumalative
                self.flowend = errcount + cumalative
            elif line.startswith('C'):
                self.chrom = line.split()[1]
        flowfile.close()
        stacks = [stacker[:] for i in range(len(initiallist))]
        colordict = {}
        lastcol = 0
        lastcol2 = 0
        lastcol3 = 0
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
                    gapped = False
                    for j in q[1]:
                        if j == '+s':
                            thedirs += '+'
                            start = True
                            if len(thedirs) == 1:
                                x1s.append((i + count) * 4 + 4)
                            else:
                                x2s.append((i + count) * 4 + 4)
                        elif j == '+':
                            thedirs += '+'
                            start = True
                            if len(thedirs) == 1:
                                x1s.append((i + count) * 4 + 3)
                            else:
                                x2s.append((i + count) * 4 + 3)
                        elif j == 'e':
                            if len(thedirs) == 1:
                                x1s[-1] = x1s[-1] - 2
                            else:
                                x2s[-1] = x2s[-1] - 2
                        elif j == '-s':
                            thedirs += '-'
                            start = True
                            if len(thedirs) == 1:
                                x1s.append((i + count) * 4 + 3)
                            else:
                                x2s.append((i + count) * 4 + 3)
                        elif j == '-':
                            thedirs += '-'
                            start = True
                            if len(thedirs) == 1:
                                x1s.append((i + count) * 4 + 3)
                            else:
                                x2s.append((i + count) * 4 + 3)
                        elif j == '_':
                            count += 1
                            gapped = True
                        elif (j == 'x' or int(j) >= self.maxvar) and len(thedirs) == 1:
                            y1s.append(stacks[i + count][-1] + q[0]/2)
                            y1s.append(stacks[i + count][-1] + q[0]/2)
                            if start:
                                start = False
                            else:
                                x1s.append((i+count+1) * 4 - 1)
                            x1s.append((i+count+1) * 4 + 1)
                            stacks[i + count][-1] += q[0]
                            count += 1
                        elif (j == 'x' or int(j) >= self.maxvar) and len(thedirs) == 2:
                            y2s.append(stacks[i + count][-1] + q[0]/2)
                            y2s.append(stacks[i + count][-1] + q[0]/2)
                            if start:
                                start = False
                            else:
                                x2s.append((i+count+1) * 4 - 1)
                            x2s.append((i+count+1) * 4 + 1)
                            stacks[i + count][-1] += q[0]
                            count += 1
                        elif len(thedirs) == 1:
                            y1s.append(stacks[i + count][int(j)] + q[0]/2)
                            y1s.append(stacks[i + count][int(j)] + q[0]/2)
                            if start:
                                start = False
                            else:
                                x1s.append((i+count + 1) * 4 - 1)
                            x1s.append((i+count+1) * 4 + 1)
                            stacks[i + count][int(j)] += q[0]
                            count += 1
                        else:
                            y2s.append(stacks[i + count][int(j)] + q[0]/2)
                            y2s.append(stacks[i + count][int(j)] + q[0]/2)
                            if start:
                                start = False
                            else:
                                x2s.append((i+count + 1) * 4 - 1)
                            x2s.append((i+count+1) * 4 + 1)
                            stacks[i + count][int(j)] += q[0]
                            count += 1
                    if q[2] in colordict:
                        h = colordict[q[2]]
                    else:
                        h = (lastcol + random.randint(50,60)) % 360
                        lastcol = h
                        colordict[q[2]] = h
                    lastcol2 += 0.1 + random.random()/4
                    lastcol2 = lastcol2 % 0.5
                    lastcol3 += 0.05 + random.random()/4
                    lastcol3 = lastcol3 % 0.5
                    s = 0.5 + lastcol2
                    l = 0.25 + lastcol3
                    color = self.hsl_to_rgb(h, s, l)
                    if thedirs == '+':
                        self.flowlist[-1].append([0, x1s, y1s, q[0], color, gapped])
                    elif thedirs == '-':
                        self.flowlist[-1].append([1, x1s, y1s, q[0], color, gapped])
                    elif thedirs == '++':
                        self.flowlist[-1].append([2, (x1s, x2s), (y1s, y2s), q[0], color, gapped])
                    elif thedirs == '+-':
                        self.flowlist[-1].append([3, (x1s, x2s), (y1s, y2s), q[0], color, gapped])
                    elif thedirs == '-+':
                        self.flowlist[-1].append([4, (x1s, x2s), (y1s, y2s), q[0], color, gapped])
                    elif thedirs == '--':
                        self.flowlist[-1].append([5, (x1s, x2s), (y1s, y2s), q[0], color, gapped])
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

    def getflow(self, samfile=None, vcffile=None, outfile=None, maxdist=None):
        try:
            import pysam
        except:
            self.queue.put('pysam not found, please install.')
            return
        if samfile is None:
            samfile = self.bamfile.get()
        else:
            self.queue = clqueue()
        if vcffile is None:
            vcffile = self.vcffile.get()
        if outfile is None:
            outfile = self.flowfile.get()
        if maxdist is None:
            maxdist = self.maxdist.get()
        snps = []
        vcf = open(vcffile)
        out = open(outfile, 'w')
        depthcount = []
        maxvar = 0
        maxeachvar = [] # maximum count for most prevalent variant, second most prevalent variant etc.
        maxerror = 0 # maximum number of reads with no called variant at a single site
        snp2count = {}
        count = 0
        self.queue.put('Reading VCF file.')
        for line in vcf:
            if not line.startswith('#'):
                chrom, pos, id, ref, alt, qual, filt, info, form, unknown = line.split()
                aninstance = variation(chrom, pos, ref, alt, qual)
                snp2count[int(pos)] = count
                count += 1
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
                #out.write('V ' + str(pos) + ' ' + ref + ',' + alt + '\n')
                snps.append(aninstance)
        out.write('C ' + aninstance.chrom + '\n')
        sam = pysam.Samfile(samfile, 'rb')
        reads = {}
        varlist = []
        lastpos = None
        consflow = {}
        flownum = 0
        newflows3 = []
        groupcon = {}
        self.queue.put('Finding flows in BAM.')
        for snp in snps: # for each variant in the vcf file
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
                if '_' in i:
                    newflows3.append(i + ',0\n')
                else:
                    flow = tuple(filter(lambda x: not x in ['+s', '+', '-s', '-', 'e'], i.split(',')[1:-1]))
                    if lastpos != int(i.split(',')[0]):
                        tempdict = {}
                        for j in consflow:
                            temp = j[1:]
                            tempnum = consflow[j]
                            if not temp == ():
                                if temp in tempdict:
                                    if tempdict[temp][1] < tempnum[1]:
                                        tempdict[temp] = tempnum
                                else:
                                    tempdict[temp] = tempnum
                        consflow = tempdict
                        lastpos = int(i.split(',')[0])
                    foundone = False
                    bestcons = 0
                    for j in consflow:
                        consensus = True
                        for k in range(min([len(j), len(flow)])):
                            if j[k] != flow[k]:
                                consensus = False
                        if consensus:
                            if consflow[j][1] >= bestcons:
                                bestcons = consflow[j][1]
                                changeflow = j
                            foundone = True
                    if foundone:
                        tempnum = consflow[changeflow][0]
                        tempnum2 = consflow[changeflow][1]
                        del consflow[changeflow]
                        if len(changeflow) > len(flow):
                            additionalbit = ()
                            tempnewflow = changeflow
                        else:
                            additionalbit = flow[len(changeflow):]
                            tempnewflow = flow
                        consflow[tempnewflow] = (tempnum, max([tempnum2, int(i.split(',')[-1])]))
                        tempgroup = groupcon[tempnum]
                        tempgroup = (tempgroup[0], tempgroup[1] + additionalbit)
                        groupcon[tempnum] = tempgroup
                    else:
                        flownum += 1
                        consflow[flow] = (flownum, int(i.split(',')[-1]))
                        groupcon[flownum] = (snp2count[int(i.split(',')[0])], flow)
                        tempnum = flownum
                    newflows3.append(i + ',' + str(tempnum) + '\n')
            newflows3.sort(key=lambda x: int(x.split(',')[-2]), reverse=True)
            newflows3.sort(key=lambda x: orderflow(x))
            newflows3.sort(key=lambda x: int(x.split(',')[0]))
            iterit = True
            while iterit and newflows3 != []:
                i = newflows3.pop(0)
                if i[-3:] == ',0\n':
                    pos = int(i.split(',')[0])
                    flow = tuple(filter(lambda x: not x in ['+s', '+', '-s', '-', 'e'], i.split(',')[1:-2]))
                    toremove = set()
                    gotit = False
                    for j in groupcon:
                        if groupcon[j][0] + len(groupcon[j][1]) < snp2count[pos]: # TODO get this working properly
                            toremove.add(j)
                        elif snp2count[pos] >= groupcon[j][0] and len(flow) + snp2count[pos] <= len(groupcon[j][1]) + groupcon[j][0]:
                            index = snp2count[pos] - groupcon[j][0]
                            getit = True
                            for q in range(len(flow)):
                                if groupcon[j][1][q+index] == flow[q] or flow[q] == '_':
                                    pass
                                else:
                                    getit = False
                            if getit:
                                i = i[:-2] + str(j) + '\n'
                                out.write('F ' + i)
                                gotit = True
                                break
                    if not gotit:
                        if pos + maxdist * 3 < snp.pos:
                            flownum += 1
                            i = i[:-2] + str(flownum) + '\n'
                            out.write('F ' + i)
                        else:
                            newflows3 = [i] + newflows3
                            iterit = False
                   # for q in toremove:
                    #    del groupcon[q]
                else:
                    out.write('F ' + i)
            for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1):
                if pileupcolumn.pos == snp.pos - 1:
                    vardict = {}
                    varorder = []
                    varcount = {}
                    currerr = 0
                    varlength = len(snp.ref)
                    for i in [snp.ref] + snp.alt.split(','):
                        varcount[i] = 0
                    for pileupread in pileupcolumn.pileups:
                        rvar = pileupread.alignment.seq[pileupread.query_position:pileupread.query_position +
                               pileupread.alignment.get_overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
                        if rvar in varcount:
                            varcount[rvar] += 1
                        else:
                            currerr += 1
                    if currerr > maxerror:
                        maxerror = currerr
                    for i in varcount:
                        varorder.append((varcount[i], i))
                    varorder.sort(reverse=True)
                    count = 0
                    varlist.append('V ' + str(snp.pos))
                    for i in varorder:
                        if i[1] == snp.ref:
                            varlist[-1] += ',*' + i[1]
                        else:
                            varlist[-1] += ',' + i[1]
                        vardict[i[1]] = str(count)
                        count += 1
                        if len(vardict) > maxvar:
                            maxeachvar += [0 for x in range(len(vardict) - maxvar)]
                            maxvar = len(vardict)
                        if i[0] > maxeachvar[count-1]:
                            maxeachvar[count-1] = i[0]
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
                            if rvar in vardict:
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
        newflows3 = []
        for i in newflows2:
            flow = tuple(filter(lambda x: not x in ['+s', '+', '-s', '-', 'e'], i.split(',')[1:-1]))
            if lastpos != int(i.split(',')[0]):
                tempdict = {}
                for j in consflow:
                    temp = j[1:]
                    tempnum = consflow[j]
                    if not temp == ():
                        tempdict[temp] = tempnum
                consflow = tempdict
                lastpos = int(i.split(',')[0])
            foundone = False
            bestcons = 0
            for j in consflow:
             #   print 'dong', j
                consensus = True
                newflow = list(j)
                for k in range(min([len(j), len(flow)])):
                    if j[k] == '_':
                        if flow[k] != '_':
                            newflow[k] = flow[k]
                    elif flow[k] == '_' or j[k] == flow[k]:
                        pass
                    else:
                        consensus = False
                if consensus:
                    if consflow[j][1] >= bestcons:
                        bestcons = consflow[j][1]
                        changeflow = j
                    foundone = True
            if foundone:
                tempnum = consflow[changeflow][0]
                tempnum2 = consflow[changeflow][1]
                del consflow[changeflow]
                if len(changeflow) > (flow):
                    newflow += list(changeflow[len(flow):])
                else:
                    newflow += list(flow[len(changeflow):])
                consflow[tuple(newflow)] = (tempnum, max([tempnum2, int(i.split(',')[-1])]))
            else:
                flownum += 1
                consflow[flow] = (flownum, int(i.split(',')[-1]))
                tempnum = flownum
            newflows3.append(i + ',' + str(tempnum) + '\n')
        newflows3.sort(key=lambda x: int(x.split(',')[-2]), reverse=True)
        newflows3.sort(key=lambda x: orderflow(x))
        newflows3.sort(key=lambda x: int(x.split(',')[0]))
        for i in newflows3:
            out.write('F ' + i)
        for i in varlist:
            out.write(i + '\n')
        out.write('I ' + str(maxvar) + ',' + str(depthcount[3*len(depthcount)/4]) + ',' + str(depthcount[-1]) + '\n')
        out.write('G ' + ','.join(map(str, maxeachvar)) + ',' + str(maxerror) + '\n')
        out.close()
        self.queue.put('Flow file successfully created.')



def orderflow(flow):
    flow.replace(',x,', ',20,')
    for i in range(21):
        for j in range(21):
            if ',' + str(i) + ',' in flow and ',' + str(i) + ',' in flow and i != j:
                if flow.find(',' + str(i) + ',') < flow.find(',' + str(j) + ','):
                    return 2
                elif flow.find(',' + str(i) + ',') < flow.find(',' + str(j) + ','):
                    return 0
    return 1




if len(sys.argv) > 1 and sys.argv[1] == '-cl':
    aninstance = App(None)
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