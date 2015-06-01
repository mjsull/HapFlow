#!/usr/bin/env python
# HapFlow   Written by: Mitchell Sullivan   mjsull@gmail.com
# Supervisor: Dr. Adam Polkinghorne
# Version 1.1.0 31.05.2015
# License: GPLv3



__author__ = 'mjsull'
from Tkinter import *
import threading
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox
import random
import os
import Queue
import platform
import webbrowser
import sys

# get sequence of read aligned to var_pos to var_pos + var-len in the reference
def get_seq_read(read, var_pos, var_len):
    a = None
    b = None
    var_pos -= 1
    for i, j  in read.get_aligned_pairs():
        if j == var_pos and a is None:
            a = i
        if j == var_pos + var_len:
            b = i
            break
    if not a is None and not b is None:
        return read.query_sequence[a:b]
    else:
        return None



# dummy queue when running in command-line mode strings pushed here will be printed to the console instead of the GUI
class clqueue:
    def put(self, theval):
        sys.stdout.write(theval + '\n')


# variation class - holds information about variants from the vcf file
class variation:
    def __init__(self, chrom, pos, ref, alt, qual):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.qual = qual


# Main app
class App:
    def __init__(self, master):
        if master is None:
            try:
                maxdist = int(sys.argv[5])
            except IndexError:
                maxdist = 1000
            try:
                import pysam
            except ImportError:
                sys.stderr.write('pysam not found, please install.')
                return
            testbam = pysam.Samfile(sys.argv[2], 'rb')
            references = testbam.references
            testbam.close()
            for i, j in enumerate(references):
                self.theref = j
                self.getflow(sys.argv[2], sys.argv[3], sys.argv[4] + '.' + str(i) + '.flw', maxdist)
            return
        self.otu = IntVar(value=0)
        self.flowlist = None
        self.otufilename = StringVar(value='')
        self.writeotu_options = StringVar(value='Assign reads to multiple OTUs')
        self.menubar = Menu(master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Create flow file", command=self.create_flow)
        self.filemenu.add_command(label="Load flow file", command=self.get_flow_name)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Change max. variants", command=self.change_max_variants)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.quit)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.toolmenu = Menu(self.menubar, tearoff=0)
        self.otumenu = Menu(self.toolmenu, tearoff=0)
        self.toolmenu.add_command(label="Goto base", command=self.goto_base)
        self.toolmenu.add_command(label="Create image", command=self.create_image)
        self.otumenu.add_radiobutton(label="OTU 1", selectcolor='black', variable=self.otu, value=0)
        self.otumenu.add_radiobutton(label="OTU 2", selectcolor='black', variable=self.otu, value=1)
        self.otumenu.add_radiobutton(label="OTU 3", selectcolor='black', variable=self.otu, value=2)
        self.otumenu.add_radiobutton(label="OTU 4", selectcolor='black', variable=self.otu, value=3)
        self.otumenu.add_radiobutton(label="OTU 5", selectcolor='black', variable=self.otu, value=4)
        self.otumenu.add_radiobutton(label="OTU 6", selectcolor='black', variable=self.otu, value=5)
        self.otumenu.add_radiobutton(label="OTU 7", selectcolor='black', variable=self.otu, value=6)
        self.otumenu.add_radiobutton(label="OTU 8", selectcolor='black', variable=self.otu, value=7)
        self.otumenu.add_radiobutton(label="OTU 9", selectcolor='black', variable=self.otu, value=8)
        self.otumenu.add_radiobutton(label="OTU 10", selectcolor='black', variable=self.otu, value=9)
        self.toolmenu.add_cascade(label="Mark OTU", menu=self.otumenu)
        self.toolmenu.add_command(label="Write OTUs", command=self.write_otus)
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
        self.fontsize2 = 20
        self.customFont = tkFont.Font(family="Courier", size=self.fontsize)
        self.customFont2 = tkFont.Font(family="Courier", size=self.fontsize2, weight='bold')
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
        self.queue = Queue.Queue()
        self.rcmenu = Menu(root, tearoff=0)
        self.rcmenu.add_command(label="Details", command=self.details)
        self.rcmenu.add_command(label="Write flow readnames", command=self.write_flow_names)
        self.rcmenu.add_command(label="Write flow to BAM", command=self.write_flow_bam)
        self.rcmenu.add_command(label="Write group readnames", command=self.write_group_names)
        self.rcmenu.add_command(label="Write group to BAM", command=self.write_group_bam)
        if platform.system() in ['Windows', 'Linux']:
            self.canvas.tag_bind('map', '<Button-3>', self.rightClick)
        else:
            self.canvas.tag_bind('map', '<Button-2>', self.rightClick)
        self.canvas.tag_bind('map', '<Button-1>', self.select_flow)
        self.canvas.bind('<Double-Button-1>', self.select_otu)
        self.canvas.bind('<Button-1>', self.remove_rc)
        root.bind('w', self.stretch_y)
        root.bind('s', self.shrink_y)
        root.bind('a', self.shrink_x)
        root.bind('d', self.stretch_x)
        self.canvas.bind('<Configure>', self.update_frame)
        self.selected = [None, None, None]
        self.lastxmod = None
        self.lastymod = None
        self.otus = [[None, [], []] for i in range(10)]
        if len(sys.argv) == 2:
            self.flowfile.set(sys.argv[1])
            self.load_flow()

    # posts menu when right click on flow - records position of click
    def rightClick(self, event):
        self.rctag = self.canvas.gettags(CURRENT)
        self.rcmenu.unpost()
        self.rcmenu.post(event.x_root, event.y_root)
        self.rcpos = (event.x_root, event.y_root)

    # color flow black to highlight entire flow
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

    # change the maximum number of variants shown per site
    def change_max_variants(self):
        x = tkSimpleDialog.askinteger('Change max. variants.', 'Please choose maximum variants per site to show.')
        self.maxvar = x

    # write defined OTUs to separate BAM files
    def write_otus(self):
        try:
            self.write_otu.destroy()
        except:
            pass
        self.write_otu = Toplevel()
        self.write_otu.grab_set()
        self.write_otu.wm_attributes("-topmost", 1)
        self.otuframe = Frame(self.write_otu)
        self.write_otu.geometry('+20+30')
        self.write_otu.title('Write OTUs to SAM')
        self.samfilenamelabel = Label(self.otuframe, text='Original BAM file:', anchor=E)
        self.samfilenamelabel.grid(column=0, row=0, sticky=E)
        self.samfilenameentry = Entry(self.otuframe, textvariable=self.bamfile)
        self.samfilenameentry.grid(column=1, row=0, sticky=EW)
        self.samfilebutton = Button(self.otuframe, text='...', command=self.load_bam)
        self.samfilebutton.grid(column=2, row=0)
        self.otufilenamelabel = Label(self.otuframe, text='Prefix for output BAM files:', anchor=E)
        self.otufilenamelabel.grid(column=0, row=1, sticky=E)
        self.otufilenameentry = Entry(self.otuframe, textvariable=self.otufilename)
        self.otufilenameentry.grid(column=1, row=1, sticky=EW)
        self.otufilebutton = Button(self.otuframe, text='...', command=self.load_otubam)
        self.otufilebutton.grid(column=2, row=1)
        self.otuopt = OptionMenu(self.otuframe, self.writeotu_options, 'Assign reads to multiple OTUs', 'Ignore reads with multiple OTUs')
        self.otuopt.grid(column=0, row=2, columnspan=2)
        self.okotu = Button(self.otuframe, text='Ok', command=self.ok_otu)
        self.okotu.grid(column=1, row=3, sticky=E)
        self.otuframe.grid(padx=10, pady=10)

    # choose original BAM file to find alignments
    def load_bam(self):
        filename = tkFileDialog.askopenfilename(title='Please select alignment file (BAM) from which flow was generated.')
        self.bamfile.set(filename)

    # select prefix for output BAM files (write OTUs)
    def load_otubam(self):
        filename = tkFileDialog.asksaveasfilename(title='Prefix for output BAM files.')
        self.otufilename.set(filename)

    # Seperates alignments into seperate BAM files
    def ok_otu(self):
        if self.bamfile.get() == '':
            tkMessageBox.showerror('BAM file not found.', 'Please select valid BAM file.')
            return
        if self.otufilename.get() == '':
            tkMessageBox.showerror('Prefix not found.', 'Please select valid prefix.')
            return
        self.write_otu.destroy()
        varnum = 0
        lastvar = ''
        outflows = [[] for i in range(10)]
        outpos = [[] for i in range(10)]
        with open(self.flowfile.get()) as f:
            for line in f:
                if line.startswith('F '):
                    outflow = line.split()[1]
                    pos = outflow.split(',')[0]
                    flow = filter(lambda x: not x in ['+s', '+', '-s', '-', 'e'], outflow.split(',')[1:-2])
                    if pos != lastvar:
                        varnum += 1
                        lastvar = pos
                    getotus = None
                    for theotu, i in enumerate(self.otus):
                        if not i[0] is None and varnum >= i[0] and varnum + len(flow) <= i[0] + len(i[1]):
                            compflow = map(str, i[1][varnum - i[0]:])
                            gotit = True
                            for j, k in enumerate(flow):
                                if k == '_' or k == compflow[j] or (k == 'x' and compflow[j] == str(self.maxvar)) or (k != 'x' and int(k) >= self.maxvar and compflow[j] == str(self.maxvar)):
                                    pass
                                else:
                                    gotit = False
                                    break
                            if gotit and self.writeotu_options.get() == 'Assign reads to multiple OTUs':
                                outflows[theotu].append(outflow)
                                outpos[theotu].append(varnum - 1)
                            elif gotit:
                                if getotus is None:
                                    getotus = theotu
                                else:
                                    getotus = None
                                    break
                    if self.writeotu_options.get() != 'Assign reads to multiple OTUs' and not getotus is None:
                        outflows[getotus].append(outflow)
                        outpos[getotus].append(varnum - 1)
        for i, j in enumerate(outflows):
            if len(j) > 0:
                if os.path.exists(self.otufilename.get() + '.' + str(i+1) + '.bam'):
                    answer = tkMessageBox.askyesno('File already exists', 'Overwrite ' + self.otufilename.get() + '.' + str(i+1) + '.bam')
                    if not answer:
                        return
                self.get_names(j, outpos[i], True, self.otufilename.get() + '.' + str(i+1) + '.bam')


    # keeps track of define OTUs and adds OTU numbers when double clicking the canvas
    def select_otu(self, event):
        absx, absy = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        xnum = int((absx+self.xmod*2) / self.xmod/4)
        xcoord = (xnum) * self.xmod * 4 - self.xmod * 2
        for i, j in enumerate(self.stacker):
            if absy <= self.ypos1 + j*self.ymod:
                ynum = i-1
                break
        if absy >= self.ypos1 + self.stacker[-1]*self.ymod:
            ynum = len(self.stacker) - 1
        if ynum != -1:
            stacknum = (xnum - 1) * (self.maxvar + 1) + ynum
            ycoord = self.ypos1 + self.stacker[ynum] * self.ymod + int(self.otustack[(xnum - 1) * (self.maxvar + 1) + ynum]) * self.fontsize2
            if self.otus[self.otu.get()][1] == []:
                self.otus[self.otu.get()][2].append(self.canvas.create_text(xcoord, ycoord, anchor=NW, text=str(self.otu.get()+1), font=self.customFont2, tags='otus'))
                self.otus[self.otu.get()][0] = xnum
                self.otus[self.otu.get()][1].append(ynum)
                self.otustack = self.otustack[:stacknum] + str(int(self.otustack[stacknum]) + 1) + self.otustack[stacknum+1:]
            elif xnum == self.otus[self.otu.get()][0] - 1:
                self.otus[self.otu.get()][2].insert(0, self.canvas.create_text(xcoord, ycoord, anchor=NW, text=str(self.otu.get()+1), font=self.customFont2, tags='otus'))
                self.otus[self.otu.get()][0] = xnum
                self.otus[self.otu.get()][1].insert(0, ynum)
                self.otustack = self.otustack[:stacknum] + str(int(self.otustack[stacknum]) + 1) + self.otustack[stacknum+1:]
            elif xnum == self.otus[self.otu.get()][0] + len(self.otus[self.otu.get()][1]):
                self.otus[self.otu.get()][2].append(self.canvas.create_text(xcoord, ycoord, anchor=NW, text=str(self.otu.get()+1), font=self.customFont2, tags='otus'))
                self.otus[self.otu.get()][1].append(ynum)
                self.otustack = self.otustack[:stacknum] + str(int(self.otustack[stacknum]) + 1) + self.otustack[stacknum+1:]
            elif xnum == self.otus[self.otu.get()][0]:
                self.canvas.delete(self.otus[self.otu.get()][2][0])
                self.otus[self.otu.get()][2].pop(0)
                self.otus[self.otu.get()][1].pop(0)
                self.otus[self.otu.get()][0] += 1
                self.otustack = self.otustack[:stacknum] + str(int(self.otustack[stacknum]) - 1) + self.otustack[stacknum+1:]
            elif xnum == self.otus[self.otu.get()][0] + len(self.otus[self.otu.get()][1]) - 1:
                self.canvas.delete(self.otus[self.otu.get()][2][-1])
                self.otus[self.otu.get()][2].pop()
                self.otus[self.otu.get()][1].pop()
                self.otustack = self.otustack[:stacknum] + str(int(self.otustack[stacknum]) - 1) + self.otustack[stacknum+1:]



    # remove the right-click menu
    def remove_rc(self, event):
        self.rcmenu.unpost()

    # find flow in the flow file
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

    # get the names of the reads associated with a flow or flows
    def get_names(self, flows, positions, BAM=True, outfile=None):
        try:
            import pysam
        except ImportError:
            tkMessageBox.showerror('Pysam not installed.', 'This functionality is only available if pysam is isntalled.')
            return
        if self.bamfile.get() == '':
            filename = tkFileDialog.askopenfilename(title='Please select alignment file (BAM) from which flow was generated.')
            self.bamfile.set(filename)
        try:
            sam = pysam.Samfile(self.bamfile.get(), 'rb')
        except IOError:
            tkMessageBox.showerror('File not found', 'Please select a valid BAM file.')
            self.bamfile.set('')
            return
        except ValueError:
            tkMessageBox.showerror('File not BAM file', 'Please make sure the file is a valid, indexed BAM file.')
            self.bamfile.set('')
            return
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
            vardict = {}
            varorder = []
            varcount = {}
            currerr = 0
            varlength = len(snp.ref)
            for i in [snp.ref] + snp.alt.split(','):
                varcount[i] = 0
            for theread in sam.fetch(snp.chrom, snp.pos, snp.pos + 1):
                rvar = get_seq_read(theread, snp.pos, varlength)
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
            for theread in sam.fetch(snp.chrom, snp.pos, snp.pos+1):
                readname = theread.query_name
                if not readname in gottenreads: # ignore the start of the second pair in dovetailed reads, there might be a better solution to this but would take a long time to implement and wouldn't provide much additional clarity.
                    rvar = get_seq_read(theread, snp.pos, varlength)
                    if readname in reads and not rvar is None:
                        if reads[readname][1] != theread.is_read1:
                            if theread.is_reverse and theread.reference_start + 1 == snp.pos:
                                reads[readname].append('-s')
                            elif theread.is_reverse:
                                reads[readname].append('-')
                            elif theread.reference_start + 1 == snp.pos:
                                reads[readname].append('+s')
                            else:
                                reads[readname].append('+')
                            reads[readname][1] = theread.is_read1
                    elif not rvar is None:
                        reads[readname] = [snp.pos, theread.is_read1]
                        if theread.is_reverse and theread.reference_start + 1 == snp.pos:
                            reads[readname].append('-s')
                        elif theread.is_reverse:
                            reads[readname].append('-')
                        elif theread.reference_start + 1 == snp.pos:
                            reads[readname].append('+s')
                        else:
                            reads[readname].append('+')
                    if rvar in vardict:
                        reads[readname].append(vardict[rvar])
                    elif not rvar is None:
                        reads[readname].append('x')
                    gottenreads.add(readname)
                    if not rvar is None and theread.reference_end + 1 == snp.pos:
                        reads[readname].append('e')
            for i in reads:
                if not i in gottenreads:
                    reads[i].append('_')
        outset = set()
        for i in reads:
            for j in flows:
                if int(j.split(',')[0]) == reads[i][0] and ','.join(reads[i][2:]).strip(',_') == ','.join(j.split(',')[1:-2]):
                    outset.add(i)
        if BAM:
            if outfile is None:
                outfile = tkFileDialog.asksaveasfilename(title='Choose path to write alignments to (BAM).')
            if outfile[-4:] != '.bam':
                outfile += '.bam'
            try:
                newsam = pysam.Samfile(outfile, 'wb', template=sam)
            except IOError:
                tkMessageBox.showerror('File not valid', 'Please select a valid output file.')
                return
            for read in sam.fetch():
                if read.qname in outset:
                     newsam.write(read)
            newsam.close()
        else:
            if outfile is None:
                outfile = tkFileDialog.asksaveasfilename(title='Choose path to write read names to.')
            try:
                out = open(outfile, 'w')
            except IOError:
                tkMessageBox.showerror('File not valid', 'Please select a valid output file.')
                return
            for i in outset:
                out.write(i + '\n')
            out.close()
        sam.close()


    # create window with details of the flow
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
        self.detail_frame = Frame(self.detail_window)
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


    # write read names of flow
    def write_flow_names(self):
        pos, num, amap, securrent = self.rctag
        posnum = int(pos[1:])
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        self.get_names([flow], [posnum], False)

    # write bam alignment of flow
    def write_flow_bam(self):
        pos, num, amap, securrent = self.rctag
        posnum = int(pos[1:])
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        self.get_names([flow], [posnum])

    # get all flows associated with group
    def get_group(self, group):
        flowfile = open(self.flowfile.get())
        count = 0
        flows, posnums = [], []
        for line in flowfile:
            if line.startswith('F'):
                if line.rstrip()[-len(group)-1:] == ',' + group:
                    flows.append(line.split()[1])
                    posnums.append(count)
                count += 1
        return flows, posnums

    # write bam file of all flows in group
    def write_group_bam(self):
        pos, num, amap, securrent = self.rctag
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        group = flow.split(',')[-1]
        flows, posnums = self.get_group(group)
        self.get_names(flows, posnums)

    # write read names of all flows in group
    def write_group_names(self):
        pos, num, amap, securrent = self.rctag
        pos = str(self.poslist[int(pos[1:])][0])
        num = int(num[1:])
        flow = self.find_flow(pos, num)
        group = flow.split(',')[-1]
        flows, posnums = self.get_group(group)
        self.get_names(flows, posnums, False)

    # convert hue/saturation/lightness to red green blue
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

    # update the frame - remove flows not near the frame of focus add flows coming near the frame of focus
    def update_frame(self, temp=None):
        if self.flowlist is None:
            return
        tilt = 0.5
        x1 = self.canvas.canvasx(0)
        x2 = self.canvas.canvasx(self.canvas.winfo_width())
        self.canvas.delete('top')
        indexa = int(x1 / self.xmod/4)
        indexb = int(x2 / self.xmod/4)
        if indexb >= len(self.flowlist):
            indexb = len(self.flowlist) - 1
        positions = []
        suppositions = []
        if self.xmod == self.lastxmod and self.ymod == self.lastymod:
            toremove = set()
            for i in self.currflows:
                if i < indexa -10 or i > indexb + 10:
                    self.canvas.delete('p' + str(i))
                    toremove.add(i)
            for i in toremove:
                self.currflows.remove(i)
        else:
            self.currflows = set()
            self.canvas.delete('top')
            self.canvas.delete('map')
            self.canvas.delete('gapped')
            self.lastxmod = self.xmod
            self.lastymod = self.ymod
        for i in self.stacker:
            self.canvas.create_line(x1+5, self.ypos1 + i * self.ymod - 3, x2-5, self.ypos1 + i * self.ymod - 3, tags='top', fill='gray')
        for i in range(max([indexa, 0]), indexb):
            if not i in self.currflows:
                count = 0
                self.currflows.add(i)
                for j in self.flowlist[i]:
                    count += 1
                    if j[0] == 0: # if single forward
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1]) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + tilt) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - tilt) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k]) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 1: # if single reverse
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1]) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + tilt) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - tilt) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k]) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 2: # if pair F F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1]) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + tilt) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - tilt) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k]) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1]) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + tilt) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - tilt) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k]) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + tilt) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - tilt) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 3: # if pair F R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1]) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + tilt) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - tilt) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k]) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1]) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + tilt) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - tilt) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k]) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + tilt) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - tilt) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 4: # if pair R F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1]) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + tilt) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - tilt) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k]) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1]) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + tilt) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - tilt) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k]) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + tilt) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - tilt) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
                    elif j[0] == 5: # if pair R R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1]) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + tilt) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - tilt) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k]) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1]) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + tilt) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - tilt) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k]) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], tags=('p' + str(i), 'f' + str(count), 'map'))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=('p' + str(i), 'f' + str(count), 'map'))
                        if j[5]:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), fill='#000000', dash=(5,2), tags=('p' + str(i), 'gapped'), state=self.gapped_state)
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + tilt) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - tilt) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=min([10, j[3] * self.ymod / 8]), fill=j[4], dash=(5,2), tags=('p' + str(i), 'f' + str(count), 'map'))
            xpos = (i + 1) * self.xmod * 4
            if i >= indexa:
                suppositions.append(xpos)
                positions.append(self.poslist[i][0])
                thetext = str(self.poslist[i][0]) + '\n' + '\n'.join(self.poslist[i][1:])
                self.canvas.create_text(xpos + 4, self.ymod * self.flowend + self.ypos1 + 5, anchor=NW, text=thetext, font=self.customFont, tags='top')
        self.canvas.create_rectangle(x1+10, self.ypossnp + 20, x2 - 10, self.ypossnp, tags='top', fill='#E1974C')
        self.canvas.tag_raise('gapped')
        self.canvas.tag_raise('otus')
        if len(positions) > 1:
            for i in range(len(positions)):
                xpos = x1+15 + int((positions[i] - positions[0]) * 1.0 / (positions[-1] - positions[0]) * (x2 - x1 - 30))
                self.canvas.create_line(suppositions[i], self.ypos1 + self.ymod * self.flowend, suppositions[i], self.ypos1 - 5, xpos, self.ypossnp + 20, xpos, self.ypossnp, tags='top', width=2)
            self.canvas.create_rectangle(x1 + 10, self.yposref + 20, x2 - 10, self.yposref, tags='top', fill='#7293CB')
            starto = x1 + 10 + int(positions[0] * 1.0 / self.reflength * (x2 - x1 - 20))
            endo = x1 + 10 + int(positions[-1] * 1.0 / self.reflength * (x2 - x1 - 20))
            self.canvas.create_rectangle(starto, self.yposref + 20, endo, self.yposref, tags='top', fill='#E1974C')
            self.canvas.create_text(x1 + 10, self.ypossnp - 2, anchor=SW, text='SNP block start..stop: ' + str(positions[0]) + '..' + str(positions[-1]), font=self.customFont, tags='top')
    
    # open a window that can initiate creating a flow file
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


    # ask for bam file name
    def loadbam(self):
        filename = tkFileDialog.askopenfilename(parent=self.create_flow_top)
        if filename == '':
            return
        self.bamfile.set(filename)

    # ask for vcf file name
    def loadvcf(self):
        filename = tkFileDialog.askopenfilename(parent=self.create_flow_top)
        if filename == '':
            return
        self.vcffile.set(filename)

    # ask for flow file name
    def loadflow(self):
        filename = tkFileDialog.asksaveasfilename(parent=self.create_flow_top)
        if filename == '':
            return
        else:
            self.flowfile.set(filename)

    # initiate the flow file creation process
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
        testsam = pysam.Samfile(self.bamfile.get(), 'rb')
        self.create_flow_top.destroy()
        if len(testsam.references) == 0:
            tkMessageBox.showerror('No reference in BAM file', 'Can this even happen?')
            return
        elif len(testsam.references) == 1:
            self.theref = testsam.references[0]
            testsam.close()
            self.run_flow()
        else:
            self.choice_top = Toplevel()
            self.choice_top.grab_set()
            self.choice_top.wm_attributes("-topmost", 1)
            self.choice_top.geometry('+20+30')
            self.choice_top.title('Choose reference')
            self.choice_frame = Frame(self.choice_top)
            self.choice_label = Label(self.choice_frame, text='Please choose reference\n from which to create flow:')
            self.choice_label.grid(row=0, column=0)
            self.choice_scroll = Scrollbar(self.choice_frame, orient=VERTICAL)
            self.choice_entry = Listbox(self.choice_frame, yscrollcommand=self.choice_scroll.set)
            self.choice_scroll.config(command=self.choice_entry.yview)
            self.choice_scroll.grid(row=1, column=1, sticky=NS)
            self.choice_entry.grid(row=1, column=0, sticky=EW)
            for i in testsam.references:
                self.choice_entry.insert(END, i)
            testsam.close()
            self.okchoice = Button(self.choice_frame, text='Ok', command=self.ok_choice)
            self.okchoice.grid(row=2, column=0, columnspan=2, sticky=E)
            self.choice_frame.grid(padx=10, pady=10)



    # ok button for choosing reference from references available in BAM
    def ok_choice(self):
        self.theref = self.choice_entry.get(ACTIVE)
        self.choice_top.destroy()
        self.run_flow()
    
    # open a window with updates about the progress of creating the flow file
    def run_flow(self):
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

    # remove the console
    def ok_console(self):
        self.run_flow_top.destroy()

    # get messages from thread and print to console
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

    # hide the paired end connectors for flows that skip a variant
    def hide_gapped(self):
        self.canvas.itemconfig('gapped', state=HIDDEN)
        self.gapped_state = HIDDEN

    # show the connecotrs for flows that skip a variant
    def show_gapped(self):
        self.canvas.itemconfig('gapped', state=NORMAL)
        self.gapped_state = NORMAL

    # stretch flows in the x direction
    def stretch_x(self, stuff=None):
        x1 = self.canvas.canvasx(0)
        indexa = int(x1 / self.xmod/4)
        self.xmod = self.xmod * 1.0526315789473684210526315789474
        self.currxscroll = self.currxscroll * 1.0526315789473684210526315789474
        self.canvas.scale('otus', 0, 0, 1.0526315789473684210526315789474, 1)
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.goto_base(indexa)

    # shrink flows in the x direction
    def shrink_x(self, stuff=None):
        x1 = self.canvas.canvasx(0)
        indexa = int(x1 / self.xmod/4)
        self.xmod = self.xmod * 0.95
        self.currxscroll = self.currxscroll * 0.95
        self.canvas.scale('otus', 0, 0, 0.95, 1)
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.goto_base(indexa)

    # stretch flows in the y direction
    def stretch_y(self, stuff=None):
        self.fontsize2 = self.fontsize2 * 1.0526315789473684210526315789474
        self.customFont2.configure(size=int(round(self.fontsize2)))
        self.ymod = self.ymod * 1.0526315789473684210526315789474
        for i in self.canvas.find_withtag('otus'):
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], (x[1] - self.ypos1) * 1.0526315789473684210526315789474 + self.ypos1)
        self.update_frame()

    # shrink flows in the y direction
    def shrink_y(self, stuff=None):
        self.fontsize2 = self.fontsize2 * 0.95
        self.customFont2.configure(size=int(round(self.fontsize2)))
        self.ymod = self.ymod * 0.95
        for i in self.canvas.find_withtag('otus'):
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], (x[1] - self.ypos1) * 0.95 + self.ypos1)
        self.update_frame()

    # quit the program
    def quit(self):
        root.quit()

    # ask the user for a base number then go ot the next variant
    def goto_base(self, fraction=None):
        if fraction is None:
            base = tkSimpleDialog.askinteger('Goto base', 'Base number')
            for i in range(len(self.poslist)):
                if base < self.poslist[i][0]:
                    fraction = (i-1) * 1.0 / len(self.poslist)
                    break
        else:
            fraction = (fraction - 1) * 1.0 / len(self.poslist)
        self.canvas.xview_moveto(fraction)
        self.update_frame()

    # create and SVG image of the canvas
    def create_image(self):
        try:
            saveas = tkFileDialog.asksaveasfilename(parent=root)
        except IOError:
            tkMessageBox.showerror('File not valid', 'Please choose another file.')
        try:
            import canvasvg
        except:
            tkMessageBox.showerror('canvasvg not found', 'Please install canvasvg.')
            return
        canvasvg.saveall(saveas, self.canvas)

    # open the hapflow README
    def help(self):
        webbrowser.open_new('https://github.com/mjsull/HapFlow/blob/master/README.md')

    # create a window with information about Hapflow
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
        self.about1label = Label(self.frame7, text='HapFlow', font='TkDefaultFont 13 bold')
        self.about1label.grid(row=0, column=0)
        self.about2label = Label(self.frame7, text='HapFlow is a Python application for visualising\n\
haplotypes present in sequencing data.\n\n\
Version 1.1\n')
        self.about2label.grid(row=1, column=0)
        self.frame7.grid(padx=10, pady=10)

    # tell people to email me with problems (not sure why I sound so eager)
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
        self.about1label1 = Label(self.frame9, text='HapFlow', font='TkDefaultFont 13 bold')
        self.about1label1.grid(row=0, column=0)
        self.supportlabel2 = Label(self.frame9, text='written by Mitchell Sullivan - mjsull@gmail.com\n\
Please do not hesitate to email with issues or bug reports.')
        self.supportlabel2.grid(row=1, column=0)
        self.frame9.grid(padx=10, pady=10)

    # create window with citation information
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
        self.about1label1 = Label(self.frame9, text='HapFlow', font='TkDefaultFont 13 bold')
        self.about1label1.grid(row=0, column=0)
        self.supportlabel2 = Label(self.frame9, text='''Sullivan MJ, Bachmann NL, Timms P, Polkinghorne A. (2015)
HapFlow: Visualising haplotypes in sequencing data
PeerJ PrePrints 3:e1105
https://dx.doi.org/10.7287/peerj.preprints.895v1
''')
        self.supportlabel2.grid(row=1, column=0)
        self.frame9.grid(padx=10, pady=10)

    # ask for flow file name
    def get_flow_name(self):
        filename = tkFileDialog.askopenfilename()
        if filename == '' or filename == ():
            return
        self.flowfile.set(filename)
        self.load_flow()

    # load a flow file
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
                if maxvar < self.maxvar:
                    self.maxvar = maxvar
            elif line.startswith('G'):
                count = 0
                cumalative = 0
                errcount = 0
                self.stacker = [0 for i in range(self.maxvar + 1)]
                for i in map(int, line[1:].rstrip().split(',')):
                    count += 1
                    if count > self.maxvar:
                        errcount += i
                    else:
                        cumalative += i
                        self.stacker[count] = cumalative
                self.flowend = errcount + cumalative
            elif line.startswith('C'):
                self.chrom = line.split()[1]
        flowfile.close()
        stacks = [self.stacker[:] for i in range(len(self.poslist))]
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
        self.ymod = 800.0 / self.flowend
        self.fontsize2 = self.stacker[1] * self.ymod / 20
        self.customFont2.configure(size=int(round(self.fontsize2)))
        self.otustack = '0' * ((self.maxvar + 1) * len(self.poslist))
        self.update_frame()


    # scroll in x direction
    def xscroll(self, stuff, stuff2, stuff3=None):
        if not stuff3 is None:
            self.canvas.xview_scroll(stuff2, 'units')
        else:
            self.canvas.xview_moveto(stuff2)
        self.update_frame()

    # create the flow file
    def getflow(self, samfile=None, vcffile=None, outfile=None, maxdist=None):
        if samfile is None:
            samfile = self.bamfile.get()
        else:
            self.queue = clqueue()
        if vcffile is None:
            vcffile = self.vcffile.get()
        try:
            import pysam
        except ImportError:
            self.queue.put('pysam not found, please install.')
            return
        if outfile is None:
            outfile = self.flowfile.get()
        if maxdist is None:
            maxdist = self.maxdist.get()
        snps = []

        out = open(outfile, 'w')
        depthcount = []
        maxvar = 0
        maxeachvar = [] # maximum count for most prevalent variant, second most prevalent variant etc.
        maxerror = 0 # maximum number of reads with no called variant at a single site
        snp2count = {}
        count = 0
        self.queue.put('Reading VCF file.') # read the vcf file
        with open(vcffile) as vcf:
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
                    if aninstance.chrom == self.theref:
                        snps.append(aninstance)
        if len(snps) == 0:
            self.queue.put('No variants found in vcf.')
            return
        out.write('C ' + snps[0].chrom + '\n')
        sam = pysam.Samfile(samfile, 'rb')
        reads = {}
        varlist = []
        lastpos = None
        consflow = {}
        flownum = 0
        newflows3 = []
        groupcon = {}
        snps.sort(key=lambda x: x.pos)
        self.queue.put('Finding flows in BAM.')
        for snp in snps: # for each variant in the vcf file
            newflows = {}
            removereads = []
            for i in reads: # move reads from dictionary to a list when we get past their first variant by maxdist bases
                if reads[i][0] < snp.pos - maxdist:
                    rstring = str(reads[i][0]) + ',' + ','.join(reads[i][2:])
                    if rstring.count('-') + rstring.count('+') > 2: # + and - indicate start of reads in flow - if there are more than two starts something has gone wrong. May need to change for future read types.
                        self.queue.put('Duplicate read names detected, please rename reads (or ensure reads only map once)')
                        return
                    if rstring in newflows: # if haplotype profile exists increment by one
                        newflows[rstring] += 1
                    else:
                        newflows[rstring] = 1
                    removereads.append(i)
            for i in removereads: # remove added reads from dictionary
                del reads[i]
            newflows2 = []
            for i in newflows:
                newflows2.append(i.strip(',_') + ',' + str(newflows[i])) # remove padding added to end of read
            newflows2.sort(key=lambda x: int(x.split(',')[-1]), reverse=True) # sort by position and then frequency (descending)
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
            newflows3.sort(key=lambda x: int(x.split(',')[-2]), reverse=True) # sort by frequency
            newflows3.sort(key=lambda x: orderflow(x)) # move flows that span multiple rows into the centre to create less criss crossing
            newflows3.sort(key=lambda x: int(x.split(',')[0])) # sort by position
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
                else:
                    out.write('F ' + i)
            # first order the variants
            vardict = {}
            varorder = []
            varcount = {}
            currerr = 0
            varlength = len(snp.ref)
            for i in [snp.ref] + snp.alt.split(','): # create counts of each called variant
                varcount[i] = 0
            for theread in sam.fetch(snp.chrom, snp.pos, snp.pos+1):
                rvar = get_seq_read(theread, snp.pos, varlength)
                if rvar in varcount:
                    varcount[rvar] += 1
                elif rvar != None:
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
            for theread in sam.fetch(snp.chrom, snp.pos, snp.pos+1):
                readname = theread.query_name
                if not readname in gottenreads: # ignore the start of the second pair in dovetailed reads, there might be a better solution to this but would take a long time to implement and wouldn't provide much additional clarity.
                    rvar = get_seq_read(theread, snp.pos, varlength)
                    if readname in reads and not rvar is None:
                        if reads[readname][1] != theread.is_read1:
                            if theread.is_reverse and theread.reference_start + 1 == snp.pos:
                                reads[readname].append('-s')
                            elif theread.is_reverse:
                                reads[readname].append('-')
                            elif theread.reference_start + 1 == snp.pos:
                                reads[readname].append('+s')
                            else:
                                reads[readname].append('+')
                            reads[readname][1] = theread.is_read1
                    elif not rvar is None:
                        reads[readname] = [snp.pos, theread.is_read1]
                        if theread.is_reverse and theread.reference_start + 1 == snp.pos:
                            reads[readname].append('-s')
                        elif theread.is_reverse:
                            reads[readname].append('-')
                        elif theread.reference_start + 1 == snp.pos:
                            reads[readname].append('+s')
                        else:
                            reads[readname].append('+')
                    if rvar in vardict:
                        reads[readname].append(vardict[rvar])
                    elif not rvar is None:
                        reads[readname].append('x')
                    gottenreads.add(readname)
                    if not rvar is None and theread.reference_end + 1 == snp.pos:
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


# order flows so that there is less clutter
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


if len(sys.argv) > 1 and sys.argv[1] in set(['-h', '-help', '--help']):
    sys.stdout.write('''HapFlow versions 1.1.0
written by Mitchell J Sullivan (mjsull@gmail.com)
This software is freely available under a gplv3 license.

To create a flow file from the command line run:

 $ python HapFlow.py -cl <bam_file> <vcf_file> <output_prefix>

Where bam_file is and indexed bam file
vcf file is a vcf file containing called variants and
output_prefix is the prefix for all output flow files (one for each flow)

For in depth instructions please consult the manual available at
mjsull.github.io/HapFlow\n\n''')
elif len(sys.argv) > 1 and sys.argv[1] == '-cl':
    aninstance = App(None)
else:
    root = Tk()
    root.title('HapFlow')
    root.option_add('*Font', 'Courier 10')
    root.option_add("*Background", "#E0E0FF")
    root.option_add("*Foreground", "#2B3856")
    root.option_add("*Listbox.Background", '#FFFFFF')
    root.option_add("*Scrollbar.Background", "#C0C0FF")
    root.option_add("*Entry.Background", "#FFFFFF")
    root.geometry("600x400")
    app = App(root)
    root.mainloop()
