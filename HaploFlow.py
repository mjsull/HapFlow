__author__ = 'mjsull'
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox


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
        self.curryscroll = 1000
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
        self.block_width = 82
        self.gap_size = 8
        self.block_height = 60
        self.flowfile = None
        try:
            self.flowfile = sys.argv[1]
            self.load_flow()
        except:
            pass

    def update_frame(self, x2=None):
        lastone = True
        x1 = self.canvas.canvasx(0)
        if x2 is None:
            x2 = self.canvas.canvasx(self.canvas.winfo_width())
        print self.canvas.winfo_width(), x1, x2
        snps = self.canvas.find_overlapping(x1, self.ypos1 - 2, x2, self.ypos1 - 4)
        self.canvas.delete(ALL)
        indexa = int(x1 / (self.block_width + self.gap_size))
        indexb = int(x2 / (self.block_width + self.gap_size))
        if indexb >= len(self.flowlist):
            indexb = len(self.flowlist) - 1
        positions = []
        suppositions = []
        for i in range(indexa, indexb + 1):
            xpos = i * (self.block_width + self.gap_size) + self.gap_size
            three = self.flowlist[i][:8]
            two = self.flowlist[i][8:12]
            one = self.flowlist[i][12:]
            if sum(three) >= self.min_flow:
                if lastone:
                    self.drawThree(map(lambda x: x * 1.0 / sum(three), three), xpos, self.ypos3)
                else:
                    self.drawThree(map(lambda x: x * 1.0 / sum(three), three), xpos, self.ypos4)
                lastone = not lastone
            else:
                lastone = True
            if sum(two) >= self.min_flow:
                self.drawTwo(map(lambda x: x * 1.0 / sum(two), two), xpos, self.ypos2)
            self.drawOne(map(lambda x: x * 1.0 / sum(one), one), xpos, self.ypos1)
            suppositions.append(xpos - self.gap_size/2)
            positions.append(self.poslist[i])
            thetext = str(self.poslist[i]) + '\nr   ' + str(self.flowlist[i][12]) + '\na   ' + str(self.flowlist[i][13]) \
             + '\nrr  ' + str(self.flowlist[i][8]) + '\nra  ' + str(self.flowlist[i][9]) + '\nar  ' + str(self.flowlist[i][10])\
             + '\naa  ' + str(self.flowlist[i][11]) + '\nrrr ' + str(self.flowlist[i][0]) + '\nrra ' + str(self.flowlist[i][1]) \
             + '\nraa ' + str(self.flowlist[i][2]) + '\nrar ' + str(self.flowlist[i][3]) + '\nara ' + str(self.flowlist[i][4]) \
             + '\narr ' + str(self.flowlist[i][5]) + '\naar ' + str(self.flowlist[i][6]) + '\naaa ' + str(self.flowlist[i][7])
            self.canvas.create_text(xpos + 4, self.ypos4 + self.block_height + 5, anchor=NW, text=thetext, font=self.customFont, tags='top')
        self.canvas.create_rectangle(x1+10, self.ypossnp + 20, x2 - 10, self.ypossnp, tags='top', fill='#E1974C')
        positions = positions[1:]
        suppositions = suppositions[1:]
        if len(positions) > 1:
            for i in range(len(positions)):
                xpos = x1+15 + int((positions[i] - positions[0]) * 1.0 / (positions[-1] - positions[0]) * (x2 - x1 - 30))
                self.canvas.create_line(suppositions[i], self.ypos4 + self.block_height + 5, suppositions[i], self.ypos1 - 5, xpos, self.ypossnp + 20, xpos, self.ypossnp, tags='top', width=2)
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
        xpos = base * (self.block_width + self.gap_size) - 20
        self.canvas.xview_scroll(xpos, 'units')

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
        for line in open(self.flowfile):
            if sum(map(int, line.split()[1:])) >= 1:
                self.poslist.append(int(line.split()[0]))
                self.flowlist.append(map(int, line.split()[1:]))
        self.currxscroll = len(self.poslist) * (self.block_width + self.gap_size) + 100
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        if self.reflength is None:
            self.reflength = self.poslist[-1]
        self.update_frame(600)



    def xscroll(self, stuff, stuff2, stuff3=None):
        if not stuff3 is None:
            self.canvas.xview_scroll(stuff2, 'units')
        else:
            self.canvas.xview_moveto(stuff2)
        self.update_frame()

    def drawOne(self, vals, x, y):
        twidth, theight = self.block_width / 2, self.block_height
        x -= twidth/2 + self.gap_size /2
        #self.canvas.create_rectangle(x, y, x+twidth, y+theight, width=1)
        curry1 = y
        curry2 = y + theight
        theight = 0.5 * theight
        if vals[0] != 0:
            lwidth = vals[0] * theight + 2
            self.canvas.create_rectangle(x, curry1, x+twidth, curry1+lwidth, fill='#7293CB')
        if vals[1] != 0:
            lwidth = vals[1] * theight + 2
            self.canvas.create_rectangle(x, curry2, x+twidth, curry2-lwidth, fill='#808585')


    def drawTwo(self, vals, x, y):
        twidth, theight = self.block_width, self.block_height
        #self.canvas.create_rectangle(x, y, x+twidth, y+theight, width=1)
        curry1 = y
        curry2 = y + theight
        theight = 0.5 * theight
        if vals[0] != 0:
            lwidth = vals[0] * theight + 2
            self.canvas.create_rectangle(x, curry1, x+twidth, curry1+lwidth, fill='#7293CB')
            curry1 += lwidth
        if vals[3] != 0:
            lwidth = vals[3] * theight + 2
            self.canvas.create_rectangle(x, curry2, x+twidth, curry2-lwidth, fill='#808585')
            curry2 -= lwidth
        if vals[1] != 0:
            lwidth = vals[1] * theight + 2
            self.canvas.create_polygon(x, curry1, x, curry1, x+twidth/2 + lwidth, curry1, x+twidth/2 + lwidth,
                                    curry2 - lwidth, x+twidth, curry2 - lwidth, x+twidth, curry2 - lwidth, x + twidth,
                                    curry2, x+twidth, curry2, x+twidth/2, curry2, x+twidth/2, curry1+lwidth,
                                    x, curry1 + lwidth, x, curry1 + lwidth,
                                    fill='#E1974C', smooth=1, outline='black')
        if vals[2] != 0:
            lwidth = vals[2] * theight + 2
            self.canvas.create_polygon(x, curry2, x, curry2, x+twidth/2, curry2, x+twidth/2,
                                    curry1 + lwidth, x+twidth, curry1+lwidth, x+twidth, curry1+lwidth, x+twidth,
                                    curry1, x+twidth, curry1, x+twidth / 2 - lwidth, curry1, x+twidth / 2 - lwidth, curry2-lwidth,
                                    x, curry2-lwidth, x, curry2-lwidth,
                                    fill='#AB6857', smooth=1, outline='black')


    def drawThree(self, vals, x, y):
        twidth, theight = self.block_width * 2 + self.gap_size, self.block_height
        #self.canvas.create_rectangle(x, y, x+twidth, y+theight, width=1)
        curry11 = y
        curry21 = y + theight
        curry12 = y
        curry22 = y + theight
        curry13 = y
        curry23 = y + theight
        before1, before2, before3, before4 = 0, 0, 0, 0
        theight = 0.5 * theight
        if vals[0] != 0:
            lwidth = vals[0] * theight + 2
            self.canvas.create_rectangle(x, curry11, x+twidth, curry13+lwidth, fill='#7293CB')
            curry11 += lwidth
            curry12 += lwidth
            curry13 += lwidth
        if vals[7] != 0:
            lwidth = vals[7] * theight + 2
            self.canvas.create_rectangle(x, curry21, x+twidth, curry23-lwidth, fill='#808585')
            curry21 -= lwidth
            curry22 -= lwidth
            curry23 -= lwidth
        if vals[1] != 0:
            lwidth = vals[1] * theight + 2
            self.canvas.create_rectangle(x, curry11, x+twidth/2, curry12+lwidth, fill='#E1974C')
            self.canvas.create_polygon(x+twidth/2, curry12, x+twidth/2, curry12, x+twidth*3/4 + lwidth, curry12, x+twidth*3/4  + lwidth,
                                    curry23 - lwidth, x+twidth, curry23-lwidth, x+twidth, curry23-lwidth, x+twidth,
                                    curry23, x+twidth, curry23, x+twidth*3/4, curry23, x+twidth*3/4, curry12+lwidth,
                                    x+twidth/2, curry12+lwidth, x+twidth/2, curry12+lwidth,
                                    fill='#E1974C', smooth=1, outline='black')
            curry11 += lwidth
            curry12 += lwidth
            curry23 -= lwidth
            before1 = lwidth
        if vals[6] != 0:
            lwidth = vals[6] * theight + 2
            self.canvas.create_rectangle(x, curry21, x+twidth/2, curry22-lwidth, fill='#AB6857')
            self.canvas.create_polygon(x+twidth/2, curry22, x+twidth/2, curry22, x+twidth*3/4, curry22, x+twidth*3/4,
                                    curry13 + lwidth, x+twidth, curry13+lwidth, x+twidth, curry13+lwidth, x+twidth,
                                    curry13, x+twidth, curry13, x+twidth*3/4 - lwidth, curry13, x+twidth*3/4 - lwidth, curry22-lwidth,
                                    x+twidth/2, curry22-lwidth, x+twidth/2, curry22-lwidth,
                                    fill='#AB6857', smooth=1, outline='black')
            curry21 -= lwidth
            curry22 -= lwidth
            curry13 += lwidth
            before2 = lwidth
        if vals[2] != 0:
            lwidth =vals[2] * theight + 2
            if curry22 < curry23:
                mod = -1
            else:
                mod = 1
            self.canvas.create_polygon(x, curry11, x, curry11, x+twidth/4+lwidth, curry11, x+twidth/4+lwidth,
                                    curry22 - lwidth, x+twidth/2, curry22-lwidth, x+twidth/2, curry22-lwidth, x+twidth/2,
                                    curry22, x+twidth/2, curry22, x+twidth/4, curry22, x+twidth/4, curry12+lwidth,
                                    x, curry11+lwidth, x, curry11+lwidth,
                                    fill='#84BA5B', smooth=1, outline='black')
            self.canvas.create_polygon(x+twidth/2, curry22, x+twidth/2, curry22, x+twidth*3/4 + lwidth/2 * mod, curry22,
                                    x+twidth*3/4 + lwidth/2 * mod, curry23, x+twidth, curry23, x+twidth, curry23, x+twidth,
                                    curry23-lwidth, x+twidth, curry23-lwidth, x+twidth*3/4 - lwidth/2 * mod, curry23-lwidth,
                                    x+twidth*3/4 - lwidth/2 * mod, curry22-lwidth, x+twidth/2, curry22-lwidth, x+twidth/2,
                                    curry22-lwidth, fill='#84BA5B', smooth=1, outline='black')
            curry11 += lwidth
            curry22 -= lwidth
            curry23 -= lwidth
            before3 = lwidth
        if vals[5] != 0:
            lwidth = vals[5] * theight + 2
            self.canvas.create_polygon(x, curry21, x, curry21, x+twidth/4, curry21, x+twidth/4,
                                    curry12 + lwidth, x+twidth/2, curry12+lwidth, x+twidth/2, curry12+lwidth, x+twidth/2,
                                    curry12, x+twidth/2, curry12, x+twidth/4 - lwidth, curry12, x+twidth/4 - lwidth, curry21-lwidth,
                                    x, curry21-lwidth, x, curry21-lwidth,
                                    fill='#CCC210', smooth=1, outline='black')
            if curry12 < curry13:
                mod = -1
            else:
                mod = 1
            self.canvas.create_polygon(x+twidth/2, curry12, x+twidth/2, curry12, x+twidth*3/4 - lwidth/2 * mod, curry12,
                                    x+twidth*3/4 - lwidth/2 * mod, curry13, x+twidth, curry13, x+twidth, curry13, x+twidth,
                                    curry13+lwidth, x+twidth, curry13+lwidth, x+twidth*3/4 + lwidth/2 * mod, curry13 + lwidth,
                                    x+twidth*3/4 + lwidth/2 * mod, curry12+lwidth, x+twidth/2, curry12+lwidth, x+twidth/2,
                                    curry12+lwidth, fill='#CCC210', smooth=1, outline='black')
            curry21 -= lwidth
            curry12 += lwidth
            curry13 += lwidth
            before4 = lwidth
        if vals[3] != 0:
            lwidth =vals[3] * theight + 2
            self.canvas.create_polygon(x, curry11, x, curry11, x+twidth/4-before4, curry11, x+twidth/4-before4,
                                    curry22 - lwidth, x+twidth/2, curry22-lwidth, x+twidth/2, curry22-lwidth, x+twidth/2,
                                    curry22, x+twidth/2, curry22, x+twidth/4 - lwidth-before4, curry22, x+twidth/4 - lwidth - before4, curry11+lwidth,
                                    x, curry11+lwidth, x, curry11+lwidth,
                                    fill='#D35E60', smooth=1, outline='black')
            self.canvas.create_polygon(x+twidth/2, curry22, x+twidth/2, curry22, x+twidth*3/4 - before2, curry22, x+twidth*3/4 - before2,
                                    curry13 + lwidth, x+twidth, curry13+lwidth, x+twidth, curry13+lwidth, x+twidth,
                                    curry13, x+twidth, curry13, x+twidth*3/4 - lwidth - before2, curry13, x+twidth*3/4 - lwidth - before2, curry22-lwidth,
                                    x+twidth/2, curry22-lwidth, x+twidth/2, curry22-lwidth,
                                    fill='#D35E60', smooth=1, outline='black')
            curry11 += lwidth
            curry22 -= lwidth
            curry13 += lwidth
        if vals[4] != 0:
            lwidth = vals[4] * theight + 2
            self.canvas.create_polygon(x, curry21, x, curry21, x+twidth/4 + lwidth + before3, curry21, x+twidth/4 + lwidth + before3,
                                    curry12 + lwidth, x+twidth/2, curry12+lwidth, x+twidth/2, curry12+lwidth, x+twidth/2,
                                    curry12, x+twidth/2, curry12, x+twidth/4 + before3, curry12, x+twidth/4 + before3, curry21-lwidth,
                                    x, curry21-lwidth, x, curry21-lwidth,
                                    fill='#9067A7', smooth=1, outline='black')
            self.canvas.create_polygon(x+twidth/2, curry12, x+twidth/2, curry12, x+twidth*3/4 + lwidth + before1, curry12, x+twidth*3/4 + lwidth + before1,
                                    curry23 - lwidth, x+twidth, curry23-lwidth, x+twidth, curry23-lwidth, x+twidth,
                                    curry23, x+twidth, curry23, x+twidth*3/4 + before1, curry23, x+twidth*3/4 + before1, curry12+lwidth,
                                    x+twidth/2, curry12+lwidth, x+twidth/2, curry12+lwidth,
                                    fill='#9067A7', smooth=1, outline='black')

def getflow(samfile, vcffile, outfile, maxdist=600, skipindels=True):
    try:
        import pysam
    except:
        return 0
    snps = []
    vcf = open(vcffile)
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
            if skipindels:
                variantlens = set()
                for i in alt.split(','):
                    variantlens.add(len(i))
                if len(variantlens) == 1 and len(ref) in variantlens:
                    snps.append(aninstance)
            else:
                snps.append(aninstance)
    sam = pysam.Samfile(samfile, 'rb')
    out = open(outfile, 'w')

    for snp in snps:
        newvars = []
        for i in thevars:
            if i[0] < snp.pos - maxdist:
                out.write('POS\t' + str(i[0]) + '\t' + '\t'.join(i[1:-1]) + '\n')
                countDict = {}
                for j in i[-1]:
                    if i[-1][j] in countDict:
                        countDict[i[-1][j]] += 1
                    else:
                        countDict[i[-1][j]] = 1
                for j in countDict:
                    out.write('HAPLO\t' + '\t'.join(map(str, j)) + '\t' + str(countDict[j]) + '\n')
            else:
                newvars.append(i)
        thevars = newvars
        for pileupcolumn in sam.pileup(snp.chrom, snp.pos, snp.pos + 1):
            if pileupcolumn.pos == snp.pos - 1:
                variants = set(snp.alt.split(','))
                varlength = len(snp.ref)
                for pileupread in pileupcolumn.pileups:
                    rvar = pileupread.alignment.seq[pileupread.qpos:pileupread.qpos +
                           pileupread.alignment.overlap(pileupcolumn.pos, pileupcolumn.pos + varlength)]
                    
    out.close()






if sys.argv[1] == '-cl':
    getflow(sys.argv[2], sys.argv[3], sys.argv[4])
else:
    root = Tk()
    root.title('strainVis')
    root.option_add('*Font', 'Courier 10')
    root.option_add("*Background", "#E0E0FF")
    root.option_add("*Foreground", "#2B3856")
    root.option_add("*Listbox.Background", '#FFFFFF')
    root.option_add("*Scrollbar.Background", "#C0C0FF")
    root.option_add("*Entry.Background", "#FFFFFF")
    root.geometry("600x400")
    app = App(root)
    root.mainloop()