__author__ = 'mjsull'
from Tkinter import *
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox
import random
import time


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
        self.currflows = set()
        if len(sys.argv) > 1:
            self.flowfile = sys.argv[1]
            self.load_flow()

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
                self.canvas.delete(str(i))
                toremove.add(i)
        thetime = time.time()
        for i in toremove:
            self.currflows.remove(i)
        thetime = time.time()
        count = 0
        for i in range(max([indexa, 0]), indexb + 1):
            if not i in self.currflows:
                count += 1
                self.currflows.add(i)
                for j in self.flowlist[i]:
                    if j[0] == 0: # if single forward
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1] - 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - 1) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k] + 1) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                    elif j[0] == 1: # if single reverse
                        self.canvas.create_line([j[1][0] * self.xmod, j[2][0] * self.ymod + self.ypos1, j[1][1] * self.xmod, j[2][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1]), 2):
                            self.canvas.create_line([(j[1][k-1] - 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1, (j[1][k-1] + 1) * self.xmod, j[2][k-1] * self.ymod + self.ypos1,
                                                     (j[1][k] - 1) * self.xmod, j[2][k] * self.ymod + self.ypos1, (j[1][k] + 1) * self.xmod, j[2][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][k] * self.xmod, j[2][k] * self.ymod + self.ypos1, j[1][k+1] * self.xmod, j[2][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                    elif j[0] == 2: # if pair F F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        if j[5]:
                            pass
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=str(i))
                    elif j[0] == 3: # if pair F R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        if j[5]:
                            pass
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=str(i))
                    elif j[0] == 4: # if pair R F
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=LAST, arrowshape=(5, 5, 0), tags=str(i))
                        if j[5]:
                            pass
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=str(i))
                    elif j[0] == 5: # if pair R R
                        self.canvas.create_line([j[1][0][0] * self.xmod, j[2][0][0] * self.ymod + self.ypos1, j[1][0][1] * self.xmod, j[2][0][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][0]), 2):
                            self.canvas.create_line([(j[1][0][k-1] - 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1, (j[1][0][k-1] + 1) * self.xmod, j[2][0][k-1] * self.ymod + self.ypos1,
                                                     (j[1][0][k] - 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1, (j[1][0][k] + 1) * self.xmod, j[2][0][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][0][k] * self.xmod, j[2][0][k] * self.ymod + self.ypos1, j[1][0][k+1] * self.xmod, j[2][0][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        self.canvas.create_line([j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][1] * self.xmod, j[2][1][1] * self.ymod + self.ypos1],
                                                width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        for k in range(2, len(j[1][1]), 2):
                            self.canvas.create_line([(j[1][1][k-1] - 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1, (j[1][1][k-1] + 1) * self.xmod, j[2][1][k-1] * self.ymod + self.ypos1,
                                                     (j[1][1][k] - 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1, (j[1][1][k] + 1) * self.xmod, j[2][1][k] * self.ymod + self.ypos1],
                                                     smooth=True, width=j[3] * self.ymod / 8, fill=j[4], tags=str(i))
                            self.canvas.create_line([j[1][1][k] * self.xmod, j[2][1][k] * self.ymod + self.ypos1, j[1][1][k+1] * self.xmod, j[2][1][k+1] * self.ymod + self.ypos1],
                                                    width=j[3] * self.ymod, fill=j[4], arrow=FIRST, arrowshape=(5, 5, 0), tags=str(i))
                        if j[5]:
                            pass
                        else:
                            self.canvas.create_line((j[1][0][-1] * self.xmod, j[2][0][-1] * self.ymod + self.ypos1, (j[1][0][-1] + 1) * self.xmod, j[2][0][-1] * self.ymod + self.ypos1,
                                                 (j[1][1][0] - 1) * self.xmod, j[2][1][0] * self.ymod + self.ypos1, j[1][1][0] * self.xmod, j[2][1][0] * self.ymod + self.ypos1),
                                                 smooth=True, width=int(j[3]/4), fill=j[4], dash=(5,2), tags=str(i))
            xpos = (i + 1) * self.xmod * 4
            if i >= indexa:
                suppositions.append(xpos)
                positions.append(self.poslist[i][0])
                thetext = str(self.poslist[i][0]) + '\n' + '\n'.join(self.poslist[i][1:])
                self.canvas.create_text(xpos + 4, self.ymod * (self.ypos1 + self.flowend) + 5, anchor=NW, text=thetext, font=self.customFont, tags='top')
        self.canvas.create_rectangle(x1+10, self.ypossnp + 20, x2 - 10, self.ypossnp, tags='top', fill='#E1974C')
        if len(positions) > 1:
            for i in range(len(positions)):
                xpos = x1+15 + int((positions[i] - positions[0]) * 1.0 / (positions[-1] - positions[0]) * (x2 - x1 - 30))
                self.canvas.create_line(suppositions[i], self.ymod * (self.ypos1 + self.flowend), suppositions[i], self.ypos1 - 5, xpos, self.ypossnp + 20, xpos, self.ypossnp, tags='top', width=2)
            self.canvas.create_rectangle(x1 + 10, self.yposref + 20, x2 - 10, self.yposref, tags='top', fill='#7293CB')
            starto = x1 + 10 + int(positions[0] * 1.0 / self.reflength * (x2 - x1 - 20))
            endo = x1 + 10 + int(positions[-1] * 1.0 / self.reflength * (x2 - x1 - 20))
            self.canvas.create_rectangle(starto, self.yposref + 20, endo, self.yposref, tags='top', fill='#E1974C')
            self.canvas.create_text(x1 + 10, self.ypossnp - 2, anchor=SW, text='SNP block start..stop: ' + str(positions[-1]) + '..' + str(positions[0]), font=self.customFont, tags='top')

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
                self.poslist.append([int(line.split()[1].split(',')[0])] + line.split()[1].split(',')[1:])
        flowfile.close()
        flowfile = open(self.flowfile)
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
        flowfile.close()
        stacks = [stacker[:] for i in range(len(initiallist))]
        colordict = {}
        lastcol = 0
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
                    s = 0.5 #random.random()
                    l = 0.5 #random.random()
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
    maxeachvar = [] # maximum count for most prevalent variant, second most prevalent variant etc.
    maxerror = 0 # maximum number of reads with no called variant at a single site
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
            #out.write('V ' + str(pos) + ' ' + ref + ',' + alt + '\n')
            snps.append(aninstance)
    sam = pysam.Samfile(samfile, 'rb')
    reads = {}
    varlist = []
    lastpos = None
    consflow = {}
    flownum = 0
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
        newflows3 = []
        for i in newflows2:
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
            print 'ding', flow, int(i.split(',')[-1]), i.split(',')[0]
            print consflow
            bestcons = 0
            for j in consflow:
                print 'dong', j, consflow[j]
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
                        tempnewflow = newflow
                        print 'fuck'
                    foundone = True
            if foundone:
                tempnum = consflow[changeflow][0]
                tempnum2 = consflow[changeflow][1]
                print 'dongas', changeflow
                del consflow[changeflow]
                if len(changeflow) > (flow):
                    tempnewflow += list(changeflow[len(flow):])
                else:
                    tempnewflow += list(flow[len(changeflow):])
                print 'dangas', tempnewflow
                consflow[tuple(tempnewflow)] = (tempnum, max([tempnum2, int(i.split(',')[-1])]))
                print 'dang', tempnum, changeflow
            else:
                flownum += 1
                consflow[flow] = (flownum, int(i.split(',')[-1]))
                tempnum = flownum
                print 'wang', flow, flownum
            newflows3.append(i + ',' + str(tempnum) + '\n')
        newflows3.sort(key=lambda x: int(x.split(',')[-2]), reverse=True)
        newflows3.sort(key=lambda x: orderflow(x))
        newflows3.sort(key=lambda x: int(x.split(',')[0]))
        for i in newflows3:
            out.write('F ' + i )
            #print i
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
        #print 'ding', flow, int(i.split(',')[-1])
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