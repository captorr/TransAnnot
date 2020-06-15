"""
Author: Zhang Chengsheng, @2020.06.15
"""

__version__ = '0.0.1'
from tkinter import *
import pickle
from matplotlib import pyplot as plt
import matplotlib.lines as lines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class TAGUI:
    def __init__(self,db):
        self.db = db['gene2transcript']
        self.boxplot = boxplot2(db)
        self.UI = Tk()
        self._initWindow()

    def _initWindow(self):
        row1,row2,row3,rowSearch,rowSearchLable=2,3,10,1,0
        col1,col2,col3=0,7,14
        self.UI.title("TransAnnot Viewer V0.0.1")
        self.UI.geometry('568x381+200+200')
        self.geneLable = Label(self.UI,text="Gene Level")
        self.geneLable.grid(row=row1,column=col1)
        self.transcriptLable = Label(self.UI,text="Transcript Level")
        self.transcriptLable.grid(row=row1,column=col2)
        self.readLable = Label(self.UI,text="Read Level")
        self.readLable.grid(row=row1,column=col3)

        # searchbox
        checkInput = self.UI.register(self._inputCheck)
        self.searchLable = Label(self.UI, text="Gene Keyword Search")
        self.searchLable.grid(row=rowSearchLable, column=col1)
        self.searchBox = Entry(self.UI, width=65, textvariable=StringVar(), validate='key',validatecommand=(checkInput,'%P','%v','%W'))
        self.searchBox.grid(row=rowSearch, column=col1, columnspan=15, padx=(20, 0), pady=(2, 2), sticky=W + N + S + E)

        # GENE BOX
        self.geneListFrame = LabelFrame(self.UI,text='')
        self.geneListFrame.grid(row=row2,column=col1,padx=(20,0),pady=(0,10))
        geneScrollX = Scrollbar(self.geneListFrame, width=15, orient=HORIZONTAL)
        geneScrollX.pack(side=BOTTOM,fill=X)
        geneScrollY = Scrollbar(self.geneListFrame, width=15, orient=VERTICAL)
        geneScrollY.pack(side=RIGHT,fill=Y)
        self.geneList = Listbox(self.geneListFrame,width=20,height=10,exportselection=False,yscrollcommand=geneScrollY.set,xscrollcommand=geneScrollX.set)
        geneScrollX.config(command=self.geneList.xview)
        geneScrollY.config(command=self.geneList.yview)
        for i in self.db:
            self.geneList.insert(END,i)
        self.geneList.pack(fill=BOTH)

        # transcript box
        self.transcriptFrame = LabelFrame(self.UI,text='')
        self.transcriptFrame.grid(row=row2,column=col2,padx=(20,0),pady=(0,10))
        transcriptScrollX = Scrollbar(self.transcriptFrame, width=15, orient=HORIZONTAL)
        transcriptScrollX.pack(side=BOTTOM, fill=X)
        transcriptScrollY = Scrollbar(self.transcriptFrame, width=15, orient=VERTICAL)
        transcriptScrollY.pack(side=RIGHT, fill=Y)
        self.transcriptList = Listbox(self.transcriptFrame,width=20,height=10,exportselection=False,yscrollcommand=transcriptScrollY.set,xscrollcommand=transcriptScrollX.set)
        transcriptScrollX.config(command=self.transcriptList.xview)
        transcriptScrollY.config(command=self.transcriptList.yview)
        self.transcriptList.pack(fill=BOTH)

        # read box
        self.readFrame = LabelFrame(self.UI, text='')
        self.readFrame.grid(row=row2, column=col3, padx=(20, 0), pady=(0, 10))
        readScrollX = Scrollbar(self.readFrame, width=15, orient=HORIZONTAL)
        readScrollX.pack(side=BOTTOM, fill=X)
        readScrollY = Scrollbar(self.readFrame, width=15, orient=VERTICAL)
        readScrollY.pack(side=RIGHT, fill=Y)
        self.readList = Listbox(self.readFrame, width=20, height=10, exportselection=False,yscrollcommand=readScrollY.set, xscrollcommand=readScrollX.set)
        readScrollX.config(command=self.readList.xview)
        readScrollY.config(command=self.readList.yview)
        self.readList.pack(fill=BOTH)

        ## 列表绑定函数
        self.geneList.bind('<<ListboxSelect>>',self._clickGeneList)
        self.transcriptList.bind('<<ListboxSelect>>', self._clickTranscriptList)
        ## 按钮
        self.drawButton = Button(self.UI,text='Show Figure',command=self._draw)
        self.drawButton.grid(row=row3,column=col1,columnspan=15,padx=(20,0),sticky=EW)

        #self.searchBox.bind('<key-input>',print(self.KEYWORD))

    def _inputCheck(self,content,reason,name):
        self._refreshGenebox(content)
        #print(content,reason,name)
        return True

    def _refreshGenebox(self,keyword):
        self.geneList.delete(0,END)
        for g in self.db:
            if keyword == '':
                self.geneList.insert(END,g)
            elif keyword.upper() in g.upper():
                self.geneList.insert(END,g)

    def _clickGeneList(self,event):
        self.transcriptList.delete(0,END)
        self.readList.delete(0, END)
        g = self.geneList.get(self.geneList.curselection())
        for t in self.db[g]:
            self.transcriptList.insert(END,t)
            for r in self.db[g][t]:
                self.readList.insert(END,r)

    def _clickTranscriptList(self,event):
        self.readList.delete(0,END)
        g = self.geneList.get(self.geneList.curselection())
        t = self.transcriptList.get(self.transcriptList.curselection())
        for r in self.db[g][t]:
            self.readList.insert(END,r)

    def _draw(self):
        g = self.geneList.get(self.geneList.curselection()) if self.geneList.curselection() else ''
        t = self.transcriptList.get(self.transcriptList.curselection()) if self.transcriptList.curselection() else ''
        r = self.readList.get(self.readList.curselection()) if self.readList.curselection() else ''
        plt = self.boxplot.draw_gene(g,t,r)
        x,y = plt.get_size_inches()
        draw = Toplevel()
        draw.geometry('1000x500+10+10')
        #draw.geometry('{}x{}+10+10'.format(int(x*100),int(y*100)))
        title = r if r else ( t if t else (g if g else ''))
        draw.title(title)
        p = FigureCanvasTkAgg(plt,draw)
        p.draw()
        p.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)

    def show(self):
        self.UI.mainloop()


class boxplot2:
    def __init__(self,db):
        self.db = db

    def draw_gene(self,gene,transcript,read):
        reads, anno, g2t = self.db['reads'], self.db['annotation'], self.db['gene2transcript']
        test_list = []
        for i in g2t[gene]:
            for j in g2t[gene][i]:
                test_list.append(j)
        test_dict = {}
        for i in test_list:
            if transcript and i not in g2t[gene][transcript]:
                continue
            if read and i != read:
                continue
            test_dict[i] = reads[i]
        anno_db = anno[gene]
        fig = self.normal_boxplot(test_dict, anno_db=anno_db, save=0, TAUI=1)
        return fig

    def normal_boxplot(self,draw_db, anno_db, save=0, fusion=0, TAUI=0):
        y_num = len(draw_db) + len(anno_db)
        # path = outputdir
        width = 0.5
        colTable = {0: 'gray', 10: '#3B3D42', 1: '#1975C4', 2: '#E68E29', 3: '#199F15', 4: '#9F3241', 5: '#18D6CB',
                    6: '#D670D2', 7: '#8BB6D6', 8: '#C6D611', 9: '#9425D6'}
        col_id_list, col_anno_list = [], []
        col_idx = 0
        colDict = {}
        all_lr = []
        for id in draw_db:
            col_id_list.append(draw_db[id]['info'][1])
            all_lr += draw_db[id]['start']
            all_lr += draw_db[id]['end']
        left, right = min(all_lr), max(all_lr)
        for transcript in anno_db:
            col_anno_list.append(transcript)
            all_lr += anno_db[transcript]['start']
            all_lr += anno_db[transcript]['end']
        for i in col_anno_list:
            if i not in colDict:
                colDict[i] = 10
        for i in col_id_list:
            if i in colDict:
                if colDict[i] == 10:
                    col_idx += 1
                    colDict[i] = min(col_idx, 10)
            else:
                colDict[i] = 0
        r = max(all_lr) - min(all_lr)
        multi = min(int(r / 150000) + 1, 6)
        fig = plt.figure(figsize=(multi * 10, (1 + (multi - 1) * 0.5) * 5 * (int(y_num / 20) + 1)))
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
        y_ticks = [idx for idx, i in enumerate(draw_db)]
        y_tickslabels = [i.split('_')[-1] for i in draw_db]
        for idx1, id in enumerate(draw_db):
            this_lr = []
            this_lr += draw_db[id]['start']
            this_lr += draw_db[id]['end']
            colBase = colTable[colDict[draw_db[id]['info'][1]]]
            for idx2, i in enumerate(draw_db[id]['start']):
                col = colBase
                if not idx2 or idx2 == len(draw_db[id]['anno']) - 1:
                    if draw_db[id]['anno'][idx2] not in ['左边缘外显子', '右边缘外显子']:
                        col = 'red'
                else:  # idx2 and idx2 != len(draw_db[id]['anno']) -1:
                    if fusion:
                        if draw_db[id]['anno'][idx2] not in ['', '左边缘外显子', '右边缘外显子']:
                            col = 'red'
                    elif draw_db[id]['anno'][idx2]:
                        col = 'red'
                length = abs(int(draw_db[id]['start'][idx2]) - int(draw_db[id]['end'][idx2]))
                rect = plt.Rectangle((int(i), idx1 - width / 2), length, width, color=col, joinstyle='miter',
                                     capstyle='butt', alpha=1)
                ax.add_patch(rect)
            ax.add_line(lines.Line2D((min(this_lr), max(this_lr)), (idx1, idx1), linewidth=0.5, solid_capstyle='butt',
                                     solid_joinstyle='miter', color='black', alpha=0.5))
        for idx1, transcript in enumerate(sorted(anno_db)):
            start = anno_db[transcript]['start']
            end = anno_db[transcript]['end']
            colBase = 'black' if fusion else colTable[colDict[transcript]]
            r = start + end
            s = idx1 + len(draw_db) + 1
            y_ticks.append(s)
            y_tickslabels.append(transcript)
            ax.add_line(
                lines.Line2D((min(r), max(r)), (s, s), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',
                             color='black', antialiased=False, alpha=0.5))
            for idx2, k in enumerate(start):
                length = abs(start[idx2] - end[idx2])
                rect = plt.Rectangle((min(start[idx2], end[idx2]), s - width / 2), length, width, color=colBase,
                                     joinstyle='miter', capstyle='butt', alpha=1)
                ax.add_patch(rect)

        ax.set_xticks(
            [min(all_lr), int(0.75 * min(all_lr) + 0.25 * max(all_lr)), int(0.5 * (min(all_lr) + max(all_lr))),
             int(0.25 * min(all_lr) + 0.75 * max(all_lr)), max(all_lr)])
        ax.set_xticklabels(
            [min(all_lr), int(0.75 * min(all_lr) + 0.25 * max(all_lr)), int(0.5 * (min(all_lr) + max(all_lr))),
             int(0.25 * min(all_lr) + 0.75 * max(all_lr)), max(all_lr)])
        ax.set_xlim([min(all_lr), max(all_lr)])
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tickslabels)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis=u'y', which=u'both', length=0)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
        ax.plot()
        # fig.show()
        if save:
            try:
                fig.savefig(save, bbox_inches='tight')
            except:
                pass
        else:
            if TAUI:
                return fig
            else:
                fig.show()


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def showGUI():
    dbfile = r'test.db.pickle'
    db = load_pickle(dbfile)
    A = TAGUI(db)
    A.show()


if __name__ == '__main__':
    showGUI()
