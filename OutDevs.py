#!/usr/bin/env python

import numpy as np
import BS_config as BS
from platform import system

import datetime
from PyQt5.QtCore import (QFile, Qt, QFileInfo, QTextStream, QDir, QSettings,QRectF)
from PyQt5.QtGui import QColor,QPalette, QPixmap, QFont, QPainter, QPen, QIcon
from PyQt5.QtWidgets import (QAction, QFileDialog, QLabel, QMessageBox, QApplication, QStyleFactory,
                             QScrollArea, QSizePolicy, QVBoxLayout, QWidget, QToolBar)

# class object that will handle output to a file
# will subclass this for output to RTF/ASCII/PDF
# noinspection PyMethodMayBeStatic
class Filedev():

    def __init__(self):
# create the reference points for instance variables that will hold all the data to be processed by this instance
        self.seqs = np.array([[' ', ' '], [' ', ' ']], dtype=str)
        self.cols = np.full(self.seqs.shape, 0, dtype=np.int32)
        self.seqnames =[]
        self.no_seqs = 0
        self.LHprenums = []
        self.RHprenums = []
        self.seqlens = []
        self.startnums = []
        self.lcs = []
        self.interlines = 0
        self.file_filter = ("All files (*);;PDF files (*.pdf);;RTF files (*.rtf);;ASCII files (*.txt)")
        self.file = None
        self.outstream = None
        self.parent = None


    def rgb(self, col):
        return col.red(), col.green(), col.blue()

    def open_output_file(self):
        QApplication.restoreOverrideCursor()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        TDialog = QFileDialog()
        fileName, _ = TDialog.getSaveFileName(self.parent,"Save file as:", BS.lastdir, self.file_filter, options=options)
        if fileName:
            self.file = QFile(fileName)
            BS.lastdir = QFileInfo(fileName).absolutePath()
            if not self.file.open(QFile.WriteOnly | QFile.Text):  # open to write
                mb = QMessageBox()
                mb.setTextFormat(Qt.RichText)
                mb.setText("<p style='font-size: 18pt'>Open File error</p>"
                           "<p style='font-size: 14pt; font-weight: normal'> Can't open file <i>{}</i> for writing.<br><br>"
                           " File error was: \"{}\".</p>".format(fileName, self.file.errorString()))
                mb.setIcon(QMessageBox.Warning)
                mb.exec()
                return False
            else:
                self.outstream = QTextStream(self.file)
                return True
        else:
            return False

    def make_lowercase(self, rulerflag):
        if self.lcs[0]:
            self.seqs[self.cols==0]=np.char.lower(self.seqs[self.cols==0])
        if self.lcs[1]:
            self.seqs[self.cols==1]=np.char.lower(self.seqs[self.cols==1])
        if self.lcs[2]:
            self.seqs[self.cols==2]=np.char.lower(self.seqs[self.cols==2])
        if self.lcs[3]:
            self.seqs[self.cols==3]=np.char.lower(self.seqs[self.cols==3])

    def exit(self):
        if self.file.isWritable():
            self.file.close()
        self.file=None

class RTFdev(Filedev):

    def __init__(self, fname):
        super(RTFdev, self).__init__()
        self.file_filter = ("RTF files (*.rtf);;All files (*)")
        self.Alignment = fname
        settings = QSettings("Boxshade", "Boxshade")
        self.bgds = settings.value("PSbgds")
        self.fgds = settings.value("PSfgds")
        self.FSize = settings.value("PSFsize", type=int)
        self.lcs = settings.value("PSLCs", type=bool)
        simflag = settings.value("simflag", type=bool)
        globalflag = settings.value("globalflag", type=bool)
        if not simflag:
            self.fgds[2] = self.fgds[0]
            self.bgds[2] = self.bgds[0]
            self.lcs[2] = self.lcs[0]
        if not globalflag:
            self.fgds[3] = self.fgds[1]
            self.bgds[3] = self.bgds[1]
            self.lcs[3] = self.lcs[1]

        dev_miny = 1.0
        dev_maxy = 15000.0
        dev_ysize = self.FSize * 20.0
        self.lines_per_page = int((dev_maxy - dev_miny) / dev_ysize)

    def graphics_init(self):
        if self.open_output_file():
            self.outstream << '{\\rtf1\\ansi\\deff0\n{\\fonttbl{\\f0\\fmodern Courier New;}}\n'
            self.outstream << '{{\\info{{\\author BOXSHADE}}{{\\title {}}}}}\n'.format(self.Alignment)
            self.outstream << '{\\colortbl\n'
            for i in range(4):
                self.outstream << '\\red{}\\green{}\\blue{};'.format(*self.rgb(self.fgds[i]))

            self.outstream << '\\red0\\green0\\blue0;' # equivalent to fgds[4]
            for i in range(4):
                self.outstream << '\\red{}\\green{}\\blue{};'.format(*self.rgb(self.bgds[i]))

            self.outstream << '\\red255\\green255\\blue255;}\n' # equivalent to bgds[4]
            self.outstream << '\\paperw11880\\paperh16820\\margl1000\\margr500\n'
            self.outstream << '\\margt910\\margb910\\sectd\\cols1\\pard\\plain\n'
            self.outstream << '\\fs{}\n\\b\n'.format(self.FSize * 2)
            self.outstream.flush()
            return True
        else:
            return False

    def set_colour(self, c):
        colno = c
        self.outstream << '\n\\chshdng0\\chcbpat{0}\\cb{0}\\cf{1} '.format(5+colno, colno)

    def char_out(self,c):
        self.outstream << c

    def string_out(self, str):
        self.outstream << str

    def newline(self):
        self.outstream << '\n\\cb{}\\cf{} \\line'.format(9, 4)

    def newpage(self):
        self.outstream << '\\page'

    def exit(self):
        self.outstream << '\\b0}\n'
        super().exit()


# noinspection PyMethodMayBeStatic
class ImageDisp(QWidget):
    def __init__(self, mw, parent=None):
        super(ImageDisp, self).__init__(parent)
        self.MW = mw
        QApplication.setStyle(QStyleFactory.create("Fusion"))
        QApplication.setPalette(QApplication.style().standardPalette())
        self.file_filter = ("PNG files (*.png);;All files (*)")
        self.setWindowFlag(Qt.Window, True)
        self.scaleFactor = 0.0
        self.imageLabel = QLabel()
        self.imageLabel.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self.imageLabel.setScaledContents(True)

        self.scrollArea = QScrollArea()
        self.scrollArea.setBackgroundRole(QPalette.Dark)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setAlignment(Qt.AlignLeft| Qt.AlignTop)
        central = QVBoxLayout()
        self.tb = QToolBar()
        central.setMenuBar(self.tb)
        central.addWidget(self.scrollArea)
        self.setLayout(central)
        self.resize(930,400)

        self.createActions()
        self.createToolBar()

    def closeEvent(self, event):
        self.MW.viewList.remove(self)
        event.accept()

    def createToolBar(self):
        self.tb.addAction(self.zoomInAct)
        self.tb.addAction(self.zoomOutAct)
        self.tb.addSeparator()
        self.tb.addAction(self.saveAct)


    def createActions(self):
        root = QFileInfo(__file__).absolutePath()
        self.zoomInAct = QAction(QIcon(root + '/images/zoomin.png'),"Zoom In (33%)", self, shortcut="Ctrl++",
                enabled=False, triggered=self.zoomIn)

        self.zoomOutAct = QAction(QIcon(root + '/images/zoomout.png'),"Zoom Out (33%)", self, shortcut="Ctrl+-",
                enabled=False, triggered=self.zoomOut)

        self.saveAct = QAction(QIcon(root + '/images/png.png'),"Save image", self, shortcut = "Ctrl+s",
                               enabled = False, triggered = self.savePNG)

    def zoomIn(self):
        self.scaleImage(1.33)

    def zoomOut(self):
        self.scaleImage(0.75)

    def updateActions(self):
        self.scaleFactor = 1.0
        self.zoomInAct.setEnabled(True)
        self.zoomOutAct.setEnabled(True)
        self.saveAct.setEnabled(True)

    def scaleImage(self, factor):
        self.scaleFactor *= factor
        self.imageLabel.resize(self.scaleFactor * self.imageLabel.pixmap().size())

        self.adjustScrollBar(self.scrollArea.horizontalScrollBar(), factor)
        self.adjustScrollBar(self.scrollArea.verticalScrollBar(), factor)
        self.scrollArea.ensureVisible(30,30)
        self.zoomInAct.setEnabled(self.scaleFactor < 3.0)
        self.zoomOutAct.setEnabled(self.scaleFactor > 0.333)

    def adjustScrollBar(self, scrollBar, factor):
        scrollBar.setValue(int(factor * scrollBar.value()
                                + ((factor - 1) * scrollBar.pageStep()/2)))

    def savePNG(self):
        self.save2output_file("PNG")

    def save2output_file(self, fmt):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        TDialog = QFileDialog()
        fileName, _ = TDialog.getSaveFileName(None,"Save file as:", BS.lastdir, self.file_filter, options=options)
        if fileName:
            self.file = QFile(fileName)
            BS.lastdir = QFileInfo(fileName).absolutePath()
            if not self.file.open(QFile.WriteOnly):  # open to write
                mb = QMessageBox(self)
                mb.setTextFormat(Qt.RichText)
                mb.setText("<p style='font-size: 18pt'>Open File error</p>"
                           "<p style='font-size: 14pt; font-weight: normal'> Can't open file <i>{}</i> for writing.<br><br>"
                           " File error was: \"{}\".</p>".format(fileName, self.file.errorString()))
                mb.setIcon(QMessageBox.Warning)
                mb.exec()

                return False
            else:
                self.imageLabel.pixmap().save(self.file, fmt)
                return True
        else:
            return False


class Paintdev(Filedev):

    def __init__(self,mw):
        super(Paintdev, self).__init__()
        self.MW=mw
        self.top_mar = 30.0
        self.left_mar = 30.0
        settings = QSettings("Boxshade", "Boxshade")
        self.bgds = settings.value("PSbgds")
        self.fgds = settings.value("PSfgds")
        self.FSize = settings.value("PSFsize", type=int)
        self.lcs = settings.value("PSLCs", type=bool)
        self.simflag = settings.value("simflag", type=bool)
        self.globalflag = settings.value("globalflag", type=bool)
        self.interlines = settings.value("interlines", type=int)
        self.outlen = settings.value("outlen", type=int)
        self.snameflag = settings.value("snameflag", type=bool)
        self.LHsnumsflag = settings.value("LHsnumsflag", type=bool)
        self.RHsnumsflag = settings.value("RHsnumsflag", type=bool)

        if not self.simflag:
            self.fgds[2] = self.fgds[0]
            self.bgds[2] = self.bgds[0]
            self.lcs[2] = self.lcs[0]
        if not self.globalflag:
            self.fgds[3] = self.fgds[1]
            self.bgds[3] = self.bgds[1]
            self.lcs[3] = self.lcs[1]
        self.fgds.append(QColor(0,0,0))
        self.bgds.append(QColor(255, 255, 255))
        self.dev_miny = self.top_mar
        self.dev_minx = self.left_mar
        if system() == "Darwin":
            self.dev_xsize = self.FSize * 0.8
            self.dev_ysize = self.FSize
        else:
            self.dev_xsize = self.FSize * 0.9
            self.dev_ysize = self.FSize * 1.2
        self.lines_per_page = 10000

    def graphics_init(self):
# For this to work, I have to calculate how big a drawing I am going to make based on data passed in from calling routine
# For this reason, I have moved gr_out.graphics_init() to the end of prep_out in the calling routine
# as this "device" needs to know what it is drawing in order to initialise itself.
#
        blocks = (self.seqs.shape[1]//self.outlen)+1
        canvas_height = int(self.top_mar+(self.dev_ysize*(blocks*(self.no_seqs+self.interlines)-self.interlines))+self.top_mar+0.5)
        nchars = self.outlen
        if self.snameflag:
            nchars += 1+len(self.seqnames[0])
        if self.LHsnumsflag:
            nchars += 1+len(self.LHprenums[0][0])
        if self.RHsnumsflag:
            nchars += 1+len(self.RHprenums[0][0])
        canvas_width = int(self.left_mar + (self.dev_xsize*(nchars))+self.left_mar+0.5)
        self.canvas = QPixmap(canvas_width,canvas_height)
        if self.canvas.isNull():
            return False
        self.canvas.fill(QColor(255,255,255))

        BS.monofont.setPointSize(self.FSize)
        BS.monofont.setWeight(QFont.Bold)
        self.pen = QPen()
        self.paint = QPainter(self.canvas)
        self.paint.setFont(BS.monofont)
        self.paint.setRenderHint(QPainter.Antialiasing, True)
        self.paint.setRenderHint(QPainter.TextAntialiasing, True)

        self.xpos = self.top_mar
        self.ypos = self.left_mar
        self.act_col = 4
        return True

    def set_colour(self, c):
        self.act_col = c

    def char_out(self, ch):
        myrect=QRectF(self.xpos,self.ypos, self.dev_xsize, self.dev_ysize)
        self.paint.fillRect(myrect,self.bgds[self.act_col])
        self.paint.setPen(self.fgds[self.act_col])
        self.paint.drawText(myrect, Qt.AlignCenter, ch)
        self.xpos += self.dev_xsize

    def string_out(self, str):
        for i in range(len(str)):
            self.char_out(str[i])

    def newline(self):
        self.xpos = self.dev_minx
        self.ypos += self.dev_ysize

    def newpage(self):
        pass

    def exit(self):
        self.paint.end()


class PSdev(Filedev):
    pscc = ['C0', 'C1', 'C2', 'C3', 'C4']
    ctypes = ['% -- different residues\n', "% -- identical residues\n", "% -- similar residues\n", "% -- conserved residues\n", "% -- normal text\n"]
    def __init__(self, fname):
        super(PSdev, self).__init__()
        self.file_filter = ("PS files (*.ps);;All files (*)")
        self.Alignment = fname
        settings = QSettings("Boxshade", "Boxshade")
        self.bgds = settings.value("PSbgds")
        self.fgds = settings.value("PSfgds")
        self.FSize = settings.value("PSFsize", type=int)
        self.lcs = settings.value("PSLCs", type=bool)
        simflag = settings.value("simflag", type=bool)
        globalflag = settings.value("globalflag", type=bool)
        self.landscapeflag = settings.value("PSlandscapeflag", type=bool)
        self.outlen = settings.value("outlen", type=int)
        self.snameflag = settings.value("snameflag", type=bool)
        self.LHsnumsflag = settings.value("LHsnumsflag", type=bool)
        self.RHsnumsflag = settings.value("RHsnumsflag", type=bool)

        if not simflag:
            self.fgds[2] = self.fgds[0]
            self.bgds[2] = self.bgds[0]
            self.lcs[2] = self.lcs[0]
        if not globalflag:
            self.fgds[3] = self.fgds[1]
            self.bgds[3] = self.bgds[1]
            self.lcs[3] = self.lcs[1]
        self.dev_minx = 30.0
        self.dev_miny = 30.0
        if self.landscapeflag:
            self.dev_maxx = 800.0
            self.dev_maxy = 545.0
        else:
            self.dev_maxx = 575.0
            self.dev_maxy = 760.0
        self.dev_xsize = self.FSize * 0.7
        self.dev_ysize = self.FSize
        self.xpos = self.dev_minx
        self.ypos = self.dev_maxy - self.dev_ysize
        self.lines_per_page = int((self.dev_maxy - self.dev_miny) / self.dev_ysize)

    def PSrgb(self, col):
        col = [x/255 for x in self.rgb(col)]
        return col

    def PageSetup(self, pn):
        self.outstream << "\n%%Page: {} {}\n%%BeginPageSetup\npsetup\n%%EndPageSetup\n".format(pn, pn)

    def psfp(self, num, dp):
        s = '{:.{}f}'.format(num,dp)
        s = s.rstrip('0').rstrip('.') if '.' in s else s
        self.outstream << s
        self.outstream << ' '
        self.count += len(s)+1

    def close_sb(self):
        if self.save_sb != []:
            sl = len(self.save_sb)
            if (self.count+sl) > 200:
                self.outstream << "\n"
                self.count = 0
            xsl = 'S' if sl > 1 else 'C'
            self.outstream << "({}){} ".format(''.join(self.save_sb), xsl)
            self.count += 4+len(self.save_sb)
            self.save_sb = []

    def add_sb(self, ch):
        self.save_sb.append(ch)
        sl = len(self.save_sb)
        if (self.count+sl) > 200:
            xsl = 'S' if sl > 1 else 'C'
            self.outstream << "({}){} ".format(''.join(self.save_sb), xsl)
            self.outstream << "\n"
            self.count = 0
            self.save_sb = []

    def graphics_init(self):
        nchars = self.outlen
        if self.snameflag:
            nchars += 1 + len(self.seqnames[0])
        if self.LHsnumsflag:
            nchars += 1 + len(self.LHprenums[0][0])
        if self.RHsnumsflag:
            nchars += 1 + len(self.RHprenums[0][0])
        line_length = self.dev_xsize * nchars
        if line_length > (self.dev_maxx-self.dev_minx):
            mb = QMessageBox()
            mb.setTextFormat(Qt.RichText)
            mb.setText("<p style='font-size: 18pt'>Picture too wide for page!</p>"
                       "<p style='font-size: 14pt; font-weight: normal'> Be aware, at your current settings, the PostScript image "
                       "will be wider than the page and will be clipped.<br>"
                       "The page width is {:n} pixels and your output would have a width of {:n} pixels.<br><br>"
                       "Do you want to continue?</p>".format((self.dev_maxx-self.dev_minx), line_length))
            mb.setIcon(QMessageBox.Information)
            mb.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
            mb.setDefaultButton(QMessageBox.No)
            ret = mb.exec()
            if ret == QMessageBox.No:
                return False

        if self.open_output_file():
            self.outstream << "%!PS-Adobe-2.0\n"
            self.outstream << "%%Creator: PyBoxshade\n"
            self.outstream << "%%Title: BOXSHADE document from: {} \n".format(self.Alignment)
            d = datetime.datetime.now()
            self.outstream << "%%CreationDate: {:%H:%M:%S %B %d, %Y}\n".format(d)
            self.outstream << "%%Pages: (atend)\n"
            if self.landscapeflag:
                self.outstream << "%%BoundingBox: {:.0f} {:.0f} {:.0f} {:.0f}\n".format(self.dev_miny-1,self.dev_minx-1,self.dev_maxy+1,self.dev_maxx+1)
                self.outstream << "%%Orientation: landscape\n"
            else:
                self.outstream << "%%BoundingBox: {:.0f} {:.0f} {:.0f} {:.0f}\n".format(self.dev_minx-1,self.dev_miny-1,self.dev_maxx+1,self.dev_maxy+1)
                self.outstream << "%%Orientation: portrait\n"

            self.outstream << "%%PaperSize: a4\n%%DocumentNeededFonts: Courier-Bold\n%%DocumentData: Clean7Bit\n"
            self.outstream << "%%LanguageLevel: 1\n%%EndComments\n%%BeginProlog\n"
            self.outstream << "/bd { bind def } bind def /xd { exch def } bd\n"
            self.outstream << "%\n% custom color selection\n%\n% grayscale:\n%\n% '<gray> setgray'\n%\n"
            self.outstream << "% where<gray> is a real number between\n%  0.0 (black) and 1.0 (white)\n%\n"
            self.outstream << "% RGB (red/green/blue) colors:\n%\n%  '<r> <g> <b> setrgbcolor'\n%\n% each color component is a real'\n"
            self.outstream << "% number between 0.0 (zero intensity) and\n%  1.0 (max intensity)\n%\n% Change the following definitions for your needs!\n%\n"
            for i in range(4):
                self.outstream << self.ctypes[i]
                self.outstream << "/bg{} {{ {:.2f} {:.2f} {:.2f} setrgbcolor }} bd % background\n".format(str(i), *self.PSrgb(self.bgds[i]))
                self.outstream << "/fg{} {{ {:.2f} {:.2f} {:.2f} setrgbcolor }} bd % foreground\n".format(str(i), *self.PSrgb(self.fgds[i]))
            self.outstream << self.ctypes[4]
            self.outstream << "/bg4 { 1 setgray } bd % background\n"
            self.outstream << "/fg4 { 0 setgray } bd % foreground\n"

            self.outstream << "%\n% end of custom color selection\n%\n"
            self.outstream << "/px 0 def /py 0 def /fg {0 setgray} bd /bg {1 setgray} bd "
            self.outstream << "/C {{px py moveto gsave {:.1f} {:.1f} rmoveto {:.1f} 0 rlineto 0 {:.1f} rlineto {:.1f} 0 rlineto closepath ".format(-0.03*self.FSize, \
                                                                                                     -0.05*self.FSize, 0.70*self.FSize, self.FSize, -0.70*self.FSize)
            self.outstream << "bg fill grestore fg 0 2 rmoveto show /px px {:.2f} add def}} bd\n".format(self.dev_xsize)
            self.outstream << "/X {/px xd} bd /Y {/py xd} bd /A {Y X} bd /S {0 1 2 index length 1 sub {2 copy 1 getinterval C pop} for pop} bd\n"
            d = " 575 0 translate 90 rotate" if self.landscapeflag else ""
            self.outstream << "/psetup {{/Courier-Bold findfont {} scalefont setfont 120 currentscreen 3 -1 roll pop setscreen {}}} bd\n".format(self.FSize, d)
            for i in range(5):
                self.outstream << "/{} {{/bg {{bg{}}} bd /fg {{fg{}}} bd}} bd ".format(self.pscc[i], (i), str(i))
            self.outstream << "\n%%EndProlog\n%%%BeginSetup\nsave initgraphics"
            self.outstream << "\n%%EndSetup"
            self.act_page = 1
            self.PageSetup(self.act_page)
            self.last_pscl = ''
            self.act_col = 4
            self.new_x, self.new_y  = True, True
            self.count = 0
            self.save_sb = []
            return True
        else:
            return False

    def set_colour(self, c):
        self.act_col = c

    def char_out(self, ch):
        if self.pscc[self.act_col] != self.last_pscl:
            self.close_sb()
            self.last_pscl = self.pscc[self.act_col]
            self.outstream << self.last_pscl+" "
            self.count += len(self.last_pscl+" ")
        if self.new_x and self.new_y:
            self.close_sb()
            self.psfp(self.xpos, 1)
            self.psfp(self.ypos, 1)
            self.outstream << "A "
            self.count += 2
            self.new_x, self.new_y = False, False
        elif self.new_y:
            self.close_sb()
            self.psfp(self.ypos, 1)
            self.outstream << "Y "
            self.count += 2
            self.new_y = False
        elif self.new_x:
            self.close_sb()
            self.psfp(self.xpos, 1)
            self.outstream << "X "
            self.count += 2
            self.new_x = False
        self.add_sb(ch)
        self.xpos +=self.dev_xsize

    def string_out(self,str):
        for i in range(len(str)):
            self.char_out(str[i])

    def newline(self):
        self.close_sb()
        self.xpos = self.dev_minx
        self.ypos -= self.dev_ysize
        self.new_x, self.new_y = True, True

    def newpage(self):
        self.close_sb()
        self.xpos = self.dev_minx
        self.outstream << "showpage erasepage "
        self.count += len("showpage erasepage ")
        self.act_page += 1
        self.PageSetup(self.act_page)
        self.ypos = self.dev_maxy - self.dev_ysize
        self.new_x, self.new_y = True, True
        self.last_pscl = ' '

    def exit(self):
        self.close_sb()
        self.outstream << "showpage erasepage\n%%Trailer\nrestore\n"
        self.outstream << "%%Pages: {}\n".format(self.act_page)
        self.outstream << "%%EOF\n"
        super().exit()

class ASCIIdev(Filedev):

    def __init__(self, fname):
        super(ASCIIdev, self).__init__()
        self.file_filter = ("Text files (*.txt);;All files (*)")
        self.Alignment = fname
        settings = QSettings("Boxshade", "Boxshade")
        self.Achars = settings.value("ASCIIchars")
        self.Achars.append('L')
        simflag = settings.value("simflag", type=bool)
        globalflag = settings.value("globalflag", type=bool)
        if not simflag:
            self.Achars[2] = self.Achars[0]
        if not globalflag:
            self.Achars[3] = self.Achars[1]
        self.lcs = [(x.isalpha()) and (x == x.lower()) for x in self.Achars]
        self.lines_per_page = 9999

    def graphics_init(self):
        if self.open_output_file():
            self.current_char = ''
            self.outstream << "Alignment file: {}\n".format(self.Alignment)
            d = datetime.datetime.now()
            self.outstream << "Created by Boxshade: {:%H:%M:%S %B %d, %Y}\n\n".format(d)
            return True
        else:
            return False

    def set_colour(self, c):
        self.current_char = self.Achars[c]

    def char_out(self, ch):
        if self.current_char.upper() == 'L':
            self.outstream << ch
        else:
            self.outstream << self.current_char

    def string_out(self, str):
        self.outstream << str

    def newline(self):
        self.outstream << '\n'

    def newpage(self):
        pass

    def exit(self):
        self.outstream << '\n'
        super().exit()




