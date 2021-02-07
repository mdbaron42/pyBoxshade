#!/usr/bin/env python


from PyQt5.QtCore import Qt, QSize, QSettings, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog,
        QGridLayout, QPushButton, QRadioButton, QColorDialog,
        QDialogButtonBox, QFrame, QGroupBox, QHeaderView, QLabel, QLineEdit,
        QTabWidget, QVBoxLayout, QHBoxLayout, QWidget, QSpinBox, QStyleFactory,
        QDoubleSpinBox, QTableWidget, QTableWidgetItem, QMessageBox)
from PyQt5.QtGui import QFont
import numpy as np

import BS_config as BS
from platform import system

# a small function to make QSpinbox objects with the defined value/max/min/step
# if I create them with the normal function, a value >99 is not accepted.
def makeQSP(value=1, minimum=0, maximum=99, step=1):
    QSP=QSpinBox()
    QSP.setRange(minimum, maximum)
    QSP.setValue(value)
    QSP.setSingleStep(step)
    return QSP

# noinspection PyMethodMayBeStatic
class prefsDialog(QDialog):
    def __init__(self, no_seqs, consensnum=0, parent=None):
        super(prefsDialog, self).__init__(parent)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.changeStyle("Fusion")
        self.setFixedWidth(580)

        pDFont = self.font()
        pDFont.setFamily("Arial")
        if system() == "Darwin":
            pDFont.setPointSize(13)
        else:
            pDFont.setPointSize(10)
        self.setFont(pDFont)

        self.tabWidget = QTabWidget()
        self.GenTab = GeneralTab(no_seqs, consensnum)
        self.PSTab = PSSettingsTab()
        self.AsciiTab = ASCIITab()
        self.SimTab = simsTab()
        self.GrpTab = grpsTab()
        self.tabWidget.addTab(self.GenTab, "General")
        self.tabWidget.addTab(self.PSTab, "RTF/PS/PNG")
        self.tabWidget.addTab(self.AsciiTab, "Text (ASCII) output")
        self.tabWidget.addTab(self.SimTab, "Similarities")
        self.tabWidget.addTab(self.GrpTab, "Groups")
        self.GenTab.simflg_changed.connect(self.PSTab.simflag_changed)
        self.GenTab.grpflg_changed.connect(self.PSTab.grpflag_changed)
        self.GenTab.simflagbox.toggle()
        self.GenTab.simflagbox.toggle()
        self.GenTab.grpflagbox.toggle()
        self.GenTab.grpflagbox.toggle()
        self.proc_flag=False

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.exit)
        buttonBox.rejected.connect(self.reject)

        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.tabWidget)
        mainLayout.addWidget(buttonBox)
        self.setLayout(mainLayout)

        self.setWindowTitle("Boxshade Preferences")

    def exit(self):
        if self.GenTab.exit():
            self.AsciiTab.exit()
            self.PSTab.exit()
            self.SimTab.exit()
            self.GrpTab.exit()
            if self.GenTab.proc_flag or self.SimTab.proc_flag or self.GrpTab.proc_flag:
                self.proc_flag = True
            self.accept()

    def changeStyle(self, styleName):
        QApplication.setStyle(QStyleFactory.create(styleName))
        QApplication.setPalette(QApplication.style().standardPalette())


class GeneralTab(QWidget):
    simflg_changed = pyqtSignal(bool)
    grpflg_changed = pyqtSignal(bool)

    def __init__(self, n, c):
        super(GeneralTab, self).__init__()
        self.no_seqs = n
        self.consensnum = c
        if c > self.no_seqs:
            self.consensnum = 0
        self.proc_flag = False # flag to indicate that need to remake consensus and colours

        self.settings = QSettings("Boxshade", "Boxshade")
        self.scflag = self.settings.value("scflag", type=bool)
        self.consflag = self.settings.value("consflag", type=bool)
        self.symbcons = self.settings.value("symbcons")
        self.snameflag = self.settings.value("snameflag", type=bool)
        self.LHsnumsflag = self.settings.value("LHsnumsflag", type=bool)
        self.RHsnumsflag = self.settings.value("RHsnumsflag", type=bool)
        self.defnumsflag=self.settings.value("defnumsflag", type=bool)
        self.simflag = self.settings.value("simflag", type=bool)
        self.globalflag = self.settings.value("globalflag", type=bool)
        self.thrfrac = self.settings.value("thrfrac", type=float)
        self.outlen = self.settings.value("outlen", type=int)
        self.interlines = self.settings.value("interlines", type=int)
        self.rulerflag = self.settings.value("rulerflag",  type=bool)
        self.pepseqsflag = self.settings.value("pepseqsflag", True, type=bool)
        self.startnums = [] # placeholder

        mainlayout = QGridLayout()
        mainlayout.setColumnMinimumWidth(1,20)
# consensus to single sequence and combobox to get that sequence
        self.scbox = QCheckBox("Consensus to a single sequence")
        self.consnum = QComboBox()
        self.consnum.setMinimumContentsLength(2)
        if self.no_seqs > 1:
            self.consnum.addItems(['{:>2}'.format(i) for i in list(map(str, range(1, self.no_seqs + 1)))])
        else:
            self.consnum.addItems([' 1'])
        consnumlayout = QHBoxLayout()
        consnumlabel = QLabel("Make consensus to sequence:")
        consnumlayout.addWidget(consnumlabel)
        consnumlayout.addWidget(self.consnum)
        if self.consensnum > 0:
            self.consnum.setCurrentIndex(self.consensnum - 1)
        if self.scflag:
            self.scbox.setCheckState(Qt.Checked)
            consnumlabel.setEnabled(True)
            self.consnum.setEnabled(True)
        else:
            consnumlabel.setEnabled(False)
            self.consnum.setEnabled(False)
        self.scbox.stateChanged.connect(consnumlabel.setEnabled)
        self.scbox.stateChanged.connect(self.consnum.setEnabled)

# simple checkbox for whether sequence names will be printed
        self.seqnamebox = QCheckBox("Print sequence names")
        if self.snameflag:
            self.seqnamebox.setCheckState(Qt.Checked)
# Double checkbox for whether sequence numbers will be printed at one end or both ends
# or not at all
        seqnumslabel = QLabel("Print sequence numbers:")
        self.LHseqnumbox = QCheckBox("at left end")
        self.RHseqnumbox = QCheckBox("at right end")
        if self.LHsnumsflag:
            self.LHseqnumbox.setCheckState(Qt.Checked)
        if self.RHsnumsflag:
            self.RHseqnumbox.setCheckState(Qt.Checked)
        seqnumslayout = QGridLayout()
        seqnumslayout.addWidget(seqnumslabel,0,0,1,2)
        seqnumslayout.addWidget(self.LHseqnumbox,1,0)
        seqnumslayout.addWidget(self.RHseqnumbox, 1,1)

        self.defnumsbox = QCheckBox("Default numbering (start at 1)")
        if self.defnumsflag:
            self.defnumsbox.setCheckState(Qt.Checked)
        self.LHseqnumbox.stateChanged.connect(lambda ch: self.defnumsbox.setEnabled(ch) if not self.RHseqnumbox.isChecked() else None)
        self.RHseqnumbox.stateChanged.connect(lambda ch: self.defnumsbox.setEnabled(ch) if not self.LHseqnumbox.isChecked() else None)
        self.LHseqnumbox.stateChanged.connect(lambda ch: self.startstable.setEnabled(ch) if (not self.defnumsbox.isChecked())\
                                                                        and (not self.RHseqnumbox.isChecked()) else None)
        self.RHseqnumbox.stateChanged.connect(lambda ch: self.startstable.setEnabled(ch) if (not self.defnumsbox.isChecked())\
                                                                        and (not self.LHseqnumbox.isChecked()) else None)

# checkbox for default sequence numbering and the table of startnumbers if required
        if self.RHsnumsflag or self.LHsnumsflag:
            self.defnumsbox.setEnabled(True)
        else:
            self.defnumsbox.setEnabled(False)
        if self.no_seqs>1:
            self.startstable = QTableWidget(self.no_seqs,2)
        else:
            self.startstable = QTableWidget(1, 2)
        self.startstable.setHorizontalHeaderLabels(["Seq. No.", "Start number"])
        self.startstable.horizontalHeader().setSectionResizeMode(0,QHeaderView.ResizeToContents)
        self.startstable.setMaximumWidth(190)
        for i in range(self.no_seqs+1):
            self.startstable.setRowHeight(i,10)
            b = QTableWidgetItem(str(i+1))
            b.setFlags(b.flags() ^ ~Qt.ItemIsSelectable )
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.startstable.setItem(i,0,b)

        if (self.LHsnumsflag or self.RHsnumsflag) and not self.defnumsflag:
            self.startstable.setEnabled(True)
        else:
            self.startstable.setEnabled(False)
        self.defnumsbox.stateChanged.connect(self.startstable.setDisabled)

#Put left hand side of dialog in place
        mainlayout.addWidget(self.scbox, 0, 0)
        mainlayout.addLayout(consnumlayout,1,0)
        mainlayout.addWidget(self.seqnamebox,2,0)
        mainlayout.addLayout(seqnumslayout, 3, 0)
        mainlayout.addWidget(self.defnumsbox,4,0)
        mainlayout.setAlignment(self.defnumsbox, Qt.AlignRight)
        mainlayout.addWidget(self.startstable, 5, 0, 4, 1)

#Commence RHS of layout
        thrlayout = QHBoxLayout()
        conslayout = QHBoxLayout()
        lenlayout = QHBoxLayout()
        interlayout = QHBoxLayout()
        thrlabel = QLabel("Threshold Fraction (0-1):")
        self.thrbox = QDoubleSpinBox(value=self.thrfrac, minimum=0.0,
                                maximum = 1.0, singleStep=0.01, decimals=2 )
        self.thrbox.setMaximumSize(QSize(60,20))
        thrlayout.addWidget(thrlabel)
        thrlayout.addWidget(self.thrbox)
        thrlayout.setContentsMargins(0,0,55,0)
        mainlayout.addLayout(thrlayout,0,2)

        self.simflagbox = QCheckBox("Special shading for similar residues")
        if self.simflag:
            self.simflagbox.setCheckState(Qt.Checked)
        self.simflagbox.stateChanged.connect(self.simbox_clicked)

        self.grpflagbox = QCheckBox("Special shading for completely\nconserved residues")
        if self.globalflag:
            self.grpflagbox.setCheckState(Qt.Checked)
        self.grpflagbox.stateChanged.connect(self.grpbox_clicked)
        self.conslinebox = QCheckBox("Print consensus line")
        self.rulerbox = QCheckBox("Print ruler line")
        mainlayout.addWidget(self.simflagbox,1,2)
        mainlayout.addWidget(self.grpflagbox,2,2)
        mainlayout.addWidget(self.rulerbox,3,2)
        mainlayout.addWidget(self.conslinebox,4,2)

        conslabel = QLabel("Symbols for consensus:")
        self.consbox = QLineEdit(self.symbcons, maxLength=3, alignment=Qt.AlignHCenter)
        self.consbox.setMaximumSize(QSize(40, 20))
        conslayout.addWidget(conslabel)
        conslayout.addWidget(self.consbox)
        conslayout.setContentsMargins(0, 0, 90, 0)
        if self.consflag:
            self.conslinebox.setCheckState(Qt.Checked)
            conslabel.setEnabled(True)
            self.consbox.setEnabled(True)
        else:
            conslabel.setEnabled(False)
            self.consbox.setEnabled(False)
        if self.rulerflag:
            self.rulerbox.setChecked(True)

        self.conslinebox.stateChanged.connect(conslabel.setEnabled)
        self.conslinebox.stateChanged.connect(self.consbox.setEnabled)
        mainlayout.addLayout(conslayout, 5, 2)

        lenlabel = QLabel("Sequence characters per line:")
        self.lenbox = makeQSP(value=self.outlen, step=1, maximum=250)
        self.lenbox.setMaximumSize(QSize(80, 20))
        lenlayout.addWidget(lenlabel)
        lenlayout.addWidget(self.lenbox)
        lenlayout.setContentsMargins(0, 0, 30, 0)
        mainlayout.addLayout(lenlayout, 6, 2)

        interlabel = QLabel("Number of lines between blocks:")
        self.interbox = makeQSP(value=self.interlines, minimum=1, step=1)
        self.interbox.setMaximumSize(QSize(40,20))
        interlayout.addWidget(interlabel)
        interlayout.addWidget(self.interbox)
        interlayout.setContentsMargins(0,0,30,0)
        mainlayout.addLayout(interlayout,7,2)
        self.ProtRB = QRadioButton("Protein")
        DNARB = QRadioButton("DNA/RNA")
        if self.pepseqsflag:
            self.ProtRB.setChecked(True)
        else:
            DNARB.setChecked(True)
        typecont = QHBoxLayout()
        typecont.setAlignment(Qt.AlignRight)
        typecont.addWidget(self.ProtRB)
        typecont.addWidget(DNARB)
        mainlayout.addLayout(typecont, 8, 2)

        self.setLayout(mainlayout)

    @pyqtSlot()
    def simbox_clicked(self):
        self.simflg_changed.emit(self.simflagbox.isChecked())

    @pyqtSlot()
    def grpbox_clicked(self):
        self.grpflg_changed.emit(self.grpflagbox.isChecked())

    def filltable(self):
        for i in range(self.no_seqs):
            b = QTableWidgetItem(str(self.startnums[i]))
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.startstable.setItem(i,1,b)

    def exit(self):
        if (self.LHseqnumbox.isChecked() or self.RHseqnumbox.isChecked()) and not self.defnumsbox.isChecked():  # wants numbers, but not the default
            if np.amax(self.startnums) == 1:
                mb = QMessageBox(self)
                mb.setTextFormat(Qt.RichText)
                mb.setText(
                    "<span style='text-align: center'><p style='font-size: 18pt'>Possible output options error</p></span>")
                mb.setInformativeText(
                    "<p style='font-size: 14pt'> You have selected non-default sequence numbering, but have not set the actual start numbers to be used for each sequence."\
                    "<br><br>Continue with saving preferences?</p>")
                mb.setIcon(QMessageBox.Information)
                mb.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
                mb.setDefaultButton(QMessageBox.No)
                ret_code = mb.exec()
                if ret_code == QMessageBox.No:
                    return False
        if (self.scflag != self.scbox.isChecked() or self.pepseqsflag != self.ProtRB.isChecked()
            or self.consflag != self.conslinebox.isChecked() or
            (self.conslinebox.isChecked() and self.symbcons != self.consbox.text()) or
                (self.scbox.isChecked() and self.consensnum != 1+self.consnum.currentIndex()) or
                (self.thrfrac != self.thrbox.value())):
            self.proc_flag =True
        self.settings.setValue("scflag", self.scbox.isChecked())
        self.settings.setValue("consflag", self.conslinebox.isChecked())
        self.settings.setValue("symbcons", self.consbox.text())
        self.settings.setValue("snameflag", self.seqnamebox.isChecked())
        self.settings.setValue("LHsnumsflag", self.LHseqnumbox.isChecked())
        self.settings.setValue("RHsnumsflag", self.RHseqnumbox.isChecked())
        self.settings.setValue("defnumsflag", self.defnumsbox.isChecked())
        self.settings.setValue("simflag", self.simflagbox.isChecked())
        self.settings.setValue("globalflag", self.grpflagbox.isChecked())
        self.settings.setValue("thrfrac", self.thrbox.value())
        self.settings.setValue("outlen", self.lenbox.value())
        self.settings.setValue("interlines", self.interbox.value())
        self.settings.setValue("rulerflag", self.rulerbox.isChecked())
        self.settings.setValue("pepseqsflag", self.ProtRB.isChecked() )
        if self.startstable.isEnabled(): #only collect these values if we exit with the table enabled
            for i in range(self.no_seqs):
                b = self.startstable.item(i,1).text()
                try:
                    a=int(b)
                except:
                    a=1
                self.startnums[i]=a

        if self.scflag:
            self.consensnum = 1+self.consnum.currentIndex()
        else:
            self.consensnum = 1
# on exit from this dialog, there are some values that are not global preferences, i.e. not stored
# in the QSettings file; those parameters are specifically related to the current alignment that called the dialog.
# The relevant parameters are: (a) consensnum, the number of the seqence in the alignment which will be taken as consensus
# if scflag is True and (b) the startnums, the number of the first residue n each sequence if sequence numbers are
# to be printed and yet defnums=False (if defnums=True then all sequences start at 1)
# I aim to get these values from the Dialog when it returns by connecting a function from the calling window to the accepted
# signal from the dialog that retrieves those values as dialog_instance.GenTab.consensnum and dialog_instance.GenTab.startnums
        return True

class PSSettingsTab(QWidget):
    ss = "background-color : {0} ; color : {1} ;"

    def __init__(self, parent=None):
        super(PSSettingsTab, self).__init__(parent)

        PSGroup1 = QGroupBox("Different from consensus")
        PSGroup2 = QGroupBox("Identical to consensus")
        self.PSGroup3 = QGroupBox("Similar to consensus")
        self.PSGroup4 = QGroupBox("All the same residue")
        PSGroup5 = QGroupBox("Other prefs")

        self.settings = QSettings("Boxshade", "Boxshade")
        self.bgds = self.settings.value("PSbgds")
        self.fgds = self.settings.value("PSfgds")
        self.FSize = self.settings.value("PSFsize", type=int)
        self.LC = self.settings.value("PSLCs", type=bool)
        self.PSlandscapeflag = self.settings.value("PSlandscapeflag", type=bool)
        if system() == "Darwin":
            BS.monofont.setPointSize(32)
        else:
            BS.monofont.setPointSize(24)
        BS.monofont.setWeight(QFont.Bold)

        PSLayout1 = QGridLayout()
        DF_button = QPushButton("Foreground")
        DB_button = QPushButton("Background")
        Diff = QLabel("ACGT")
        Diff.setFrameStyle(QFrame.NoFrame)
        Diff.setAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        Diff.setAutoFillBackground(True)
        Diff.setFont(BS.monofont)
        Diff.setFixedSize(80, 32)
        Diff_UC = QRadioButton("Uppercase")
        self.Diff_LC = QRadioButton("Lowercase")
        if self.LC[0]:
            Diff.setText("acgt")
            self.Diff_LC.setChecked(True)
        else:
            Diff_UC.setChecked(True)
        Diff_UC.toggled.connect(lambda ch: Diff.setText("ACGT") if ch else Diff.setText("acgt"))
        Diffcont = QVBoxLayout()
        Diffcont.addWidget(Diff_UC)
        Diffcont.addWidget(self.Diff_LC)
        PSLayout1.addLayout(Diffcont, 1, 0)
        DF_button.clicked.connect(lambda ch: self.get_fg(Diff, 0))
        DB_button.clicked.connect(lambda ch: self.get_bg(Diff, 0))

        PSLayout1.addWidget(DF_button, 0, 0)
        PSLayout1.addWidget(DB_button, 0, 1)
        PSLayout1.addWidget(Diff, 1, 1, Qt.AlignCenter)
        PSGroup1.setLayout(PSLayout1)
        Diff.setStyleSheet(self.ss.format(self.bgds[0].name(), self.fgds[0].name()))

        PSLayout2 = QGridLayout()
        IDF_button = QPushButton("Foreground")
        IDB_button = QPushButton("Background")
        ID = QLabel("ACGT")
        ID.setFrameStyle(QFrame.NoFrame)
        ID.setAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        ID.setAutoFillBackground(True)
        ID.setFont(BS.monofont)
        ID.setFixedSize(80, 32)
        ID_UC = QRadioButton("Uppercase")
        self.ID_LC = QRadioButton("Lowercase")
        if self.LC[1]:
            ID.setText("acgt")
            self.Diff_LC.setChecked(True)
        else:
            ID_UC.setChecked(True)
        ID_UC.toggled.connect(lambda ch: ID.setText("ACGT") if ch else ID.setText("acgt"))
        IDcont = QVBoxLayout()
        IDcont.addWidget(ID_UC)
        IDcont.addWidget(self.ID_LC)
        PSLayout2.addLayout(IDcont, 1, 0)
        IDF_button.clicked.connect(lambda ch: self.get_fg(ID, 1))
        IDB_button.clicked.connect(lambda ch: self.get_bg(ID, 1))

        PSLayout2.addWidget(IDF_button, 0, 0)
        PSLayout2.addWidget(IDB_button, 0, 1)
        PSLayout2.addWidget(ID, 1, 1, Qt.AlignCenter)
        PSGroup2.setLayout(PSLayout2)
        ID.setStyleSheet(self.ss.format(self.bgds[1].name(), self.fgds[1].name()))

        PSLayout3 = QGridLayout()
        SimF_button = QPushButton("Foreground")
        SimB_button = QPushButton("Background")
        Sim = QLabel("ACGT")
        Sim.setFrameStyle(QFrame.NoFrame)
        Sim.setAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        Sim.setAutoFillBackground(True)
        Sim.setFont(BS.monofont)
        Sim.setFixedSize(80, 32)
        Sim_UC = QRadioButton("Uppercase")
        self.Sim_LC = QRadioButton("Lowercase")
        if self.LC[2]:
            Sim.setText("acgt")
            self.Sim_LC.setChecked(True)
        else:
            Sim_UC.setChecked(True)
        Sim_UC.toggled.connect(lambda ch: Sim.setText("ACGT") if ch else Sim.setText("acgt"))
        Simcont = QVBoxLayout()
        Simcont.addWidget(Sim_UC)
        Simcont.addWidget(self.Sim_LC)
        PSLayout3.addLayout(Simcont, 1, 0)
        SimF_button.clicked.connect(lambda ch: self.get_fg(Sim, 2))
        SimB_button.clicked.connect(lambda ch: self.get_bg(Sim, 2))

        PSLayout3.addWidget(SimF_button, 0, 0)
        PSLayout3.addWidget(SimB_button, 0, 1)
        PSLayout3.addWidget(Sim, 1, 1, Qt.AlignCenter)
        self.PSGroup3.setLayout(PSLayout3)
        Sim.setStyleSheet(self.ss.format(self.bgds[2].name(), self.fgds[2].name()))

        PSLayout4 = QGridLayout()
        GlobF_button = QPushButton("Foreground")
        GlobB_button = QPushButton("Background")
        Glob = QLabel("ACGT")
        Glob.setFrameStyle(QFrame.NoFrame)
        Glob.setAlignment(Qt.AlignVCenter | Qt.AlignHCenter)
        Glob.setAutoFillBackground(True)
        Glob.setFont(BS.monofont)
        Glob.setFixedSize(80, 32)
        Glob_UC = QRadioButton("Uppercase")
        self.Glob_LC = QRadioButton("Lowercase")
        if self.LC[3]:
            Glob.setText("acgt")
            self.Glob_LC.setChecked(True)
        else:
            Glob_UC.setChecked(True)
        Glob_UC.toggled.connect(lambda ch: Glob.setText("ACGT") if ch else Glob.setText("acgt"))
        Globcont = QVBoxLayout()
        Globcont.addWidget(Glob_UC)
        Globcont.addWidget(self.Glob_LC)
        PSLayout4.addLayout(Globcont, 1, 0)
        GlobF_button.clicked.connect(lambda ch: self.get_fg(Glob, 3))
        GlobB_button.clicked.connect(lambda ch: self.get_bg(Glob, 3))

        PSLayout4.addWidget(GlobF_button, 0, 0)
        PSLayout4.addWidget(GlobB_button, 0, 1)
        PSLayout4.addWidget(Glob, 1, 1, Qt.AlignCenter)
        self.PSGroup4.setLayout(PSLayout4)
        Glob.setStyleSheet(self.ss.format(self.bgds[3].name(), self.fgds[3].name()))

        PSLayout5 = QGridLayout()
        Fsizelayout = QHBoxLayout()
        Fsizelabel = QLabel("Font size:")
        self.Fsizebox = makeQSP(value=self.FSize, step=1, maximum=48, minimum=6)
        self.Fsizebox.setMaximumSize(QSize(60, 20))
        Fsizelayout.addWidget(Fsizelabel)
        Fsizelayout.addWidget(self.Fsizebox)
        Fsizelayout.setContentsMargins(0, 0, 0, 50)
        PSLayout5.addLayout(Fsizelayout, 0, 0)
        PortraitRB = QRadioButton("Portrait")
        self.LandscapeRB = QRadioButton("Landscape")
        if self.PSlandscapeflag:
            self.LandscapeRB.setChecked(True)
        else:
            PortraitRB.setChecked(True)
        Rotatecont = QVBoxLayout()
        Rotatecont.addWidget(PortraitRB)
        Rotatecont.addWidget(self.LandscapeRB)
        PSLayout5.addLayout(Rotatecont, 1, 0)

        PSGroup5.setLayout(PSLayout5)

        mainLayout = QGridLayout()
        mainLayout.addWidget(PSGroup1, 0, 0)
        mainLayout.addWidget(PSGroup2, 0, 1)
        mainLayout.addWidget(self.PSGroup3, 1, 0)
        mainLayout.addWidget(self.PSGroup4, 1, 1)
        mainLayout.addWidget(PSGroup5, 0, 2)
        self.setLayout(mainLayout)

    @pyqtSlot(bool)
    def simflag_changed(self, message):
        self.PSGroup3.setEnabled(message)
        self.raise_()

    @pyqtSlot(bool)
    def grpflag_changed(self, message):
        self.PSGroup4.setEnabled(message)
        self.raise_()

    def get_fg(self, widget, i):
        color = QColorDialog.getColor()
        if color.isValid():
            self.fgds[i] = color
            widget.setStyleSheet(self.ss.format(self.bgds[i].name(), self.fgds[i].name()))
        return

    def get_bg(self, widget, i):
        color = QColorDialog.getColor()
        if color.isValid():
            self.bgds[i] = color
            widget.setStyleSheet(self.ss.format(self.bgds[i].name(), self.fgds[i].name()))
        return

    def exit(self):
        self.settings.setValue("PSbgds", self.bgds)
        self.settings.setValue("PSfgds", self.fgds)
        self.settings.setValue("PSFsize", int(self.Fsizebox.value()))
        self.LC[0] = self.Diff_LC.isChecked()
        self.LC[1] = self.ID_LC.isChecked()
        self.LC[2] = self.Sim_LC.isChecked()
        self.LC[3] = self.Glob_LC.isChecked()
        self.settings.setValue("PSLCs", self.LC)
        self.settings.setValue("PSlandscapeflag", self.LandscapeRB.isChecked())
        self.settings.sync()


class ASCIITab(QWidget):
    def __init__(self, parent=None):
        super(ASCIITab, self).__init__(parent)

        topLabel = QLabel("<big><b>Text (ASCII) output is only for comparison to a specific sequence.</b></big><br><br>"
                           "In the boxes below, enter the single character to be used to represent residues that are<br>"
                           "Different from, Identical to, or Similar to the consensus, or to be used when all sequences<br>"
                           "have the same residue (Conserved). Use upper or lowercase 'L' (L or l) to indicate the<br>"
                           "character itself or the lower case version of that character<br>")

        layout = QGridLayout()

        self.settings = QSettings("Boxshade", "Boxshade")
        self.Achars = self.settings.value("ASCIIchars")
        self.inp=['']*4
        layout.addWidget(topLabel, 0, 0, 1, 4)
        rowlabel1 = QLabel("Different: ")
        rowlabel2 = QLabel("Identical: ")
        rowlabel3 = QLabel("Similar:  ")
        rowlabel4 = QLabel("Conserved:")
        rowlabel1.setAlignment(Qt.AlignHCenter)
        rowlabel2.setAlignment(Qt.AlignHCenter)
        rowlabel3.setAlignment(Qt.AlignHCenter)
        rowlabel4.setAlignment(Qt.AlignHCenter)
        layout.addWidget(rowlabel1, 1, 0)
        layout.addWidget(rowlabel2, 1, 1)
        layout.addWidget(rowlabel3, 1, 2)
        layout.addWidget(rowlabel4, 1, 3)
        self.inp[0] = QLineEdit(self.Achars[0], maxLength=1, alignment=Qt.AlignHCenter)
        self.inp[0].setMaximumSize(QSize(25, 25))
        layout.addWidget(self.inp[0], 2, 0, alignment=Qt.AlignHCenter|Qt.AlignTop)
        self.inp[1] = QLineEdit(self.Achars[1], maxLength=1, alignment=Qt.AlignHCenter)
        self.inp[1].setMaximumSize(QSize(25, 25))
        layout.addWidget(self.inp[1], 2, 1, alignment=Qt.AlignHCenter|Qt.AlignTop)
        self.inp[2] = QLineEdit(self.Achars[2], maxLength=1, alignment=Qt.AlignHCenter)
        self.inp[2].setMaximumSize(QSize(25, 25))
        layout.addWidget(self.inp[2], 2, 2, alignment=Qt.AlignHCenter|Qt.AlignTop)
        self.inp[3] = QLineEdit(self.Achars[3], maxLength=1, alignment=Qt.AlignHCenter)
        self.inp[3].setMaximumSize(QSize(25, 25))
        layout.addWidget(self.inp[3], 2, 3, alignment=Qt.AlignHCenter|Qt.AlignTop)
        pad=QLabel("")
        pad.setMinimumHeight(150)
        layout.addWidget(pad,3,0,1,4)

        self.setLayout(layout)

    def exit(self):
        tChars=[self.inp[i].text() for i in range(4)]
        if tChars != self.Achars:
            self.settings.setValue("ASCIIchars", tChars)

aaset = 'ACDEFGHIKLMNPQRSTVWY'
naset = 'ACGTKMRSWYBDHV'

class simsTab(QWidget):
    
    def __init__(self, parent=None):
        super(simsTab, self).__init__(parent)
        self.proc_flag=False
        self.settings = QSettings("Boxshade", "Boxshade")
        self.simsline = self.settings.value("simsline")
        self.simsline = self.simsline.split(":")
        self.simsline = self.simsline[1:len(self.simsline)-1]
        self.DNAline = self.settings.value("DNAsimsline")
        self.DNAline = self.DNAline.split(":")
        self.DNAline = self.DNAline[1:len(self.DNAline) - 1]

        leftLabel = QLabel('<b>Residues to be considered similar for shading purposes</b><br><br>'
                          'The "Sims" table allows one to set, for each amino acid (or nucleic acid), which residues are '
                          'considered to be similar to it if it is the consensus residue. Each amino acid/base can have '
                          'multiple amino acids/bases that are to be considered similar, or none, and amino acids/bases can '
                          'appear multiple times on the "similar to" side. The relationship is one way, '
                          'e.g. saying that H is similar to a consensus K does not make K similar to a consensus H.')
        leftLabel.setWordWrap(True)
        layout = QGridLayout()
        leftLabel.setAlignment(Qt.AlignLeft | Qt.AlignTop)
#        pad = QLabel("")
#        pad.setMinimumWidth(10)
#        layout.addWidget(pad, 0, 1)
        layout.addWidget(leftLabel, 0, 0)
        self.simstable = QTableWidget(len(aaset), 2)
        self.simstable.setHorizontalHeaderLabels(["Amino Acid", "Similar to:"])
#        self.simstable.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.simstable.setColumnWidth(0, 65)
        self.simstable.setColumnWidth(1, 65)
        self.simstable.setFixedWidth(165)
        self.simstable.setFixedHeight(308)
        for i in range(len(aaset)):
            self.simstable.setRowHeight(i, 10)
            b = QTableWidgetItem(aaset[i])
            b.setFlags(b.flags() ^ ~Qt.ItemIsSelectable)
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.simstable.setItem(i, 0, b)
            ii = [j for j in self.simsline if j.startswith(aaset[i])]
            if len(ii)>0:
                b = QTableWidgetItem(ii[0][2:])
            else:
                b = QTableWidgetItem('')
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.simstable.setItem(i, 1, b)
        layout.addWidget(self.simstable,0,2, Qt.AlignTop)
        self.DNAtable = QTableWidget(len(naset), 2)
        self.DNAtable.setHorizontalHeaderLabels(["Nucleic Acid", "Similar to:"])
#        self.DNAtable.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.DNAtable.setColumnWidth(0,75)
        self.DNAtable.setColumnWidth(1,65)
        self.DNAtable.setFixedWidth(160)
        self.DNAtable.setFixedHeight(257)
        for i in range(len(naset)):
            self.DNAtable.setRowHeight(i, 10)
            b = QTableWidgetItem(naset[i])
            b.setFlags(b.flags() ^ ~Qt.ItemIsSelectable)
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.DNAtable.setItem(i, 0, b)
            ii = [j for j in self.DNAline if j.startswith(naset[i])]
            if len(ii)>0:
                b = QTableWidgetItem(ii[0][2:])
            else:
                b = QTableWidgetItem('')
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.DNAtable.setItem(i, 1, b)
        layout.addWidget(self.DNAtable,0,3, Qt.AlignTop|Qt.AlignRight)
        self.setLayout(layout)

    def exit(self):
        simsline = []
        for i in range(len(aaset)):
            b = self.simstable.item(i, 1).text().strip()
            if len(b) > 0 and b.isalpha():
                simsline.append(aaset[i] + ' ' + b)
        simsline = 'SIMS:'+ ':'.join(simsline)+':END'
        if simsline != self.settings.value("simsline"):
            self.proc_flag=True
            self.settings.setValue("simsline", simsline)
        DNAline = []
        for i in range(len(naset)):
            b = self.DNAtable.item(i, 1).text().strip()
            if len(b) > 0 and b.isalpha():
                DNAline.append(naset[i] + ' ' + b)
        DNAline = 'SIMS:' + ':'.join(DNAline) + ':END'
        if DNAline != self.settings.value("DNAsimsline"):
            self.proc_flag=True
            self.settings.setValue("DNAsimsline", DNAline)


class grpsTab(QWidget):
    def __init__(self, parent=None):
        super(grpsTab, self).__init__(parent)
        self.proc_flag=False
        self.settings = QSettings("Boxshade", "Boxshade")
        self.grpsline = self.settings.value("grpsline")
        self.grpsline = self.grpsline.split(":")
        self.grpsline = self.grpsline[1:len(self.grpsline) - 1]
        self.DNAgline = self.settings.value("DNAgrpsline")
        self.DNAgline = self.DNAgline.split(":")
        self.DNAgline = self.DNAgline[1:len(self.DNAgline) - 1]
        self.no_lines = 13

        leftLabel = QLabel('<b>Residues to be considered as groups for shading purposes</b><br><br>'
                          'The "Groups" table defines groups of amino acids (or nucleic acids) which are to be considered similar when '
                          'establishing a group consensus (e.g. all positive charged, all large hydrophobic, all purines). These '
                          'relationships are multiway, e.g. making a group "LMFY" establishes relationships between '
                            'all four amino acids.')
        leftLabel.setWordWrap(True)
        layout = QGridLayout()
        leftLabel.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        pad = QLabel("")
        pad.setMinimumWidth(40)
        layout.addWidget(pad, 0,1)
        layout.addWidget(leftLabel, 0, 0)
        self.grpstable = QTableWidget(self.no_lines, 1)
        self.grpstable.setHorizontalHeaderLabels(["Amino Acid Groups"])
        self.grpstable.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.grpstable.setFixedWidth(131)
        self.grpstable.setFixedHeight(279)
        for i in range(len(self.grpsline)):
            self.grpstable.setRowHeight(i, 20)
            b = QTableWidgetItem(self.grpsline[i])
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.grpstable.setItem(i, 0, b)
        for i in range(len(self.grpsline),self.no_lines):
            self.grpstable.setRowHeight(i, 20)
            b = QTableWidgetItem('')
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.grpstable.setItem(i, 0, b)
        layout.addWidget(self.grpstable,0,2)
        self.DNAgtable = QTableWidget(self.no_lines, 1)
        self.DNAgtable.setHorizontalHeaderLabels(["Nucleic Acid Groups"])
        self.DNAgtable.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.DNAgtable.setFixedWidth(135)
        self.DNAgtable.setFixedHeight(279)
        for i in range(len(self.DNAgline)):
            self.DNAgtable.setRowHeight(i, 20)
            b = QTableWidgetItem(self.DNAgline[i])
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.DNAgtable.setItem(i, 0, b)
        for i in range(len(self.DNAgline), self.no_lines):
            self.DNAgtable.setRowHeight(i, 20)
            b = QTableWidgetItem('')
            b.setTextAlignment(int(Qt.AlignHCenter | Qt.AlignVCenter))
            self.DNAgtable.setItem(i, 0, b)
        layout.addWidget(self.DNAgtable, 0, 3)
        self.setLayout(layout)

    def exit(self):
        grpsline = []
        for i in range(self.no_lines):
            b = self.grpstable.item(i, 0).text().strip()
            if len(b) > 1 and b.isalpha():
                grpsline.append(b)
        grpsline = 'GRPS:'+ ':'.join(grpsline)+':END'
        if grpsline != self.settings.value('grpsline'):
            self.proc_flag=True
            self.settings.setValue("grpsline", grpsline)
        grpsline = []
        for i in range(self.no_lines):
            b = self.DNAgtable.item(i, 0).text().strip()
            if len(b) > 1 and b.isalpha():
                grpsline.append(b)
        grpsline = 'GRPS:' + ':'.join(grpsline) + ':END'
        if grpsline != self.settings.value('DNAgrpsline'):
            self.proc_flag = True
            self.settings.setValue("DNAgrpsline", grpsline)

if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    no_seqs=8
    consensnum=3
    tabdialog = prefsDialog(no_seqs,consensnum)
    tabdialog.GenTab.startnums = [1]*no_seqs
    tabdialog.GenTab.filltable()

    if 1 == tabdialog.exec():
        consensnum = tabdialog.GenTab.consensnum
        startnums = tabdialog.GenTab.startnums
#    sys.exit(app.exec_())
