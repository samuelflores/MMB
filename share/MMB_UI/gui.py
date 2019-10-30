import os
import chimera
import Tkinter as Tk
import Tix
import ttk
import tkMessageBox

from utils import *

from Tkinter import W,E,S,N

from chimera.baseDialog import ModelessDialog
from chimera import tkoptions

from VolumeViewer import Volume_Menu

from MMB_UI import *

####################################################################


## Try to execute a pyMMB command and display an error dialog in case of exception
#  @param func the function to execute
def tryMMB(func, *args):
    if not hasattr(command, '__call__'):
        return
    try:
        func(*args)
    except pyMMB.MMBError as e:
        tkMessageBox.showerror("MMB Error", e.msg)

class MMBDialog(ModelessDialog):
    """
    MMB main dialog description
    """

    name = 'MMB'
    buttons = ("Close","Save","ApplyViz","ResetMMB")

    title = "MMB Control GUI"

    ## GUI code
    def fillInUI(self, parent):
        """
        Implementation of the MMB_UI
        """

        # We hide the window while all the widgets are not created
        # Prevents showing something that looks like a glitch
        parent.winfo_toplevel().withdraw()

        # style = ttk.Style()
        # style.theme_use("default")

        ## Chimera's selection's trigger handler
        self.selHandlerID = None

        ## Flag to avoid continuous refresh
        self.stopRefreshSelection = False

        ## Tooltip handler
        self.balloon = Tix.Balloon()

        ## Tabbed frames
        self.tabs = ttk.Notebook(parent)
        self.tabs.pack(side='top',fill='both', expand=True)

        self.createExamplesTab()
        self.createInputTab()
        self.createCmdsTab()
        self.createBaseInteractionsTab()
        self.createMobilizersTab()
        self.createConstraintsTab()
        self.createContactsTab()
        self.createThreadingTab()
        self.createDensitiesTab()
        self.createPhysicsTab()
        self.createSimulationTab()

        self.refreshAll()

        self.tabs.select(1)

        self.enter()

    ## Bind an handler to "selection changed" trigger only when the dialog is active
    def map(self, *ignore):
        self.selHandlerID = None
        #self.selHandlerID = chimera.triggers.addHandler("selection changed", self.synchronizeBondSelection, None)

    ## Unbind the "selection changed" trigger
    def unmap(self, *ingore):
        if self.selHandlerID:
            chimera.triggers.deleteHandler("selection changed", self.selHandlerID)
            self.selHandlerID = None

    ## Create an empty tab named "name"
    #  @param name Name of the tab
    def createTab(self, name):
        dummyF = Tk.Frame(self.tabs)
        dummyF.pack(fill='both')
        self.tabs.add(dummyF,text=name)

###################################################################################################
    ## Create the Examples Tab
    def createExamplesTab(self):
        # Frame
        self.examplesFrame = Tk.Frame(self.tabs)
        self.examplesFrame.pack(fill='both')

        rnaBut = Tix.Button(self.examplesFrame, text="1ARJ", command=runRNAExample)
        rnaBut.pack(side='left')
        trnaBut = Tix.Button(self.examplesFrame, text="tRNA", command=runtRNAExample)
        trnaBut.pack(side='left')
        peptideBut = Tix.Button(self.examplesFrame, text="AAAPAFYAA", command=runPeptideExample)
        peptideBut.pack(side='left')
        morphBut = Tix.Button(self.examplesFrame, text="Morph", command=runMorph)
        morphBut.pack(side='left')

        self.tabs.add(self.examplesFrame, text="Examples")

###################################################################################################
    ## Create the tab for PDB input and Polymers editing
    def createInputTab(self):
        #print "howdy"
        # Frame
        inputFrame = Tk.Frame(self.tabs)
        inputFrame.pack(fill='both')

        # Polymers frame 
        self.polymersFrame = ScrolledFrame(inputFrame, tkTags=("PolymersFrameTag",))

        self.drawPolymersFrame(self.polymersFrame)

        self.tabs.add(inputFrame, text='Input', compound='left')

        # Button bar
        buttons = Tk.Frame(inputFrame, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh/Undo', command=self.refreshSequences)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh sequences or undo unvalidated modifications")

        self.deleteButton = Tix.Button(buttons, text='Delete', command=self.deleteSequences)
        self.deleteButton.pack(side='left')
        self.balloon.bind_widget(self.deleteButton, msg="Delete selected sequences")

        self.loadButton = Tix.Button(buttons, text="START", command=self.loadPolymers)
        self.loadButton.pack(side='right')
        self.balloon.bind_widget(self.loadButton, msg="Generate structure in MMB. You can't modify the sequences after this operation.")

        separator = Tk.Label(buttons, text=" | ")
        separator.pack(side='right')

        self.openButton = Tix.Button(buttons, text="Open PDB", 
                                     command=lambda ipd=InputPdbDialog: ipd(command=self.extractFromPdb, 
                                                                            multiple=False).enter()
                                     )
        self.openButton.pack(side='right')
        self.balloon.bind_widget(self.openButton, msg="Extract sequences from a PDB file which will also be used to generate the structure (by default).")

        self.openMMBButton = Tix.Button(buttons, text="Open MMB file", 
                                     command=lambda ipd=InputMMBFileDialog: ipd(command=self.openMMBFile, 
                                                                            multiple=False).enter()
                                     )
        self.openMMBButton.pack(side='right')
        self.balloon.bind_widget(self.openMMBButton, msg="Load an MMB commands file.")


        self.importButton = Tix.Button(buttons, text="Import from Chimera", 
                                     command=self.importFromChimera)
        self.importButton.pack(side='right')
        self.selOnlyCheck = TracedCheckButton(buttons, text="Curent selection only")
        self.selOnlyCheck.pack(side='right') 
        self.balloon.bind_widget(self.importButton, msg="Extract sequences from a PDB file which will also be used to generate the structure (by default).")

    ## Draw the frame containing the polymers widgets
    def drawPolymersFrame(self, frame):
        # Clear the frame first
        for w in frame.children.values():
            w.destroy()

        # Headers
        polyHeaders = ("Type", "ChainID", "First Res. #", "Sequence", "Topology", "","Select")
        for i, header in enumerate(polyHeaders):
            headerLabel = Tk.Label(frame,text=header)
            headerLabel.grid(row=0,column=i)
            # print header, headerLabel.grid_info()
            # frame.grid_columnconfigure(i, minsize=len(header))
            # print frame.grid_bbox(column=i, row=0)

        self.polymersWidgetsLines = []
        # Polymers data
        for i,poly in enumerate(sortedPolymers()):
            row = i+1
            widgetsLine = PolymerWidgetsLine(frame, row, poly)
            # widgetsLine.grid(row=row,column=0, columnspan=len(polyHeaders), sticky="we")
            self.polymersWidgetsLines.append(widgetsLine)
            widgetsLine.display()

        if MMB_UI.polymersInitialized == False:
            m = BioPolymer(new=True)
            widgetsLine = PolymerWidgetsLine(frame, len(self.polymersWidgetsLines)+1, m, 
                                               addCommand= lambda x=m: self.addPolymer(m))
            # widgetsLine.grid(row=len(self.polymersWidgetsLines)+1,column=0, columnspan=len(polyHeaders))
            widgetsLine.display()
            self.newPolymerWidgetsLine = widgetsLine

    # Callbacks
    ## Remove the polymer with chain x and refresh the tabs
    #  @param x chain ID
    def removePolymer(self, x):
        """ param: x is the polymer's chainID """
        """ TODO: remove the interactions/forces associated with the polymer"""
        del MMB_UI.polymers[x]
        self.drawPolymersFrame(self.polymersFrame)
        self.refreshAll()

    ## Extract the polymers from the PDB file
    #  @param filepath. string path of the pdb file
    def extractFromPdb(self, filepath):
        """ 
        Use MMB pdb reader to extract polymers from the file.
        Clear current polymers
        """
        if(os.path.isfile(filepath)):
            try:
                if filepath != MMB_UI.currentPdb:
                    MMB_UI.currentPdb = filepath
                    MMB_UI.loadSequencesFromPdb(filepath) # add polymers' sequences
                    self.drawPolymersFrame(self.polymersFrame)
            except pyMMB.MMBError as e:
                tkMessageBox.showerror("MMB Error", e.msg)
                return False
            return True
        return False

    def openMMBFile(self, filepath):
        """ 
        Open an MMB command file.
        """
        if(os.path.isfile(filepath)):
            try:
                MMB_UI.setWorkDir(os.path.dirname(filepath))
                fileContent = open(filepath,'r').readlines()
                MMB_UI.sendMMBCmds(fileContent)
                self.refreshSequences()
            except pyMMB.MMBError as e:
                tkMessageBox.showerror("MMB Error", e.msg)
                return False
            return True
        return False      

    def importFromChimera(self):
        try:
            selectedAtoms = chimera.selection.currentAtoms()
            if not self.selOnlyCheck.get():
                selectedAtoms = []
            MMB_UI.importSequencesFromChimera(selectedAtoms)
            self.drawPolymersFrame(self.polymersFrame)
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        return True

    ## Clear all polymers and refresh the tabs
    def clearPolymers(self):
        """
        Clear current polymers and pdb file
        """
        MMB_UI.polymers = {}
        MMB_UI.currentPdb = ""
        self.inputPdbOption.set("")
        self.drawPolymersFrame(self.polymersFrame)
        self.refreshAll()

    ## Remove the selected base pair interactions from MMB
    def deleteSequences(self):
        try:
            # MMB IDs change after each deletion    
            # So we have to make sure we use the good ID
            # First: extract and sort selected MMBObjects via WidgetsLines
            lines = [l for l in self.polymersWidgetsLines if l.varSelect.get()]

            # Second: delete taking in to account the IDs decrement (-i)
            for line in lines:
                line.data.mmbDelete()
            self.refreshSequences()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        return True

    def refreshSequences(self):
        try:
            MMB_UI.refreshSequences()
            self.drawPolymersFrame(self.polymersFrame)
            if MMB_UI.polymersInitialized == True:
                self.loadButton.config(state="disabled")
                self.openButton.config(state="disabled")
                self.deleteButton.config(state="disabled")
                self.importButton.config(state="disabled")
                self.selOnlyCheck.config(state="disabled")
            else:
                self.loadButton.config(state="normal")
                self.openButton.config(state="normal")
                self.importButton.config(state="normal")
                self.selOnlyCheck.config(state="normal")
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        return True

    ## Add a new polymer to MMB
    def addPolymer(self, poly):
        try:
            poly.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshSequences()
        return True


    ## Load the polymers to MMB and Chimera
    def loadPolymers(self):
        if not MMB_UI.polymers:
            tkMessageBox.showwarning("MMB Error", "Nothing to initialize, please add BioPolymers.")
            return False
        try:
            if not MMB_UI.checkSequences():
                message = "One or several Sequences are not yet validated. Do you want to validate them and continue ?"
                if tkMessageBox.askyesno(title="MMB - Validate Sequences ?", message=message) == 0:
                    return False
            MMB_UI.validateSequences()

            # Let the user choose a working directory, where MMB results will be
            # WIP - this breaks the simulation for now...
            # import OpenSave
            # self._selectWDDialog = OpenSave.OpenModal(
            #             title="Select Working Directory",
            #             dirsOnly=1,
            #             multiple=False)
            # self._selectWDDialog.run(self.polymersFrame)
            # dirList = self._selectWDDialog.getPaths()
            # if not dirList:
            #     tkMessageBox.showerror("MMB Error", "Unable to open the selected directory, sorry.")
            #     return
            # MMB_UI.setWorkDir(dirList[0])
            # 
            # Why load before init?
            # initPolymers load the structures in the real MMB system but they are not OK until we init the simulation
            # To counter this, loadPolymers loads structures generated from a dummy MMB system initialized entirely that we discard after that. I know, we load twice, but it is the best solution right now (oct. 2013) compare to the (re)development cost of making a proper system.
            loadPolymers()
            initPolymers(MMB_UI.polymers, MMB_UI.currentPdb)

            self.tabs.select(1)
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshAll()
        return True


###################################################################################################
    ## Create the tab to manage MMB commands 
    def createCmdsTab(self):
        cmdsFrame = Tk.Frame(self.tabs)
        cmdsFrame.pack(fill='both')

        # File Input 
        inputFrame = Tk.Frame(cmdsFrame)
        from OpenSave import OpenModeless
        class InputMMBDialog(OpenModeless):
            default = 'Open MMB input file'
            title = 'MMB -- Input File'

        ofo = tkoptions.InputFileOption(inputFrame,
                                         0, 'Open an MMB input file ',
                                         MMB_UI.currentPdb,
                                         self.openCmdsFile,
                                         balloon='Open an MMB input file'
                                         )
        ofo.dialogType = InputMMBDialog
        self.inputMMBfile = ofo
        inputFrame.pack()

        self.cmdGlobalEditor = createCustomTextEditors(cmdsFrame, "Global MMB parameters")
        self.cmdGlobalEditor.pack(fill='both', expand=True)
        
        # self.cmdPromptVar = Tk.StringVar()
        # self.cmdPrompt = Tk.Entry(cmdsFrame, textvariable=self.cmdPromptVar)
        self.cmdPrompt = CmdPrompt(cmdsFrame, command=self.applyCmd, label="Command: ")
        self.cmdPrompt.pack(fill='x')
        self.balloon.bind_widget(self.cmdPrompt, msg="MMB prompt. Enter to send the command. UP and DOWN to navigate in the history.")

        # Button bar
        buttons = Tk.Frame(cmdsFrame, bd=2, relief="flat")
        buttons.pack(side='bottom', fill='both')

        cleanBut = Tix.Button(buttons, text="Clean", command=self.cleanCmds)
        cleanBut.pack(side='left')

        applyBut = Tix.Button(buttons, text="Apply", command=self.applyCmds)
        applyBut.pack(side='left')

        self.tabs.add(cmdsFrame, text="Commands")

    ## Reads and puts the content of a file into cmdGlobalEditor
    def openCmdsFile(self, x):
        if os.path.isfile(x.get()):
            self.cmdGlobalEditor.delete(1.0,'end')
            filename = x.get()
            MMB_UI.setWorkDir(os.path.dirname(filename))
            fileContent = open(filename,'r').read()
            self.cmdGlobalEditor.insert(1.0, fileContent)
    ## Send commands in the command editor to MMB.
    def applyCmds(self):
        cmds = self.cmdGlobalEditor.get(1.0,'end')
        try:
            error = sendMMBCmds(cmds.splitlines())
            if isinstance(error, basestring):
                msg = "Sorry, MMB GUI doesn't support the following command at the moment:\n"
                msg += error
                tkMessageBox.showerror("MMB Error", msg)
                self.cleanCmds()
            elif error:
                msg = "Ignored commands:\n"
                for c in error:
                    msg += c+"\n"
                tkMessageBox.showerror("MMB Warning", msg)
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            self.cleanCmds()
            return False
        self.refreshAll()
        return True

    ## Send a command to MMB.
    #  @param mmbCmd string of an MMB command
    def applyCmd(self, mmbCmd):
        try:
            error = sendMMBCmds([mmbCmd])
            if error != "":
                msg = "Sorry, MMB GUI doesn't support the following command at the moment:\n"
                tkMessageBox.showerror("MMB Error", msg+error)
                self.cleanCmds()
            else:
                self.cmdGlobalEditor.insert(Tk.END, "\n"+mmbCmd+"\n")
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshAll()
        return True

    ## Clear the Forces and Constraints in MMB
    def cleanCmds(self):
        try:
            clearForcesAndConstraints()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshAll()
        return True



###################################################################################################
    ## Create the tabs for the mobilizers
    def createMobilizersTab(self):
        self.mobiTab = Tk.Frame(self.tabs)
        self.mobiTab.pack(fill='both')

        mobilizersFrame = Tk.Frame(self.mobiTab)
        mobilizersFrame.pack(fill='both', expand=True)

        self.mobiFrame = ScrolledFrame(mobilizersFrame, tkTags=("MobilizersFrameTag",))
        self.drawMobilizersFrame(self.mobiFrame)

        mobiWithinFrame = Tk.Frame(self.mobiTab, bd=2)
        mobiWithinFrame.pack(fill='both', expand=True)
        self.mobiWithinFrame = ScrolledFrame(mobiWithinFrame, tkTags=("MobilizersWithinFrameTag",))
        self.drawMobilizersWithinFrame(self.mobiWithinFrame)

        rootMobiFrame = Tk.Frame(self.mobiTab, bd=2)
        rootMobiFrame.pack(fill='both', expand=True)
        self.rootMobiFrame = ScrolledFrame(rootMobiFrame, tkTags=("RootMobilizersFrameTag",))
        self.drawRootMobilizerFrame(self.rootMobiFrame)

        # Button bar
        buttons = Tk.Frame(self.mobiTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshMobilizers)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh mobilizers")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.deleteMobilizers)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected mobilizer")

        fromSelButton = Tix.Button(buttons, text='Mobilizer From Selection', command=self.setNewMobilizerFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new Mobilizer chain and residues according to selection")

        fromSelWithinButton = Tix.Button(buttons, text='MobilizerWithin From Selection', command=self.setNewMobilizerWithinFromSelection)
        fromSelWithinButton.pack(side='left')
        self.balloon.bind_widget(fromSelWithinButton, msg="Set the new MobilizerWhithin chain and residue according to selection")

        self.tabs.add(self.mobiTab, text="Mobilizers")

    ## Draw the frame containing the mobilizers widgets
    def drawMobilizersFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        mobiWidgetsLine = []

        # Headers
        for i, header in enumerate(mobilizersHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)

        # Widgets lines
        for i,p in enumerate(MMB_UI.mobilizers):
            widgetsLine = MobilizerWidgetsLine(frame, i+1, p)
            mobiWidgetsLine.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            m = Mobilizer()
            widgetsLine = MobilizerWidgetsLine(frame, len(MMB_UI.mobilizers)+1, m, 
                                               addCommand= lambda x=m: self.addMobilizer(m))
            widgetsLine.display()
            self.newMobilizerWidgetsLine = widgetsLine

        return mobiWidgetsLine

    ## Synchronize the mobilizers with MMB
    def refreshMobilizers(self):
        try:
            MMB_UI.refreshMobilizerStretches()
            MMB_UI.refreshMobilizerWithin()
            MMB_UI.refreshRootMobilizers()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.mobiWidgets = self.drawMobilizersFrame(self.mobiFrame)
        self.mobiWithinWidgets = self.drawMobilizersWithinFrame(self.mobiWithinFrame)
        self.rootMobiWidgets = self.drawRootMobilizerFrame(self.rootMobiFrame)
        return True

    ## Remove the selected base pair interactions from MMB
    def deleteMobilizers(self):
        self.deleteSelectedObjects(self.mobiWidgets + self.mobiWithinWidgets)
        self.refreshMobilizers()
        return True

    ## Add a new mobilizer stretch to MMB
    #  @param mobi a MMB_UI.Mobilizer
    def addMobilizer(self, mobi):
        try:
            mobi.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshMobilizers()
        return True

    def drawMobilizersWithinFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        mobiWithinWidgets = []
        # Headers
        for i, header in enumerate(mobilizersWithinHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)
            
        # Widgets lines
        for i,p in enumerate(MMB_UI.mobilizersWithin):
            wl = MobilizerWithinWidgetsLine(frame, i+1, p)
            mobiWithinWidgets.append(wl)
            wl.display()

        if len(MMB_UI.polymers) > 0:
            m = MobilizerWithin(model=MMB_UI.currentModel)
            wl = MobilizerWithinWidgetsLine(frame, len(MMB_UI.mobilizersWithin)+1, m, 
                                            addCommand=lambda x=m:self.addMobilizer(m))
            wl.display()
            self.newMobilizerWithinWidgetsLine = wl

        return mobiWithinWidgets

    def drawRootMobilizerFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        rootMobiWidgets = []
        #Headers
        rootMobiHeaders = ["Mobilizer", "Chain", "Root","Valid","Select"]
        for i, header in enumerate(rootMobiHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)

        # Widgets lines
        for i,p in enumerate(MMB_UI.rootMobilizers):
            wl = RootMobilizerWidgetsLine(frame, i+1, p)
            rootMobiWidgets.append(wl)
            wl.display()
        return rootMobiWidgets

    ## Sets the residues and chain for a new mobilizer from Chimera's current selection.
    #  Exactlt two residues must be selected.
    def setNewMobilizerFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) == 0:
            tkMessageBox.showerror("MMB Error", "Current selection contains no residue")
            return
        chainIds = [r.id.chainId for r in res]
        if len(set(chainIds)) > 1:
            tkMessageBox.showerror("MMB Error", "Residues must be part of only one chain")
            return
        #else:
        resIds = [r.id.position for r in res]
        resId1 = min(resIds)
        resId2 = max(resIds)
        poly = MMB_UI.polymers[chainIds[0]]
        widgetLine = self.newMobilizerWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widStartRes.set(poly.getResNameId(resId1))
        widgetLine.widEndRes.set(poly.getResNameId(resId2))

    ## Sets the residue and chain for a new mobilizerWithin from Chimera's current selection.
    #  Exactlt two residues must be selected.
    def setNewMobilizerWithinFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) != 1:
            tkMessageBox.showerror("MMB Error", "Exactly one residue must be selected")
            return
        res = res[0]
        chainId = res.id.chainId

        #else:
        resId = res.id.position
        poly = MMB_UI.polymers[chainId]
        widgetLine = self.newMobilizerWithinWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widRes.set(poly.getResNameId(resId))



###################################################################################################
    ## Create the tab for the base interactions
    def createBaseInteractionsTab(self):
        fullFrame = Tk.Frame(self.tabs)
        fullFrame.pack(fill='both')

        pairsFrame = Tk.Frame(fullFrame)
        pairsFrame.pack(fill='both', expand=True)

        self.bpiFrame = ScrolledFrame(pairsFrame, tkTags=("BaseInteractionFrameTag",))

        self.basePairWidgets = self.drawBaseInteractionsFrame(self.bpiFrame)


        # Button bar
        buttons = Tk.Frame(fullFrame, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Undo', command=self.refreshBasePairInteractions)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Undo non-validated modifications")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.delBasePairInteractions)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected base-pairs")

        fromSelButton = Tix.Button(buttons, text='From Selection', command=self.setNewPairResFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set new pair's residues from current Chimera's selection.\n"+ \
                                                     "Exactly 2 residues must be selected.")

        self.tabs.add(fullFrame, text="Base Interactions")

    ## Draw the frame containing the base interactions widgets
    def drawBaseInteractionsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        basePairWidgets = []

        # Headers names for base interactions
        pairsHeaders= ("Id","Chain1","Res1", "Edge1","Chain2","Res2","Edge2","Bond Orient.","Valid", "Select")
        for i, header in enumerate(pairsHeaders):
            pairLabel = Tk.Label(frame,text=header,bd=2, relief='raised')
            pairLabel.grid(row=0,column=i,sticky=W+E)

        for i,p in enumerate(MMB_UI.basePairs):
            widgetsLine = BasePairWidgetsLine(frame, i+1, p)
            widgetsLine.display()
            basePairWidgets.append(widgetsLine)

        if MMB_UI.polymersInitialized and MMB_UI.sortedNucleicAcidsIDs():
            p = BasePairInteraction.emptyPair()
            widgetsLine = BasePairWidgetsLine(frame, len(MMB_UI.basePairs)+1, p, 
                                              addCommand=lambda x=p:self.addBasePairInteraction(x))
            widgetsLine.display()
            self.newBasePairWidgetsLine = widgetsLine

        duplexFrame = Tk.Frame(frame, relief='groove',bd=2)
        duplexFrame.grid(row=len(MMB_UI.basePairs)+2, columnspan=10,sticky=W+E)
        duplexTitle = Tk.Label(duplexFrame, text="Nucleic Acid Duplex")
        duplexTitle.grid(row=0, column=0, columnspan=8)
        duplexHeaders = ("Chain1","Start1","End1", "Chain2", "Start2", "End2", " ", "Select")
        for i, header in enumerate(duplexHeaders):
            pairLabel = Tk.Label(duplexFrame,text=header,bd=2, relief='raised')
            pairLabel.grid(row=1,column=i,sticky=W+E)
        if MMB_UI.polymersInitialized and MMB_UI.sortedNucleicAcidsIDs():
            nad = NucleicAcidDuplex()
            widgetsLine = NucleicAcidDuplexWidgetsLine(duplexFrame, 2, nad, 
                                              addCommand=lambda x=nad:self.addBasePairInteraction(x))
            widgetsLine.display()


        return basePairWidgets

    ## Add the base pair interaction to MMB
    #  @param x a MMB_UI.BasePairInteraction
    def addBasePairInteraction(self, bp):
        try:
            bp.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshBasePairInteractions()
        return True

    ## Remove the selected base pair interactions from MMB
    def delBasePairInteractions(self):
        self.deleteSelectedObjects(self.basePairWidgets)
        self.refreshBasePairInteractions()
        return True

    ## Synchronize the base pair interactions with MMB
    def refreshBasePairInteractions(self):
        try:
            MMB_UI.refreshBasePairInteractions()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.basePairWidgets = self.drawBaseInteractionsFrame(self.bpiFrame)
        return True

    ## Sets the residues for a new base pair interaction from Chimera's current selection.
    #  Exactly two residues must be selected.
    def setNewPairResFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) != 2:
            tkMessageBox.showerror("MMB Error", "You must select exactly two residues.")
            return
        #else
        res1 = res[0]
        res2 = res[1]
        poly1 = MMB_UI.polymers[res1.id.chainId]
        poly2 = MMB_UI.polymers[res2.id.chainId]
        widgetLine = self.newBasePairWidgetsLine
        widgetLine.widPoly1.set(poly1.chainID)
        widgetLine.widPoly2.set(poly2.chainID)
        widgetLine.widRes1.set(poly1.getResNameId(res1.id.position))
        widgetLine.widRes2.set(poly2.getResNameId(res2.id.position))


###################################################################################################
    ## Create the tabs for the constraints
    def createConstraintsTab(self):
        self.constTab = Tk.Frame(self.tabs)
        self.constTab.pack(fill='both')

        constraintsFrame = Tk.Frame(self.constTab)
        constraintsFrame.pack(fill='both', expand=True)

        self.constFrame = ScrolledFrame(constraintsFrame, tkTags=("ConstraintsFrameTag",))
        self.drawConstraintsFrame(self.constFrame)

        springsFrame = Tk.Frame(self.constTab)
        springsFrame.pack(fill='both', expand=True)

        self.springFrame = ScrolledFrame(springsFrame, tkTags=("SpringsFrameTag",))
        self.drawSpringsFrame(self.springFrame)

        # Button bar
        buttons = Tk.Frame(self.constTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshConstraints)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh constraints")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.delConstraints)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected mobilizer")

        fromSelButton = Tix.Button(buttons, text='Constraint From Selection', 
                                   command=self.setNewConstraintFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new Constraint according to selection")

        springFromSelButton = Tix.Button(buttons, text='Spring From Selection', 
                                   command=self.setNewSpringFromSelection)
        springFromSelButton.pack(side='left')
        self.balloon.bind_widget(springFromSelButton, msg="Set the new Spring according to selection")

        self.tabs.add(self.constTab, text="Constraints")

    ## Draw the frame containing the base interactions widgets
    def drawConstraintsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        constraintsWidgets = []

        # Headers
        for i, header in enumerate(constraintsHeaders):
            constraintLabel = Tk.Label(frame,text=header,bd=2, relief='raised')
            constraintLabel.grid(row=0,column=i,sticky=W+E)

        for i,c in enumerate(MMB_UI.constraints):
            widgetsLine = ConstraintWidgetsLine(frame, i+1, c)
            constraintsWidgets.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            c = MMB_UI.Constraint(model=MMB_UI.currentModel)
            widgetsLine = ConstraintWidgetsLine(frame, len(MMB_UI.constraints)+1, c, 
                                              addCommand=lambda x=c:self.addConstraint(x))
            widgetsLine.display()
            self.newConstraintWidgetsLine = widgetsLine

        return constraintsWidgets

    ## Draw the frame containing the base interactions widgets
    def drawSpringsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        springsWidgets = []

        # Headers
        ## Headers names for constraints
        springsHeaders=("Id","Chain1","Res1","Atom1"," ", "Chain2", "Res2", "Atom2","Tether","Constant","Length", "Valid", "Select")
        for i, header in enumerate(springsHeaders):
            springLabel = Tk.Label(frame,text=header,bd=2, relief='raised')
            springLabel.grid(row=0,column=i,sticky=W+E)

        for i,c in enumerate(MMB_UI.atomSprings):
            widgetsLine = AtomSpringWidgetsLine(frame, i+1, c)
            springsWidgets.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            c = MMB_UI.AtomSpring(model=MMB_UI.currentModel)
            widgetsLine = AtomSpringWidgetsLine(frame, len(MMB_UI.atomSprings)+1, c, 
                                              addCommand=lambda x=c:self.addConstraint(x))
            widgetsLine.display()
            self.newSpringWidgetsLine = widgetsLine

        return springsWidgets

    ## Synchronize the constraints with MMB
    def refreshConstraints(self):
        try:
            MMB_UI.refreshConstraints()
            MMB_UI.refreshAtomSprings()
            pass
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.constraintsWidgets = self.drawConstraintsFrame(self.constFrame)
        self.springsWidgets = self.drawSpringsFrame(self.springFrame)
        return True

    ## Add a new Constraint stretch to MMB
    #  @param const a MMB_UI.Constraint
    def addConstraint(self, const):
        try:
            const.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshConstraints()
        return True

    ## Remove the selected constraints from MMB
    def delConstraints(self):
        self.deleteSelectedObjects(self.constraintsWidgets + self.springsWidgets)
        self.refreshConstraints()
        return True

    ## Sets the chains, residues and atoms for a new constraint from Chimera's current selection.
    def setNewConstraintFromSelection(self):
        wl = self.newConstraintWidgetsLine

        atoms = chimera.selection.currentAtoms()
        if len(atoms) < 1:
            tkMessageBox.showerror("MMB Error", "At least one atom must be selected")
            return
        atomNames = [a.name for a in atoms]

        res = []
        [res.extend(a.oslParents()) for a in atoms]
        if len(res) > 2:
            tkMessageBox.showerror("MMB Error", "No more than two residues must be selected")
            return
        resIds = [r.id.position for r in res]
        chainIds = [r.id.chainId for r in res]

        poly1 = MMB_UI.polymers[chainIds[0]]
        wl.widChain1.set(poly1.chainID)
        wl.widRes1.set(poly1.getResNameId(resIds[0]))

        if len(atoms) <= 2:
            wl.widAtom1.set(atomNames[0])

        if len(res) > 1:
            poly2 = MMB_UI.polymers[chainIds[1]]
            wl.widChain2.set(poly2.chainID)
            wl.widRes2.set(poly2.getResNameId(resIds[1]))
            if len(atoms) == 2:
                wl.widAtom2.set(atomNames[1])
        else:
            wl.widChain2.set("Ground")

    ## Sets the chains, residues and atoms for a new constraint from Chimera's current selection.
    def setNewSpringFromSelection(self):
        wl = self.newSpringWidgetsLine

        atoms = chimera.selection.currentAtoms()
        if len(atoms) != 2:
            tkMessageBox.showerror("MMB Error", "Exactly two atoms must be selected")
            return
        atomNames = [a.name for a in atoms]
        print atoms[0].oslParents()
        print atoms[1].oslParents()

        res = atoms[0].oslParents() + atoms[1].oslParents()
        print res
        resIds = [r.id.position for r in res]
        chainIds = [r.id.chainId for r in res]

        poly1 = MMB_UI.polymers[chainIds[0]]
        wl.widChain1.set(poly1.chainID)
        wl.widRes1.set(poly1.getResNameId(resIds[0]))
        wl.widAtom1.set(atomNames[0])

        poly2 = MMB_UI.polymers[chainIds[1]]
        wl.widChain2.set(poly2.chainID)
        wl.widRes2.set(poly2.getResNameId(resIds[1]))
        wl.widAtom2.set(atomNames[1])

###################################################################################################
    ## Create the tabs for the contacts
    def createContactsTab(self):
        self.contactTab = Tk.Frame(self.tabs)
        self.contactTab.pack(fill='both')

        contactsFrame = Tk.Frame(self.contactTab)
        contactsFrame.pack(fill='both', expand=True)

        self.contactFrame = ScrolledFrame(contactsFrame, tkTags=("ContactsFrameTag",))
        self.drawContactsFrame(self.contactFrame)

        contactWithinFrame = Tk.Frame(self.contactTab, bd=2)
        contactWithinFrame.pack(fill='both', expand=True)
        self.contactWithinFrame = ScrolledFrame(contactWithinFrame, tkTags=("ContactWithinFrameTag",))
        self.drawContactsWithinFrame(self.contactWithinFrame)

        # Button bar
        buttons = Tk.Frame(self.contactTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshContacts)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh contacts")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.deleteContacts)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected contact")

        fromSelButton = Tix.Button(buttons, text='Contact From Selection', command=self.setNewContactFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new Contact chain and residues according to selection")

        fromSelWithinButton = Tix.Button(buttons, text='ContactWithin From Selection', command=self.setNewContactWithinFromSelection)
        fromSelWithinButton.pack(side='left')
        self.balloon.bind_widget(fromSelWithinButton, msg="Set the new ContactWhithin chain and residue according to selection")

        self.tabs.add(self.contactTab, text="Contacts")

    ## Draw the frame containing the contacts widgets
    def drawContactsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        contactWidgetsLine = []

        contactsHeaders=("Id","Contact","Contact Scheme","Chain","Start Res","End Res","Valid","Select")
        # Headers
        for i, header in enumerate(contactsHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)

        # Widgets lines
        for i,p in enumerate(MMB_UI.contacts):
            widgetsLine = ContactWidgetsLine(frame, i+1, p)
            contactWidgetsLine.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            m = Contact()
            widgetsLine = ContactWidgetsLine(frame, len(MMB_UI.contacts)+1, m, 
                                               addCommand= lambda x=m: self.addContact(m))
            widgetsLine.display()
            self.newContactWidgetsLine = widgetsLine

        return contactWidgetsLine

    ## Synchronize the contacts with MMB
    def refreshContacts(self):
        try:
            MMB_UI.refreshContactStretches()
            MMB_UI.refreshContactWithin()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.contactWidgets = self.drawContactsFrame(self.contactFrame)
        self.contactWithinWidgets = self.drawContactsWithinFrame(self.contactWithinFrame)
        return True

    ## Remove the selected base pair interactions from MMB
    def deleteContacts(self):
        self.deleteSelectedObjects(self.contactWidgets + self.contactWithinWidgets)
        self.refreshContacts()
        return True

    ## Add a new contact stretch to MMB
    #  @param contact a MMB_UI.Contact
    def addContact(self, contact):
        try:
            contact.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshContacts()
        return True

    def drawContactsWithinFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        contactWithinWidgets = []
        # Headers
        contactsWithinHeaders=("Id","Contact","ContactScheme","Chain","Res","Radius","Valid","Select")
        for i, header in enumerate(contactsWithinHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)
            
        # Widgets lines
        for i,p in enumerate(MMB_UI.contactsWithin):
            wl = ContactWithinWidgetsLine(frame, i+1, p)
            contactWithinWidgets.append(wl)
            wl.display()

        if len(MMB_UI.polymers) > 0:
            m = ContactWithin(model=MMB_UI.currentModel)
            wl = ContactWithinWidgetsLine(frame, len(MMB_UI.contactsWithin)+1, m, 
                                            addCommand=lambda x=m:self.addContact(m))
            wl.display()
            self.newContactWithinWidgetsLine = wl

        return contactWithinWidgets


    ## Sets the residues and chain for a new contact from Chimera's current selection.
    #  Exactlt two residues must be selected.
    def setNewContactFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) == 0:
            tkMessageBox.showerror("MMB Error", "Current selection contains no residue")
            return
        chainIds = [r.id.chainId for r in res]
        if len(set(chainIds)) > 1:
            tkMessageBox.showerror("MMB Error", "Residues must be part of only one chain")
            return
        #else:
        resIds = [r.id.position for r in res]
        resId1 = min(resIds)
        resId2 = max(resIds)
        poly = MMB_UI.polymers[chainIds[0]]
        widgetLine = self.newContactWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widStartRes.set(poly.getResNameId(resId1))
        widgetLine.widEndRes.set(poly.getResNameId(resId2))

    ## Sets the residue and chain for a new contactWithin from Chimera's current selection.
    #  Exactlt two residues must be selected.
    def setNewContactWithinFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) != 1:
            tkMessageBox.showerror("MMB Error", "Exactly one residue must be selected")
            return
        res = res[0]
        chainId = res.id.chainId

        #else:
        resId = res.id.position
        poly = MMB_UI.polymers[chainId]
        widgetLine = self.newContactWithinWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widRes.set(poly.getResNameId(resId))


###################################################################################################
    ## Create the tabs for the mobilizers
    def createPhysicsTab(self):
        self.physicsTab = Tk.Frame(self.tabs)
        self.physicsTab.pack(fill='both')

        physicsFrame = Tk.Frame(self.physicsTab)
        physicsFrame.pack(fill='both', expand=True)

        self.deactivatePhysicsFrame = ScrolledFrame(physicsFrame, tkTags=("DeactivatePhysicsFrameTag",))
        self.drawDeactivatePhysicsFrame(self.deactivatePhysicsFrame)

        self.allResiduesWithinFrame = ScrolledFrame(physicsFrame, tkTags=("PhysicsFrameTag",))
        self.drawAllResiduesWithinFrame(self.allResiduesWithinFrame)

        frame2 = Tk.Frame(self.physicsTab, bd=2)
        frame2.pack(fill='both', expand=True)
        self.includeAllAtomsInResiduesFrame = ScrolledFrame(frame2, tkTags=("IncludeAllAtomsInResiduesFrameTag",))
        self.drawIncludeAllNonBondAtomsInResiduesFrame(self.includeAllAtomsInResiduesFrame)

        # Button bar
        buttons = Tk.Frame(self.physicsTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshPhysics)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh physics")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.deletePhysics)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected physics")

        fromSelButton = Tix.Button(buttons, text='AllResiduesWithin From Selection', command=self.setNewAllResiduesWithinFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new AllResiduesWithin chain and residues according to selection")

        fromSelIANBAIRButton = Tix.Button(buttons, text='IncludeAllNonBondAtomsInResidues From Selection', command=self.setIncludeAllNonBondAtomsInResiduesFromSelection)
        fromSelIANBAIRButton.pack(side='left')
        self.balloon.bind_widget(fromSelIANBAIRButton, msg="Set the new IncludeAllNonBondAtomsInResidues chain and residue according to selection")

        self.tabs.add(self.physicsTab, text="Physics")

    def drawDeactivatePhysicsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        physicsParamsFrame = Tk.Frame(frame)
        physicsParamsFrame.grid(row=0, columnspan=5)
        mdParamsButton = Tix.Button(physicsParamsFrame, text="Default MD Parameters", command=self.setDefaultMDParameters)
        mdParamsButton.pack(side='left')
        physicsRadiusEntry = ParameterEntry(physicsParamsFrame, "physicsRadius")
        physicsRadiusEntry.pack(side='left')
        self.balloon.bind_widget(physicsRadiusEntry, msg="Set default MD parameters for Physics zones")

        deactivatePhysicsWidgets = []
        #Headers
        deactivatePhysicsHeaders = ["", "Chain","Activated","Valid","Select"]
        for i, header in enumerate(deactivatePhysicsHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=1,column=i,sticky=W+E)

        # Widgets lines
        for i,p in enumerate(MMB_UI.polymers.values()):
            wl = DeactivatePhysicsWidgetsLine(frame, i+2, p)
            deactivatePhysicsWidgets.append(wl)
            wl.display()
        return deactivatePhysicsWidgets

    allResiduesWithinHeaders = ["id", "", "Radius", "Chain", "ResID","Valid","Select"]
    def drawAllResiduesWithinFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        allResiduesWithinWidgets = []
        # Headers
        for i, header in enumerate(MMBDialog.allResiduesWithinHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)
            
        # Widgets lines
        for i,p in enumerate(MMB_UI.allResiduesWithins):
            wl = AllResiduesWithinWidgetsLine(frame, i+1, p)
            allResiduesWithinWidgets.append(wl)
            wl.display()

        if len(MMB_UI.polymers) > 0:
            m = AllResiduesWithin(model=MMB_UI.currentModel,new=True)
            wl = AllResiduesWithinWidgetsLine(frame, len(MMB_UI.allResiduesWithins)+1, m, 
                                            addCommand= lambda x=m:self.addPhysics(x))
            wl.display()
            self.newAllResiduesWithinWidgetsLine = wl

        return allResiduesWithinWidgets

    includeAllNonBondAtomsInResiduesHeaders = ["id", "", "Chain", "ResID", "", "Valid","Select"]
    def drawIncludeAllNonBondAtomsInResiduesFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        includeAllNonBondAtomsInResiduesWidgets = []
        # Headers
        for i, header in enumerate(MMBDialog.includeAllNonBondAtomsInResiduesHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)
            
        # Widgets lines
        for i,p in enumerate(MMB_UI.includeAllNonBondAtomsInResidues):
            wl = IncludeNonBondAllAtomsInResidueWidgetsLine(frame, i+1, p)
            includeAllNonBondAtomsInResiduesWidgets.append(wl)
            wl.display()

        if len(MMB_UI.polymers) > 0:
            m = IncludeAllNonBondAtomsInResidue(model=MMB_UI.currentModel,new=True)
            wl = IncludeNonBondAllAtomsInResidueWidgetsLine(frame, len(MMB_UI.includeAllNonBondAtomsInResidues)+1, m, 
                                                            addCommand= lambda x=m:self.addPhysics(x))
            wl.display()
            self.newIncludeNonBondAllAtomsInResidueWidgetsLine = wl

        return includeAllNonBondAtomsInResiduesWidgets

    def addPhysics(self, p):
        try:
            p.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshPhysics()
        return True

    def refreshPhysics(self):
        try:
            MMB_UI.refreshAllResiduesWithins()
            MMB_UI.refreshIncludeAllNonBondAtomsInResidues()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.deactivatePhysicsWidgets = self.drawDeactivatePhysicsFrame(self.deactivatePhysicsFrame)
        self.allResiduesWithinWidgets = self.drawAllResiduesWithinFrame(self.allResiduesWithinFrame)
        self.includeAllNonBondAtomsInResiduesWidgets = self.drawIncludeAllNonBondAtomsInResiduesFrame(self.includeAllAtomsInResiduesFrame)
        return True

    def deletePhysics(self):
        self.deleteSelectedObjects(self.allResiduesWithinWidgets + self.includeAllNonBondAtomsInResiduesWidgets)
        self.refreshPhysics()

    def setNewAllResiduesWithinFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) != 1:
            tkMessageBox.showerror("MMB Error", "Exactly one residue must be selected")
            return
        res = res[0]
        chainId = res.id.chainId
        resId = res.id.position
        poly = MMB_UI.polymers[chainId]
        widgetLine = self.newAllResiduesWithinWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widRes.set(poly.getResNameId(resId))

    def setIncludeAllNonBondAtomsInResiduesFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) == 0:
            tkMessageBox.showerror("MMB Error", "Current selection contains no residue")
            return
        chainIds = [r.id.chainId for r in res]
        if len(set(chainIds)) > 1:
            tkMessageBox.showerror("MMB Error", "Residues must be part of only one chain")
            return
        #else:
        resIds = [r.id.position for r in res]
        resId1 = min(resIds)
        resId2 = max(resIds)
        poly = MMB_UI.polymers[chainIds[0]]
        widgetLine = self.newIncludeNonBondAllAtomsInResidueWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widRes.set(poly.getResNameId(resId1))
        widgetLine.widResEnd.set(poly.getResNameId(resId2))

###################################################################################################
    ## Create the tabs for the threading commands
    def createThreadingTab(self):
        self.threadingTab = Tk.Frame(self.tabs)
        self.threadingTab.pack(fill='both')

        threadingsFrame = Tk.Frame(self.threadingTab)
        threadingsFrame.config(relief="ridge",borderwidth=1)
        threadingsFrame.pack(fill='both')

        self.threadingFrame = ScrolledFrame(threadingsFrame, tkTags=("ThreadingsFrameTag",))
        self.drawThreadingsFrame(self.threadingFrame)

        gappedThreadingsFrame = Tk.Frame(self.threadingTab)
        gappedThreadingsFrame.config(relief="ridge",borderwidth=1)
        gappedThreadingsFrame.pack(fill='both', expand=True)

        self.gappedThreadingFrame = ScrolledFrame(gappedThreadingsFrame, tkTags=("GappedThreadingsFrameTag",))
        self.drawGappedThreadingsFrame(self.gappedThreadingFrame)

        # Button bar
        buttons = Tk.Frame(self.threadingTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshThreadings)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh threadings")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.delThreadings)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected mobilizer")

        fromSelButton = Tix.Button(buttons, text='Threading From Selection', 
                                   command=self.setNewThreadingFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new Threading according to selection")

        self.tabs.add(self.threadingTab, text="Threadings")

    ## Draw the frame containing the base interactions widgets
    def drawThreadingsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        threadingsWidgets = []

        threadingsTitle = Tk.Label(frame, text="Threading")
        threadingsTitle.grid(row=0, column=0, columnspan=8)
        # Headers
        threadingsHeaders=("Id","Chain1","Start1","End1", "Chain2", "Start2", "End2", "Constant", "Backbone Only","    ", "    ","Select")
        for i, header in enumerate(threadingsHeaders):
            threadingLabel = Tk.Label(frame,text=header,bd=2, relief='raised')
            threadingLabel.grid(row=1,column=i,sticky=W+E)

        for i,c in enumerate(MMB_UI.threadings):
            widgetsLine = ThreadingWidgetsLine(frame, i+2, c)
            threadingsWidgets.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            c = MMB_UI.Threading(model=MMB_UI.currentModel)
            widgetsLine = ThreadingWidgetsLine(frame, len(MMB_UI.threadings)+2, c, 
                                              addCommand=lambda x=c:self.addThreading(x))
            widgetsLine.display()
            self.newThreadingWidgetsLine = widgetsLine

        return threadingsWidgets

    ## Draw the frame containing the base interactions widgets
    def drawGappedThreadingsFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        gappedThreadingsWidgets = []

        gappedThreadingsTitle = Tk.Label(frame, text="Gapped Threading")
        gappedThreadingsTitle.grid(row=0, column=0, columnspan=8)
        # Headers
        gappedThreadingsHeaders=("Id","Chain1", "Chain2", "Constant", "Backbone Only", "    ", "    ", "Select")
        for i, header in enumerate(gappedThreadingsHeaders):
            gappedThreadingLabel = Tk.Label(frame,text=header,bd=2, relief='raised')
            gappedThreadingLabel.grid(row=1,column=i,sticky=W+E)

        for i,c in enumerate(MMB_UI.gappedThreadings):
            widgetsLine = GappedThreadingWidgetsLine(frame, i+2, c)
            gappedThreadingsWidgets.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            c = MMB_UI.GappedThreading(model=MMB_UI.currentModel)
            widgetsLine = GappedThreadingWidgetsLine(frame, len(MMB_UI.gappedThreadings)+2, c, 
                                              addCommand=lambda x=c:self.addThreading(x))
            widgetsLine.display()
            self.newGappedThreadingWidgetsLine = widgetsLine

        return gappedThreadingsWidgets

    ## Synchronize the threadings with MMB
    def refreshThreadings(self):
        try:
            MMB_UI.refreshThreadings()
            MMB_UI.refreshGappedThreadings()
            pass
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.threadingsWidgets = self.drawThreadingsFrame(self.threadingFrame)
        self.gappedThreadingsWidgets = self.drawGappedThreadingsFrame(self.gappedThreadingFrame)
        return True

    ## Add a new Threading to MMB
    #  @param thread a MMB_UI.Threading
    def addThreading(self, thread):
        try:
            thread.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshThreadings()
        return True

    ## Remove the selected threadings from MMB
    def delThreadings(self):
        self.deleteSelectedObjects(self.threadingsWidgets + self.gappedThreadingsWidgets)
        self.refreshThreadings()
        return True

    # ## Sets the chains, residues and atoms for a new threading from Chimera's current selection.
    def setNewThreadingFromSelection(self):
        wl = self.newThreadingWidgetsLine

        res = chimera.selection.currentResidues()
        if len(res) == 0:
            tkMessageBox.showerror("MMB Error", "Selection must contain at least 2 residues on one chain and at least 1 on another chain")
            return
        chainIds = list(set([r.id.chainId for r in res]))
        if len(chainIds) != 2:
            tkMessageBox.showerror("MMB Error", "Selection must contain at least 2 residues on one chain and at least 1 on another chain")
            return
        #else:
        resIds1 = [r.id.position for r in res if r.id.chainId == chainIds[0]]
        start1 = min(resIds1)
        end1 = max(resIds1)

        resIds2 = [r.id.position for r in res if r.id.chainId == chainIds[1]]
        start2 = min(resIds2)
        end2 = max(resIds2)

        # Complete the shortest stretch to have the same length
        len1 = end1 - start1
        len2 = end2 - start2
        if len1 > len2:
            end2 = start2 + len1
        else:
            end1 = start1 + len2

        poly = MMB_UI.polymers[chainIds[0]]
        if not poly.getResName(end1):
            tkMessageBox.showerror("MMB Error", "Cannot complete the Threading stretch in chain "+poly.chainID)
            return
        wl.widChain1.set(poly.chainID)
        wl.widStartRes1.set(poly.getResNameId(start1))
        wl.widEndRes1.set(poly.getResNameId(end1))

        poly = MMB_UI.polymers[chainIds[1]]
        if not poly.getResName(end2):
            tkMessageBox.showerror("MMB Error", "Cannot complete the Threading stretch in chain "+poly.chainID)
            return
        wl.widChain2.set(poly.chainID)
        wl.widStartRes2.set(poly.getResNameId(start2))
        wl.widEndRes2.set(poly.getResNameId(end2))

###################################################################################################
    ## Create the tabs for the densities
    def createDensitiesTab(self):
        self.densityTab = Tk.Frame(self.tabs)
        self.densityTab.pack(fill='both')

        densityParametersFrame = Tk.Frame(self.densityTab)
        densityParametersFrame.pack(side="top")
        self.openMapButton = LabeledButton(densityParametersFrame, labelText=ParametersTraces.get("densityFileName"), buttonText="Open Map",
                                           command=lambda ipd=OpenXplorMapDialog: ipd(command=self.openXplorMap, 
                                                                                      multiple=False).enter()
                                           )
        self.openMapButton.pack(side='left')

        densityForceEntry = ParameterEntry(densityParametersFrame, "densityForceConstant")
        densityForceEntry.pack(side='left')

        densitiesFrame = Tk.Frame(self.densityTab)
        densitiesFrame.pack(fill='both', expand=True)

        self.densityFrame = ScrolledFrame(densitiesFrame, tkTags=("DensitiesFrameTag",))
        self.drawDensitiesFrame(self.densityFrame)

        # Button bar
        buttons = Tk.Frame(self.densityTab, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refreshButton = Tix.Button(buttons, text='Refresh', command=self.refreshDensities)
        refreshButton.pack(side='left')
        self.balloon.bind_widget(refreshButton, msg="Refresh densities")

        deleteButton = Tix.Button(buttons, text='Delete', command=self.deleteDensities)
        deleteButton.pack(side='left')
        self.balloon.bind_widget(deleteButton, msg="Delete selected density")

        fromSelButton = Tix.Button(buttons, text='Density From Selection', command=self.setNewDensityFromSelection)
        fromSelButton.pack(side='left')
        self.balloon.bind_widget(fromSelButton, msg="Set the new Density chain and residues according to selection")

        self.tabs.add(self.densityTab, text="Densities")

    ## Draw the frame containing the densities widgets
    def drawDensitiesFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        densityWidgetsLine = []

        densitiesHeaders=("Id","Density","Chain","Start Res","End Res","Valid","Select")
        # Headers
        for i, header in enumerate(densitiesHeaders):
            label = Tk.Label(frame,text=header,bd=2, relief='raised')
            label.grid(row=0,column=i,sticky=W+E)

        # Widgets lines
        for i,p in enumerate(MMB_UI.densities):
            widgetsLine = DensityWidgetsLine(frame, i+1, p)
            densityWidgetsLine.append(widgetsLine)
            widgetsLine.display()

        if len(MMB_UI.polymers) > 0:
            m = Density()
            widgetsLine = DensityWidgetsLine(frame, len(MMB_UI.densities)+1, m, 
                                               addCommand= lambda x=m: self.addDensity(m))
            widgetsLine.display()
            self.newDensityWidgetsLine = widgetsLine

        return densityWidgetsLine

    ## Synchronize the densities with MMB
    def refreshDensities(self):
        try:
            MMB_UI.refreshDensityStretches()
            self.openMapButton.setLabelText(ParametersTraces.get("densityFileName"))
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.densityWidgets = self.drawDensitiesFrame(self.densityFrame)
        return True

    ## Remove the selected base pair interactions from MMB
    def deleteDensities(self):
        self.deleteSelectedObjects(self.densityWidgets)
        self.refreshDensities()
        return True

    ## Add a new density stretch to MMB
    #  @param density a MMB_UI.Density
    def addDensity(self, density):
        try:
            density.mmbAdd()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        self.refreshDensities()
        return True

    def setNewDensityFromSelection(self):
        res = chimera.selection.currentResidues()
        if len(res) == 0:
            tkMessageBox.showerror("MMB Error", "Current selection contains no residue")
            return
        chainIds = [r.id.chainId for r in res]
        if len(set(chainIds)) > 1:
            tkMessageBox.showerror("MMB Error", "Residues must be part of only one chain")
            return
        #else:
        resIds = [r.id.position for r in res]
        resId1 = min(resIds)
        resId2 = max(resIds)
        poly = MMB_UI.polymers[chainIds[0]]
        widgetLine = self.newDensityWidgetsLine
        widgetLine.widChain.set(poly.chainID)
        widgetLine.widStartRes.set(poly.getResNameId(resId1))
        widgetLine.widEndRes.set(poly.getResNameId(resId2))

    def openXplorMap(self, filepath):
        if(os.path.isfile(filepath)):
            try:
                MMB_UI.openDensityMap(filepath)
                ParametersTraces.set("densityFileName", filepath)
            except pyMMB.MMBError as e:
                tkMessageBox.showerror("MMB Error", e.msg)
                return False
            self.openMapButton.setLabelText(ParametersTraces.get("densityFileName"))
            return True
        return False

###################################################################################################
    ## Create the tab for the control of the simulation
    def createSimulationTab(self):
        simFrame = Tk.Frame(self.tabs)
        simFrame.pack(fill='both')

        self.paramsFrame = Tk.Frame(simFrame)
        self.drawParamsFrame(self.paramsFrame)
        self.paramsFrame.pack(fill='both')

        self.simControlFrame = Tk.Frame(simFrame, bd=2, relief="sunken")
        self.simControlFrame.pack(fill='both')
        self.drawSimulationControlFrame(self.simControlFrame)

        self.simInfoFrame = Tk.Frame(simFrame)
        self.simInfoFrame.pack(fill='both', expand=True)
        self.drawSimulationInfoFrame(self.simInfoFrame)


        # Button bar
        buttons = Tk.Frame(simFrame, bd=2, relief="sunken")
        buttons.pack(side='bottom', fill='both')

        refresh = Tix.Button(buttons, text='Refresh', command=self.refreshSimulationTab)
        refresh.pack(side='left')

        self.tabs.add(simFrame, text="Simulation")

    def refreshSimulationTab(self):
        self.drawParamsFrame(self.paramsFrame)
        self.drawSimulationControlFrame(self.simControlFrame)
        self.drawSimulationInfoFrame(self.simInfoFrame)

    def drawParamsFrame(self, frame):
        [w.destroy() for w in frame.children.values()]

        reportingIntervalEntry = ParameterEntry(frame, "reportingInterval")
        reportingIntervalEntry.grid(row=0, column=0, sticky=Tk.S+Tk.N+Tk.E)

        # numReportingIntervalsEntry = ParameterEntry(frame, "numReportingIntervals")
        # numReportingIntervalsEntry.grid(row=1, column=0, sticky=Tk.S+Tk.N+Tk.E)

        setInitialVelocities = ParameterEntry(frame, "setInitialVelocities")
        setInitialVelocities.grid(row=2, column=0, sticky=Tk.S+Tk.N+Tk.E)

        # integratorAccuracy = ParameterEntry(frame, "integratorAccuracy")
        # integratorAccuracy.grid(row=3, column=0, sticky=Tk.S+Tk.N+Tk.E)

        temperatureEntry = ParameterEntry(frame, "temperature")
        temperatureEntry.grid(row=0, column=1, sticky=Tk.S+Tk.N+Tk.E)

        removeRigidBodyMomentum = ParameterEntry(frame, "removeRigidBodyMomentum")
        removeRigidBodyMomentum.grid(row=1, column=1, sticky=Tk.S+Tk.N+Tk.E)

        smallGroupInertiaMultiplier = ParameterEntry(frame, "smallGroupInertiaMultiplier")
        smallGroupInertiaMultiplier.grid(row=2, column=1, sticky=Tk.S+Tk.N+Tk.E)

        dutyCycleEntry = ParameterEntry(frame, "dutyCycle")
        dutyCycleEntry.grid(row=0, column=2, sticky=Tk.S+Tk.N+Tk.E)
        self.balloon.bind_widget(dutyCycleEntry, msg="For a fraction of the time (1 - dutyCycle) all forces will be turned off")

        scrubberPeriodEntry = ParameterEntry(frame, "scrubberPeriod")
        scrubberPeriodEntry.grid(row=1, column=2, sticky=Tk.S+Tk.N+Tk.E)
        self.balloon.bind_widget(scrubberPeriodEntry, msg="potential rescaling period, in ps")

        forceMultiplierEntry = ParameterEntry(frame, "twoTransformForceMultiplier", guitext="Force Multiplier")
        forceMultiplierEntry.grid(row=3, column=0, sticky=Tk.S+Tk.N+Tk.E)

        allparamsButton = Tix.Button(frame, text="All Parameters", command=self.showAllParametersDialog)
        allparamsButton.grid(row=4,column=0)
        self.balloon.bind_widget(allparamsButton, msg="For experimented users !!!")

        convergeCheck = ParameterEntry(frame, "detectConvergence")
        convergeCheck.grid(row=0, column=3, sticky=Tk.S+Tk.N+Tk.E)

        convergeTimeout = ParameterEntry(frame,"convergenceTimeout")
        convergeTimeout.grid(row=1, column=3, sticky=Tk.S+Tk.N+Tk.E)

        convergeEpsilon = ParameterEntry(frame,"convergenceEpsilon")
        convergeEpsilon.grid(row=2, column=3, sticky=Tk.S+Tk.N+Tk.E)

    def drawSimulationControlFrame(self, frame):
        [w.destroy() for w in frame.children.values()]

        self.initButton = Tix.Button(frame, text='Init', command=self.initSimulation)
        self.initButton.pack(side='left')
        changeWidgetState(self.initButton, True)

        self.numIntervalEntry = IntegerEntry(frame, value=20, width=5)
        self.numIntervalEntry.pack(side='left')
        self.balloon.bind_widget(self.numIntervalEntry, msg="Number of Reporting Interval for this run")

        self.runButton = Tix.Button(frame, text='Run', command=self.runIntervals)
        self.runButton.pack(side='left')
        changeWidgetState(self.runButton, False)

        # self.pauseButton = Tix.Button(frame, text='Pause', command=self.pauseSimulation)
        # self.pauseButton.pack(side='left')
        # changeWidgetState(self.pauseButton, False)

        # self.continueButton = Tix.Button(frame, text='Continue', command=self.continueSimulation)
        # self.continueButton.pack(side='left')
        # changeWidgetState(self.continueButton, False)

        self.stopButton = Tix.Button(frame, text='Stop', command=self.stopSimulation)
        self.stopButton.pack(side='left')
        changeWidgetState(self.stopButton, False)

    def drawSimulationInfoFrame(self, frame):
        # Clear the frame first
        [w.destroy() for w in frame.children.values()]

        textLabel = Tk.Label(frame, text="Current Step")
        textLabel.grid(row=0, column=0)
        stepsLabel = Tk.Label(frame, text=MMB_UI.getCurrentFrameNumber())
        stepsLabel.grid(row=0, column=1)

        oks,notoks = MMB_UI.getSatisfiedBasePairs()

        textLabel = Tk.Label(frame, text="# Satisfied base pairs:")
        textLabel.grid(row=1, column=0)
        satisfiedBPLabel = Tk.Label(frame, text=oks)
        satisfiedBPLabel.grid(row=1, column=1)
        textLabel = Tk.Label(frame, text="# Unsatisfied base pairs:")
        textLabel.grid(row=2, column=0)
        unsatisfiedBPLabel = Tk.Label(frame, text=notoks)
        unsatisfiedBPLabel.grid(row=2, column=1)

        potEnergyEntry = ParameterEntry(frame,"potentialEnergy", guitext="Potential Energy");
        potEnergyEntry.grid(row=3,column=0)
        changeWidgetState(potEnergyEntry.entry, False)

        potKineticEntry = ParameterEntry(frame,"kineticEnergy", guitext="Kinetic Energy");
        potKineticEntry.grid(row=4,column=0)
        changeWidgetState(potKineticEntry.entry, False)

        potTotalEntry = ParameterEntry(frame,"totalEnergy", guitext="Total Energy");
        potTotalEntry.grid(row=5,column=0)
        changeWidgetState(potTotalEntry.entry, False)

        convergedCheck = ParameterEntry(frame, "converged")
        convergedCheck.grid(row=6, column=0)
        changeWidgetState(convergedCheck.entry, False)

    ## Init MMB's simulation parameters and MDMovie dialog
    def initSimulation(self):
        try:
            MMB_UI.initSimulation()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        changeWidgetState(self.initButton, False)
        changeWidgetState(self.runButton, True)
        [changeWidgetState(child, False) for child in self.paramsFrame.children.values()]
        return True

    ## Execute n intervals until the number of intervals specified in MMB.
    #  @param n int
    def runIntervals(self):
        try:
            print "runIntervals"
            MMB_UI.runIntervals(n=self.numIntervalEntry.get(), 
                                callback=lambda f=self.simInfoFrame: self.drawSimulationInfoFrame(f),
                                endCallback=self.stopSimulation)
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        changeWidgetState(self.runButton, False)
        # changeWidgetState(self.pauseButton, True)
        changeWidgetState(self.stopButton, True)
        [changeWidgetState(child, False) for child in self.paramsFrame.children.values()]
        return True

    ## Pause the current running simulation
    def pauseSimulation(self):
        try:
            MMB_UI.pauseSimulation()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        # changeWidgetState(self.pauseButton, False)
        # changeWidgetState(self.continueButton, True)
        return True

    ## Continue a currently paused simulation
    def continueSimulation(self):
        try:
            MMB_UI.continueSimulation()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        # changeWidgetState(self.pauseButton, True)
        # changeWidgetState(self.continueButton, False)
        return True

    ## Stop the current simulation
    def stopSimulation(self):
        try:
            MMB_UI.stopSimulation()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False

        changeWidgetState(self.runButton, True)
        # changeWidgetState(self.pauseButton, False)
        # changeWidgetState(self.continueButton, False)
        changeWidgetState(self.stopButton, False)
        [changeWidgetState(child, True) for child in self.paramsFrame.children.values()]
        return True

    ## Display a dialog with widgets for all MMB parameters
    def showAllParametersDialog(self):
        self.allParamsDialog = chimera.dialogs.display(AllParametersDialog.name)

    def setDefaultMDParameters(self):
        try:
            MMB_UI.setDefaultMDParameters()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False

###################################################################################################

    # Callbacks ------------------------------------------------
 
    ## Callback for previous tab buttons
    def previousTab(self):
        curTab = self.tabs.index(self.tabs.select())
        self.tabs.select(curTab - 1)

    ## Callback for next tab buttons. If command return None, the tab is unchanged.
    #  @param command the function to exectute before changing the tab.
    def nextTab(self, command=None):
        """
        param: command is the function to execute before changing tab
        """
        if hasattr(command, '__call__'):
            # if command return False, don't go to the next tab
            if not command() : return
        curTab = self.tabs.index(self.tabs.select())
        self.tabs.select(curTab + 1)

 

    ## Delete selected objects among a list of widgetsLines
    def deleteSelectedObjects(self, widgetsLines):
        try:
            # MMB IDs change after each deletion    
            # So we have to make sure we use the good ID
            # First: extract and sort selected MMBObjects via WidgetsLines
            lines = [l for l in widgetsLines if l.varSelect.get()]
            lines = sorted(lines, key=lambda wl:wl.data.mmbID)

            # Second: delete taking in to account the IDs decrement (-i)
            for i,line in enumerate(lines):
                line.data.mmbID = line.data.mmbID-i
                line.data.mmbDelete()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)
            return False
        return True

    ## Refresh all data
    def refreshAll(self):
        self.refreshSequences()
        self.refreshMobilizers()
        self.refreshContacts()
        self.refreshDensities()
        self.refreshConstraints()
        self.refreshThreadings()
        self.refreshPhysics()
        self.refreshBasePairInteractions()
        self.refreshSimulationTab()

        # If MMB is not initialized, disable tabs except input
        state = "normal"
        if not MMB_UI.polymersInitialized:
            state = "hidden"
            #state = "disabled" # tabs are inactive but not grayed... (on Mac at least) 
        for i in range(2, len(self.tabs.tabs())):
            self.tabs.tab(i, state=state)

    ## Synchronize the base pair interaction form selections with the current Chimera's selection.
    #  Called each time the selection is changed.
    def synchronizeBondSelection(self, trigger, additional, currentSel):
        if not self.stopRefreshSelection:
            self.stopRefreshSelection = True
            selBondsChimera = currentSel.pseudobonds()
            if len(selBondsChimera) == 0:
                self.stopRefreshSelection = False
            widgets = self.basePairWidgets + self.constraintsWidgets
            selBondsGUI = [p.data.chimeraBond for p in widgets if p.varSelect.get()]

            # Deselect Chimera's unselected bonds
            notInter = [b for b in selBondsGUI if b not in selBondsChimera]
            pairs = [p for p in widgets if p.data.chimeraBond in notInter]
            [p.select(False) for p in pairs]

            # Select Chimera's selected bonds
            pairs = [p for p in widgets if p.data.chimeraBond in selBondsChimera]
            [p.select(True) for p in pairs]
        else:
            self.stopRefreshSelection = False

    ## Callback for the Save button
    #  Open a Save dialog to create a dump file of MMB commands according to the forms
    def Save(self):
        from OpenSave import SaveModeless
        SaveModeless(command=self.createDumpFile)

    ## Callback for the Save dialog
    def createDumpFile(self, okayed, dialog):
        if not okayed:
            return
        filename = dialog.getPaths()[0]
        comms = MMB_UI.dumpCommands(filename)
        r = Tk.Tk()
        r.withdraw()
        r.clipboard_clear()
        r.clipboard_append(comms)
        r.destroy()

    def ApplyViz(self):
        MMB_UI.applyMMBVisualization()

    def ResetMMB(self):
        commands = MMB_UI.commandsToStr()
        previous_model = MMB_UI.currentModel
        print "previous_model:",previous_model
        MMB_UI.initMMB()
        self.refreshAll()

        result = tkMessageBox.askquestion("Reset MMB","Would you like to reload the last structure and parameters?")
        if result == "yes":
            MMB_UI.importSequencesFromChimera(models=[previous_model])
            mmbDiag.refreshSequences()
            self.cmdGlobalEditor.delete(1.0,'end')
            self.cmdGlobalEditor.insert(1.0, commands)
        
        self.tabs.select(1)
        
###############################################################################
# Auto-fill the GUI
class FakeWidget:
    def __init__(self,x):
        self.x = x
    def get(self):
        return self.x

mmbDiag = None
def runMMBUI():
    global mmbDiag
    mmbDiag = chimera.dialogs.display(MMBDialog.name)
    # pdbname = "/Users/alex/Development/MMB_gui/tests/1ARJ.short.pdb"
    #pdbname = "/Users/alex/Dropbox/MMB/Paper/tRNA_mRNA.pdb"
    # pdbname = "/Users/alex/Development/MMB_gui/tests/tRNA.pdb"
    # pdbname = "/Users/alex/Projects/5-DEBS/hinges/2QO3/2QO3_AB_full.pdb"
    # MMB_UI.loadSequencesFromPdb(pdbname)
    # pyMMB.cmd("protein H 522 RSPGVGCVPAAEHRLREEILAKFLHWLMSVYVVELLRSFFYVTETTFQKNRLFFYRKSV")
    # mmbDiag.refreshSequences()
    # mmbDiag.loadPolymers()
    # mmbDiag.openCmdsFile(FakeWidget("/Users/alex/Development/MMB_gui/tests/commands.TAR.gui.dat"))
    # mmbDiag.openCmdsFile(FakeWidget("/Users/alex/Development/MMB_gui/tests/commands.tRNA-fitting.gui.dat"))
    # mmbDiag.openCmdsFile(FakeWidget("/Users/alex/Projects/5-DEBS/hinges/2QO3/commands.gui.tests.dat"))
    # mmbDiag.openCmdsFile(FakeWidget("/Users/alex/Dropbox/MMB/Paper/tRNA_mRNA.dat"))
    # mmbDiag.applyCmds()
    # mmbDiag.tabs.select(9)

def runtRNAExample():
    mmbDiag.ResetMMB()
    pdbname = os.path.join(MMB_UI.mmbUIPath,"examples/tRNA_mRNA.pdb")
    MMB_UI.loadSequencesFromPdb(pdbname)
    mmbDiag.refreshSequences()
    mmbDiag.loadPolymers()
    mmbDiag.openCmdsFile(FakeWidget(os.path.join(MMB_UI.mmbUIPath,"examples/tRNA_mRNA.dat")))
    mmbDiag.applyCmds()
    MMB_UI.applyMMBVisualization()

def runRNAExample():
    mmbDiag.ResetMMB()
    pdbname = os.path.join(MMB_UI.mmbUIPath,"examples/1ARJ.short.pdb")
    MMB_UI.loadSequencesFromPdb(pdbname)
    mmbDiag.refreshSequences()
    mmbDiag.loadPolymers()
    mmbDiag.openCmdsFile(FakeWidget(os.path.join(MMB_UI.mmbUIPath,"examples/commands.TAR.gui.dat")))
    mmbDiag.applyCmds()
    MMB_UI.applyMMBVisualization()

def runPeptideExample():
    mmbDiag.ResetMMB()
    pyMMB.cmd("protein A 1 AAAPAFWAA")
    mmbDiag.refreshSequences()
    mmbDiag.loadPolymers()

def runMorph():
    mmbDiag.ResetMMB()
    mmbDiag.importFromChimera()
    mmbDiag.refreshSequences()
    mmbDiag.loadPolymers()
    mmbDiag.openCmdsFile(FakeWidget("/Users/alex/Projects/6-Paper/Nucleotidase/commands.morph_1OID.gui.dat"))
    mmbDiag.applyCmds()

###############################################################################
#   Now we register the above dialog with Chimera, so that it may be 
#   invoked via the 'display(name)' method of the chimera.dialogs module.
#   Here the class itself is registered, but since it is a named dialog
#   deriving from ModalDialog/ModelessDialog, the instance will automatically
#   reregister itself when first created.  This allows the 'dialogs.find()'
#   function to return the instance instead of the class.
chimera.dialogs.register(MMBDialog.name, MMBDialog)
chimera.dialogs.register(AllParametersDialog.name, AllParametersDialog)

#    Create the Chimera toolbar button that displays the dialog when
#    pressed.  Note that since the package is not normally searched for
#    icons, we have to prepend the path of this package to the icon's
#    file name.
dir, file = os.path.split(__file__)
icon = os.path.join(dir, 'MMB_icon.tiff')
# chimera.tkgui.app.toolbar.add(icon, lambda d=chimera.dialogs.display, n=MMBDialog.name: d(n), 'MMB configuration form', None)



