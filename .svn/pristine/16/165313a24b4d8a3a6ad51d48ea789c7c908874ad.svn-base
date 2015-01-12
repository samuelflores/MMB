import chimera.extension

class MMB_UIEMO(chimera.extension.EMO):
    def name(self):
        return 'MMB User Interface'
    def description(self):
        return 'GUI for MacroMoleculeBuilder (MMB).'
    def categories(self):
        return ['MMB']
    def icon(self):
        return self.path('MMB_icon.tiff')
    def activate(self):
        self.module('gui').runMMBUI()
        return None

chimera.extension.manager.registerExtension(MMB_UIEMO(__file__))

# def cmdRotamers(cmdName, args):
#     from Midas.midas_text import doExtensionFunc
#     from Rotamers import useBestRotamers
#     doExtensionFunc(useBestRotamers, args, specInfo=[
#         ("spec", "targets", "residues"),
#         ("densitySpec", "density", "models")
#     ])
# from Midas.midas_text import addCommand
# addCommand("swapaa", cmdRotamers, help=True)