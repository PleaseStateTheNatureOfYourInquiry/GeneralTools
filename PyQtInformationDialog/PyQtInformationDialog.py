# Author: Maarten Roos-Serote
# ORCID author: 0000 0001 5001 1347

# Version: 20240527


# Qt imports.
from PyQt6.QtWidgets import QMessageBox, QLabel


# A handy general class to create an information dialog window for PyQt6.
class informationDialog ():
    '''
    :param informationType: type of information window: 'warning' (default) or 'information', determines its layout / look.
    :type informationType: str

    :param windowTitle: title printed at the top bar of the window.
    :type windowTitle: str

    :param messageText: message text printed in the window.
    :type messageText: str

    :param yesAndCancel: if ``True`` and  numberOfChoices == 2, then create **Yes** and **Cancel** buttons, with **Cancel** button set to default, else create **Yes** and **No** button, with **Yes** the default.
    :type yesAndCancel: bool
    
    :param numberOfChoices: When set to 1, create one **OK** button. When set to 2, create **Yes* and **Cancel** or **No** buttons (see yesAndCancel). When set to 3, create **Yes**, **No**, **Cancel** buttons, with **Cancel** as default.
    :type numberOfChoices: int

    **Description:##
    Create an information dialog window, with one (OK), two (Yes, Cancel) or three (Yes, No, Cancel) options.
    '''

    def __init__ (self, informationType = 'warning', windowTitle = 'Attention', messageText = '', yesAndCancel = False, numberOfChoices = 3):
        
        self.informationDialogBox = QMessageBox ()
        
        self.informationDialogBox.setWindowTitle (windowTitle)
        
        self.informationDialogBox.setText (messageText)
        
        if informationType == 'warning':
        
            self.informationDialogBox.setIcon (QMessageBox.Icon.Warning)
            
        elif informationType == 'information':

            self.informationDialogBox.setIcon (QMessageBox.Icon.Information)

        if numberOfChoices == 3:
        
            self.informationDialogBox.setStandardButtons (QMessageBox.StandardButton.Cancel | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Yes)
            self.informationDialogBox.setDefaultButton (QMessageBox.StandardButton.Cancel)
 
        elif numberOfChoices == 2: 
        
            if yesAndCancel:

                self.informationDialogBox.setStandardButtons (QMessageBox.StandardButton.Cancel | QMessageBox.StandardButton.Yes)
                self.informationDialogBox.setDefaultButton (QMessageBox.StandardButton.Yes)
            
            else:
            
                self.informationDialogBox.setStandardButtons (QMessageBox.StandardButton.No | QMessageBox.StandardButton.Yes)
                self.informationDialogBox.setDefaultButton (QMessageBox.StandardButton.No)


        elif numberOfChoices == 1:
        
            self.informationDialogBox.setStandardButtons (QMessageBox.StandardButton.Ok)
            self.informationDialogBox.setDefaultButton (QMessageBox.StandardButton.Ok)

