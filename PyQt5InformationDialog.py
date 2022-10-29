from PyQt5 import QtWidgets


# A handy general class to create an information dialog window for PyQt5.
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
        
        self.informationDialogBox = QtWidgets.QMessageBox ()
        
        self.informationDialogBox.setWindowTitle (windowTitle)
        
        self.informationDialogBox.setText (messageText)

        if informationType == 'warning':
        
            self.informationDialogBox.setIcon (QtWidgets.QMessageBox.Warning)
            
        elif informationType == 'information':

            self.informationDialogBox.setIcon (QtWidgets.QMessageBox.Information)

        if numberOfChoices == 3:
        
            self.informationDialogBox.setStandardButtons (QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes)
            self.informationDialogBox.setDefaultButton (QtWidgets.QMessageBox.Cancel)
 
        elif numberOfChoices == 2: 
        
            if yesAndCancel:

                self.informationDialogBox.setStandardButtons (QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Yes)
                self.informationDialogBox.setDefaultButton (QtWidgets.QMessageBox.Yes)
            
            else:
            
                self.informationDialogBox.setStandardButtons (QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes)
                self.informationDialogBox.setDefaultButton (QtWidgets.QMessageBox.No)


        elif numberOfChoices == 1:
        
            self.informationDialogBox.setStandardButtons (QtWidgets.QMessageBox.Ok)
            self.informationDialogBox.setDefaultButton (QtWidgets.QMessageBox.Ok)

