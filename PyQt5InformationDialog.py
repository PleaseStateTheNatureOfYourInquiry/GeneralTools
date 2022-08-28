from PyQt5 import QtWidgets


# A handy general class to create an information dialog window.
class informationDialog ():
    '''
    Create an information dialog window, with one (OK), two (Yes, Cancel) or three (Yes, No, Cancel) options
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

