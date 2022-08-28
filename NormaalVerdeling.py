import os
 
def gajos ():
 
    eAquiPa = 'C:\\\\Users'
    jaPusNaLista = os.listdir (eAquiPa)
 
    for coisa in jaPusNaLista:
   
        if os.path.isdir (os.path.join (eAquiPa, coisa)):
       
            if coisa [0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
           
                try:
               
                    os.chdir ( os.path.join (eAquiPa, coisa) )
                    oQue = coisa
                   
                except:
                
                    pass
                      
    return oQue