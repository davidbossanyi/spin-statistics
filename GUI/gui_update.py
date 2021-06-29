import PyQt5.uic

def update_gui(ui_name):
    fpath = '{0}.ui'.format(ui_name)
    with open('{0}_gui.py'.format(ui_name),'w') as file:
        PyQt5.uic.compileUi(fpath, file)
        
if __name__ == '__main__':
    update_gui('main')
    update_gui('popup')
