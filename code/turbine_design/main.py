from tkinter import *

from gui import initial_screen
import turbine_design


def main():
    # Creating a user interface app
    print('roncha')
    myGui = Tk()
    myGui.title("Axial Turbine Design")
    myGui_width = 1600
    myGui_height = 900
    myGui_posx = (myGui.winfo_screenwidth() // 2) - (myGui_width // 2)
    myGui_posy = (myGui.winfo_screenheight() // 2) - (myGui_height // 2) - 50
    myGui.geometry('{}x{}+{}+{}'.format(myGui_width,
                                        myGui_height, myGui_posx, myGui_posy))
    window = initial_screen(myGui)
    myGui.mainloop()


if __name__ == '__main__':
    main()
