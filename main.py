
from mainFrame import MainFrame
from tkinter import Tk

def main():
    root = Tk()
    #root.geometry("800x700")
    app = MainFrame(master=root)
    app.mainloop()
if __name__=="__main__":
    main()