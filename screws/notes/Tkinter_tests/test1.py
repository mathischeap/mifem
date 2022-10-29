from tkinter import *

window = Tk()
window.title("First Window")
window.geometry("350x200")
lbl1 = Label(window, text="Hello")
lbl1.grid(column=0, row=0)


lbl2 = Label(window, text="world")
lbl2.grid(column=300, row=200)

window.mainloop()