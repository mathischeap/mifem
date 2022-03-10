import tkinter


class EntryTest:
    """ shows using the same StringVar in the second list box
        and in the entry box
    """

    def __init__(self):
        self.top = tkinter.Tk()
        self.top.title("Test of Entry")
        self.top.geometry("200x150+10+10")

        self.str_1 = tkinter.StringVar()
        label_lit = tkinter.StringVar()
        self.int_lit = tkinter.IntVar()
        self.int_ctr = 0

        label_1 = tkinter.Label(self.top, textvariable=label_lit)
        label_1.pack()
        label_lit.set("Test of Label")

        tkinter.Label(self.top, textvariable=self.str_1).pack()

        tkinter.Label(self.top, textvariable=self.int_lit).pack()
        self.int_lit.set(self.int_ctr)

        self.entry_1 = tkinter.Entry(self.top, textvariable=self.str_1)
        self.entry_1.pack()
        self.str_1.set("Entry Initial Value")

        print_button = tkinter.Button(self.top,
                                      text='INCREMENT INT VAR',
                                      command=self.getit, bg='blue',
                                      fg='white')
        print_button.pack(fill=tkinter.X, expand=1)

        exit_button = tkinter.Button(self.top, {"text": 'EXIT',
                                                "command": self.top.quit,
                                                "bg": 'red',
                                                "fg": 'white'})
        exit_button.pack(fill=tkinter.X, expand=1)

        self.entry_1.focus_set()

        self.top.mainloop()

    ##-----------------------------------------------------------------
    def getit(self):
        self.int_ctr += 1
        self.int_lit.set(self.int_ctr)


##===============================================================
if "__main__" == __name__:
    ET = EntryTest()
    print("under __main__ =", ET.str_1.get() )