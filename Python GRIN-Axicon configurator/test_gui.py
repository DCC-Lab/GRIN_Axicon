import tkinter as tk
from tkinter import N, S, E, W, NE, NW, SE, SW

from tkinter import ttk
import tkinter.font as font
import unittest
import re 

gRoot = None

class TestTk(unittest.TestCase):
    def testInit(self):
        self.assertTrue(True)

    def testLoadTk(self):
        root = tk.Tk()
        self.assertIsNotNone(root)

    def testShowTkWindow(self):
        root = tk.Tk()
        root.title = "Test"

    """
    Reading this https://www.tutorialspoint.com/how-to-get-the-screen-size-in-tkinter#:~:text=Tkinter%20components%20adjust%20the%20window,of%20the%20screen%20in%20pixels.
    I come up with the following test to confirm I understand what is going on.
    """
    def testWindowGeometry(self):
        root = tk.Tk()
        root.title = "Test"
        geometry = root.geometry() # Something like that: '774x396+34+59' widthxheight+offsetx+offsety
        match = re.search(r"(\d+)x(\d+)\D(\d+)\D(\d+)", geometry)
        if match is not None:
            groups = match.groups()
            # Window is 1x1+0+0 if not showed yet
            self.assertTrue(int(groups[0]) == 1)
            self.assertTrue(int(groups[1]) == 1)
            self.assertTrue(int(groups[2]) == 0)
            self.assertTrue(int(groups[3]) == 0)

    def testScreenGeometry(self):
        root = tk.Tk()
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()

        # Set this to the resolution of your screen: I have no way of knowing at this point
        # and certainly not across plateforms

        actual_width = 1792
        actual_height = 1120
        self.assertEqual(screen_width , actual_width)
        self.assertEqual(screen_height , actual_height)

    def testDPI(self):
        root = tk.Tk()
        dpi = root.winfo_fpixels('1i') # Found here: https://stackoverflow.com/questions/42961810/detect-dpi-scaling-factor-in-python-tkinter-application
        screen_width_in_inches = root.winfo_screenwidth() / dpi
        screen_height_in_inches = root.winfo_screenheight() / dpi

        self.assertTrue(screen_width_in_inches > 10)
        self.assertTrue(screen_width_in_inches < 30)
        self.assertTrue(screen_height_in_inches > 10)
        self.assertTrue(screen_height_in_inches < 30)


    def testFirstButton(self):
        root = tk.Tk()
        ttk.Button(root, text="Hello World").grid()

    @unittest.skip("import has been modified at the top")
    def testNorthSouthEastWestUndefined(self):
        root = tk.Tk()
        mainframe = ttk.Frame(root, padding="3 3 12 12")

        with self.assertRaises(NameError):
            mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

    def testNorthSouthEastWestImported(self):
        root = tk.Tk()
        mainframe = ttk.Frame(root, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

    def testNorthSouthEastWestShortcuts(self):
        root = tk.Tk()
        # There are shortcuts NW, SE, etc... We can give them together or separately
        mainframe = ttk.Frame(root, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        mainframe.grid(column=0, row=0, sticky=(NW, SE))
        mainframe.grid(column=0, row=0, sticky=(N,W, SE))

    @unittest.expectedFailure
    def testStickiness(self):
        # I don't understand why this fails often
        root = tk.Tk()
        root.title = "Test resize"
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        root.geometry("750x250+400+300")
        mainframe = ttk.Frame(root, padding="0 0 0 0")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

        geometry = root.geometry()
        match = re.search(r"(\d+)x(\d+)\D(\d+)\D(\d+)", geometry)
        if match is not None:
            groups = match.groups()
            self.assertEqual(int(groups[0]), mainframe.winfo_width())
            self.assertEqual(int(groups[1]), mainframe.winfo_height())


if __name__ == '__main__':
    unittest.main()
