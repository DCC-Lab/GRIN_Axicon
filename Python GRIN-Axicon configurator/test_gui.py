import tkinter as tk
from tkinter import ttk
import tkinter.font as font
import unittest
import re 

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


if __name__ == '__main__':
    unittest.main()
