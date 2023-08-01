import tkinter as tk
from tkinter import ttk
import tkinter.font as font
import unittest

class TestTk(unittest.TestCase):
    def testInit(self):
        self.assertTrue(True)

    def testLoadTk(self):
        tk_root = tk.Tk()
        self.assertIsNotNone(tk_root)
        


if __name__ == '__main__':
    unittest.main()
