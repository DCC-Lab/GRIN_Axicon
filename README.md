# GRIN_Axicon
CODE V macro to model a GRIN-axicon and python scripts to fit the GRIN-axicon parameters to an existing axicon's parameters. Make sure to read the intership overview (```IntershipRport.pdf```)

# **GRIN-axicon configurator setup**
**Intallation process:**
1. Install the ```GUI_GRIN-axicon.py``` file
2. Run it with any python IDLE (like VS Code)

**Configuring a GRIN-axicon**
1. When running the ```GUI_GRIN-axicon.py``` file, an interface should pop up, start by entering the axicon parameters in the Axicon configuration pannel.
2. Play with the GRIN-axicon sliders to match the input beam diameter to the output beam diameter (if you don't necessarily want an exact copy of the axicon, you can play around until you find your desired GRIN-axcion setup)

# **Macro setup**
**Intallation process:**
1. Install the ```GRIN-axicon_macro.seq``` file
2. Copy the sequence file into the working directory
3. For ease of use, when first opening CODE V and, using the Macro Manager, assign the macro to a desired menu
4. Install the modified version of ```usergrn.seq```
5. In the local repository where the built-in CODE V macro are saved, delete the original usergrn.seq and replace it by the modified one (in my case, C:\CODEV112_SR1\macro).

**Creating a GRIN**
1. Open CODE V
2. Run the GRIN-axicon_macro.seq either by a menu shortcut, the Macro Manager or by typing ```in GRIN-axicon_macro```
3. Type the desired grin glass name
2. Select the desired mode (theo vs exp:  when using the experimental mode, make sure that the file is in the working directory)
4. Run macro with the desired parameters

Also, when using the experimental mode, make sure to format the refractive index profile correctly using this Python progam (if you are not in the DCC Lab team, copy the code to your own Google Collab file and change the file names):
[Google Collab]([https://github.com](https://colab.research.google.com/drive/1K1_1UN_SLIXMgnj2FgyKjVjrnuUiGcWT?usp=sharing)https://colab.research.google.com/drive/1K1_1UN_SLIXMgnj2FgyKjVjrnuUiGcWT?usp=sharing)


