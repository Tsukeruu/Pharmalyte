"""
The format for naming css elements in textual (my way)
i use "customname_widgetname", for example i have a collapsible that is an about me for my program
i do "about_collapsible"
"""
from textual import on
import math
from textual.app import App, ComposeResult
from textual.widgets import Footer, Label, Header, Collapsible, Input, Button, Tabs, Tab
from textual.css.query import NoMatches
from textual.containers import Container, Horizontal, Vertical
from textual.screen import Screen, ModalScreen
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors, Draw 
import pubchempy
from pubchempy import PubChemHTTPError

OPTION_LIST = [
    "Molecular mass calculator",
    "PH Calculation",
    "About"
]

class PHscreen(ModalScreen[None]):
    BINDINGS = [('escape','app.pop_screen','escape')]
    def compose(self) -> ComposeResult:
        self.PHContainer = Container(
            Label("""This simple mode will take in concentrations and returns the PH aswell as explains the buffer capacity, this mode is essentially useful for pharmacists as pharmacists usually almost always deal with PH, buffers and solutions, aswell as the fact that pharmacists use buffers to make solutions. 
            """,id="PHLabel_container_PHscreen"),
            Input(placeholder="Is the solution an acid, base or a buffer?",type="text",id="PH_Intro"),
            id="PHscreen_container"
        )
        yield self.PHContainer
    
    def return_buffer_capacity(self):
        acid_conc = self.firstvalue
        base_conc = self.secondvalue
        if acid_conc <= 0 or base_conc <= 0:
            return "No buffer present."
        ratio = acid_conc / base_conc
        if 0.5 <= ratio <= 2:
            return "Strong buffer capacity. Solution resists pH changes well."
        else:
            return "Weak buffer capacity. Solution resists pH changes, but less effectively"

    @on(Input.Submitted)
    def accept_input(self,event: Input.Submitted):
        PH_INPUT = self.query_one("#PH_Intro")
        input_o = event.input
        value = input_o.value 
        input_id = input_o.id
        try:
            oldphcontainerlabel = self.query("#PHCONTAINER_LABEL").first()
        except NoMatches as e:
            oldphcontainerlabel = None
        if input_id == "PH_Intro":
            if value.lower().strip() == "acid":
                input_o.display = False
                self.PHContainer.mount(Input(placeholder="Enter acid concentration (M)",id="phcontainer_input_acid"))
                self.PHContainer.mount(Label("",id="PHCONTAINER_LABEL"))
            elif value.lower().strip() == "base":
                input_o.display = False
                self.PHContainer.mount(Input(placeholder="Enter base concentration (M)",id="phcontainer_input_base"))
                self.PHContainer.mount(Label("",id="PHCONTAINER_LABEL"))
            elif value.lower().strip() == "buffer":
                input_o.display = False
                self.PHContainer.mount(Input(placeholder="Enter weak acid concentration (M)",id="phcontainer_input_buffer1"))
                self.PHContainer.mount(Input(placeholder="Enter conjugate base concentration (M)",id="phcontainer_input_buffer2"))
                self.PHContainer.mount(Label("",id="PHCONTAINER_LABEL"))
            else:
                app.pop_screen()
        elif input_id == "phcontainer_input_acid":
            if int(value) <= 0:
                app.pop_screen()
            else:
                self.query_one("#PHCONTAINER_LABEL").update(f" Ph of strong acid: {round(-math.log10(int(value)))}")
        elif input_id == "phcontainer_input_base":
            if int(value) <= 0:
                app.pop_screen()
            else:
                calculated = 1e-14 / int(value)
                self.query_one("#PHCONTAINER_LABEL").update(f" PH of strong base: {round(-math.log10(calculated))}")

        elif input_id == "phcontainer_input_buffer1":
            self.firstvalue = int(value)
        elif input_id == "phcontainer_input_buffer2":
            self.secondvalue = int(value)
            try:
                if self.firstvalue > 0 and self.secondvalue > 0:
                    self.pka = 4.76

                    calculated_ph = round(self.pka+math.log10(self.secondvalue/self.firstvalue))
                    self.query_one("#PHCONTAINER_LABEL").update(f"""
 PH of the buffer is: {calculated_ph}
 Buffer capacity: {self.return_buffer_capacity()}
""")
            except (AttributeError, ValueError) as e:
                app.pop_screen()
    def on_mount(self):
        pass

class SmilesScreen(ModalScreen[None]):
    BINDINGS = [("escape","app.pop_screen","escape")]
    def compose(self) -> ComposeResult:
        self.smilesContainer = Container(
            Label("""
To get started, simply input molecules in SMILES format, note that SMILES Is a shorthand way to write
molecules in text format, smiles is known across multiple pharmacytical programs including many famous ones! an example would be CC0, which codes for ethanol, type it down to test it out!
            """,id="smilesabout_container_label1"),
            Input(placeholder="Enter smiles: ",id="smilesabout_container_input1"),
            Label("",id="smilesplaceholder_container_label2"),
            id="smilesmain_container1"
        )
        yield self.smilesContainer
            
    @on(Input.Submitted)
    def accept_molecules(self):
        smilesinput = self.query_one("#smilesabout_container_input1")
        smiles = smilesinput.value
        molecularsmiles = Chem.MolFromSmiles(smiles, sanitize=True)
        try:
            pcpcompounds = pubchempy.get_compounds(smiles, 'smiles')
            pcpcompounds = pcpcompounds[0]
            Draw.MolToFile(molecularsmiles, f"smile_images/{smiles}.png")
            commonname = pcpcompounds.synonyms[0]
        except (PubChemHTTPError) as e:
            app.pop_screen()
            return
        except (IndexError, TypeError) as c:
            commonname = "NO COMMON NAME HAS BEEN FOUND!" #setting common name is because of an index error where the synonym at index 0 does not exist
        if molecularsmiles is None:
            app.pop_screen()
            return
        else:
            molecularweight = Descriptors.MolWt(molecularsmiles)
            numofatoms = molecularsmiles.GetNumAtoms()
            numofheavyatoms = Descriptors.HeavyAtomCount(molecularsmiles)
            numofrings = Descriptors.RingCount(molecularsmiles)
            numofvalence = Descriptors.NumValenceElectrons(molecularsmiles) 
            placeholder = self.query_one("#smilesplaceholder_container_label2")
            placeholder.update(f"""
SMILES: {smiles}
MOLECULAR WEIGHT: {molecularweight}
COMMON NAME: {commonname}
NUMBER OF VALENCE ELECTRONS: {numofvalence}
NUMBER OF ATOMS: {numofatoms}
NUMBER OF HEAVY ATOMS: {numofheavyatoms}
NUMBER OF RINGS: {numofrings}
            """)

    def on_mount(self):
        RDLogger.DisableLog('rdApp.*')


class pharmalyte_Main(App):
    CSS_PATH = "main.tcss"
    SCREENS = {"smile": SmilesScreen,"ph": PHscreen}

    def compose(self) -> ComposeResult:
        yield Tabs(
            Tab(OPTION_LIST[0],id="smiles"),
            Tab(OPTION_LIST[1],id="PH"),
            Tab(OPTION_LIST[2],id="ABOUT"),
        )
        yield Collapsible(
            Label("""A simple program that looks forward to replicating and simulating the jobs of pharmacists in regards of chemistry
this program will allow you to convert SMILES (molecules) into real molecular data aswell as calculate ph dosage.
            """,
            id="about_label"),
            title="About this program",
            id="about_collapsible"
        )
        yield Container( 
            Label("""
Smiles is a way to represent molecules as text, an example of this would be CC0, which represents ethanol
smiles are the building blocks of many pharmacytical programs made by worldwide programmers!
Clicking the button here will transport you into a new screen to input your smiles conversion.
            """,id="smiles_container_label1"),
            Button("Convert smiles to molecular data!",id="smiles_container_button1"),
            id="smiles_container1"   
        )
        yield Container(
            Label("""
This system calculates PH and explains the buffer capacity, which is quite strong for pharmacy because buffers are core in pharmacy!
            """,id="PH_container_label1"),
            Button("Calculate PH",id="PH_container_button1"),
            id="PH_container1"
        )
        yield Container(
            Label("""
This project was built by Laith aljundi, Adam zuraik and omar arafah, all three worked together to provide you
with this wonderful program that pharmacists use all together for specific molecules and medicines.
Special credits to ameeer al shobaki for the ideas.
            """,id="ABOUT_label1"),
            id="ABOUT_container1"
        )
        

    def on_button_pressed(self, event: Button.Pressed):
        pressed_button = event.button
        if pressed_button.id == "smiles_container_button1":
            self.push_screen(SmilesScreen())
        elif pressed_button.id == "PH_container_button1":
            self.push_screen(PHscreen())
        else:
            pass
   
    def on_tabs_tab_activated(self, event: Tabs.TabActivated) -> None:
        for content in ["PH_container1","smiles_container1","ABOUT_container1"]:
            self.query_one(f"#{content}").display = False

        self.query_one(f"#{event.tab.id}_container1").display = True

    def on_mount(self) -> None:
        self.theme = "catppuccin-mocha"
        smiles_horizontal_container1 = self.query_one("#smiles_container1")
        smiles_horizontal_container1.border_title = "Moles Calculator"
        PH_horizontal_container1 = self.query_one("#PH_container1")
        PH_horizontal_container1.border_title = "PH Dosage"
        about_container1 = self.query_one("#ABOUT_container1")
        about_container1.border_title = "ABOUT & CREDITS"


app = pharmalyte_Main()
app.run()

    
