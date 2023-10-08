#  @auther: Brian Kyanjo
#  @date: August 8th, 2023
#  @description: a GUI with terminal access to run the geoflood model by 
#                collecting user input and feeding it into a python script 
#                and later runs and dispalys the output in the terminal of the GUI

# import necessary libraries
import tkinter as tk
import tkinter.scrolledtext as st
import subprocess
import threading
import sys

# GUI layout
# create a class for the GUI; include interactive elements such as buttons, input fields, and drop-down menus
#  for user interaction, divide the GUI into sections for script configuration, script execution, script output,
#  and terminal-like output
class GeofloodScriptGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Geoflood Script GUI with Terminal")

        # Create GUI elements
        self.input_label = tk.Label(root, text="Enter Input:")
        self.input_entry = tk.Entry(root)
        self.run_button = tk.Button(root, text="Run Script", command=self.run_script)
        self.output_text = st.ScrolledText(root, height=10, width=40)
        
        # Place GUI elements using grid layout
        self.input_label.grid(row=0, column=0, sticky=tk.W)
        self.input_entry.grid(row=0, column=1, padx=5, pady=5)
        self.run_button.grid(row=0, column=2, padx=5, pady=5)
        self.output_text.grid(row=1, column=0, columnspan=3, padx=5, pady=5)

        # Initialize terminal-like interface
        self.terminal_text = st.ScrolledText(root, height=10, width=40, state=tk.DISABLED)
        self.terminal_text.grid(row=2, column=0, columnspan=3, padx=5, pady=5)

    def run_script(self):
        input_data = self.input_entry.get()
        self.output_text.delete("1.0", tk.END)
        
        # Start a separate thread to run the script and capture terminal-like output
        thread = threading.Thread(target=self.run_script_thread, args=(input_data,))
        thread.start()

    def run_script_thread(self, input_data):
        # Redirect terminal-like output to the GUI
        self.redirect_output_to_gui()
        
        try:
            # Run the Python script using subprocess
            script_output = subprocess.check_output(["python", "setrun.py", input_data], stderr=subprocess.STDOUT)
            self.output_text.insert(tk.END, script_output.decode())
        except subprocess.CalledProcessError as e:
            self.output_text.insert(tk.END, "Error: " + e.output.decode())

        # Restore normal stdout and stderr
        self.restore_stdout_stderr()

    def redirect_output_to_gui(self):
        sys.stdout = self.TerminalWriter(self.terminal_text)

    def restore_stdout_stderr(self):
        sys.stdout = sys.__stdout__

    class TerminalWriter:
        def __init__(self, terminal_text):
            self.terminal_text = terminal_text

        def write(self, text):
            self.terminal_text.configure(state=tk.NORMAL)
            self.terminal_text.insert(tk.END, text)
            self.terminal_text.configure(state=tk.DISABLED)
            self.terminal_text.see(tk.END)

        def flush(self):
            pass

# main function
if __name__ == "__main__":
    root = tk.Tk()
    app = GeofloodScriptGUI(root)
    root.mainloop()



