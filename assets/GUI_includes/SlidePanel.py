import customtkinter as ctk
# Slide Panel ----------------------------------------------

class SlidePanel(ctk.CTkFrame):
    def __init__(self, parent, start_pos, end_pos):
        super().__init__(master=parent)
        # general attributes
        self.start_pos = start_pos + 0.05
        self.end_pos = end_pos - 0.01
        self.width = abs(start_pos - end_pos)

        # animation logic
        self.pos = self.start_pos
        self.in_start_pos = True

        # layout
        self.place(relx=self.start_pos, rely=0.05, relwidth=0.3, relheight=0.9)

    def animate(self):
        if self.in_start_pos:
            self.animate_forward()
        else:
            self.animate_backwards()

    def animate_forward(self):
        if self.pos > self.end_pos:
            self.pos -= 0.008
            self.place(relx=self.pos, rely=0.05, relwidth=0.3, relheight=0.9)
            self.lift()  # Cette ligne ramène le panel au premier plan

            self.after(10, self.animate_forward)
        else:
            self.in_start_pos = False

    def animate_backwards(self):
        if self.pos < self.start_pos:
            self.pos += 0.008
            self.place(relx=self.pos, rely=0.05, relwidth=0.3, relheight=0.9)
            self.lift()  # Cette ligne ramène le panel au premier plan

            self.after(10, self.animate_backwards)
        else:
            self.in_start_pos = True

    def set_color(self, new_color):
        self.configure(fg_color=new_color)

    def set_border_color(self, new_border_color):
        self.configure(border_color=new_border_color)