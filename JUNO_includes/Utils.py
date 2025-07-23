import turtle

"""
This file contains some purely utility functions that would not make up for a full class/file 
"""

def trueValue(datas, param, hiera):
        """
        Easy way to access pdg and article values (see article [2]) to compare sensitivity results
        """ 
        data = datas[hiera]
        if param == "m21":
             trueValue = data.m_21
             articleValues = {100:0.074e-5, 6*365.25:0.024e-5, 20*365.25:0.017e-5}
             ep_pdg = 0.18e-5
        elif param == "s12":
             trueValue = data.s_12
             articleValues = {100:0.0058, 6*365.25:0.0016, 20*365.25:0.0010}
             ep_pdg = 0.013
        elif param == "s13":
             trueValue = data.s_13
             articleValues = {100:0.010, 6*365.25:0.0026, 20*365.25:0.0016}
             ep_pdg = 0.07e-2      
        elif param=="m31_NH" and hiera==0:
             trueValue = 2.5303e-3
             articleValues = {100:0.021e-3, 6*365.25:0.0047e-3, 20*365.25:0.0029e-3}
             ep_pdg = 2.8058e-5
        elif param=="m13_IH" and hiera==1:
             trueValue = 2.4537e-3 
             articleValues = None   # no article values
             ep_pdg = 2.8944e-5
        else:
             print(f"Parameter {param} incompatible with hierarchy {hiera}")
        
        return trueValue,articleValues,ep_pdg


def b_intensifies(GUIturtle=None):
    """
    Little artifact we left as a joke showed to Mr. João Pedro Athaydes Marcondes de André
    It is only accessible from the options of the GUI (or CTRL+B)
    It is purely aestethic and doesn't affect the other programs and results.

    code taken from      https://github.com/Chando0185/brazil_argentina_flag/tree/main
    """
    if GUIturtle == None:
        obj=turtle.Turtle()
    else:
         obj = GUIturtle
    obj.speed(3)

    if GUIturtle == None:
        win=turtle.Screen()
        win.bgcolor("white")

    obj.pencolor("seagreen")

    obj.penup()
    obj.goto(-100,200)
    obj.pendown()
    obj.begin_fill()
    obj.fillcolor("seagreen")
    obj.setheading(0)
    obj.forward(300)
    obj.setheading(270)
    obj.forward(200)
    obj.setheading(180)
    obj.forward(300)
    obj.setheading(90)
    obj.forward(200)
    obj.end_fill()

    obj.setheading(270)
    obj.forward(180)
    obj.setheading(360)

    obj.penup()
    obj.forward(150)
    obj.pendown()

    obj.begin_fill()
    obj.fillcolor("yellow")
    obj.setheading(142)
    obj.forward(120)
    obj.setheading(40)
    obj.forward(120)

    obj.setheading(-40)
    obj.forward(120)
    obj.setheading(-141)
    obj.forward(117)
    obj.end_fill()


    obj.penup()
    obj.setheading(141)
    obj.forward(70)


    obj.setheading(90)
    obj.forward(30)
    obj.pendown()

    obj.begin_fill()
    obj.fillcolor("#002776")
    obj.circle(-50)
    obj.end_fill()
    obj.penup()
    obj.setheading(345)
    obj.forward(95)
    obj.pendown()
    obj.setheading(120)
    obj.pencolor("white")
    obj.pensize(10)
    obj.circle(70,85)
    obj.penup()
    obj.pensize(8)
    obj.forward(105)
    obj.pendown()

    obj.begin_fill()
    obj.fillcolor("seagreen")

    obj.setheading(90)
    obj.pencolor("seagreen")
    obj.forward(200)
    obj.setheading(180)
    obj.forward(10)
    obj.setheading(270)
    obj.forward(500)


    obj.setheading(180)
    obj.forward(80)
    obj.setheading(270)
    obj.forward(20)
    obj.setheading(360)
    obj.forward(180)

    obj.setheading(90)
    obj.forward(20)
    obj.setheading(180)
    obj.forward(90)
    obj.setheading(90)
    obj.forward(500)

    obj.end_fill()
    obj.hideturtle()
    
    if GUIturtle == None:
        turtle.done()
    turtle.TurtleScreen._RUNNING = True


def B_mode(GUIturtle=None):
     """
     Little artifact we left as a joke showed to Mr. João Pedro Athaydes Marcondes de André
     It is only accessible from the options of the GUI (or CTRL+B)
     It is purely aestethic and doesn't affect the other programs and results.
     """
     try:
          b_intensifies(GUIturtle=GUIturtle)
     finally:
          turtle.TurtleScreen._RUNNING = True

def str2list(res):      
     """
     Utility function used for the saving of graphs from the GUI
     """      
     k=0
     L=[]
     el = ""
     while res[k] not in [")","&"]:
          if res[k]=="(":
               k+=1
               continue
          if res[k] == ",":
               try:
                    L.append(float(el))
                    el=""
               except:
                    L.append(el)
                    el=""
          else:
               el+=res[k]
          k+=1
     return L

def str2nbPoints(nbstr,other=None,default=20,min=3):
     """
     Utility function used in the GUI to handle user entries for choosing the number of points for an analysis
     """
     if nbstr == "":
          nbpoints = default if other==None  else other
     else:
          nbpoints = max(int(float(nbstr)),min)
     return nbpoints

def convert_path_linux(path):
    """
    Utility function used in the GUI for compatibility with Linux
    """
    for i in range(len(path)):
        if path[i:i+1]=="\\":
            path = path[:i] +"/"+ path[i+1:]
    return(path)

def incompability_params(param,hiera):
     """
     Warning for incompatible parameters used in the main
     """
     if ((param=="m31_NH" and hiera==1) or (param=="m13_IH" and hiera==0)):
          print(f"\n /!\\ Parameter {param} incompatible with hierarchy {hiera} !\n")
          return True