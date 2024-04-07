# Python program to illustrate a stop watch 
# using Tkinter 
#importing the required libraries 
import tkinter as Tkinter 
from datetime import datetime
counter = 30*60
running = False
slide_times = [0.5, 0.5, 2, 2.5, 1.5, 1.5, 2.5, 1, 2, 1, 2, 1, 1.5, 0.5]
slide_seconds = [int(x*60) for x in slide_times]
current_slide = 0

def counter_label(label, actual_slide, time_left): 
    def count(): 
        if running: 
            global counter
            global slide_seconds
            global current_slide

            # To manage the initial delay. 
            if counter==30*60:			 
                display="Starting..."
            elif counter <= 0:
                display="Time's up!"
            else:
                tt = datetime.fromtimestamp(counter)
                string = tt.strftime("%M:%S")
                display=string 

            label['text']=display # Or label.config(text=display)
            if 30*60 > counter > 10*60:
                label['fg'] = 'black'
            elif 5*60 < counter <= 10*60:
                label['fg'] = 'orange'
            elif counter <= 5*60:
                label['fg'] = 'red'
                
            if slide_seconds[current_slide] > 0:
                actual_slide['text'] = 'Actual slide: ' + str(current_slide + 1)
                ss = datetime.fromtimestamp(slide_seconds[current_slide])
                string_ss = ss.strftime("%M:%S")
                time_left['text'] = 'Time left: ' + string_ss
                if slide_seconds[current_slide] <= 20:
                    time_left['fg'] = 'green'
            elif current_slide < len(slide_seconds)-1:
                current_slide += 1
                actual_slide['text'] = 'Actual slide: ' + str(current_slide + 1)
                ss = datetime.fromtimestamp(slide_seconds[current_slide])
                string_ss = ss.strftime("%M:%S")
                time_left['text'] = 'Time left: ' + string_ss
                time_left['fg'] = 'black'
            else:
                actual_slide['text'] = 'Questions'
                time_left['text'] = 'Finished!'
                time_left['fg'] = 'red'

            # label.after(arg1, arg2) delays by 
            # first argument given in milliseconds 
            # and then calls the function given as second argument. 
            # Generally like here we need to call the 
            # function in which it is present repeatedly. 
            # Delays by 1000ms=1 seconds and call count again. 
            label.after(1000, count) 
            counter -= 1
            slide_seconds[current_slide] -= 1

    # Triggering the start of the counter. 
    count()	 

# start function of the stopwatch 
def Start(label, actual_slide, time_left): 
    global running 
    running=True
    counter_label(label, actual_slide, time_left) 
    start['state']='disabled'
    stop['state']='normal'
    reset['state']='normal'

# Stop function of the stopwatch 
def Stop(): 
    global running 
    start['state']='normal'
    stop['state']='disabled'
    reset['state']='normal'
    running = False

# Reset function of the stopwatch 
def Reset(label, actual_slide, time_left): 
    global counter
    global current_slide
    global slide_seconds
    global slide_times
    counter=30*60
    current_slide = 0
    slide_seconds = [int(x*60) for x in slide_times]

    # If rest is pressed after pressing stop. 
    if running==False:	 
        reset['state']='disabled'
        label['text']='Welcome!'
        label['fg'] = 'black'
        actual_slide['text'] = 'Actual slide: 1'
        time_left['text'] = 'Time left: 00:00'
        time_left['fg'] = 'black'

    # If reset is pressed while the stopwatch is running. 
    else:			 
        label['text']='Starting...'
        label['fg'] = 'black'
        actual_slide['text'] = 'Actual slide: 1'
        time_left['text'] = 'Time left: 00:00'
        time_left['fg'] = 'black'

root = Tkinter.Tk() 
root.title("Stopwatch") 

# Fixing the window size. 
root.minsize(width=250, height=210) 
label = Tkinter.Label(root, text="Welcome!", fg="black", font="Verdana 30 bold") 
label.pack() 
actual_slide = Tkinter.Label(root, text="Actual slide: 1", fg="black", font="Verdana 20 bold")
actual_slide.pack()
time_left = Tkinter.Label(root, text="Time left: 00:00", fg="black", font="Verdana 20 bold")
time_left.pack()
f = Tkinter.Frame(root)
start = Tkinter.Button(f, text='Start', width=6, command=lambda:Start(label, actual_slide, time_left)) 
stop = Tkinter.Button(f, text='Stop',width=6,state='disabled', command=Stop) 
reset = Tkinter.Button(f, text='Reset',width=6, state='disabled', command=lambda:Reset(label, actual_slide, time_left)) 
f.pack(anchor = 'center',pady=5)
start.pack(side="left") 
stop.pack(side ="left") 
reset.pack(side="left") 
root.mainloop()
