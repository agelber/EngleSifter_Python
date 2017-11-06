import Tkinter as tk
import TF
import tkMessageBox
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import FileDialog
import DE


'''
expression = pd.read_csv("BSq_Ordered.csv")
expression = expression.set_index("Unnamed: 0")
expression.index.name = None
expression = expression.astype(float)
exp_ind = list(expression.index)
expression.to_pickle("BSq_Ordered.pkl")
print (expression.head())
'''




expression = pd.read_pickle("BSq_Ordered.pkl")
exp_ind = list(expression.index)
print (expression.head())


def accept():
    text = gene_ent.get()

    if text == "exit":
        win.quit()
        return
    else:
        smps = get_samples()

        try:
            plt.close()
            plt.figure(figsize=(12.6, 7))
            plt.bar(range(len(smps)), expression[text].loc[smps], align='center', alpha=0.5)
            plt.ylabel("TPM")
            plt.xticks(range(len(smps)), rename(smps), rotation=70, fontsize=8, ha='right')
            plt.subplots_adjust(left=0.07, bottom=0.33, right=0.96, top=0.99, wspace=0, hspace=0)
            plt.show()
        except KeyError:
                tkMessageBox.showerror("Error", "Could not find gene: %s" % text)


def get_samples():
    ret = []
    for j in range(len(samples)):
        if samples[j].get() == 1:
            ret.append(exp_ind[j])
    return ret


def select(g, x):
    if g == 'on':
        y = 1
    if g == 'off':
        y = 0
    if x == 'all':
        for k in range(len(samples)):
            samples[k].set(y)
    if x == 'SC':
        for k in range(0, 13):
            samples[k].set(y)

    if x == 'cn12':
        for k in range(13, 25):
            samples[k].set(y)

    if x == 'cn7':
        for k in range(25, 33):
            samples[k].set(y)

    if x == 'cn6':
        for k in range(33, 39):
            samples[k].set(y)

    if x == 'cn5':
        for k in range(39,48):
            samples[k].set(y)

    if x == 'cn4':
        for k in range(48, 58):
            samples[k].set(y)

    if x == 'cn3':
        for k in range(58, 68):
            samples[k].set(y)


def change_frame(f1, f2):
    f1.grid_remove()
    f2.grid(row=0, column=0)


def rename(x):
    meep = []
    for i in range(len(x)):
        meep.append(x[i].replace("_", " "))
    return meep


win = tk.Tk()

win.title("EngleSifter")

# Frames
menu = tk.Frame(win)
viz = tk.Frame(win)
de = tk.Frame(win)
de_ent = tk.Frame(de)
viz_ent = tk.Frame(viz)
viz_smps = tk.Frame(viz)

# Menu Formatting
font=("Helvetica", 23),

menu_title = tk.Label(menu, text="Welcome to EngleSifter\n", font=("Helvetica", 30))
fb1 = tk.Button(menu, text="Visualization", font=("Helvetica", 23), command=lambda: change_frame(menu, viz))
fb2 = tk.Button(menu, text="Differential Expression", font=("Helvetica", 23), command=lambda: change_frame(menu, de))
spacer = tk.Label(menu, text='\n')

menu_title.grid(row=0, column=1)
# menu_title.grid_columnconfigure(5, minsize=6000)
fb1.grid(row=3, column=0)
fb2.grid(row=3, column=2)
# fb1.grid_rowconfigure(4, minsize=400000)

spacer.grid(row=4)

# Viz Formatting

viz_title = tk.Label(viz_ent, text="Choose samples and enter a gene: ")
gne = tk.StringVar()
gene_ent = tk.Entry(viz_ent, bd=5, textvariable=gne)
go = tk.Button(viz_ent, text="Go!", command=accept)
sel_all = tk.Button(viz_ent, text="Select All", command=lambda: select('on', 'all'))
unsell_all = tk.Button(viz_ent, text="Unselect All", command=lambda: select('off', 'all'))
all_sc = tk.Button(viz_smps, text="All SC", command=lambda: select('on', 'SC'))
all_cn12 = tk.Button(viz_smps, text="All CN12", command=lambda: select('on', 'cn12'))
all_cn7 = tk.Button(viz_smps, text="All CN7", command=lambda: select('on', 'cn7'))
all_cn6 = tk.Button(viz_smps, text="All CN6", command=lambda: select('on', 'cn6'))
all_cn5 = tk.Button(viz_smps, text="All CN5", command=lambda: select('on', 'cn5'))
all_cn4 = tk.Button(viz_smps, text="All CN4", command=lambda: select('on', 'cn4'))
all_cn3 = tk.Button(viz_smps, text="All CN3", command=lambda: select('on', 'cn3'))

selection_list = [all_sc, all_cn12, all_cn7, all_cn6, all_cn5, all_cn4, all_cn3]

back1 = tk.Button(viz_ent, text="Back to menu", command=lambda: change_frame(viz, menu))

viz_title.grid(row=0, column=0, sticky=tk.W)
gene_ent.grid(row=0, column=1, sticky=tk.W)
go.grid(row=0, column=4, sticky=tk.W)
back1.grid(row=0, column=5, sticky=tk.E)
sel_all.grid(row=1, column=0, sticky=tk.W)
unsell_all.grid(row=1, column=0, sticky=tk.E)

samples = []
sc = TF.ToggledFrame(viz_smps, text='Spinal Cord', relief="raised", borderwidth=1)
c12 = TF.ToggledFrame(viz_smps, text='CN12          ', relief="raised", borderwidth=1)
c7 = TF.ToggledFrame(viz_smps, text='CN7            ', relief="raised", borderwidth=1)
c6 = TF.ToggledFrame(viz_smps, text='CN6           ', relief="raised", borderwidth=1)
c5 = TF.ToggledFrame(viz_smps, text='CN5            ', relief="raised", borderwidth=1)
c4 = TF.ToggledFrame(viz_smps, text='CN4            ', relief="raised", borderwidth=1)
c3 = TF.ToggledFrame(viz_smps, text='CN3            ', relief="raised", borderwidth=1)


sample_names = ["Spinal Cord", "CN12", "CN7", "CN6", "CN5", "CN4", "CN3"]

sampleFrames = [sc, c12, c7, c6, c5, c4, c3]

for i in range(len(expression.index)):
    samples.append(tk.IntVar())

sampleButtons = []
for i in range(len(expression.index)):
    if i in range(0, 13):
        s1 = sc
    if i in range(13, 25):
        s1 = c12
    if i in range(25, 33):
        s1 = c7
    if i in range(33, 39):
        s1 = c6
    if i in range(39, 48):
        s1 = c5
    if i in range(48, 58):
        s1 = c4
    if i in range(58, 68):
        s1 = c3
    sampleButtons.append(
        tk.Checkbutton(s1.sub_frame, text=list(expression.index)[i].replace("_", " "), variable=samples[i], onvalue=1,
                       offvalue=0))

for i in range(len(expression.index)):
    sampleButtons[i].grid(row=(i / 3) + 5, column=i % 3, sticky=tk.W)


for i in range(7):
    sampleFrames[i].grid(row=(i % 4)+1, column=int(2*(i/4)), sticky=tk.W)
    selection_list[i].grid(row=(i % 4)+1, column=int(2*(i/4)+1), sticky=tk.W)

viz_ent.grid(row=0, sticky=tk.W)
viz_smps.grid(row=1, sticky=tk.W)




# DE formatting

cond1 = tk.Frame(de)
cond2 = tk.Frame(de)

de_title = tk.Label(de, text="Differential Expression", font=("Helvetica", 28))
cond_one = tk.Label(cond1, text="Condition 1", font=("Helvetica", 23))
cond_two = tk.Label(cond2, text="Condition 2", font=("Helvetica", 23))

cond1_vars = [tk.IntVar() for i in range(7)]
cond2_vars = [tk.IntVar() for i in range(7)]

cond1_buttons = [tk.Checkbutton(cond1, text=sample_names[i], variable=cond1_vars[i], onvalue=1,
                                offvalue=0) for i in range(7)]
cond2_buttons = [tk.Checkbutton(cond2, text=sample_names[i], variable=cond2_vars[i], onvalue=1,
                                offvalue=0) for i in range(7)]

for i in range(7):
    cond1_buttons[i].grid(row=i+3, column=0, sticky=tk.W)
    cond2_buttons[i].grid(row=i+3, column=1, sticky=tk.W)


cond_one.grid(row=0, column=0)
cond_two.grid(row=0, column=1, sticky=tk.W)
back2 = tk.Button(de, text="Back to menu", command=lambda: change_frame(de, menu))
test1 = tk.Button(de, text="DE Test 1", command=lambda: change_frame(de, menu))
test2 = tk.Button(de, text="DE Test 2", command=lambda: change_frame(de, menu))
test3 = tk.Button(de, text="DE Test 3", command=lambda: change_frame(de, menu))

de_ent.grid(row=0, column=0)
de_title.grid(row=0,column=1)
test1.grid(row=3, column=0, sticky=tk.E)
test2.grid(row=3, column=1)
test3.grid(row=3, column=2, sticky=tk.W)
back2.grid(row=4, column=3)
cond1.grid(row=1, column=0)
cond2.grid(row=1, column=3, sticky=tk.E)


# main loop


menu.grid(row=0, column=0)

if __name__ == '__main__':
    win.mainloop()
