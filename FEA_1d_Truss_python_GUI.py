# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 18:09:11 2021

@author: Abhijeet Marekar
"""

from tkinter import *
import tkinter
from tkinter import messagebox
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tkinter as tk
from tkinter import ttk
from matplotlib import style
import numpy as np
import math
def main():    
    '''' Initializing the root window '''
    root = Tk()
    root.geometry("+580+160")
    root.title(" FEA_1D_Truss_Structure (Abhijeet_Marekar) ")
    root.iconbitmap('truss.ico')
    #messagebox.showinfo("Information", "Do Not Close Graph Window !")
    
    '''' Defining Global Variables '''
    global node_cord, ele_conn, ele_x_cord, ele_y_cord, plot_elements
    global DEFINE_BCS, SHOW_FORCE, ele_nod, ele_const, ele_sin, ele_cos
    global ele_len, ele_mat, zero_disp_index, displacement_vector
    global bc_nodes_dict, force_node_dict, rm_node_no, strain, stress
    node_cord = []
    ele_conn = []
    ele_x_cord = []
    ele_y_cord = []
    #x_list = []
    #y_list = []
    plot_elements = False
    DEFINE_BCS = False
    SHOW_FORCE = False
    ele_nod = []
    ele_const= []
    ele_sin = []
    ele_cos = []
    ele_len = []
    ele_mat = []
    zero_disp_index = []  
    displacement_vector = []
    bc_nodes_dict = {}
    force_node_dict = {}
    rm_node_no = 0
    
    #num_of_ele = []
    i = 1
    '''' Base Figure '''
    global ax1
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.grid()
    plt.xlabel("x-axis")
    fig = plt.gcf()
    fig.canvas.set_window_title('Graphical Window !!! DO NOT CLOSE IT !!!')
    '''' Defining functions '''
    def cal_angle(point1, point2):
        x0= point1[0]
        y0= point1[1]
        
        x1=point2[0]
        y1=point2[1]
        
        if x0>=0 and x1>=0 and y0>=0 and y1>=0:
            
            if x0>= x1 and y0>=y1:
                x0,x1 = x1,x0
                y0,y1 = y1,y0
            angle= math.atan2((y1-y0),(x1-x0))
            
            #print("First Quadrant")
            #print(f"The angle with +x axis is {math.degrees(angle)}")
            return angle
        if x0<=0 and x1<=0 and y0>=0 and y1>=0:
            
            if x0<= x1 and y0>=y1:
                x0,x1 = x1,x0
                y0,y1 = y1,y0
            angle= math.atan2((y1-y0),(x1-x0))
            
            #print("Second Quadrant")
            #print(f"The angle with +x axis is {math.degrees(angle)}")
            return angle
        
        if x0<=0 and x1<=0 and y0<=0 and y1<=0:
            
            if x0<= x1 and y0<=y1:
                x0,x1 = x1,x0
                y0,y1 = y1,y0
            
            angle= math.atan2((y1-y0),(x1-x0))
            
            #print("Third Quadrant")
            #print(f"The angle with +x axis is {math.degrees(angle)+360}")
            return angle + (2*math.pi)
        if x0>=0 and x1>=0 and y0<=0 and y1<=0:
            
            if x0>= x1 and y0<=y1:
                x0,x1 = x1,x0
                y0,y1 = y1,y0
            
            angle= math.atan2((y1-y0),(x1-x0))
            # print(angle)
            
            #print("Fourth Quadrant")
            #print(f"The angle with +x axis is {math.degrees(angle+(2*math.pi))}")
            return angle + (2*math.pi)
    
    def add_area():
        global A
        A = float(area_entry.get())
        area_entry.delete(0,END)
        if A == 0:
            temp = messagebox.showerror("Invalid Area", f" Area should have non-zero value!")
            if temp.lower() == 'ok':
                reset()
        area_value_label = Label(frame2, text= f"A = {A} (mm^2) ").grid(row= 1,column = 5)
        area_button = Button(frame2, text= "Define", command = add_area, state = DISABLED).grid(row= 1,column = 4)
        E_button = Button(frame2, text= "Define", command = add_E).grid(row= 2,column = 4)
        
    def add_E():
        global E
        E = float(E_entry.get())
        E_entry.delete(0,END)
        if E == 0:
            temp = messagebox.showerror("Invalid modulus", f" modulus of Elasticity should have non-zero value!")
            if temp.lower() == 'ok':
                reset()
        area_value_label = Label(frame2, text= f"E = {E:.2e} (N/mm^2) ").grid(row= 2,column = 5)
        E_button = Button(frame2, text= "Define", command = E, state = DISABLED).grid(row= 2,column = 4)
        add_node_button = Button(frame2, text= "Add Node", command = add_node).grid(row=5, column= 3)
        
    def add_node():
        global x_list, y_list
        x = float(node_x_entry.get())
        y = float(node_y_entry.get())
        node_cord.append((x,y))
        node_x_entry.delete(0,END)
        node_y_entry.delete(0,END)
        # for x,y in node_cord:
        #     x_list.append(x)
        #     y_list.append(y)
        
        ani.event_source.start()
        del_node_button = Button(frame2, text = "Del Node", command = del_node).grid(row = 5, column = 4)
        if len(node_cord) <= 1:
            add_ele_button = Button(frame2, text= "Add Element", command = add_element, state = DISABLED).grid(row=9, column= 3)
        else:
            add_ele_button = Button(frame2, text= "Add Element", command = add_element).grid(row=9, column= 3)
    def finish_geometry_definition():
        
        global global_mat,displacement_vector,force_vector
        add_node_button = Button(frame2, text= "Add Node", command = add_node, state = DISABLED).grid(row=5, column= 3)
        
        number_of_nodes = len(node_cord)
        global_mat = np.zeros((number_of_nodes*2, number_of_nodes*2))
        ii = 1
        for n in range(len(node_cord)):
            
            displacement_vector.append('x'+str(i))
            displacement_vector.append('y'+str(i))
            ii+=1
        force_vector = np.zeros(((len(node_cord))*2,1))
        
        add_ele_button = Button(frame2, text= "Add Element", command = add_element, state = DISABLED).grid(row=9, column= 3)
        del_ele_button = Button(frame2, text = "Del Element", command = del_element, state = DISABLED).grid(row = 9, column = 4)
        add_node_button = Button(frame2, text= "Add Node", command = add_node, state = DISABLED).grid(row=5, column= 3)
        del_node_button = Button(frame2, text = "Del Node", command = del_node, state = DISABLED).grid(row = 5, column = 4)
        fix_bc_button = Button(frame2, text = "Fix", command = fix_bc).grid(row = 13, column = 3)
        fx0_bc_button = Button(frame2, text = "Fx = 0", command = fx_0_bc).grid(row = 13, column = 4)
        fy0_bc_button = Button(frame2, text = "Fy = 0", command = fy_0_bc).grid(row = 13, column = 5)
        force_add_buttn = Button(frame2, text = "Add", command = add_force).grid(row = 16, column = 4)
        
    def add_element():
         global plot_elements 
         # global i
         # i += 1
         # if i > num_of_ele:
         #     add_ele_button = Button(frame2, text= "Add Element", command = add_element, state= DISABLED).grid(row=9, column= 3)
         # if i <= num_of_ele:    
         #     ele_entry = Label(frame2, text = f"Nodes for Element {i} ").grid(row=9,column=0)
         
         node1 = int(node_1.get())
         node2 = int(node_2.get())
         ele_conn.append((node1,node2))
         node_1.delete(0,END)
         node_2.delete(0,END)
         
         plot_elements = True
         ani.event_source.start()
         del_ele_button = Button(frame2, text = "Del Element", command = del_element).grid(row = 9, column = 4)
         define_geo_button = Button(frame2, text = "Define Geometry", command = finish_geometry_definition).grid(row = 10, column = 2)
         
        
    def add_num_ele():
        
        global num_of_ele, node_1, node_2
        num_of_ele = int(num_ele_entry.get())
        #num_ele_entry = Entry(frame2, width= 12, borderwidth = 2, state = DISABLED).grid(row = 3, column= 2)
        add_num_ele_button = Button(frame2, text= "Confirm", command = add_num_ele, state = DISABLED).grid(row=7, column= 3)
        num_ele_entry.delete(0, END)
        
    
        
    def element_calculation(node1,node2):
        
        l = math.sqrt((node_cord[node2-1][0]-node_cord[node1-1][0])**2+(node_cord[node2-1][1]-node_cord[node1-1][1])**2)
        constant = A*E /l
        #print(node_cord[node1-1], node_cord[node2-1])
        ele_angle = cal_angle(node_cord[node1-1], node_cord[node2-1])
        # cos = (node_cord[node2-1][0]-node_cord[node1-1][0])/l
        # sin = (node_cord[node2-1][1]-node_cord[node1-1][1])/l
        #print(f"element angle with +x axis is {ele_angle}")
        cos = (math.cos(ele_angle))
        sin = (math.sin(ele_angle))
        #print(f"sin: {sin},  cos: {cos}")
        e_mat_temp = constant*np.array([[cos**2, cos*sin, -cos**2, -cos*sin],
                                        [cos*sin, sin**2, -cos*sin, -sin**2],
                                        [-cos**2, -cos*sin, cos**2, cos*sin],
                                        [-cos*sin, -sin**2, cos*sin, sin**2]])
        #print(e_mat_temp)
        ele_mat.append(e_mat_temp)
        ele_nod.append((node1,node2))
        ele_const.append(constant)
        ele_sin.append(sin)
        ele_cos.append(cos)
        ele_len.append(l)
        global_mat_index = [node1*2-2, node1*2-1,node2*2-2, node2*2-1]
        #print(global_mat_index)
    
        ele_mat_index = [0,1,2,3]
        k,o=0,0
        for i in global_mat_index:
    
            o=0
            for j in global_mat_index:
                #print(i,j,k,o)
                global_mat[i][j] = global_mat[i][j]+ e_mat_temp[k][o]
                o = o+1
            k = k + 1
    '''    
    def animate(i):
        x_list = []
        y_list = []
        for x,y in node_cord:
            x_list.append(x)
            y_list.append(y)
            
        a.clear()
        a.scatter(x_list, y_list)
        a.plot([1,2,3],[4,5,6])
        
        #f.annotation('Try',(5,5))
        
    '''
    def delete_connected_elements():
        rm_node_no = int(len(node_cord) + 1)
        if len(ele_conn) == 0:
            return
        # index = 0
        # index_to_remove = []
        # print("Before : ",ele_conn)
        for n in list(ele_conn):
            
            print('Remove Node No: ',rm_node_no)
            print(f"{n}")
            if n[0] == rm_node_no:
                index_to_remove.append(index)
                ele_conn.remove(n)
                index += 1
                continue
                
            if n[1] == rm_node_no:
                ele_conn.remove(n)
                index_to_remove.append(index)
            # index += 1
        # print("After : ",ele_conn)
    def animate(i):
        # x_list = [] 
        # y_list = []
        global node_cord, ax1, ele_conn
        # for x,y in node_cord:
        #     x_list.append(x)
        #     y_list.append(y)
     
        ax1.clear()
        for x,y in node_cord:
            ax1.scatter(x, y, c = 'blue')
        k= 1
        for i,j in node_cord:
            plt.text(i+0.15,j+0.15,f'{k}')
            k +=1
            
        if plot_elements:
            
            for nod1, nod2 in ele_conn:
                x = []
                y = []
                x.append(node_cord[nod1-1][0])
                y.append(node_cord[nod1-1][1])
                x.append(node_cord[nod2-1][0])
                y.append(node_cord[nod2-1][1])
                # print(x)
                # print(y)
                ax1.plot(x,y,'c')
        if DEFINE_BCS:
            
            for bc_node in bc_nodes_dict.keys():
                plt.text((node_cord[bc_node-1][0])-0.15,(node_cord[bc_node-1][1])-1, bc_nodes_dict[bc_node], c = 'r')
        if SHOW_FORCE:
            
            for f_node in force_node_dict.keys():
                plt.text((node_cord[f_node-1][0])-0.15,(node_cord[f_node-1][1])-1, force_node_dict[f_node] + 'N', c = 'g')
                
        ani.event_source.stop()
        ax1.grid()
        if len(node_cord) > 0:
            if (np.min(node_cord)) and (np.max(node_cord)) != 0:
                ax1.set_xlim(np.min(node_cord)-5, np.max(node_cord)+5)
    
    def fix_bc():
        global displacement_vector,bc_nodes_dict,DEFINE_BCS
        node_num = int(bc_node_entry.get())
        displacement_vector[node_num*2-2] = 0
        displacement_vector[node_num*2-1] = 0
        zero_disp_index.append(node_num*2-2)
        zero_disp_index.append(node_num*2-1)
        bc_nodes_dict[node_num] = "(1,1)"
        bc_node_entry.delete(0,END)
        DEFINE_BCS = True
        ani.event_source.start()
        
    def fx_0_bc():
        global displacement_vector,bc_nodes_dict,DEFINE_BCS
        node_num = int(bc_node_entry.get())
        displacement_vector[node_num*2-2] = 0
        zero_disp_index.append(node_num*2-2)
        bc_nodes_dict[node_num] = "(1,0)"
        bc_node_entry.delete(0,END)
        DEFINE_BCS = True
        ani.event_source.start()
        
    def fy_0_bc():
        global displacement_vector,bc_nodes_dict,DEFINE_BCS
        node_num = int(bc_node_entry.get())
        displacement_vector[node_num*2-1] = 0
        zero_disp_index.append(node_num*2-1)
        bc_nodes_dict[node_num] = "(0,1)"
        bc_node_entry.delete(0,END)
        DEFINE_BCS = True
        ani.event_source.start()
        
    
    def del_element():
        global ele_conn
        if ele_conn is not None:
            ele_conn.pop(-1)
        if len(ele_conn) == 0:
            define_geo_button = Button(frame2, text = "Define Geometry", command = finish_geometry_definition, state = DISABLED).grid(row = 10, column = 2)
        ani.event_source.start()
    
    def del_node():
        global node_cord, ele_conn, plot_elements
        if node_cord is not None:
            node_cord.pop(-1)
            #x_list.pop(-1)
            #y_list.pop(-1)
        if len(node_cord) <= 1:
            add_ele_button = Button(frame2, text= "Add Element", command = add_element, state = DISABLED).grid(row=9, column= 3)
            define_geo_button = Button(frame2, text = "Define Geometry", command = finish_geometry_definition, state = DISABLED).grid(row = 10, column = 2)
            plot_elements = False
        else:
            add_ele_button = Button(frame2, text= "Add Element", command = add_element).grid(row=9, column= 3)
        global rm_node_no
        delete_connected_elements()
            
        ani.event_source.start()
    
    def add_force():
        global force_vector,SHOW_FORCE
        node_num = int(force_node_entry.get())
        fx = float(fx_entry.get())
        fy = float(fy_entry.get())
        force_vector[node_num*2-2][0] = fx
        force_vector[node_num*2-1][0] = fy
        force_node_dict [node_num] = f"({fx},{fy})"
        force_node_entry.delete(0,END)
        fx_entry.delete(0,END)
        fy_entry.delete(0,END)
        SHOW_FORCE = True
        ani.event_source.start()
        solve_button = Button(frame2, text = "Solve", command = solve, bg = 'green', width = 10).grid(row = 18, column = 4) 
    
    def cal_stress_strain(initial_cord, displacements):
        
        displaced_nodal_cord = []
        changed_ele_len = []
        for n in range(len(initial_cord)):
            new_cord  = (initial_cord[n][0]+displacements[2*n], initial_cord[n][1]+displacements[2*n+1])
            displaced_nodal_cord.append(new_cord)
        for node1,node2 in ele_conn:
            l = math.sqrt((displaced_nodal_cord[node2-1][0]-displaced_nodal_cord[node1-1][0])**2+
                          (displaced_nodal_cord[node2-1][1]-displaced_nodal_cord[node1-1][1])**2)
            changed_ele_len.append(l)
        ele_strain = (np.array(changed_ele_len) - np.array(ele_len)) / np.array(ele_len)  
        ele_stress = ele_strain * E
        
        return ele_strain, ele_stress
            
    def solve():
        global ele_conn,reduced_GSM,reduced_force_vector,displacement_vector
        global inv_rGSM,disp_solutions, final_force_vector,strain, stress
        for node1 , node2 in ele_conn:
            element_calculation(node1, node2)
        reduced_GSM = global_mat          
        reduced_force_vector = force_vector     
        reduced_GSM = np.delete(reduced_GSM,zero_disp_index,0)
        reduced_GSM = np.delete(reduced_GSM,zero_disp_index,1)
        reduced_force_vector = np.delete(reduced_force_vector, zero_disp_index,0)    
        try:
            inv_rGSM = np.linalg.inv(reduced_GSM)
        except np.linalg.LinAlgError:
            note = "\n\n(If have single bar element then please ensure\nthat length of bar is along y-axis.)"
            temp = messagebox.showerror("Sigular Matrix", f"Error occured while inverting matrix !\nPlease check that BC's are well define{note}")
            
            if temp.lower() == 'ok':
                reset()
        disp_solutions = np.matmul(inv_rGSM,reduced_force_vector)
        ind = 0
        for ele in range(len(displacement_vector)):
            
            if displacement_vector[ele] != 0:
                displacement_vector[ele] = disp_solutions[ind]
                ind +=1
            if displacement_vector[ele] == 0:
                displacement_vector[ele] = [0.0]
                
            
        displacement_vector = (np.array(displacement_vector))
        final_force_vector = np.matmul(global_mat,displacement_vector)
        
        #Resluts_label = Label(result_window, text = "Results").grid(row = 1, column = 1)
        for x in displacement_vector:
            print(f"{x[0]:.3e}")
        strain, stress = cal_stress_strain(node_cord,displacement_vector)
        show_results()
        show_results_button = Button(frame2, text = "Show Results", command = show_results).grid(row = 18, column = 5)
        solve_button = Button(frame2, text = "Solve", command = solve, state = DISABLED).grid(row = 18, column = 4)   
        
    def show_results():
        result_window = Toplevel()
        result_window.title(" Results ")
        empty_label0 = Label(result_window, text = " ").grid(row = 0, column = 1)
        Results_label = Label(result_window, text = "Nodal Results", bd = 1, relief = SUNKEN, bg = 'white'). grid(row = 1, column = 0, columnspan = 8, sticky = 'W'+'E')
        empty_label = Label(result_window, text = " ").grid(row = 2, column = 1)
        Column_1_label = Label(result_window, text = "Node", width = 12).grid(row = 4, column = 0)
        Column_10_label = Label(result_window, text = "|").grid(row = 4, column = 1)
        Column_2_label = Label(result_window, text = "Component", width = 12).grid(row = 4, column = 2)
        Column_20_label = Label(result_window, text = "|").grid(row = 4, column = 3)
        Column_3_label = Label(result_window, text = "Force (N)", width = 12).grid(row = 4, column = 4)
        Column_30_label = Label(result_window, text = "|").grid(row = 4, column = 5)
        Column_4_label = Label(result_window, text = "Disp (mm)", width = 12).grid(row = 4, column = 6)
        Column_40_label = Label(result_window, text = "|").grid(row = 4, column = 7)
        node_num = 1
        for item in range(int(len(displacement_vector))):
            if item%2 ==0:
                node_nub_label = Label(result_window, text = f"{node_num}", relief = SUNKEN).grid(row = 8+item, column = 0, rowspan = 2 , sticky = 'W'+'E'+'S'+'N')
                comp_x = Label(result_window, text= " x ",bg = 'lightgray').grid(row = 8+item, column = 2,sticky = 'W'+'E')
                force = Label(result_window, text = f"{final_force_vector[item,0]:.3e}",bg = 'lightgray').grid(row = 8+item, column = 4, sticky = 'W'+'E')
                disp = Label(result_window, text = f"{displacement_vector[item,0]:.3e}", bg = 'lightgray').grid(row = 8+item, column = 6, sticky = 'W'+'E')
                Column_100_label = Label(result_window, text = "|",bg = 'lightgray').grid(row = 8+item, column = 1, sticky = 'W'+'E')
                Column_200_label = Label(result_window, text = "|",bg = 'lightgray').grid(row = 8+item, column = 3, sticky = 'W'+'E')
                Column_300_label = Label(result_window, text = "|",bg = 'lightgray').grid(row = 8+item, column = 5, sticky = 'W'+'E')
                Column_400_label = Label(result_window, text = "|",bg = 'lightgray').grid(row = 8+item, column = 7, sticky = 'W'+'E')
            else:
                comp_y = Label(result_window, text= " y ").grid(row = 8+item, column = 2)
                force = Label(result_window, text = f"{final_force_vector[item,0]:.3e}").grid(row = 8+item, column = 4)
                disp = Label(result_window, text = f"{displacement_vector[item,0]:.3e}").grid(row = 8+item, column = 6)
                Column_100_label = Label(result_window, text = "|").grid(row = 8+item, column = 1)
                Column_200_label = Label(result_window, text = "|").grid(row = 8+item, column = 3)
                Column_300_label = Label(result_window, text = "|").grid(row = 8+item, column = 5)
                Column_400_label = Label(result_window, text = "|").grid(row = 8+item, column = 7)
            
            if item%2 ==0:
                node_num +=1
        pos = len(displacement_vector) + 8 + 2
        ele_empty_label0 = Label(result_window, text = " ").grid(row = pos, column = 1)
        ele_Results_label = Label(result_window, text = "Element Results", bd = 1, relief = SUNKEN, bg = 'white'). grid(row = pos +1, column = 0, columnspan = 8, sticky = 'W'+'E')
        ele_empty_label = Label(result_window, text = " ").grid(row = pos +2, column = 1)
        ele_Column_1_label = Label(result_window, text = "Element", width = 12).grid(row = pos +4, column = 0)
        ele_Column_10_label = Label(result_window, text = "|").grid(row = pos +4, column = 1)
        ele_Column_2_label = Label(result_window, text = "Nodes", width = 12).grid(row = pos +4, column = 2)
        ele_Column_20_label = Label(result_window, text = "|").grid(row = pos +4, column = 3)
        ele_Column_3_label = Label(result_window, text = "Stress (N/mm^2)", width = 12).grid(row = pos +4, column = 4)
        ele_Column_30_label = Label(result_window, text = "|").grid(row = pos +4, column = 5)
        ele_Column_4_label = Label(result_window, text = "Strain", width = 12).grid(row = pos +4, column = 6)
        ele_Column_40_label = Label(result_window, text = "|").grid(row = pos +4, column = 7)
        
        for ele in range(len(ele_conn)):
            node_nub_label = Label(result_window, text = f"{ele +1}", relief = SUNKEN).grid(row = pos+5+ele, column = 0, rowspan = 1 , sticky = 'W'+'E'+'S'+'N')
            if ele%2 == 0:
                node_conn_label = Label(result_window, text= f"{ele_conn[ele]}",bg = 'lightgray').grid(row = pos+5+ele, column = 2,sticky = 'W'+'E')
                stress_label = Label(result_window, text= f"{stress[ele,]:.3e}",bg = 'lightgray').grid(row = pos+5+ele, column = 4,sticky = 'W'+'E')
                strain_label = Label(result_window, text = f"{strain[ele,]:.3e}",bg = 'lightgray').grid(row = pos+5+ele, column = 6, sticky = 'W'+'E')
            else:
                node_conn_label = Label(result_window, text= f"{ele_conn[ele]}").grid(row = pos+5+ele, column = 2,sticky = 'W'+'E')
                stress_label = Label(result_window, text= f" {stress[ele,]:.3e} ").grid(row = pos+5+ele, column = 4,sticky = 'W'+'E')
                strain_label = Label(result_window, text = f"{strain[ele,]:.3e}").grid(row = pos+5+ele, column = 6, sticky = 'W'+'E')
    
            
    def reset():
        plt.close()
        root.destroy()
        plt.close()
        main()
    '''' Definition of frames'''
    # frame1= LabelFrame(root, text="Frame1",width = 200, height = 200, padx=50,pady=50)
    # frame1.grid(row=0,column=0)
    
    frame2= LabelFrame(root, text="Input Data",width = 200, height = 200,padx=50,pady=50)
    frame2.grid(row=0,column=2)
    
    '''' Embeding of Grpah into tkinter'''
    
    # graph= FigureCanvasTkAgg(f, master = frame1)
    # graph.draw()
    # graph.get_tk_widget().grid(row= 2, column = 1)
    
    '''' Taking Element Constants '''
    
    area_label = Label(frame2, text= "Enter cross-section area of element (mm^2) ").grid(row= 1,column = 0, columnspan = 2)
    area_entry = Entry(frame2, width = 8, borderwidth = 2)
    area_entry.grid(row= 1, column= 3)
    #while True:
    #    try:
    area_button = Button(frame2, text= "Define", command = add_area).grid(row= 1,column = 4)
    
    #        break
    #    except ValueError:
    #        messagebox.showerror('Invalid input','Entered Number is not a valid value.\nEnter vaild interger/float.')
    E_label = Label(frame2, text= "Enter modulus of Elasticity (N/mm^2) ").grid(row= 2,column = 0, columnspan = 2)
    E_entry = Entry(frame2, width = 8, borderwidth = 2)
    E_entry.grid(row= 2, column= 3)
    E_button = Button(frame2, text= "Define", command = add_E, state = DISABLED).grid(row= 2,column = 4)
    
    
    '''' Taking nodal Cordinates'''
    
    
    node_entry = Label(frame2, text = "Node Cordinate ").grid(row=5,column=0)
    x = Label(frame2, text= "x (mm)").grid(row=4, column = 1)
    y = Label(frame2, text= "y (mm)").grid(row=4, column = 2)
    node_x_entry = Entry(frame2, width= 12, borderwidth = 2)
    node_y_entry = Entry(frame2, width= 12, borderwidth = 2)
    add_node_button = Button(frame2, text= "Add Node", command = add_node, state = DISABLED).grid(row=5, column= 3)
    # done_node = Button(frame2, text= "Finish", command = finsh_node_addition, state = DISABLED).grid(row=5, column= 5)
    del_node_button = Button(frame2, text = "Del Node", command = del_node, state = DISABLED).grid(row = 5, column = 4)
    node_x_entry.grid(row = 5, column= 1)
    node_y_entry.grid(row = 5, column= 2)
    
    
    '''' Taking Element Connectivity '''
    
    
    #num_ele_label = Label(frame2, text = "Enter Number of Elements ").grid(row=7,column=0)
    #num_ele_entry = Entry(frame2, width= 12, borderwidth = 2)
    #num_ele_entry.grid(row = 7, column= 2)
    #add_num_ele_button = Button(frame2, text= " Confirm ", command = add_num_ele, state = DISABLED).grid(row=7, column= 3)
    
    # add_ele_button = Button(frame2, text= "Add Element").grid(row=9, column= 3)
    x = Label(frame2, text= "Start Node").grid(row=8, column = 1)
    y = Label(frame2, text= "End Node").grid(row=8, column = 2)
    
    ele_entry = Label(frame2, text = "Define Element ").grid(row=9,column=0)
    
    node_1 = Entry(frame2, width= 12, borderwidth = 2)
    node_1.grid(row = 9, column= 1)
    node_2 = Entry(frame2, width= 12, borderwidth = 2)
    node_2.grid(row = 9, column= 2)
    add_ele_button = Button(frame2, text= "Add Element", command = add_element, state = DISABLED).grid(row=9, column= 3)
    del_ele_button = Button(frame2, text = "Del Element", command = del_element, state = DISABLED).grid(row = 9, column = 4)
    
    ''' Finalize the Geometry '''
    
    define_geo_button = Button(frame2, text = "Define Geometry", command = finish_geometry_definition, state = DISABLED).grid(row = 10, column = 2)
    
    
    ''' Displacement Boundary Conditions  '''
    
    
    bc_entry_label = Label(frame2, text= "Displacement Constrain for Node : ").grid(row = 13, column = 0)
    bc_node_lable = Label(frame2, text = "Node Number").grid(row = 12, column = 1)
    bc_node_entry = Entry(frame2, width = 10, borderwidth = 2)
    bc_node_entry.grid(row= 13, column = 1)
    fix_bc_button = Button(frame2, text = "Fix", command = fix_bc, state = DISABLED).grid(row = 13, column = 3)
    fx0_bc_button = Button(frame2, text = "Fx = 0", command = fx_0_bc, state = DISABLED).grid(row = 13, column = 4)
    fy0_bc_button = Button(frame2, text = "Fy = 0", command = fy_0_bc, state = DISABLED).grid(row = 13, column = 5)
    
    ''' Force Definition '''
    
    force_entry_label = Label(frame2, text = "Enter Force value on Node (N) :").grid(row = 16, column = 0)
    force_node_label = Label(frame2, text = "Node Number").grid(row = 15, column = 1)
    force_node_entry = Entry(frame2, width = 10, borderwidth = 2)
    force_node_entry.grid(row = 16, column = 1)
    fx_entry_label = Label(frame2, text = "Fx").grid(row = 15, column = 2)
    fy_entry_label = Label(frame2, text = "Fy").grid(row = 15, column = 3)
    fx_entry = Entry(frame2, width = 10, borderwidth = 2)
    fx_entry.insert(0, '0')
    fx_entry.grid(row = 16, column=2)
    fy_entry = Entry(frame2, width = 10, borderwidth = 2)
    fy_entry.insert(0, '0')
    fy_entry.grid(row = 16, column=3)
    
    force_add_buttn = Button(frame2, text = "Add", command = add_force, state = DISABLED).grid(row = 16, column = 4)
    
    
    ''' Solve and Show Results'''
    solve_button = Button(frame2, text = "Solve", command = solve, state = DISABLED, bg = 'red', width = 10).grid(row = 18, column = 4) 
    ani = animation.FuncAnimation(fig, animate)
    
    
    ''' Show Results '''
    show_results_button = Button(frame2, text = "Show Results", command = show_results, state = DISABLED).grid(row = 18, column = 5)
    ''' Reset_Button '''
    
    reset_button = Button(frame2, text = "Reset", width = 10, command = reset).grid(row = 20, column = 4)
    
    plt.show()
    
    root.mainloop()
main()