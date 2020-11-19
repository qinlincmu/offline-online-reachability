"""

simulator for AA project

author:Qin Lin (qinlin@andrew.cmu.edu)

"""
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import matplotlib.animation as animation
import math
import numpy as np
from random import randint
import matplotlib.patches as patches


class State:
    """
    vehicle state class
    """

    def __init__(self, x=0, x_dot=0.0, y=0.0, y_dot=0.0, vx=0.0, vy=0.0, yaw=0.0, yaw_dot=0.0, v=0.0):
        self.x = x
        self.x_dot = x_dot
        self.y = y
        self.y_dot = y_dot
        self.yaw = yaw
        self.yaw_dot = yaw_dot
        self.v = v
        self.vx = vx
        self.vy = vy
def update_state_slip(state, v, delta):

    L_r = 1.364
    L = 2.728
    DT = 0.025
    beta = math.atan(math.tan(delta)*L_r/L)

    
    state.x_dot = v * math.cos(beta+state.yaw) 
    state.x = state.x + state.x_dot* DT
    state.y_dot = v * math.sin(beta+state.yaw)
    state.y = state.y + state.y_dot* DT


    
    state.vx = float(v*math.cos(beta))
    state.vy = float(v*math.sin(beta))

    state.yaw_dot = v*math.tan(delta)*math.cos(beta) / L
    state.yaw = state.yaw + state.yaw_dot * DT

    return state


def do_simulation(state, delta):
    x,y,yaw,vx,vy,yawdot=[],[],[],[],[],[]
    for i in range(len(delta)):
        state = update_state_slip(state, 10, delta[i])
        x.append(state.x)
        y.append(state.y)
        yaw.append(state.yaw)
        vx.append(state.vx)
        vy.append(state.vy)
        yawdot.append(state.yaw_dot)
    return x,y,yaw,vx,vy,yawdot

def main():
    print(__file__ + " start!!")
    f = open("/home/cfu1/Documents/FFAST/catkin_ws/src/reach_flow/src/bicycle_slip/lookup.txt", "r")
    data = []
    for line in f.readlines():
        temp = line.strip().split(" ")
        temp1 = [float(item) for item in temp]
        data.append(temp1)
    num = len(data)
    r = randint(0, num-1)
    delta = [data[r][2] for i in range(3)]
    initial_state = State(x=0.0, x_dot=0.0, y=0.0, y_dot=0.0, vx=0.0, vy=0.0, yaw=(data[r][1]+data[r][0])/2, yaw_dot=0.0, v=0.0)
    x,y,yaw,vx,vy,yawdot = do_simulation(initial_state, delta)
    xg = data[r][3:9]
    yg = data[r][9:15]
    yawg = data[r][15:21]
    vxg = data[r][21:27]
    vyg = data[r][27:33]
    yawdotg = data[r][33:39]
    #print("x", x, xg)
    #print("y", y, yg)
    #print("yaw", yaw, yawg)
    print("vx", vx, vxg)
    print("vy", vy, vyg)
    print("yawdot", yawdot, yawdotg)
    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(x[i],y[i], "r*")
        ax.add_patch(
        patches.Rectangle(
            (xg[i*2], yg[i*2]),
            xg[i*2+1]-xg[i*2],
            yg[i*2+1]-yg[i*2],
            edgecolor = 'blue',
            facecolor = 'red',
            fill=False
        ))
    ax.set_xlim(min(xg)-1, max(xg)+1)
    ax.set_ylim(min(yg)-1, max(yg)+1)
    plt.show()

if __name__ == '__main__':
    main()
