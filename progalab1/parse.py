import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

L = 2

data = np.loadtxt('out.txt')
dt = data[1, 0]
t = data[:, 0]
x1 = data[:, 1]
y1 = data[:, 2]
x2 = data[:, 3]
y2 = data[:, 4]

l1 = (x1**2+y1**2)**(1/2)
l2 = ( (x2-x1)**2 + (y2-y1) ** 2 ) ** (1/2)

fig1, ax1 = plt.subplots()
ax1.plot(t, x1, label='x1')
ax1.plot(t, y1, label='y1')
ax1.plot(t, l1, label='l1(t)')
ax1.set_title('x1(t), y1(t)')
ax1.grid()
ax1.legend(loc='upper right')

# fig4, ax4 = plt.subplots()
# ax4.plot(x1, y1, label='y1(x1)')
# ax4.plot(x2, y2, label='y2(x2)')

fig2, ax2 = plt.subplots()
ax2.plot(t, x2, label='x2')
ax2.plot(t, y2, label='y2')
ax2.plot(t, l2, label='l2(t)')
ax2.set_title('x2(t), y2(t)')
ax2.grid()
ax2.legend(loc='upper right')

# fig3, ax3 = plt.subplots()
# # t1 = [el % np.pi for el in t1]
# # t2 = [el % np.pi for el in t2]
# ax3.plot(y[:, 0], y[:, 2])

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, 1.))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    history_x = x2[:i]
    history_y = y2[:i]

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i * dt))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(y1), interval=1, blit=True)
plt.show()
