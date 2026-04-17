import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("output/curve.csv")

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(df['x_present'], df['x_delayed'], lw = .5, alpha = .5)
ax.set_xlabel('x(t)')
ax.set_ylabel('x(t-tau)')

# We save the figure instead of showing it
plt.savefig("output/traj_plot.png", dpi=300)

with open('output/poincare.txt') as f:
    for line in f:
        line = line.strip().strip('{}') 
        values = [float(v) for v in line.split(',')]
        plt.scatter(values[0], values[-1], color='purple', linewidths = .5) 

# We save the figure instead of showing it
plt.savefig("output/poincare_traj_plot.png", dpi=300)


df = pd.read_csv("output/newton.csv")

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(df['x_present'], df['x_delayed'], lw = .5, alpha = .5)
ax.set_xlabel('x(t)')
ax.set_ylabel('x(t-tau)')

# We save the figure instead of showing it
plt.savefig("output/newton_traj.png", dpi=300)