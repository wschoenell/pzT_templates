import matplotlib.pyplot as plt

figure = plt.figure(1)
figure.clf()

ax_hist_upp_left = figure.add_axes([0.25, 0.75, 0.325, 0.15])
ax_hist_upp_right = figure.add_axes([0.575, 0.75, 0.325, 0.15])
ax_hist_low = figure.add_axes([0.1, 0.1, 0.15, 0.65])
c = plt.scatter(range(10), range(10), c=range(10))
ax_left_low = figure.add_axes([0.25, 0.1, 0.325, 0.65], sharey=figure.get_axes()[-1])
figure.get_axes()[-1].yaxis.set_visible(False)
ax_right_low = figure.add_axes([0.575, 0.1, 0.325, 0.65], sharey=figure.get_axes()[-1])
figure.get_axes()[-1].yaxis.set_visible(False)
c = plt.scatter(range(10), range(10), c=range(10))

