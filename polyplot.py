import numpy as np
import matplotlib.pyplot as plt


def polyplot(vertices, savepath=None):
	vertices = np.array(vertices)
	x = vertices[:, 0]
	y = vertices[:, 1]

	plt.figure()

	# Fill polygon
	plt.fill(x, y, alpha=0.4)

	# Draw polygon boundary
	plt.plot(np.append(x, x[0]), np.append(y, y[0]), marker='o')

	# Force the origin to be visible
	xmin = min(x.min(), 0)
	xmax = max(x.max(), 0)
	ymin = min(y.min(), 0)
	ymax = max(y.max(), 0)

	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)

	plt.axis('equal')
	plt.grid(True)

	# Draw axes
	plt.axhline(0, linewidth=3, color='black', zorder=0)
	plt.axvline(0, linewidth=3, color='black', zorder=0)

	if savepath is not None:
		plt.savefig(savepath)
	else:
		plt.show()