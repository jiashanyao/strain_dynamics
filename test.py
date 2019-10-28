import matplotlib.pyplot as plt
import json

with open('result_matrix_beta=0.5.json') as infile:
     H = json.load(infile)
fig = plt.figure(figsize=(6, 3.2))

ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(H, extent=[10, 30, 10, 3])
ax.set_aspect(2.857)

plt.show()
