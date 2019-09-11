import matplotlib.pyplot as plt
import json

with open('ran_result.json') as infile:
    ran_result = json.load(infile)
with open('re_result.json') as infile:
    re_result = json.load(infile)

color_levels = 2
parameter = 'mu'
para_values = [d[parameter] for d in ran_result]
min_value, max_value = min(para_values), max(para_values)
interval = (max_value - min_value) / color_levels
intervals = []
for i in range(color_levels):
    intervals.append((min_value + i * interval, min_value + (i + 1) * interval))
print(para_values)
print(intervals)

plt.subplot(121)
for i in range(color_levels):
    if i < color_levels - 1:
        ran_div = [d['mean_diversity'] for d in ran_result if intervals[i][0] <= d[parameter] < intervals[i][1]]
        re_div = [d['mean_diversity'] for d in re_result if intervals[i][0] <= d[parameter] < intervals[i][1]]
        plt.scatter(ran_div, re_div,
                    label='{:.3f}'.format(intervals[i][0]) + '<=' + parameter + '<' + '{:.3f}'.format(intervals[i][1]))
    else:
        ran_div = [d['mean_diversity'] for d in ran_result if intervals[i][0] <= d[parameter] <= intervals[i][1]]
        re_div = [d['mean_diversity'] for d in re_result if intervals[i][0] <= d[parameter] <= intervals[i][1]]
        plt.scatter(ran_div, re_div,
                    label='{:.3f}'.format(intervals[i][0]) + '<=' + parameter + '<=' + '{:.3f}'.format(intervals[i][1]))
    plt.legend()
plt.plot([0, 1], [0, 1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('Diversity in random model')
plt.ylabel('Diversity in regular model')
plt.gca().set_aspect('equal', adjustable='box')
plt.subplot(122)
for i in range(color_levels):
    if i < color_levels - 1:
        ran_disc = [d['mean_discordance'] for d in ran_result if intervals[i][0] <= d[parameter] < intervals[i][1]]
        re_disc = [d['mean_discordance'] for d in re_result if intervals[i][0] <= d[parameter] < intervals[i][1]]
        plt.scatter(ran_disc, re_disc,
                    label='{:.3f}'.format(intervals[i][0]) + '<=' + parameter + '<' + '{:.3f}'.format(intervals[i][1]))
    else:
        ran_disc = [d['mean_discordance'] for d in ran_result if intervals[i][0] <= d[parameter] <= intervals[i][1]]
        re_disc = [d['mean_discordance'] for d in re_result if intervals[i][0] <= d[parameter] <= intervals[i][1]]
        plt.scatter(ran_disc, re_disc,
                    label='{:.3f}'.format(intervals[i][0]) + '<=' + parameter + '<=' + '{:.3f}'.format(intervals[i][1]))
    plt.legend()
plt.plot([0.4, 1], [0.4, 1])
plt.xlim(0.4, 1)
plt.ylim(0.4, 1)
plt.xlabel('Discordance in random model')
plt.ylabel('Discordance in regular model')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.show()
