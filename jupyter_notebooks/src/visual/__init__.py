#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 24/01/2020
# Last Modified Date: 27/01/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
# -*- coding:utf-8 -*-
from scipy import stats
from collections import OrderedDict
import matplotlib as mpl


def set_mpl_style(name):
    paper = {
        "font.weight": "normal",
        "font.size": 17,
        "axes.titleweight": "normal",
        "axes.labelweight": "normal",
        "figure.titleweight": "normal",
        "axes.grid": False,
        #"grid.color": "lightgrey",
        #"grid.alpha": 0.7,
        "axes.axisbelow": True,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 2
    }
    if name == 'paper':
        mpl.rcParams.update(paper)


def calStats(grouped_val, groups, test='t-test', return_z_score=False, **kwargs):
    func_map = {
        't-test': stats.ttest_ind,
        'wilcoxon': stats.wilcoxon
    }

    result = OrderedDict()
    for i, first in enumerate(groups):
        for second in groups[i+1:]:
            if return_z_score:
                result[(first, second)] = func_map[test](grouped_val.get_group(first),
                                                         grouped_val.get_group(
                                                             second),
                                                         **kwargs)
            else:
                result[(first, second)] = func_map[test](grouped_val.get_group(first),
                                                         grouped_val.get_group(
                                                             second),
                                                         **kwargs
                                                         )[1]
    return result


def fancy_scientific(x):
    tmp = '{:.2e}'.format(x).replace('e', ' x 10^').replace(
        '^+', '^').replace('^0', '^').replace('^-0', '^-').replace('x 10^', '\\times 10^{')+'}'
    return tmp.replace('\\times 10^{0}', '')


# colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
# # Sort colors by hue, saturation, value and name.
# by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
#                 for name, color in colors.items())
# COLORS = [name for hsv, name in by_hsv]
COLORS = [
    'black', 'gray', 'silver', 'rosybrown', 'firebrick', 'red', 'darksalmon', 'sienna', 'sandybrown', 'bisque',
    'tan', 'gold', 'darkkhaki', 'olivedrab', 'chartreuse', 'darkgreen', 'seagreen', 'mediumspringgreen', 'lightseagreen', 'paleturquoise',
    'darkcyan', 'deepskyblue', 'slategray', 'royalblue', 'navy', 'mediumpurple', 'darkorchid', 'plum', 'm', 'palevioletred',
    'indianred', 'khaki', 'olive'
]
