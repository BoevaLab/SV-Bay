def BasicStats(arr_name, arr):
    print 'Basic stats for ' + arr_name + ':'
    sorted_arr = sorted(arr)
    print   'Len: ' + str(len(sorted_arr)) + \
            ' min: ' + str(sorted_arr[0]) + \
            ' max: ' + str(sorted_arr[-1]) + \
            ' avg: ' + str(sum(sorted_arr) / len(sorted_arr)) + \
            ' mediana: ' + str(sorted_arr[len(sorted_arr) / 2]) + \
            ' non-zero: ' + str(sum([x != 0 for x in sorted_arr]))+\
            ' non zero regions: '+str(sum([x>0 for x in sorted_arr]))

bp1_region_sizes = []
bp2_region_sizes = []
bp_difs = []
bp_exp_difs = []
decisions_nums = []

f = open('bp_stats.txt')
for line in f:
    (bp1_size, bp2_size, bp_dist, exp_bp_dist, decisions) = line.split()
    bp1_region_sizes.append(int(bp1_size))
    bp2_region_sizes.append(int(bp2_size))
    bp_difs.append(int(bp_dist))
    bp_exp_difs.append(float(exp_bp_dist))
    decisions_nums.append(int(decisions))
f.close()

BasicStats('bp1_region_sizes', bp1_region_sizes)
BasicStats('bp2_region_sizes', bp2_region_sizes)
BasicStats('bp_difs', bp_difs)
BasicStats('bp_exp_difs', bp_exp_difs)
BasicStats('decisions_nums', decisions_nums)