def FindNeighbors2(current_sub_links, all_sub_links, list_centr, chrom, curr_num):
    logger.debug('Eto nastoyawii poslednii ' + all_sub_links[-1].name)
    logger.info('Number of current sub links: ' + str(len(current_sub_links)))
    candidate_sub_links = [sub_link for sub_link in all_sub_links if sub_link.link.gamma_alelles == -1 and sub_link.link.num_elements >= curr_num]
    logger.info('Number of candidate sub links: ' + str(len(candidate_sub_links)))
    ##
    candidate_sub_links_sorted_by_start = sorted(candidate_sub_links, key = sortSubLinksBeg)
    candidate_sub_links_sorted_by_end = sorted(candidate_sub_links, key = sortSubLinksEnd)
    sorted_starts = [sub_link.safe_start for sub_link in candidate_sub_links_sorted_by_start]
    sorted_ends = [sub_link.safe_end for sub_link in candidate_sub_links_sorted_by_end]
    logger.info('Generated sorted arrays')
    ##
    current_sub_link = current_sub_links[0]
    print 'Current:', current_sub_link.safe_start, current_sub_link.safe_end, current_sub_link.name
    # Nearest left neighboor without intersection
    left_neighbor_idx = bisect.bisect(sorted_ends, current_sub_link.safe_start)
    left_neighbor = candidate_sub_links_sorted_by_end[left_neighbor_idx - 1]
    current_sub_link.left_neighbor_end = left_neighbor.safe_end
    current_sub_link.left_neighbor_name = left_neighbor.name
    print 'Left neighbor:', left_neighbor.safe_start, left_neighbor.safe_end, left_neighbor.name
    # Nearest right neighboor without intersection
    right_neighbor_idx = bisect.bisect(sorted_starts, current_sub_link.safe_end)
    right_neighbor = candidate_sub_links_sorted_by_start[right_neighbor_idx]
    current_sub_link.right_neighbor_begin = right_neighbor.safe_start
    current_sub_link.right_neighbor_name = right_neighbor.name
    print 'Right neighbor:', right_neighbor.safe_start, right_neighbor.safe_end, right_neighbor.name