import pstats

with open('cumprofiled.txt', 'w') as stream:
    stats = pstats.Stats('profiled.csv', stream=stream)
    stats.sort_stats("cumulative").print_stats()
