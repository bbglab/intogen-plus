[genome]
build = option('hg18', 'hg19', 'hg38', 'mm10', 'c3h', default='hg19')


[statistic]

sampling = integer(default=100000)
sampling_chunk = integer(default=100)



[settings]
cores = integer(default=None)

[logging]
version = integer(default=1)