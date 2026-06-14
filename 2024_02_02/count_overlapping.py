import gzip
import collections
import os
for i in 0, 1, 2, 3:
	print("------")
	counters = collections.defaultdict(int)
	result_set = set()
	with (gzip.open(os.path.expanduser("~/code/str-truth-set-v2/combined.51_samples.variants.bed.gz"), "rt") as f,
		  open(f"combined.51_samples.overlapping_variants.within_{i}bp.bed", "w") as out):
		# count total and adjacent repeats
		previous = None
		for line in f:
			fields = line.strip().split("\t")
			fields[1] = int(fields[1])
			fields[2] = int(fields[2])

			counters["total"] += 1
			if previous and fields[0] == previous[0] and int(fields[1]) - int(previous[2]) <= i:
				counters["adjacent"] += 1
				if tuple(previous[0:3]) not in result_set:
					out.write("\t".join(map(str, previous)) + "\n")
					result_set.add(tuple(previous[0:3]))
				out.write(line)
				result_set.add(tuple(fields[0:3]))

			previous = fields

	print(f"Variants within {i}bp: {counters['adjacent']:,d} ({counters['adjacent']/counters['total']:.1%})")
	for key, value in counters.items():
		# print with %
		print(f"{key}: {value:,d} ({value/counters['total']:.1%})")