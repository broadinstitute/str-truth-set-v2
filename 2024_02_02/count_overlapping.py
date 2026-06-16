import gzip
import collections
import os
for i in 0, 1, 2, 3:
	print("------")
	counters = collections.defaultdict(int)
	result_set = set()
	with (gzip.open(os.path.expanduser("~/code/str-truth-set-v2/combined.51_samples.variants.bed.gz"), "rt") as f,
		  open(f"combined.51_samples.overlapping_variants.within_{i}bp.bed", "w") as out):
		# count total and adjacent repeats. The bed is position-sorted, but a variant can overlap (within i bp) an
		# earlier wide interval even when the immediately-preceding row does not, so track the running maximum end
		# seen so far on the current chromosome (and the row that produced it) rather than only the previous row.
		prev_chrom = None
		max_end = None
		max_end_row = None
		for line in f:
			fields = line.strip().split("\t")
			fields[1] = int(fields[1])
			fields[2] = int(fields[2])

			counters["total"] += 1
			if fields[0] == prev_chrom and fields[1] - max_end <= i:
				counters["adjacent"] += 1
				if tuple(max_end_row[0:3]) not in result_set:
					out.write("\t".join(map(str, max_end_row)) + "\n")
					result_set.add(tuple(max_end_row[0:3]))
				out.write(line)
				result_set.add(tuple(fields[0:3]))

			if fields[0] != prev_chrom or fields[2] > max_end:
				prev_chrom = fields[0]
				max_end = fields[2]
				max_end_row = fields

	print(f"Variants within {i}bp: {counters['adjacent']:,d} ({counters['adjacent']/counters['total']:.1%})")
	for key, value in counters.items():
		# print with %
		print(f"{key}: {value:,d} ({value/counters['total']:.1%})")