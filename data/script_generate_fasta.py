import random

# Output file
filename = "test_1gb.fasta"

# Approximate size in bytes
target_size = 1000 * 1024 * 1024  # 10 MB

# DNA bases
bases = ['A', 'C', 'G', 'T']

# Open file for writing
with open(filename, "w") as f:
    seq_size = 0
    seq_id = 1
    while seq_size < target_size:
        # Write a FASTA header
        f.write(f">seq{seq_id}\n")
        # Generate a random sequence of ~80 bases per line
        line_len = 80
        seq = ''.join(random.choices(bases, k=line_len))
        f.write(seq + "\n")
        seq_size += len(seq) + len(f">seq{seq_id}\n") + 1  # approximate line size
        seq_id += 1

print(f"Generated {filename}, ~{seq_size/1024/1024:.2f} MB")
