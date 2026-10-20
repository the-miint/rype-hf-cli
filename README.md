# rype-hf-cli

`rype-hf` is streaming host filtration: [rype](https://github.com/the-miint/rype)
removes reads that match a host index, and an optional minimap2 pass removes
any host reads rype missed. Reads flow through one pipe from gzipped input to
compressed survivors, with no intermediate files.

```
rype classify filter  ->  [minimap2 -a | samtools fastq (unmapped only)]  ->  bgzip -l1
```

## Usage

```bash
# Short read: rype, then a minimap2 rescue pass; survivors as R1/R2 files
rype-hf -i human.ryxdi -1 R1.fq.gz -2 R2.fq.gz -t 0.0001 \
        -m human-t2t.mmi -x sr -p 8 -o clean_R1.fastq.gz -O clean_R2.fastq.gz

# Long read: rype alone
rype-hf -i human.ryxdi -1 hifi.fq.gz -t 0.05 -p 8 -o clean.fastq.gz
```

| Option | Meaning |
|---|---|
| `-i` | rype host index (`.ryxdi`) |
| `-1`, `-2` | reads (FASTA/FASTQ, optionally gzipped); `-2` for paired input |
| `-t` | rype score threshold: a read or pair scoring at or above it in any bucket is host |
| `-o` | output: `-` for stdout, `*.gz` for BGZF, otherwise plain FASTQ |
| `-O` | paired input only: R2 survivors go here and R1 survivors to `-o` |
| `-m`, `-x` | minimap2 host index and preset for the rescue pass (default preset `sr`) |
| `-p` | threads (default: `$SLURM_CPUS_PER_TASK`, else all CPUs) |
| `-M` | rype `--max-memory`, e.g. `48G` (default: auto, see Memory) |
| `-s` | write a one-row TSV of read counts per stage |
| `-V` | print the version |

`rype-hf -h` prints the same list.

## What it keeps

- **rype stage.** A read or pair is host when it scores at or above `-t` in
  any bucket of the index. Pairs are scored as one query and kept or dropped
  together.
- **Rescue stage (`-m`).** rype's survivors are aligned to the minimap2 index
  (`--secondary=no`). A single-end read is kept if it stays unmapped. A pair is
  kept only if both mates stay unmapped, so a pair is dropped if either mate
  aligns.
- **Headers.** Survivors are written as `<id>/1` and `<id>/2`. Header comments
  are not carried through.

## Output

- **Paired input.** With `-O`, R1 and R2 survivors go to separate files, in the
  same order. Without it they are interleaved into `-o`.
- **Compression.** `*.gz` outputs are written as BGZF with `bgzip -l 1`, using
  `-p` threads. With `-O` each file gets its own bgzip with half the threads.
  When most reads survive, level 1 is several times faster than gzip's default
  level 6, for output about 8% larger.

## Memory

By default rype gets 90% of the memory limit minus 512 MiB. Long reads can run
slightly past rype's budget, and the kernel also charges page cache to the
limit, so a budget equal to the limit can run out of memory.

The memory limit is the first of these that is found:

1. the tightest cgroup v2 `memory.max` from the process's cgroup up to the
   root (under Slurm the limit sits on the job's cgroup, not the task's);
2. `SLURM_MEM_PER_NODE`;
3. `MemAvailable` (Linux);
4. total RAM (macOS, `sysctl hw.memsize`).

With `-m`, minimap2's footprint is subtracted first: 1.25 × the `.mmi` size
plus 1 GiB, because minimap2 holds its whole index in memory for the run.
rype-hf always passes rype an explicit budget, because rype's own `auto` only
reads its own cgroup and does not know minimap2 is running alongside.

rype has no hard memory floor: a smaller budget means more passes over the
index, which costs time. rype-hf allows budgets down to 256 MiB; set `-M` to
override the automatic value.

## Requirements

- **Linux or macOS**, bash 3.2 or later.
- **rype with `classify filter`.** Until it is merged into
  [the-miint/rype](https://github.com/the-miint/rype), build branch
  `michaeliter/classify-filter` of
  [michaeliter/RYpe](https://github.com/michaeliter/RYpe):
  ```bash
  git clone -b michaeliter/classify-filter https://github.com/michaeliter/RYpe.git
  cd RYpe && cargo build --release      # binary: target/release/rype
  ```
- **bgzip** (htslib) for `*.gz` output.
- **minimap2 and samtools** for the rescue pass (`-m`).

Tools are found on `PATH`, or set explicitly with `RYPE_BIN`, `MINIMAP2_BIN`,
`SAMTOOLS_BIN` and `BGZIP_BIN`.

## Tests

```bash
tests/run-tests.sh
```

The tests cover options, output layouts, the rescue pass, failure handling and
the stats file. rype is replaced by a stand-in that keeps every read
(`tests/fake-rype`), so no index is needed. Tests that need bgzip, minimap2 or
samtools are skipped when those are not installed.

## License

BSD 3-Clause, as rype. See [LICENSE](LICENSE).
