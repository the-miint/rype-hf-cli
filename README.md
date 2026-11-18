# rype-hf-cli

`rype-hf` is streaming host filtration: [rype](https://github.com/the-miint/rype)
removes reads that match a host index, and an optional minimap2 pass removes
any host reads rype missed. Reads flow through one pipe from gzipped input to
compressed survivors, with no intermediate files.

```
rype classify filter  ->  [minimap2 -a | samtools fastq (unmapped only)]  ->  bgzip -l1
```

`rype-hf` is a single bash script. It needs these tools on `PATH`:

| Tool | Needed for | Install from |
|---|---|---|
| [rype](https://github.com/the-miint/rype) with `classify filter` | always | source (Rust), step 2 |
| `bgzip` (htslib) | `*.gz` output | conda, Homebrew or apt, step 3 |
| `minimap2` and `samtools` | the rescue pass (`-m`) | conda, Homebrew or apt, step 3 |

Runs on Linux and macOS with bash 3.2 or later, so the macOS system bash works.

## Installation

Only the `apt` commands need root. With conda for step 3, everything installs
into your home directory, which also works on an HPC cluster.

### 1. Install the Rust toolchain and a C compiler

rype is built from source with Cargo. Skip this step if `cargo --version`
already works and `rustup update` leaves you on a recent stable Rust.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh   # accept the defaults
source "$HOME/.cargo/env"        # or open a new shell
```

Some of rype's dependencies compile C code, so you also need a C compiler:

- **macOS:** `xcode-select --install`
- **Debian/Ubuntu:** `sudo apt install build-essential`
- **HPC cluster:** usually already present; otherwise `module load gcc`

### 2. Install rype

`classify filter` is in rype's `main` branch but not yet in a crates.io
release: `cargo install rype` installs 1.0.0-rc.1, which predates it. Install
from GitHub instead:

```bash
cargo install --git https://github.com/the-miint/rype rype --bin rype --force
```

This builds a release binary and installs it as `~/.cargo/bin/rype`. The build
takes a few minutes. `--force` replaces an older rype if you have one.

To pin a commit for reproducibility, add `--rev <commit>`.

### 3. Install bgzip, minimap2 and samtools

Pick one of these.

**conda / mamba (recommended, Linux and macOS)**

```bash
conda create -n rype-hf -c conda-forge -c bioconda htslib minimap2 samtools
conda activate rype-hf
```

Activate the environment in every shell or batch job that runs `rype-hf`.

**Homebrew (macOS)**

```bash
brew install htslib minimap2 samtools
```

**apt (Debian/Ubuntu)**

```bash
sudo apt install tabix minimap2 samtools    # bgzip is in the tabix package
```

If you only run rype, with no `-m`, and write plain (not `*.gz`) output, you
can skip this step.

### 4. Install rype-hf

```bash
git clone https://github.com/the-miint/rype-hf-cli.git
mkdir -p ~/.local/bin
ln -sf "$PWD/rype-hf-cli/bin/rype-hf" ~/.local/bin/rype-hf
```

Make sure `~/.local/bin` is on your `PATH`. If it isn't, add this line to
`~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Because it is a symlink, `git -C rype-hf-cli pull` updates `rype-hf` in place.

### 5. Check the installation

```bash
rype-hf -V                              # rype-hf 0.1.0
rype classify filter --help | head -1   # fails if your rype is too old
bgzip --version | head -1
minimap2 --version
samtools --version | head -1
rype-hf-cli/tests/run-tests.sh          # should end with "N passed, 0 failed"
```

The tests need no index. They report `skip` for any test whose tool
(bgzip, minimap2, samtools) is missing.

## Building host indexes

`rype-hf` takes a rype index (`-i`) and, for the rescue pass, a minimap2
index (`-m`). Build both once from the same host genome, for example human
T2T-CHM13:

```bash
curl -LO https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

# rype index (defaults: k=64, w=50)
rype index create -o human.ryxdi -r chm13v2.0.fa.gz

# minimap2 index for the rescue pass
minimap2 -x sr -d human-t2t.mmi chm13v2.0.fa.gz
```

Build the `.mmi` with the same preset you will pass to `rype-hf -x`. minimap2
fixes its k-mer and window sizes when it builds the index, so a `-x sr` index
is for short reads. For long reads, build a separate index with, for example,
`-x map-hifi` or `-x map-ont`.

`rype index stats -i human.ryxdi` summarizes a rype index.

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

## Choosing tool binaries

Tools are found on `PATH`. To use a specific binary instead, set `RYPE_BIN`,
`MINIMAP2_BIN`, `SAMTOOLS_BIN` or `BGZIP_BIN`:

```bash
RYPE_BIN=/path/to/rype/target/release/rype rype-hf -i human.ryxdi ...
```

## Troubleshooting

| Error | Fix |
|---|---|
| `rype not found (set RYPE_BIN)` | Add `~/.cargo/bin` to `PATH` (`source "$HOME/.cargo/env"`), or set `RYPE_BIN`. |
| `... has no 'classify filter' subcommand` | Your rype is older than `classify filter`. Reinstall it (step 2), and check `which rype` is the new one. |
| `*.gz output needs bgzip` | Install htslib (step 3), or activate your conda environment. |
| `-m needs minimap2 and samtools` | Install them (step 3), or activate your conda environment. |
| `cannot detect a memory limit; set -M` | Pass rype's memory budget yourself, e.g. `-M 32G`. |
| `only N MiB left for rype` | The job's memory limit is too small once minimap2's index is reserved. Request more memory, or set `-M`. |

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
