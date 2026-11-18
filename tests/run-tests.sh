#!/usr/bin/env bash
# Plumbing tests for rype-hf: options, output layouts (interleaved, split
# R1/R2, gz and plain), the minimap2 rescue, failure handling and the stats
# file. rype itself is replaced by tests/fake-rype (every read survives), so no
# rype index is needed; classification is rype's own test suite's job.
#
# Needs gzip, awk and paste. bgzip, minimap2 and samtools are optional: tests
# that need them are skipped when they are not found (on PATH, or via
# BGZIP_BIN / MINIMAP2_BIN / SAMTOOLS_BIN).
#
# Usage: tests/run-tests.sh
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
HF=${HERE}/../bin/rype-hf
export RYPE_BIN=${HERE}/fake-rype
export BGZIP_BIN=${BGZIP_BIN:-$(command -v bgzip || true)}
export MINIMAP2_BIN=${MINIMAP2_BIN:-$(command -v minimap2 || true)}
export SAMTOOLS_BIN=${SAMTOOLS_BIN:-$(command -v samtools || true)}
W=$(mktemp -d "${TMPDIR:-/tmp}/rype-hf-tests.XXXXXX")
trap 'rm -rf "${W}"' EXIT
cd "${W}"

pass=0 fail=0 skip=0
ok ()   { pass=$((pass + 1)); echo "ok    $1"; }
bad ()  { fail=$((fail + 1)); echo "FAIL  $1"; }
skip () { skip=$((skip + 1)); echo "skip  $1 ($2)"; }
check () { local name=$1; shift; if "$@"; then ok "${name}"; else bad "${name}"; fi; }

# Run rype-hf with a watchdog (no `timeout` on stock macOS): a hang is a failure.
run () {
    bash "${HF}" "$@" 2> last.err & local pid=$! n=0
    while kill -0 "${pid}" 2>/dev/null; do
        sleep 0.2; n=$((n + 1))
        (( n > 300 )) && { kill "${pid}"; echo "rype-hf hung" >> last.err; return 124; }
    done
    wait "${pid}"
}
records () { case $1 in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac | awk 'END { print NR / 4 }'; }
names ()   { case $1 in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac | awk 'NR % 4 == 1'; }
mates_match () { diff <(names "$1" | sed 's,/1$,,') <(names "$2" | sed 's,/2$,,') > /dev/null; }
deinterleave () {  # <interleaved> <1|2>
    case $1 in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac | paste - - - - - - - - \
        | if [[ $2 == 1 ]]; then cut -f1-4; else cut -f5-8; fi | tr '\t' '\n'
}
suffixes_kept () { ! names "$1" | grep -qv '/1$' && ! names "$2" | grep -qv '/2$'; }
same_content () { cmp -s <(case $1 in (*.gz) gzip -dc "$1" ;; (*) cat "$1" ;; esac) "$2"; }

# --- test data: a 50 kb "host" genome, 100 pairs drawn from it and 100 random
# pairs (named host::pN and microbe::pN), plus a single-end set.
awk 'BEGIN { srand(1); s = ""; for (i = 0; i < 50000; i++) s = s substr("ACGT", int(rand() * 4) + 1, 1)
             print ">host"; print s }' > host.fa
awk -v host="$(sed -n 2p host.fa)" 'BEGIN {
    srand(2); q = ""; for (i = 0; i < 150; i++) q = q "I"
    comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
    for (i = 0; i < 200; i++) {
        if (i % 2 == 0) { frag = substr(host, int(rand() * 49000) + 1, 400); tag = "host" }
        else { frag = ""; for (j = 0; j < 400; j++) frag = frag substr("ACGT", int(rand() * 4) + 1, 1); tag = "microbe" }
        r2 = ""; for (j = 400; j > 250; j--) r2 = r2 comp[substr(frag, j, 1)]
        printf "@%s::p%d\n%s\n+\n%s\n", tag, i, substr(frag, 1, 150), q > "r1.fq"
        printf "@%s::p%d\n%s\n+\n%s\n", tag, i, r2, q > "r2.fq"
    } }'
gzip -k r1.fq r2.fq
touch idx.ryxdi
COMMON=(-i idx.ryxdi -t 0.1 -M 1G -p 2)

# --- options
check "-V prints the version" bash -c "bash '${HF}' -V | grep -q '^rype-hf [0-9]'"
run "${COMMON[@]}" -1 r1.fq.gz -o x.fq -O y.fq; check "-O without -2 is refused" test $? -eq 2
run "${COMMON[@]}" -1 r1.fq.gz -2 r2.fq.gz -o - -O y.fq; check "-O with -o - is refused" test $? -eq 2
run -i idx.ryxdi -1 r1.fq.gz -o x.fq; check "missing -t prints usage" test $? -eq 2

# --- rype only
run "${COMMON[@]}" -1 r1.fq.gz -o se.fq -s se.tsv
check "single-end: all 200 reads, plain FASTQ" test "$(records se.fq)" = 200
check "stats file: reads_in and rype_kept" test "$(sed -n 2p se.tsv)" = "$(printf '200\t200\tNA')"
run "${COMMON[@]}" -1 r1.fq.gz -2 r2.fq.gz -o il.fq
check "paired, interleaved: 400 records" test "$(records il.fq)" = 400
check "paired, interleaved: /1 then /2" bash -c "diff <(awk 'NR%8==1' il.fq | sed 's,/1\$,,') <(awk 'NR%8==5' il.fq | sed 's,/2\$,,') >/dev/null"
run "${COMMON[@]}" -1 r1.fq.gz -2 r2.fq.gz -o s_R1.fq -O s_R2.fq
check "paired, split: 200 + 200 records" test "$(records s_R1.fq) $(records s_R2.fq)" = "200 200"
check "paired, split: mates in the same order" mates_match s_R1.fq s_R2.fq
check "paired, split R1 == interleaved R1" same_content s_R1.fq <(deinterleave il.fq 1)
check "paired, split R2 == interleaved R2" same_content s_R2.fq <(deinterleave il.fq 2)
if [[ -x ${BGZIP_BIN} ]]; then
    run "${COMMON[@]}" -1 r1.fq.gz -2 r2.fq.gz -o g_R1.fq.gz -O g_R2.fq.gz
    check "split .gz: valid BGZF" bash -c "'${BGZIP_BIN}' -t g_R1.fq.gz && '${BGZIP_BIN}' -t g_R2.fq.gz"
    check "split .gz == plain split" bash -c "cmp -s <(gzip -dc g_R1.fq.gz) s_R1.fq && cmp -s <(gzip -dc g_R2.fq.gz) s_R2.fq"
else
    skip "split .gz output" "no bgzip"
fi

# --- minimap2 rescue: host pairs map and are dropped, microbe pairs survive
if [[ -x ${MINIMAP2_BIN} && -x ${SAMTOOLS_BIN} ]]; then
    "${MINIMAP2_BIN}" -x sr -d host.mmi host.fa 2> /dev/null
    RES=("${COMMON[@]}" -m host.mmi -x sr -1 r1.fq.gz -2 r2.fq.gz)
    run "${RES[@]}" -o m_il.fq -s m.tsv
    check "rescue, interleaved: only the 100 microbe pairs" \
        test "$(names m_il.fq | grep -c '^@microbe::') $(records m_il.fq)" = "200 200"
    run "${RES[@]}" -o m_R1.fq -O m_R2.fq
    check "rescue, split: 100 + 100 records" test "$(records m_R1.fq) $(records m_R2.fq)" = "100 100"
    check "rescue, split: /1 and /2 suffixes kept" suffixes_kept m_R1.fq m_R2.fq
    check "rescue, split R1 == interleaved R1" same_content m_R1.fq <(deinterleave m_il.fq 1)
    check "rescue, split R2 == interleaved R2" same_content m_R2.fq <(deinterleave m_il.fq 2)
    FAKE_FAIL=1 run "${RES[@]}" -o f_R1.fq -O f_R2.fq
    check "rescue, split: failing rype exits non-zero, no hang" test $? -eq 3
else
    skip "minimap2 rescue" "no minimap2/samtools"
fi
FAKE_FAIL=1 run "${COMMON[@]}" -1 r1.fq.gz -2 r2.fq.gz -o f_R1.fq -O f_R2.fq
check "split: failing rype exits non-zero, no hang" test $? -eq 3
FAKE_FAIL=1 run "${COMMON[@]}" -1 r1.fq.gz -o f.fq
check "single-end: failing rype exits non-zero" test $? -eq 3

# --- memory auto-detection gives a positive budget on this machine
run -i idx.ryxdi -t 0.1 -p 2 -1 r1.fq.gz -o auto.fq
check "-M auto: memory limit detected" grep -q 'rype --max-memory [0-9]* MiB' last.err

echo "${pass} passed, ${fail} failed, ${skip} skipped"
(( fail == 0 ))
