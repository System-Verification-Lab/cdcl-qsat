# Search-driven clause learning solver for PROD-QSAT 

Prototype implementation of a solver for Product-State Quantum SAT (under a promise gap).

## Quick start

With Rust and Cargo installed, run benchmarks on random instances from the repository root:

```bash
cargo run --release
```

The default run generates 128 samples for each qubit count from 3 through 9, with one more constraint than qubits. Each constraint acts on 3 qubits, and the promise gap is `0.125 / 512.0`. Samples are processed using 64 worker threads.

To save the results for later analysis:

```bash
cargo run --release > benchmark.csv
```

The output contains a CSV header followed by one row per sample, including the qubit and constraint counts, result, CPU time in seconds, solver statistics, and sample ID. Rows may arrive in any order because samples run concurrently.

## Customising a run

Edit `benchmark()` in [`src/utils.rs`](src/utils.rs) to change:

- `k` and `pgap`: the number of qubits per constraint and the promise gap.
- `Orchestrator::new(3, 9, 128)`: the first and last qubit counts (inclusive), and the number of samples per count. Keep the first count at least `k`.
- `threads`: the number of worker threads.

For a single random instance, call
```rust
sample_case::<ComplexFloat>(num_qubits, k, num_constraints, promise_gap, seed)`
```
from `src/main.rs` using `Float::singleton(...)` for the promise gap. It returns a `QsatDiagnosis` object containing the result and timing/statistics; the seed controls instance generation.
