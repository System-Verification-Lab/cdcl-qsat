# Search-driven clause learning for PRODSAT-QSAT 

Prototype implementation of a sound refuter for product-state satisfiability in quantum SAT.

## Quick start

Run the example in `src/main.rs`:

```bash
cargo run
```

The `main.rs` example shows how to:

- create a `QsatInstance` with a given number of qubits and search depth,
- add constraints as `(support, projector)` pairs,
- call `solve` and read the `QsatResult`.

