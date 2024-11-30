# RPQ-bench benchmark utility

This is an utility tool we've used for benchmarking linear algebra-based RPQ algorithm inside the [LAGraph infrastructre](https://github.com/GraphBLAS/LAGraph).

It just runs LAGraph algorithm and loads dataset in specific custom format. For more information on how to convert conventional [N-triples](https://www.w3.org/TR/n-triples/) into that format you may refer to the [la-rpq repository](https://github.com/SparseLinearAlgebra/la-rpq) containing the information on how to process the datasets and benchmark various different databases.

## Example usage

That's how you might evaluate benchmarks on the `dataset-mm` graph that contains 15 labels over 10 different queries from `queries-mm`.

```bash
./rpq-bench ./dataset-mm 15 ./queries-mm 10
```

