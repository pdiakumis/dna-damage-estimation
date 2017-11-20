Estimation of DNA damage
========================


<!-- vim-markdown-toc GFM -->

* [Snakemake Stuff](#snakemake-stuff)
    * [Set up](#set-up)
        * [example runs](#example-runs)
        * [expand](#expand)
    * [resources](#resources)

<!-- vim-markdown-toc -->

Here we'll use the DNA Damage Estimator scripts by Laurence Ettwiller
([GitHub repo](https://github.com/Ettwiller/Damage-estimator)).
Since these scripts use the same steps for three different modules (essentially
the code has been copied over from one script to the other), I cleaned it out
and simplified some of the steps.
I've created a [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipeline that helps to run the individual steps in a cluster environment.

Snakemake Stuff
===============

Set up
------
* Create conda environment to run snakemake
  * `conda env create --name snakepipe --file environment.yaml`
  * Activate with `source activate snakepipe`

* `environment.yaml`: contains tools you want to use in the pipeline.
* `rule.snakefile`: contains one or more pipeline steps (called 'rules'). Each
  rule can run for multiple samples, thus creating multiple 'jobs'.
* `cluster.json`: cluster configuration file to specify partition, time etc.

### example runs

* Dry run and print shell command used:

```bash
snakemake -s rule.snakefile -np
```

* Plot DAG:

```bash
snakemake -s rule.snakefile --dag | dot -Tsvg > dag.svg
```

* Run on multiple cores:

```bash
snakemake -s rule.snakefile -j 4
```

* Run on PBS/TORQUE cluster:
    - Each job will be compiled into a shell script that is submitted with the
      given command (e.g. `qsub`). The `--jobs 100` flag limits the number of
      concurrently submitted jobs to 100.
    - You can decorate the specified submission command with parameters from the
      submitted job (e.g. `threads`).

```bash
snakemake -s rule.snakefile --cluster qsub --jobs 100
snakemake -s rule.snakefile --cluster "qsub -pe threaded {threads}" --jobs 100
```

* Run on SLURM cluster:

```bash
snakemake -j 999 \
          --cluster-config cluster.json \
          --cluster "sbatch -A {cluster.account} \
          -p {cluster.partition} \
          -n {cluster.n} \
          -t {cluster.time}"
```

where `cluster.json` is:

```
{
    "__default__" :
    {
        "account" : "my account",
        "time" : "00:15:00",
        "n" : 1,
        "partition" : "core"
    },
    "compute1" :
    {
        "time" : "00:20:00"
    }
}
```

### expand
Snakemake has a function in `snakemake.io` called `expand` which helps you
'expand' different
file lists as in the following:

```
PLOTFORMATS = ["pdf", "png"]
DATASETS = ["ds1", "ds2"]
```
* `expand("{dataset}/file.txt", dataset = DATASETS)`

Output:
```
["ds1/file.txt", "ds2/file.txt"]
```

* `expand("{dataset}/file.{ext}", dataset=DATASETS, ext=PLOTFORMATS)`

Output:
```
["ds1/file.pdf", "ds2/file.pdf", "ds1/file.png", "ds2/file.pdf"]
```

* `expand(["{dataset}/plot1.{ext}", "{dataset}/plot2.{ext}"], dataset=DATASETS, ext=PLOTFORMATS)`

Output:
```
["ds1/plot1.pdf", "ds1/plot2.pdf", "ds2/plot1.pdf", "ds2/plot2.pdf",
"ds1/plot1.png", "ds1/plot2.png", "ds2/plot1.png", "ds2/plot2.png"]
```

* `expand("{dataset}/plot1.{ext} {dataset}/plot2.{ext}".split(), zip, dataset=DATASETS, ext=PLOTFORMATS)`

Output:
```
["ds1/plot1.pdf", "ds1/plot2.pdf", "ds2/plot1.png", "ds2/plot2.png"]
```

* `expand("{{dataset}}/plot1.{ext}", ext=PLOTFORMATS)`

Output:
```
['{dataset}/plot1.pdf', '{dataset}/plot1.png']
```

## resources

* Snakemake docs: <https://snakemake.readthedocs.io/en/stable/>
* Snakemake code examples: <https://codegists.com/code/snakefile/>
* Useful rules repo: <https://github.com/percyfal/snakemake-rules>
* Snakemake on SLURM: <https://hpc.nih.gov/apps/snakemake.html>
