MDT uses the [pytest](https://docs.pytest.org/en/latest/) framework.


### Running tests

To run all tests in order, just run:
```bash
cd [...]/molecular-design-toolkit/moldesign/_tests
py.test
```

For more speed, run `py.test -n N` where _N_ is the number of concurrent tests to run; it's reasonable to set _N_ to twice the number of available physical CPUs.

Note that, with Docker for Mac/Windows, Docker may not have access to all CPUs and/or RAM - see the Docker preferences to adjust this.


### Hints


##### docker

1. Make sure your development docker containers are up-to-date by running
```bash
cd molecular-design-toolkit/DockerMakefiles
docker-make --all --tag dev
```

2. To debug exceptions from inside a docker container, you'll need to know the docker container ID (this will usually be printed to the terminal). You can use the `debug-job` script in this directory to launch a debugger in the crashed container. Run `./debug-job --help` to see usage options.

3. Run `docker system prune` to remove all stopped containers and most unused images.



##### pytest

See `py.test --help` for many, many useful options. Here are some useful flags:

- `-n N`: run N tests concurrently (see above)
- `--lf`: re-run only the tests that failed last time
- `-x`: stop testing after a single failure
- `--pdb`: launch python's PDB debugger after each failure (this has no relationship to the Protein DataBank)
 
It's often useful to run
`py.test --lf -x --pdb` to launch a debugger then exit after any test failures.

