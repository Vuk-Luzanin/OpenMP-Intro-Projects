# MUPS
Homework for our [Multiprocessor Systems class](http://mups.etf.bg.ac.rs/) at [School of Electrical Engineering, University of Belgrade](https://www.etf.bg.ac.rs/) using [OpenMP](https://en.wikipedia.org/wiki/OpenMP). The programming ones involved programs that needed to be sped up using the mentioned technology. 

## Homework
All homework tasks are in Serbian.

- [Homework task (OpenMP)](https://web.archive.org/web/20230710231639im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ1_2022-2023.pdf)
    - [Initial code for the first task](https://web.archive.org/web/20230710231703im_/http://mups.etf.bg.ac.rs/dz/2022-2023/MPS_DZ1_OpenMP.zip)

## Structure
- `run.py` is our test runner that plots our speedup graphs and logs the output
- `Makefile` builds all tasks for us
- `src` contains homework files named as `name.c`.
- Tests are executed with:

```bash
python3 run.py <test_name>
```

Where <test_name> is optional (prime/feynman/moldyn) — if omitted, all tests will be executed.

- Compiled binaries, logs, and generated SVG result charts will be placed in the gen/ directory.

Results and logs will appear inside the gen/ directory.