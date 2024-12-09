```bash
bsub -P jhq -J parse -n 1 -q interactive -R "rusage[mem=4001]" -W 10:00 -Is "bash"

cd /home/jqu/project/parse_biosiences/09-cell_num_express_OLAH
```







