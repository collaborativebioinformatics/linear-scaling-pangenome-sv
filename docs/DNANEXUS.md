# DNAnexus Integration

## Setup

```bash
dx login
dx select <project>

# Start cloud workstation
dx run --instance-type mem3_ssd1_v2_x16 --ssh app-cloud_workstation

# Inside workstation
git clone <repo>
cd <repo>
bash scripts/check_environment.sh
```

## Storage Structure

- /data/hprc/
- /data/reference/
- /data/prepared/
- /graphs/baseline/
- /graphs/chunks/
- /graphs/merged/
- /variants/
- /benchmark/
- /web/
- /logs/
