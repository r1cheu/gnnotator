import os

db_name = f"{snakemake.wildcards.species}_DB_pasa.sqlite"
script_params = snakemake.params.conf["params"]

os.makedirs(os.path.dirname(snakemake.output.pasa_config), exist_ok=True)
with open(snakemake.output.pasa_config, "w") as f:
    f.write(f"DATABASE=./{db_name}\n")
    for script_name, params_dict in script_params.items():
        for param_key, param_val in params_dict.items():
            f.write(f"{script_name}:{param_key}={param_val}\n")
