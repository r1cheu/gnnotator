import os

species = snakemake.wildcards.species

evm_weight = snakemake.config["evm_weights"]

output_file = snakemake.output.weight

os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for category, items in evm_weight.items():
        for name, weight in items.items():
            if category == "TRANSCRIPT":
                dynamic_name = f"{name}-{species}_DB_pasa"
                f.write(f"{category:<20} {dynamic_name:<25} {weight}\n")
            else:
                f.write(f"{category:<20} {name:<25} {weight}\n")
