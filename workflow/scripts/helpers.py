

def get_samples_of_species(species):

    reads = config.get("species", {}).get(species, {}).get("reads")

    if not reads:
        raise ValueError(f"No reads found for species {species} in config.")

    samples = []

    for r in reads:
        sample_id = get_sample_id_for_sample_path(r) # example: Dfun01
        samples.append(sample_id)
    
    if not samples:
        raise ValueError(f"No samples could be determined for species {species}.")

    return samples

def get_path_of_sample(species, sample):
    reads = config.get("species", {}).get(species, {}).get("reads")

    if not reads:
        raise ValueError(f"No reads found for species {species} in config.")

    for r in reads:

        sample_from_filename = get_sample_id_for_sample_path(r) # example: Dfun01
        if sample_from_filename == sample:
            return r
    
    raise ValueError(f"No read could be determined for sample {sample} of species {species}.")

def get_sample_id_for_sample_path(sample_path):
    filename = os.path.basename(sample_path) # example: Dfun01.fastq.gz
    sample_id = filename.replace(".fastq.gz", "") # example: Dfun01

    if not sample_id:
        raise ValueError(f"No sample ID could be determined for sample path {sample_path}.")

    return sample_id

def get_reference_of_species(species):
    reference = config.get("species", {}).get(species, {}).get("reference")

    if not reference:
        raise ValueError(f"No reference found for species {species} in config.")

    return reference

def get_reference_id_of_species(species):
    reference_path = get_reference_of_species(species)
    reference_id = os.path.basename(reference_path).split(".")[0] # example: Dfun01_ref
    if not reference_id:
        raise ValueError(f"No reference ID could be determined for species {species}.")
    return reference_id