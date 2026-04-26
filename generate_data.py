import json
import random

diseases = {
    "Anthrax": ["Bacillus anthracis"],
    "Bacterial Vaginosis": ["Gardnerella vaginalis"],
    "Botulism": ["Clostridium botulinum"],
    "Brucellosis": ["Brucella abortus", "Brucella melitensis", "Brucella suis"],
    "Campylobacteriosis": ["Campylobacter jejuni"],
    "Cat Scratch Disease": ["Bartonella henselae"],
    "Cellulitis / Erysipelas": ["Streptococcus pyogenes", "Staphylococcus aureus"],
    "Chlamydia": ["Chlamydia trachomatis"],
    "Cholera": ["Vibrio cholerae"],
    "Diphtheria": ["Corynebacterium diphtheriae"],
    "Endocarditis": ["Staphylococcus aureus", "Streptococcus mutans", "Enterococcus faecalis"],
    "Food Poisoning": ["Salmonella enterica", "Campylobacter jejuni", "Escherichia coli", "Listeria monocytogenes", "Bacillus cereus", "Staphylococcus aureus", "Clostridium perfringens"],
    "Gas Gangrene": ["Clostridium perfringens"],
    "Gonorrhea": ["Neisseria gonorrhoeae"],
    "Helicobacter Pylori Infection": ["Helicobacter pylori"],
    "Intra-abdominal Infection": ["Escherichia coli", "Bacteroides fragilis", "Enterococcus faecalis"],
    "Legionnaires' Disease": ["Legionella pneumophila"],
    "Leprosy": ["Mycobacterium leprae"],
    "Leptospirosis": ["Leptospira interrogans"],
    "Lyme Disease": ["Borrelia burgdorferi"],
    "Meningitis": ["Neisseria meningitidis", "Streptococcus pneumoniae", "Haemophilus influenzae", "Listeria monocytogenes", "Streptococcus agalactiae"],
    "MRSA Infection": ["Staphylococcus aureus"],
    "Osteomyelitis": ["Staphylococcus aureus", "Pseudomonas aeruginosa"],
    "Otitis Media": ["Streptococcus pneumoniae", "Haemophilus influenzae", "Moraxella catarrhalis"],
    "Pertussis (Whooping Cough)": ["Bordetella pertussis"],
    "Plague": ["Yersinia pestis"],
    "Pneumonia": ["Klebsiella pneumoniae", "Streptococcus pneumoniae", "Haemophilus influenzae", "Mycoplasma pneumoniae", "Chlamydophila pneumoniae", "Legionella pneumophila", "Pseudomonas aeruginosa"],
    "Q Fever": ["Coxiella burnetii"],
    "Rocky Mountain Spotted Fever": ["Rickettsia rickettsii"],
    "Salmonellosis": ["Salmonella enterica"],
    "Sepsis / Bloodstream Infection": ["Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae", "Acinetobacter baumannii", "Enterococcus faecium", "Pseudomonas aeruginosa"],
    "Shigellosis": ["Shigella dysenteriae", "Shigella sonnei"],
    "Skin and Soft Tissue Infection": ["Staphylococcus aureus", "Streptococcus pyogenes", "Pseudomonas aeruginosa", "Acinetobacter baumannii"],
    "Syphilis": ["Treponema pallidum"],
    "Tetanus": ["Clostridium tetani"],
    "Tuberculosis": ["Mycobacterium tuberculosis"],
    "Tularemia": ["Francisella tularensis"],
    "Typhoid Fever": ["Salmonella typhi"],
    "Typhus": ["Rickettsia prowazekii"],
    "Urinary Tract Infection (UTI)": ["Escherichia coli", "Klebsiella pneumoniae", "Proteus mirabilis", "Staphylococcus saprophyticus", "Enterococcus faecalis", "Pseudomonas aeruginosa"],
    "Vibriosis": ["Vibrio parahaemolyticus", "Vibrio vulnificus"]
}

unique_bacteria = set()
for bugs in diseases.values():
    for bug in bugs:
        unique_bacteria.add(bug)

bacteria_list = sorted(list(unique_bacteria))

bacteria_genomes = {}
bases = ['A', 'C', 'G', 'T']
for bug in bacteria_list:
    abbrev = f"{bug.split()[0][0]}_{bug.split()[1]}_simulated_genome"
    lines = []
    # Gen 3 lines of 79 chars like original file
    for _ in range(3):
        lines.append("".join(random.choices(bases, k=79)))
    genome_seq = "\n".join(lines)
    header = f">{abbrev} | MicrobeNet"
    bacteria_genomes[bug] = f"{header}\n{genome_seq}"

# To ensure order is preserved, dictionaries in Python 3.7+ preserve insertion order.
# Just need to make sure we construct it in order.
ordered_bacteria_genomes = {k: bacteria_genomes[k] for k in bacteria_list}

out_data = {
    "diseases": diseases,
    "bacteria": ordered_bacteria_genomes
}

with open(r"c:\Users\tbaii\Downloads\DayDream (1)\data\bacteria_disease_variants.json", "w") as f:
    json.dump(out_data, f, indent=2)

print("Created updated dataset with", len(diseases), "diseases and", len(bacteria_list), "bacteria.")
