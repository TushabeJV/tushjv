from bs4 import BeautifulSoup
from Bio import Entrez
import csv

# Set your email address for Entrez
Entrez.email = "tushabejohnvianney@gmail.com"

def extract_taxon_ids(krona_html_file):
    with open(krona_html_file, 'r') as file:
        soup = BeautifulSoup(file, 'lxml')

    taxon_ids = set()

    # Find all 'val' elements within 'taxon' tags
    for taxon in soup.find_all('taxon'):
        val_tag = taxon.find('val')
        if val_tag and 'href' in val_tag.attrs:
            taxon_id = val_tag['href']
            taxon_ids.add(taxon_id)

    return list(taxon_ids)

def fetch_species_name(taxon_id):
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            return records[0]["ScientificName"]
        else:
            return None
    except Exception as e:
        print(f"Error fetching species name for taxon ID {taxon_id}: {e}")
        return None

def create_taxon_to_species_mapping(taxon_ids):
    taxon_to_species = {}
    for taxon_id in taxon_ids:
        species_name = fetch_species_name(taxon_id)
        if species_name:
            taxon_to_species[taxon_id] = species_name
        else:
            taxon_to_species[taxon_id] = "Unknown"
    return taxon_to_species

def save_mapping_to_file(taxon_to_species, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Taxon_ID', 'Species_Name']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for taxon_id, species_name in taxon_to_species.items():
            writer.writerow({'Taxon_ID': taxon_id, 'Species_Name': species_name})

# Example usage
krona_html_file = 'krona_output.html'
output_file = 'taxon_to_species_mapping.txt'

taxon_ids = extract_taxon_ids(krona_html_file)
print("Extracted Taxon IDs:", taxon_ids)

taxon_to_species_mapping = create_taxon_to_species_mapping(taxon_ids)
print("Taxon ID to Species Mapping:", taxon_to_species_mapping)

# Save the mapping to a file
save_mapping_to_file(taxon_to_species_mapping, output_file)
print(f"Mapping saved to {output_file}")