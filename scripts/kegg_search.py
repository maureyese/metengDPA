from Bio.KEGG import REST
import pandas as pd
from typing import DefaultDict, Tuple, List, Optional
from dataclasses import dataclass, field
import re
import time

# -------------------
# Classes
# -------------------

@dataclass
class GeneralKegg:
    """Class to store KEGG pathway and organism data"""
    pathway_df: pd.DataFrame
    organism_df: pd.DataFrame
    
    def __repr__(self):
        return f"GeneralKegg(pathway_df shape: {self.pathway_df.shape}, organism_df shape: {self.organism_df.shape})"

@dataclass
class ReactionData:
    """Class to store reaction-related information"""
    reaction_id: str
    reaction_equation: str
    substrates: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    substrate_compounds: List[str] = field(default_factory=list)
    product_compounds: List[str] = field(default_factory=list)
    
    def __repr__(self):
        return f"ReactionData({self.reaction_id}: {self.reaction_equation})"

@dataclass
class EnzymeData:
    """Class to store enzyme-related information"""
    enzyme_name: str
    ec_numbers: List[str] = field(default_factory=list)
    reaction: str = ""
    gene: str = ""
    kegg_gene_id: str = ""
    ko_numbers: List[str] = field(default_factory=list)
    reaction_data: List[ReactionData] = field(default_factory=list)
    substrates: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    
    def __repr__(self):
        return f"EnzymeData(gene: {self.gene}, EC: {self.ec_numbers}, reactions: {len(self.reaction_data)})"

@dataclass
class PathwayData:
    """Class to store pathway-related information"""
    pathway_name: str
    kegg_id: str
    global_path: bool
    enzymes: List[EnzymeData] = field(default_factory=list)
    
    def __repr__(self):
        return f"PathwayData(name: {self.pathway_name}, enzymes: {len(self.enzymes)})"

@dataclass
class OrganismData:
    """Class to store associated pathways to a specific organism"""
    organism_name: str
    abbreviation: str
    pathways_info: pd.DataFrame

    def __repr__(self):
        return f"OrganismData(organism: {self.organism_name}, pathways_info: {len(self.pathways_info.shape)})"

# -------------------
# 1. Collect KEGG complete data
# -------------------

def collect_kegg() -> GeneralKegg:
    '''Collect complete KEGG data from API and return GeneralKegg object'''
    # Print general information
    gen_info = REST.kegg_info("kegg").read()
    print("="*50)
    print("General KEGG information:\n\n")
    print(gen_info)
    print("="*50)

    print("Retrieving pathways from KEGG...")
    pathways = REST.kegg_list("pathway").read()
    path_df = pd.DataFrame(columns=["kegg_id", "pathway_name"])
    for line in pathways.split("\n"):
        if len(line.split("\t")) != 2:
            continue
        path_id, path_name = line.split("\t")
        path_df.loc[len(path_df)] = [path_id, path_name]
    print(f"Total pathways retrieved: {path_df.shape}\n")

    print("Retrieving organisms from KEGG...")
    organisms = REST.kegg_list("organism").read()
    parts = organisms.split("\n")
    # Filter out empty lines and only process non-empty ones
    non_empty_parts = [line for line in parts if line.strip() and len(line.split("\t")) == 4]

    org_df = pd.DataFrame(columns=["kegg_id", "abbreviation", "organism", "taxonomy"])
    for line in non_empty_parts:
        kegg_id, abbr, org, tax = line.split("\t")
        org_df.loc[len(org_df)] = [kegg_id, abbr, org, tax]
    print(f"Total organisms retrieved: {org_df.shape}")
    
    # Create and return the class instance
    kegg_data = GeneralKegg(pathway_df=path_df, organism_df=org_df)
    return kegg_data

# -------------------
# 2. Identify associated pathways to organism
# -------------------

def retrieve_organism_pathways(
    abbr: str, 
    general_data: GeneralKegg
) -> OrganismData:
    """Retrieve associated pathways from a specific organism"""

    print(f"Retrieving pathways for organism: {abbr}")

    # Get organism name from general data
    org_row = general_data.organism_df[general_data.organism_df['abbreviation'] == abbr]
    if org_row.empty:
        organism_name = abbr
    else:
        organism_name = org_row.iloc[0]['organism']
    
    pathways = REST.kegg_list(f"pathway/{abbr}").read()

    pathway_df = pd.DataFrame(columns=["kegg_id", "pathway_name"])

    for line in pathways.strip().split("\n"):
        if line and "\t" in line:
            pathway_id, pathway_name = line.split("\t")
            pathway_df.loc[len(pathway_df)] = [pathway_id, pathway_name]
    print(f"Total pathways retrieved from {abbr}: {pathway_df.shape}")

    return OrganismData(
        organism_name=organism_name,
        abbreviation=abbr,
        pathways_info=pathway_df
    )

# -------------------
# 3. Helper functions for parsing data from pathways of organism
# -------------------

def parse_ec_numbers(text: str) -> List[str]:
    """Extract EC numbers from text"""
    # Pattern for EC numbers: [EC:1.2.3.4] or [EC:1.2.3.4 1.2.3.5]
    ec_pattern = r'\[EC:([\d\.\s]+)\]'
    matches = re.findall(ec_pattern, text)
    if matches:
        # Split multiple EC numbers by space and clean
        ec_numbers = []
        for match in matches:
            ec_numbers.extend([ec.strip() for ec in match.split() if ec.strip()])
        return ec_numbers
    return []

def parse_ko_numbers(text: str) -> List[str]:
    """Extract KO numbers from text"""
    # Pattern for KO numbers: [KO:K12345]
    ko_pattern = r'\[KO:([^\]]+)\]'
    matches = re.findall(ko_pattern, text)
    if matches:
        ko_numbers = []
        for match in matches:
            ko_numbers.extend([ko.strip() for ko in match.split() if ko.strip()])
        return ko_numbers
    return []

def parse_gene_line(line: str, organism_abbr: str) -> Optional[EnzymeData]:
    """Parse a single GENE line from KEGG pathway data"""
    line = line.strip()
    if not line:
        return None
    
    # KEGG GENE line format: "gene_id  gene_symbol; description [KO:...] [EC:...]"
    # Example: "122622  ADSS1; adenylosuccinate synthase 1 [KO:K01939] [EC:6.3.4.4]"
    
    # Split by multiple spaces to separate gene ID from the rest
    parts = re.split(r'\s{2,}', line)  # Split by 2 or more spaces
    if len(parts) < 2:
        return None
    
    gene_id = parts[0].strip()
    rest_of_line = parts[1].strip()
    
    # Extract gene symbol and enzyme name
    if ';' in rest_of_line:
        gene_symbol = rest_of_line.split(';')[0].strip()
        enzyme_name = rest_of_line.split(';')[1].strip()
        
        # Remove everything after [ to get clean enzyme name
        if '[' in enzyme_name:
            enzyme_name = enzyme_name.split('[')[0].strip()
    else:
        gene_symbol = rest_of_line
        enzyme_name = rest_of_line
        if '[' in enzyme_name:
            enzyme_name = enzyme_name.split('[')[0].strip()
    
    # Extract EC numbers
    ec_numbers = parse_ec_numbers(rest_of_line)
    
    # Extract KO numbers
    ko_numbers = parse_ko_numbers(rest_of_line)
    
    # Create full KEGG gene ID
    kegg_gene_id = f"{organism_abbr}:{gene_id}"
    
    return EnzymeData(
        enzyme_name=enzyme_name,      # "adenylosuccinate synthase 1"
        ec_numbers=ec_numbers,        # ["6.3.4.4"]
        gene=gene_symbol,             # "ADSS1"
        kegg_gene_id=kegg_gene_id,    # "hsa:122622"
        ko_numbers=ko_numbers         # ["K01939"]
    )

# -------------------
# 3. Helper functions for parsing data from enzymes
# -------------------

def parse_all_reac_line(all_reac_line: str) -> List[ReactionData]:
    """Parse ALL_REAC line to get reaction IDs"""
    # Example: "R00256; (other) R01579 R06134"
    reaction_data = []
    
    # Extract all reaction IDs
    reaction_ids = re.findall(r'R\d+', all_reac_line)
    
    for reaction_id in reaction_ids:
        try:
            # Get detailed reaction information
            reaction_info = REST.kegg_get(f"rn:{reaction_id}").read()
            
            # Parse EQUATION section directly
            equation = ""
            if "EQUATION" in reaction_info:
                eq_start = reaction_info.find("EQUATION")
                if eq_start != -1:
                    remaining_text = reaction_info[eq_start + len("EQUATION"):]
                    next_section_match = re.search(r'\n[A-Z]+\s', remaining_text)
                    if next_section_match:
                        eq_section = remaining_text[:next_section_match.start()]
                    else:
                        eq_section = remaining_text
                    
                    # Get the first line of EQUATION section
                    equation_lines = eq_section.strip().split("\n")
                    if equation_lines:
                        equation = equation_lines[0].strip()
            
            # Parse substrates and products from equation
            substrates = []
            products = []
            substrate_compounds = []
            product_compounds = []
            
            if ' = ' in equation:
                substrates_str, products_str = equation.split(' = ', 1)
                substrates = [s.strip() for s in substrates_str.split(' + ')]
                products = [p.strip() for p in products_str.split(' + ')]
            
            # Extract compound information from SUBSTRATE and PRODUCT sections
            if "SUBSTRATE" in reaction_info:
                substrate_start = reaction_info.find("SUBSTRATE")
                if substrate_start != -1:
                    remaining_text = reaction_info[substrate_start + len("SUBSTRATE"):]
                    next_section_match = re.search(r'\n[A-Z]+\s', remaining_text)
                    if next_section_match:
                        substrate_section = remaining_text[:next_section_match.start()]
                    else:
                        substrate_section = remaining_text
                    
                    for line in substrate_section.split("\n"):
                        line = line.strip()
                        if line:
                            compound_info = parse_compound_line(line)
                            if compound_info:
                                substrate_compounds.append(compound_info)
            
            if "PRODUCT" in reaction_info:
                product_start = reaction_info.find("PRODUCT")
                if product_start != -1:
                    remaining_text = reaction_info[product_start + len("PRODUCT"):]
                    next_section_match = re.search(r'\n[A-Z]+\s', remaining_text)
                    if next_section_match:
                        product_section = remaining_text[:next_section_match.start()]
                    else:
                        product_section = remaining_text
                    
                    for line in product_section.split("\n"):
                        line = line.strip()
                        if line:
                            compound_info = parse_compound_line(line)
                            if compound_info:
                                product_compounds.append(compound_info)
                
            reaction_data.append(ReactionData(
                reaction_id=reaction_id,
                reaction_equation=equation,
                substrates=substrates,
                products=products,
                substrate_compounds=substrate_compounds,
                product_compounds=product_compounds
            ))
                
        except Exception as e:
            print(f"Error retrieving reaction {reaction_id}: {e}")
    
    return reaction_data

def parse_compound_line(compound_line: str) -> str:
    """Parse compound line to extract compound information"""
    # Example: "L-glutamine [CPD:C00064]"
    compound_match = re.match(r'(.+?)\s+\[CPD:([^\]]+)\]', compound_line)
    if compound_match:
        compound_name = compound_match.group(1).strip()
        compound_id = compound_match.group(2).strip()
        return f"{compound_name} [{compound_id}]"
    return compound_line.strip()

def parse_reaction_line(reaction_line: str) -> Optional[ReactionData]:
    """Parse a REACTION line from KEGG EC entry"""
    # Example: "L-glutamine + H2O = L-glutamate + NH3 [RN:R00256]"
    reaction_match = re.match(r'(.+)\s+\[RN:([^\]]+)\]', reaction_line)
    if reaction_match:
        equation = reaction_match.group(1).strip()
        reaction_id = reaction_match.group(2).strip()
        
        # Parse equation to get substrates and products
        if ' = ' in equation:
            substrates_str, products_str = equation.split(' = ', 1)
            substrates = [s.strip() for s in substrates_str.split(' + ')]
            products = [p.strip() for p in products_str.split(' + ')]
        else:
            substrates = []
            products = []
        
        return ReactionData(
            reaction_id=reaction_id,
            reaction_equation=equation,
            substrates=substrates,
            products=products
        )
    return None

# -------------------
# 3. Get detailed information regarding an enzyme (from EC Number)
# -------------------

def retrieve_ec_information(ec_numbers: List[str]) -> List[ReactionData]:
    """Retrieve reaction information for EC numbers"""
    reaction_data = []

    for ec_number in ec_numbers:
        try:
            # Get EC entry from KEGG
            time.sleep(0.5)
            ec_info = REST.kegg_get(f"ec:{ec_number}").read()

            # Parse REACTION section directly
            if "REACTION" in ec_info:
                # Find REACTION section
                reaction_start = ec_info.find("REACTION")
                if reaction_start != -1:
                    # Get the section after REACTION
                    remaining_text = ec_info[reaction_start + len("REACTION"):]
                    # Find the next section or end of file
                    next_section_match = re.search(r'\n[A-Z]+\s', remaining_text)
                    if next_section_match:
                        reaction_section = remaining_text[:next_section_match.start()]
                    else:
                        reaction_section = remaining_text
                    
                    # Parse each reaction line
                    for line in reaction_section.split("\n"):
                        line = line.strip()
                        if line and not line.startswith(' '):
                            reaction_info = parse_reaction_line(line)
                            if reaction_info:
                                reaction_data.append(reaction_info)
            
            # Also check for ALL_REAC if no reactions found in REACTION section
            if not reaction_data and "ALL_REAC" in ec_info:
                # Find ALL_REAC section
                all_reac_start = ec_info.find("ALL_REAC")
                if all_reac_start != -1:
                    remaining_text = ec_info[all_reac_start + len("ALL_REAC"):]
                    next_section_match = re.search(r'\n[A-Z]+\s', remaining_text)
                    if next_section_match:
                        all_reac_section = remaining_text[:next_section_match.start()]
                    else:
                        all_reac_section = remaining_text
                    
                    # Parse the ALL_REAC line(s)
                    for line in all_reac_section.split("\n"):
                        line = line.strip()
                        if line:
                            reactions_from_all_reac = parse_all_reac_line(line)
                            if reactions_from_all_reac:
                                reaction_data.extend(reactions_from_all_reac)
                        
        except Exception as e:
            print(f"Error retrieving information for EC {ec_number}: {e}")
    
    return reaction_data

# -------------------
# 3. Get detailed information regarding a pathway
# (including reactions information)
# -------------------

def retrieve_pathway_info(pathway_id:str, abbr:str) -> PathwayData:
    # Retrieve information from pathway
    try:
        pathway_info = REST.kegg_get(pathway_id).read()
        
        # Extract pathway name properly
        pathway_name = "Unknown"
        for line in pathway_info.split("\n"):
            if line.startswith("NAME"):
                pathway_name = line.replace("NAME", "").strip()
                # Remove trailing semicolon if present
                pathway_name = pathway_name.rstrip(';')
                break
        
        # Check if it's a global pathway
        global_path = False
        first_line = pathway_info.split("\n")[0] if pathway_info else ""
        
        if "Global" in first_line or "Overview" in first_line:
            global_path = True
            return PathwayData(
                pathway_name=pathway_name, 
                kegg_id=pathway_id, 
                global_path=global_path, 
                enzymes=[]
            )
        
        # Check if pathway contains EC numbers
        if "[EC:" not in pathway_info:
            global_path = True
            return PathwayData(
                pathway_name=pathway_name, 
                kegg_id=pathway_id, 
                global_path=global_path, 
                enzymes=[]
            )

        # Parse GENE section
        enzymes = []
        if "GENE" in pathway_info:
            # Split by GENE section and get the part after GENE
            gene_section = pathway_info.split("GENE")[1]
            
            # Split by next section header or end of file
            next_section_match = re.search(r'\n[A-Z]+\s', gene_section)
            if next_section_match:
                gene_lines = gene_section[:next_section_match.start()].split("\n")
            else:
                gene_lines = gene_section.split("\n")
            
            # Parse each gene line
            for gene_line in gene_lines:
                if gene_line.strip() and not gene_line.strip().startswith((' ', '\t')):
                    # This is a new gene entry
                    enzyme = parse_gene_line(gene_line, abbr)
                    if enzyme and enzyme.ec_numbers:
                        # Retrieve reaction information for EC numbers
                        enzyme.reaction_data = retrieve_ec_information(enzyme.ec_numbers)
                        enzymes.append(enzyme)
                elif gene_line.strip() and enzymes:
                    # This is a continuation line for the last gene
                    pass
        
        return PathwayData(
            pathway_name=pathway_name, 
            kegg_id=pathway_id, 
            global_path=global_path, 
            enzymes=enzymes
        )
    except Exception as e:
        print(f"    Error processing pathway {pathway_id}: {e}")
        # Add pathway without enzyme data
        return None

# -------------------
# Main execution
# -------------------

if __name__ == "__main__":
    # ----------------------------------
    # 1. Retrieve overall information from KEGG
    # ----------------------------------

    kegg_result = collect_kegg()

    if not kegg_result.organism_df.empty:
        print("="*50)
        print(f"\nStored data: {kegg_result}")
        print(f"Pathways DataFrame: {kegg_result.pathway_df.shape}")
        print(f"Organisms DataFrame: {kegg_result.organism_df.shape}")
        print("="*50)

        # Debugging
        # print("\nFirst few pathways:")
        # print(kegg_result.pathway_df.head())
        # print("\nFirst few organisms:")
        # print(kegg_result.organism_df.head())
        # print(kegg_result.organism_df["organism"][0])

        # ----------------------------------
        # 2. Test with a specific organism and specific pathways
        # ----------------------------------

        test_organism = kegg_result.organism_df.iloc[0]['abbreviation']
        print(f"\nTesting with organism: {test_organism}")

        organism_data = retrieve_organism_pathways(test_organism, kegg_result)
        print(f"Retrieved data: {organism_data}")
        #print(organism_data.pathways_info.head)
        for index, row in organism_data.pathways_info.iterrows():
            print(f"{row['kegg_id']} - {row['pathway_name']}")

        #print(REST.kegg_get("hsa01100").read())

        # ----------------------------------
        # 3. Test with specific pathways
        # ----------------------------------
        test_specific_pathways = [
            "hsa01100",  # Metabolic pathways
            "hsa00250",  # Alanine, aspartate and glutamate metabolism
            "hsa00010"   # Glycolysis / Gluconeogenesis
        ]

        for pathway_id in test_specific_pathways:
            print(f"\nTesting pathway: {pathway_id}")
            pathway_data = retrieve_pathway_info(pathway_id, "hsa")
            
            if pathway_data:
                print(f"Successfully retrieved {pathway_data.pathway_name}")
                print(f"\t-Global pathway: {pathway_data.global_path}")
                print(f"\t-Number of enzymes: {len(pathway_data.enzymes)}")
                
                # Count enzymes with EC numbers
                enzymes_with_ec = sum(1 for enzyme in pathway_data.enzymes if enzyme.ec_numbers)
                print(f"\t-Enzymes with EC numbers: {enzymes_with_ec}")
                
                # Show reactions if available
                total_reactions = sum(len(enzyme.reaction_data) for enzyme in pathway_data.enzymes)
                print(f"\t-Total reactions: {total_reactions}")
                
                # Print enzyme summary
                print(f"\n\t-Enzyme Summary:")
                for i, enzyme in enumerate(pathway_data.enzymes, 1):
                    reaction_count = len(enzyme.reaction_data)
                    print(f"\t  {i}. {enzyme.gene} ({enzyme.kegg_gene_id})")
                    print(f"\t     Name: {enzyme.enzyme_name}")
                    print(f"\t     EC: {enzyme.ec_numbers if enzyme.ec_numbers else 'None'}")
                    print(f"\t     KO: {enzyme.ko_numbers if enzyme.ko_numbers else 'None'}")
                    print(f"\t     Reactions: {reaction_count}")
                    
                    # Show first reaction equation if available
                    if enzyme.reaction_data:
                        first_reaction = enzyme.reaction_data[0]
                        print(f"\t     Example Reaction: {first_reaction.reaction_equation}")
                        
            else:
                print(f"\tFailed to retrieve pathway {pathway_id}")