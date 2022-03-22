"""## The genb_taxa submodule
This module describes classes and methods used to validate taxonomic data against
the NCBI database and convert it to a form compatible with the GBIF backbone database.

The NCBITaxon dataclass is used to store data about a taxon entry in NCBI (GenBank)
form. It is initialised with user data and then the taxon Validator class can
be used to update a NCBITaxon object with the result of NCBI validation. At present
only online validation is possible through the `RemoteNCBIValidator` class.

FURTHER DETAILS OF STUFF DEFINED.

THIS IS STILL VERY MUCH A WORK IN PROGRESS SO MORE INFORMATION IS GOING TO BE
ENTERED HERE
"""

from typing import Union, Optional
import dataclasses
from Bio import Entrez

import requests
import urllib
from enforce_typing import enforce_types

from safedata_validator.logger import (COUNTER_HANDLER, FORMATTER, LOGGER,
                                       loggerinfo_push_pop)
from safedata_validator.validators import (GetDataFrame, HasDuplicates,
                                           IsLower, blank_value)

# ADD TO RESOURCE FILE, AS A CHECK THAT USER HAS PROVIDED ONE
# Key should also be added to the resource file
Entrez.email = "jacobcook1995@gmail.com"
# Hard coding api key in for now
user_key = "1738fe86eba2d8fc287ff0d1dcbfeda44a0a"


# TODO - Correctly read in xslx data
# Lot of this has already been written by David, main thing is that I will have
# to decide on the formatting, like what data are we taking in.

# TODO - Think about data output
# Have to make sure that the indexing is compatibale with David's database
# Probably have to add in the first databasing steps as well

# POTENTIAL - Take steps to increase the speed
# Make local copy by downloading the relevant part of NCBI's database and
# running the validation locally

# TODO - Modify the resource file to ask the user to provide an email address
# This should only be done if the user actually wants to use this module as it
# isn't need elsewhere (as far as I know)
# Also should ask for api key

# QUESTIONS FOR DAVID
# SHOULD A YOU ARE NOT CONNCTED TO THE INTERNET ERROR BE SETUP?
# IS THERE A SENSIBLE WAY TO DUMMY CONNECTION ERRORS?

# Extended version of backbone ranks to capture superkingdoms
BACKBONE_RANKS_EX = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                    'family', 'genus', 'species', 'subspecies']

class NCBIError(Exception):
    """Exception class for remote NCBI errors

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message="NCBI taxa ID not found"):
        self.message = message
        super().__init__(self.message)

# New function to remove k__ notation, and to check that rank and
def taxa_strip(name: str, rank: str):
    """A function to remove k__ type notation from taxa names. The function also
    checks that the provided rank is consistent with the rank implied by prefix
    """
    if name == None:
        return (None, True)
    elif "__" in name:
        # Strip name down
        ind = name.rfind("_")
        s_name = name[ind+1:]
        # Check if ranks match
        match = (name[0].lower() == rank[0].lower())
        return (s_name, match)
    else:
        return (name, True)

# New function to generate species binomial
def species_binomial(genus: str, species: str):
    # First check if species is a single name
    if len(species.split()) == 1 and len(genus.split()) == 1:
        return genus.strip() + " " + species.strip()
    # Look for Candidatus formats
    elif "candidatus" in species.lower() or "candidatus" in genus.lower():
        if "candidatus" in species.lower():
            # Construct binomal with first word of species name removed
            bi = genus.strip()
            for i in range(1,len(species.split())):
                bi = bi + " " + species.split()[i]
            return bi
        else:
            return genus.strip() + " " + species.strip()
    # Then check that species name is more words than the genus name
    elif len(species.split()) > len(genus.split()):
        if genus in species:
            return species
        else:
            LOGGER.error(f'Species name ({species}) appears to be binomal but '
                         f'does not contain genus name ({genus})')
            return None
    else:
        LOGGER.error(f'Genus name ({genus}) appears to be too long')
        return None

# Equivalent function to generate subspecies trinominal
def subspecies_trinomial(species: str, subspecies: str):
    # First check if subspecies is a single name
    if len(subspecies.split()) == 1 and len(species.split()) == 2:
        return species.strip() + " " + subspecies.strip()
    # Look for Candidatus formats
    elif "candidatus" in subspecies.lower() or "candidatus" in species.lower():
        if "candidatus" in subspecies.lower():
            # Construct trinominal with first word of subspecies name removed
            tri = species.strip()
            for i in range(1,len(subspecies.split())):
                tri = tri + " " + subspecies.split()[i]
            return tri
        else:
            return species.strip() + " " + subspecies.strip()
    # Catch too shot species name case
    elif len(species.split()) == 1:
        LOGGER.error(f'Species name ({species}) too short')
        return None
    # Then check that species name is more words than the genus name
    elif len(subspecies.split()) > len(species.split()):
        if species in subspecies:
            return subspecies
        else:
            LOGGER.error(f'Subspecies name ({subspecies}) appears to be trinomal'
                         f'but does not contain species name ({species})')
            return None
    else:
        LOGGER.error(f'Species name ({species}) too long')
        return None


@enforce_types
@dataclasses.dataclass
class NCBITaxon:
    """Holds taxonomic information from a user on a microbial taxa. This can be
    populated using NCBI validation.

    There are 3 class properties that can be used to create an instance:
        * name
        * taxa_hier: Dictionary of valid taxonomic hierarchy (with ID's)
        * ncbi_id: NCBI ID for full taxa (i.e. including non-backbone ranks)
    The remaining properties are populated by processing functions not when an
    instance is created.
        * superseed: is supplied taxon name/ID still the accepted usage
        * orig: Records rank of the orginal taxa if it was not found
    """

    # Init properties
    name: str
    taxa_hier: dict
    ncbi_id: Optional[Union[int, float]] = None
    superseed: bool = dataclasses.field(init=False)
    orig: str = dataclasses.field(init=False)

    def __post_init__(self):
        """Sets the defaults for the post-init properties and checks inputs
        """

        if self.ncbi_id is not None:
            if isinstance(self.ncbi_id, float) and not self.ncbi_id.is_integer():
                raise TypeError('NCBI Id is not an integer')
            self.ncbi_id = int(self.ncbi_id)

        if self.taxa_hier is not None:
            if len(self.taxa_hier) == 0:
                raise ValueError('Taxa hierarchy dictonary empty')
            elif all(isinstance(x,str) for x in self.taxa_hier.keys()) == False:
                raise ValueError('Not all taxa dictionary keys are strings')
            elif all(isinstance(x,tuple) for x in self.taxa_hier.values()) == False:
                raise ValueError('Not all taxa dictionary values are tuples')
            elif all(list(map(type,x)) == [str, int] for x in
                 self.taxa_hier.values()) == False:
                 raise ValueError('Taxa tuples not all in [string integer] form')

        self.superseed = False
        self.orig = None

@enforce_types
class RemoteNCBIValidator:
    """This provides a validate method for a MicrobeTaxon using various online
    NCBI APIs. This doesn't need an __init__ method and just contains methods.
    """
    # Functionality to find taxa information from genbank ID
    def id_lookup(self, nnme: str, ncbi_id: int):
        """Method to return full taxonomic information from a NCBI ID. This
        includes details of any potential synonymus names.

        Params:
            nnme: A nickname to identify the taxon
            ncbi_id: An integer

        Returns:
            A NCBITaxon object
        """

        if not isinstance(ncbi_id, int):
            raise ValueError('Non-integer NCBI taxonomy ID')

        if not ncbi_id > 0:
            raise ValueError('Negative NCBI taxonomy ID')

        # Use efetch to find taxonomy details based on the index
        try:
            handle = Entrez.efetch(db="taxonomy",id=f"{ncbi_id}",
                                   retmode="xml",api_key=user_key)
            tax_dic = Entrez.read(handle)[0]
            handle.close()
        # And case where a connection can't be made to the remote server
        except urllib.error.HTTPError:
            raise NCBIError('Connection error to remote server')
        # Catch case where the index isn't found
        except IndexError:
            raise NCBIError("No entry found from ID, probably bad ID")

        # Then extract the lineage
        linx = tax_dic["LineageEx"]
        # Find number of taxonomic ranks
        tx_len = len(linx)
        # Extract parent taxa ranks
        rnks = [linx[i]["Rank"] for i in range(tx_len)]

        # Check that the taxon rank provided is a backbone rank
        if tax_dic["Rank"] in BACKBONE_RANKS_EX:
            # In this case use provided rank
            rnk = BACKBONE_RANKS_EX.index(tax_dic["Rank"])

        # Filter out ID's with non-backbone ranks (e.g. strains)
        else:
            # Set as not a valid taxa
            vld_tax = False
            # Find lowest index
            r_ID = len(BACKBONE_RANKS_EX) - 1

            # While loop that runs until valid taxa is found
            while vld_tax == False:
                # Check if taxa id is found
                if any([rnks[i] == f"{BACKBONE_RANKS_EX[r_ID]}" for i in range(tx_len)]):
                    # Close loop and store rank number
                    vld_tax = True
                    # Add 1 to the rank as only including lineage in this case
                    rnk = r_ID + 1
                    # Warn user that non-backbone rank has been supplied
                    LOGGER.warning(f'{nnme} of non-backbone rank: {tax_dic["Rank"]}')
                # Raise error once backbone ranks have been exhausted
                elif r_ID < 1:
                    LOGGER.error(f'Taxon hierarchy for {nnme} contains no backbone ranks')
                    return
                else:
                    r_ID -= 1

        # Make list of backbone ranks we are looking for
        actual_bb_rnks = BACKBONE_RANKS_EX[0:rnk]

        # Number of missing ranks initialised to zero
        m_rnk = 0

        # Check that all desired backbone ranks appear in the lineage
        if all(item in rnks for item in actual_bb_rnks) == False:
            # Find all missing ranks
            miss = list(set(actual_bb_rnks).difference(rnks))
            # Count missing ranks
            m_rnk = len(miss)
            # Remove missing ranks from our list of desired ranks
            for i in range(0,m_rnk):
                actual_bb_rnks.remove(miss[i])

        # Find valid indices (e.g. those corresponding to backbone ranks)
        vinds = [idx for idx, element in enumerate(rnks) if element in actual_bb_rnks]

        # Create dictonary of valid taxa lineage using a list
        red_taxa = {f"{actual_bb_rnks[i]}":(str(linx[vinds[i]]["ScientificName"]),
                    int(linx[vinds[i]]["TaxId"])) for i in range(0,rnk-m_rnk)}

        # Then add taxa information as a final entry
        red_taxa[f"{tax_dic['Rank']}"] = (str(tax_dic["ScientificName"]),
                                                 int(tax_dic["TaxId"]))

        # Create and populate microbial taxon
        mtaxon = NCBITaxon(name=nnme,ncbi_id=ncbi_id,taxa_hier=red_taxa)

        # Check if AkaTaxIds exists in taxonomic information
        if 'AkaTaxIds' in tax_dic.keys():
            # Warn user that they've given a superseeded taxa ID
            LOGGER.warning(f"NCBI ID {(tax_dic['AkaTaxIds'])[0]} has been "
                            f"superseeded by ID {tax_dic['TaxId']}")
            # Record that a superseeded NCBI ID has been provided
            mtaxon.superseed = True

        return mtaxon

    # New function to read in taxa information
    def taxa_search(self, nnme: str, taxah: dict):
        """Method that takes in taxonomic information, and finds the corresponding
        genbank ID. This genbank ID is then used to generate a NCBITaxon object,
        which is returned. This function also makes use of parent taxa information
        to distinguish between ambigious taxa names.

        Params:
            taxah: A dictonary containing taxonomic information

        Returns:
            A NCBITaxon object
        """

        if isinstance(taxah, dict) == False:
            raise TypeError('Taxa hierarchy should be a dictonary')
        elif all(isinstance(x,str) for x in taxah.keys()) == False:
            raise ValueError('Not all taxa dictionary keys are strings')
        elif all(isinstance(x,str) for x in taxah.values()) == False:
            raise ValueError('Not all taxa dictionary values are strings')

        # Warn user if no backbone ranks have been provided
        if bool(set(taxah.keys()) & set(BACKBONE_RANKS_EX)) == False:
            LOGGER.warning(f"No backbone ranks provided in {nnme}'s taxa hierarchy")

        # Find last dictonary key
        f_key = list(taxah.keys())[-1]

        # Then find corresponding entry as a search term
        s_term = taxah[f_key]

        # Search the online database
        try:
            handle = Entrez.esearch(db="taxonomy", term=s_term, api_key=user_key)
        except urllib.error.HTTPError:
            NCBIError('Connection error to remote server')

        # If it works then save the response
        record = Entrez.read(handle)
        handle.close()

        # Store count of the number of records found
        c = int(record['Count'])

        # Check that a singular record has been provided before proceeding
        if c == 1:
            # Find taxa ID as single entry in the list
            tID = int(record['IdList'][0])
            # Use ID lookup function to generate as a NCBITaxon object
            mtaxon = self.id_lookup(nnme,tID)
        # Catch cases where no record is found
        elif c == 0:
            # Check if there actually is any higher taxonomy provided
            if len(taxah.keys()) == 1:
                LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                             f'taxonomic hierarchy is absent')
                return

            # If there is then set up a loop over it
            fnshd = False
            cnt = 1
            while fnshd == False:
                # Increment counter
                cnt += 1
                # Use to find higher taxonomic level to use in search
                f_key = list(taxah.keys())[-cnt]
                new_s_term = taxah[f_key]

                # Search the online database
                try:
                    handle = Entrez.esearch(db="taxonomy", term=new_s_term, api_key=user_key)
                except urllib.error.HTTPError:
                    NCBIError('Connection error to remote server')

                # Process the response
                record = Entrez.read(handle)
                handle.close()
                c = int(record['Count'])

                # Check how many records have been found
                if c == 1:
                    # Find taxa ID as single entry in the list
                    tID = int(record['IdList'][0])
                    # Use ID lookup function to generate as a NCBITaxon object
                    mtaxon = self.id_lookup(nnme,tID)
                    # Store orginally supplied rank
                    mtaxon.orig = list(taxah.keys())[-1]
                    # Warn the user that a higher taxonomic rank is being used
                    LOGGER.warning(f'{s_term} not registered with NCBI, using '
                                   f'higher level taxon {new_s_term} instead')
                    fnshd = True
                elif c > 1: # Not going to handle ambiguities in this case
                    LOGGER.error(f'Taxa {nnme} cannot be found and its higher '
                                 f'taxonomic hierarchy is ambigious')
                    return
                # Catch when all the provided hierachy has been exhausted
                elif cnt == len(taxah.keys()):
                    fnshd = True

                # If valid higher taxon not found print error
                if 'mtaxon' not in locals():
                    LOGGER.error(f'Taxa {nnme} cannot be found and neither can '
                                 f'its higher taxonomic hierarchy')
                    return

        # Case where only one rank has been provided
        elif len(taxah) == 1:
            # Preallocate container for the ranks
            t_ranks = []
            # First check if multiple taxa have the same rank
            for i in range(c):
                # Find taxa ID as single entry in the list
                tID = int(record['IdList'][i])
                # Use ID lookup function to generate a tempoary taxon
                temp_taxon = self.id_lookup(nnme,tID)
                # Add rank to list
                t_ranks.append(list(temp_taxon.taxa_hier.keys())[-1])

            # Record whether ranks match expected rank
            if f_key == "kingdom":
                mtch = [True if x == f_key or x == "superkingdom" else False
                        for x in t_ranks]
            else:
                mtch = [True if x == f_key else False for x in t_ranks]

            # If there's only one match, we have determined the correct entry
            if sum(mtch) == 1:
                # Find relevant index
                ind = ([i for i, x in enumerate(mtch) if x])[0]
                # Find taxa ID as single entry in the list
                tID = int(record['IdList'][ind])
                # Use ID lookup function to generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)
            # Then check whether multiple taxonomic levels have been provided
            else:
                # If not raise an error
                LOGGER.error(f'Taxa {nnme} cannot be found using only one '
                             f'taxonomic level, more should be provided')
                return
        # Higher ranks provided
        else:
            # Find second from last dictonary key
            f_key = list(taxah.keys())[-2]
            # Then find corresponding entry as a search term
            s_term = taxah[f_key]

            # Then find details of this parent record
            try:
                handle = Entrez.esearch(db="taxonomy", term=s_term,
                                        api_key=user_key)
            except urllib.error.HTTPError:
                NCBIError('Connection error to remote server')

            p_record = Entrez.read(handle)
            handle.close()

            # Store count of the number of records found
            pc = int(p_record['Count'])

            # Check that single parent taxa exists in records
            if pc == 0:
                LOGGER.error(f'Provided parent taxa for {nnme} not found')
                return
            elif pc > 1:
                LOGGER.error(f'More than one possible parent taxa for {nnme} found')
                return
            else:
                # Find parent taxa ID as single entry in the list
                pID = int(p_record['IdList'][0])
                # Then use ID lookup function to find generate as a NCBITaxon object
                ptaxon = self.id_lookup("parent",pID)

            # Save parent taxa rank and name
            p_key = list(ptaxon.taxa_hier.keys())[-1]
            p_val = ptaxon.taxa_hier[p_key]

            # Use list comprehension to make list of potential taxa
            potents = [self.id_lookup(f"c{i}",int(record['IdList'][i]))
                        for i in range(0,c)]

            # Check if relevant rank exists in the child taxa
            child = [p_key in potents[i].taxa_hier for i in range(0,c)]

            # Then look for matching rank name pairs
            for i in range(0,c):
                # Only check cases where rank exists in the child taxa
                if child[i] == True:
                    child[i] = (potents[i].taxa_hier[f"{p_key}"] == p_val)

            # Check for errors relating to finding too many or few child taxa
            if sum(child) == 0:
                LOGGER.error(f'Parent taxa not actually a valid parent of {nnme}')
                return
            elif sum(child) > 1:
                LOGGER.error(f'Parent taxa for {nnme} refers to multiple '
                             f'possible child taxa')
                return
            else:
                # Find index corresponding to correct child taxa
                tID = int(record['IdList'][child.index(True)])
                # Use ID lookup function to find generate as a NCBITaxon object
                mtaxon = self.id_lookup(nnme,tID)

        # Find last dictonary key
        f_key = list(mtaxon.taxa_hier.keys())[-1]

        # Check if taxonomic rank supplied is used
        if mtaxon.orig == None and f_key != list(taxah.keys())[-1] and f_key != "superkingdom":
            # If not raise an error
            LOGGER.error(f'{list(taxah.values())[-1]} is a {f_key}'
                         f' not a {list(taxah.keys())[-1]}')
            return
        elif list(taxah.keys())[-1] == "kingdom" and f_key == "superkingdom":
            # If not print a warning
            LOGGER.warning(f'NCBI records {(mtaxon.taxa_hier[f_key])[0]} as '
                           f'a superkingdom rather than a kingdom')
        # Then check whether orginally supplied name is still used
        elif mtaxon.orig == None and taxah[f_key] != (mtaxon.taxa_hier[f_key])[0]:
            # If not print a warning
            LOGGER.warning(f'{taxah[f_key]} not accepted usage should be '
                           f'{(mtaxon.taxa_hier[f_key])[0]} instead')
            # And record the fact that usage is superseeded
            mtaxon.superseed = True

        return mtaxon

class NCBITaxa:
    # REWRITE THIS HEAVILY
    """A class to hold a list of taxon names and a validated taxonomic
    index for those taxa and their taxonomic hierarchy. The validate_taxon
    method checks that taxon details and their optional parent taxon can be
    matched into the the GBIF backbone and populates two things:

    i)  the taxon_names attribute of the dataset, which is just a set of
        names used as a validation list for taxon names used in data worksheets.
    ii) the taxon_index attribute of the dataset, which contains a set
        of lists structured as:

            [worksheet_name (str),
            gbif_id (int),
            gbif_parent_id (int),
            canonical_name (str),
            taxonomic_rank (str),
            status (str)]

        Where a taxon is not accepted or doubtful on GBIF, two entries are
        inserted for the taxon, one under the canon name and one under the
        provided name. They will share the same worksheet name and so can
        be paired back up for description generation. The worksheet name
        for parent taxa and deeper taxonomic hierarchy is set to None.

    The index_higher_taxa method can be used to extend the taxon_index to
    include all of the higher taxa linking the validated taxa.

    The index can then be used:

    a) to generate the taxonomic coverage section of the dataset description, and
    b) as the basis of a SQL dataset_taxa table to index the taxonomic coverage
    of datasets.
    """

    def __init__(self):
        # FILL OUT THIS PROPERLY ONCE I KNOW WHAT IT IS SUPPOSED TO BE
        """Sets the initial properties of a NCBITaxa object
        """

        self.taxon_index = []
        self.taxon_names = set()
        self.ncbi_t = dict()
        self.hierarchy = set()
        self.n_errors = None
        self.taxon_names_used = set()

        # At the moment we only have one validator defined
        self.validator = RemoteNCBIValidator()

    def load(self, worksheet):
        # THIS PROBABLY NEEDS A BIT OF REWRITING
        """Loads a set of taxa from the rows of a SAFE formatted NCBITaxa worksheet
        and then adds the higher taxa for those rows.

        Args:
            worksheet: An openpyxl worksheet instance following the NCBITaxa formatting

        Returns:
            Updates the taxon_names and taxon_index attributes of the class instance
            using the data in the worksheet.
        """

        start_errors = COUNTER_HANDLER.counters['ERROR']

        # Get the data read in.
        LOGGER.info("Reading NCBI taxa data")
        FORMATTER.push()
        dframe = GetDataFrame(worksheet)

        if not dframe.data_columns:
            LOGGER.error('No data or only headers in Taxa worksheet')
            return

        # Dupe headers likely cause serious issues, so stop
        if 'duplicated' in dframe.bad_headers:
            LOGGER.error('Cannot parse taxa with duplicated headers')
            return

        # Get the headers
        headers = IsLower(dframe.headers).values

        # Field cleaning
        core_fields = {'name', 'ncbi id'}
        missing_core = core_fields.difference(headers)

        if missing_core:
            # core names are not found so can't continue
            LOGGER.error('Missing core fields: ', extra={'join': missing_core})
            return

        # Check that at least two backbone taxa have been provided
        fnd_rnks = set(BACKBONE_RANKS_EX).intersection(headers)

        if len(fnd_rnks) < 2:
            # can't continue if less than two backbone ranks are provided
            LOGGER.error('Less than two backbone taxonomic ranks are provided')
            return

        # Find core field indices and use to isolate non core (i.e. taxonomic) headers
        core_inds = [headers.index(item) for item in core_fields]
        non_core = [element for i, element in enumerate(headers) if i not in core_inds]

        # Check that backbone ranks are in the correct order
        l1 = [v for v in BACKBONE_RANKS_EX if v in fnd_rnks]
        l2 = [v for v in non_core if v in fnd_rnks]

        if l1 != l2:
            LOGGER.error('Backbone taxonomic ranks not provided in the correct order')

        # Check that if subspecies has been provided species is provided
        if 'subspecies' in headers:
            if 'species' not in headers:
                LOGGER.error('If subspecies is provided so must species')
                return

        # Same check
        if 'species' in headers:
            if 'genus' not in headers:
                LOGGER.error('If species is provided so must genus')
                return

        # Any duplication in names
        dupl_taxon_names = HasDuplicates(dframe.data_columns[headers.index('name')])

        if dupl_taxon_names:
            LOGGER.error('Duplicated names found: ',
                         extra={'join': dupl_taxon_names.duplicated})

        # get dictionaries of the taxa
        taxa = [dict(zip(headers, rw)) for rw in zip(*dframe.data_columns)]
        FORMATTER.pop()

        # check number of taxa found
        if len(taxa) == 0:
            LOGGER.info('No taxon rows found')
            return

        # Change data format such that it is appropriate to be used to search
        # the NCBI database, i.e. convert from k__ notation, generate species
        # binomials and subspecies trinominals, this is then all collected in a
        # dictonary.

        for idx, row in enumerate(taxa):

            # Strip k__ notation so that the text is appropriate for the search
            for rnk in non_core:
                row[rnk], match = taxa_strip(row[rnk],rnk)
                if match == False:
                    print(f"Implied rank of {row[rnk]} in row {idx + 1} does not"
                          f" match rank it is assigned")

            # Standardise blank values to None
            row = {ky: None if blank_value(vl) else vl for ky, vl in row.items()}
            # Replace any NA values with None
            row = {ky: None if vl == "NA" else vl for ky, vl in row.items()}

            # Start with empty dictonary for taxonomic hierarchy
            taxa_hier = {}

            # Loop over all ranks to populate the dictonary
            for rnk in non_core:
                if row[rnk] != None:
                    if rnk == "species":
                        taxa_hier[rnk] = species_binomial(row["genus"], row[rnk])
                    elif rnk == "subspecies":
                        taxa_hier[rnk] = subspecies_trinomial(taxa_hier["species"], row[rnk])
                    else:
                        taxa_hier[rnk] = row[rnk]

            self.taxon_names.update([row['name']])
            LOGGER.info(f"Validating row {idx + 1}: {row['name']}")
            FORMATTER.push()
            self.validate_and_add_taxon((row['name'], taxa_hier, row['ncbi id']))
            FORMATTER.pop()

        # Add the higher taxa
        # self.index_higher_taxa()

        # summary of processing
        self.n_errors = COUNTER_HANDLER.counters['ERROR'] - start_errors
        if self.n_errors > 0:
            LOGGER.info('Taxa contains {} errors'.format(self.n_errors))
        else:
            LOGGER.info('{} taxa loaded correctly'.format(len(self.taxon_names)))

        FORMATTER.pop()


    def validate_and_add_taxon(self, gb_taxon_input):
        """ User information is provided that names a taxon and (optionally) gives
        a NCBI taxonomy ID. This information is then used to find the closest GBIF
        backbone compatible entry in the NCBI database. This information along with
        the NCBI ID for the orginal entry is then saved, along with a list of possible
        synoyms. All this information is then used to find the closest taxa match
        in the GBIF database. This information is then used to update the relevant
        NCBITaxa instance.

        The gb_taxon_input has the form:

        ['m_name', 'taxon_hier', 'ncbi_id']

        If there is no NCBI ID, the structure is:
        ['m_name', 'taxon_hier', None]

        Args:
            taxon_input: Taxon information in standard form as above

        Returns:
            IS THE BELOW CORRECT FOR MY CASE???
            Updates the taxon_names and taxon_index attributes of the class instance.
        """

        m_name, taxon_hier, ncbi_id = gb_taxon_input

        # Sanitise worksheet names for taxa - only keep unpadded strings.
        if m_name is None or not isinstance(m_name, str) or m_name.isspace() or not m_name:
            LOGGER.error('Worksheet name missing, whitespace only or not text')
            return
        elif m_name != m_name.strip():
            LOGGER.error(f"Worksheet name has whitespace padding: {repr(m_name)}")
            m_name = m_name.strip()
            self.taxon_names.add(m_name)
        else:
            self.taxon_names.add(m_name)

        # Check that ID provided is reasonable (if provided)
        i_fail = False

        # ID can be None or an integer (openpyxl loads all values as float)
        if not(ncbi_id is None or
               (isinstance(ncbi_id, float) and ncbi_id.is_integer()) or
               isinstance(ncbi_id, int)) :
            LOGGER.error('NCBI ID contains value that is not an integer')
            i_fail = True

        # Check the main taxon details
        h_fail = False

        # Check that a dictonary with at least one entry has been provided
        if not isinstance(taxon_hier,dict) or len(taxon_hier) == 0:
            LOGGER.error('Taxa hierarchy should be a (not empty) dictonary')
            h_fail = True
        # Otherwise check for padding of dictonary keys and values
        else:
            # Make a translation table
            translate = {}
            # Loop over all dictonary keys
            for idx in taxon_hier.keys():
                if not isinstance(idx,str):
                    LOGGER.error(f"Non-string dictonary key used: {repr(idx)}")
                    h_fail = True
                elif idx == "" or idx.isspace():
                    LOGGER.error('Empty dictonary key used')
                    h_fail = True
                elif idx != idx.strip():
                    LOGGER.error(f"Dictonary key has whitespace padding: {repr(idx)}")
                    # Save keys to swap to new translation table
                    translate[idx] = idx.strip()

                # Extract corresponding dictonary value
                val = taxon_hier[idx]
                # Then perform similar checks on dictonary value
                if not isinstance(val,str):
                    LOGGER.error(f"Non-string dictonary value used: {repr(val)}")
                    h_fail = True
                elif val == "" or val.isspace():
                    LOGGER.error('Empty dictonary value used')
                    h_fail = True
                elif val != val.strip():
                    LOGGER.error(f"Dictonary value has whitespace padding: {repr(val)}")
                    taxon_hier[idx] = val.strip()

            # Use translation table to replace whitespaced dictonary keys
            for old, new in translate.items():
                taxon_hier[new] = taxon_hier.pop(old)

        # Now check that the taxa hierarchy is correctly ordered
        if i_fail:
            LOGGER.error(f'Improper NCBI ID provided, cannot be validated')

        if h_fail:
            LOGGER.error(f'Taxon details not properly formatted, cannot validate')

        if h_fail or i_fail:
            return

        # Now check that dictonary containing taxa hierarchy is properly ordered
        if len(taxon_hier) > 1: # Only matters if it contains multiple entries
            # Find all keys in taxa hierarchy
            t_ord = list(taxon_hier.keys())
            # Find corresponding keys from backbone
            b_ord = [x for x in BACKBONE_RANKS_EX if x in t_ord]
            # Remove non-backbone keys from taxa order
            t_ord = [x for x in t_ord if x in b_ord]
            # Then catch cases where orders don't match
            if b_ord != t_ord:
                LOGGER.error(f'Taxon hierarchy not in correct order')
                return

        # Now that inputs are sanitised, continue with checking...
        # First gather the info needed to index the entry
        if ncbi_id != None:
            ncbi_info = [list(taxon_hier.values())[-1], list(taxon_hier.keys())[-1],
                         int(ncbi_id)]
        else:
            ncbi_info = [list(taxon_hier.values())[-1], list(taxon_hier.keys())[-1],
                         None]

        # Then go straight ahead and search for the taxon
        hr_taxon = self.validator.taxa_search(m_name, taxon_hier)
        # Catch case where errors are returned rather than a taxon
        if hr_taxon == None:
            LOGGER.error(f'Search based on taxon hierarchy failed')
            return

        # Then check if a genbank ID number has been provided
        if ncbi_id != None:
            id_taxon = self.validator.id_lookup(m_name, int(ncbi_id))
            # Check whether this matches what was found using hierarchy
            if id_taxon != hr_taxon:
                # Check if taxonomy hierarchy superseeded
                if hr_taxon.superseed == True:
                    LOGGER.warning(f'Taxonomic classification superseeded for '
                    f'{m_name}, using new taxonomic classification')
                elif id_taxon.superseed == True:
                    LOGGER.warning(f'NCBI taxa ID superseeded for {m_name}'
                    f', using new taxa ID')
                else:
                    LOGGER.error(f"The NCBI ID supplied for {m_name} does "
                    f"not match hierarchy: expected {hr_taxon.ncbi_id}"
                    f" got {ncbi_id}")
                    return
        else:
            # Warn user that superseeded taxonomy won't be used
            if hr_taxon.superseed:
                LOGGER.warning(f'Taxonomic classification superseeded for '
                f'{m_name}, using new taxonomic classification')

        # Check that the hierachy found matches
        match = self.compare_hier(hr_taxon,taxon_hier)

        # Store the NCBI taxon keyed by NCBI information (needs tuple)
        f_key = list(hr_taxon.taxa_hier.keys())[-1]
        # THIS INDEX SHOULD INCLUDE SUPERSEEDED TAXA AS A SPECIAL CASE
        # PARENT ID SHOULD BE INCLUDED TO ALLOW FOR PROPER INDEXING
        # ALSO ADD IF IT DIVERGES FROM ORGINAL TAXA
        self.taxon_index.append([hr_taxon.name, hr_taxon.ncbi_id,
                                    hr_taxon.taxa_hier[f_key], f_key,
                                    hr_taxon.superseed])
        self.hierarchy.update(list(hr_taxon.taxa_hier.items()))

        # Check if this has succeded without warnings or errors
        if hr_taxon.superseed == False:
            # Straight forward when there isn't a genbank id, or previously processed
            if ncbi_id == None:
                # If so inform the user of this
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')
            # Otherwise need to check for superseeded ID's
            elif id_taxon.superseed == False:
                LOGGER.info(f'Taxon ({m_name}) found in NCBI database')

    @loggerinfo_push_pop('Indexing taxonomic hierarchy')
    def index_higher_taxa(self):

        # Use the taxon hierarchy entries to add higher taxa
        # - drop taxa with a GBIF ID already in the index

        # If usage is superseed or the taxon is not included in known
        known = [tx[1] for tx in self.taxon_index if tx[5] != True]
        to_add = [tx for tx in self.hierarchy if tx[1] not in known]
        to_add.sort(key=lambda val: BACKBONE_RANKS_EX.index(val[0]))

        # Look up the taxonomic hierarchy
        for tx_lev, tx_id in to_add:
            # REMOVE LOOK UP WHEN I REWRITE THIS FUNCTION
            higher_taxon = self.validator.id_lookup(tx_id)
            self.taxon_index.append([None,
                                     higher_taxon.gbif_id,
                                     higher_taxon.parent_id,
                                     higher_taxon.name,
                                     higher_taxon.rank,
                                     higher_taxon.taxon_status])
            LOGGER.info(f'Added {tx_lev} {higher_taxon}')

    def compare_hier(self, mtaxon: NCBITaxon, taxon_hier: dict):
        """ Function to compare the hierachy of a taxon with the hierachy that was
        initially supplied. This function only checks that provided information
        matches, missing levels or entries into levels are ignored"""
        # Not worth checking hierachy in superseeded case
        if mtaxon.superseed == True:
            return

        # Find all common ranks between the hierarchies
        rnks = list(set(taxon_hier.keys()) & set(mtaxon.taxa_hier.keys()))

        # Find if all taxa names match
        t_match = [taxon_hier[r] == (mtaxon.taxa_hier[r])[0] for r in rnks]
        if all(t_match) == False:
            # Find indices of non-matching values
            inds = [i for i in range(len(t_match)) if t_match[i] == False]
            # Then loop over these indices giving appropriate warnings
            for ind in inds:
                LOGGER.warning(f'Hierarchy mismatch for {mtaxon.name} its {rnks[ind]}'
                               f' should be {(mtaxon.taxa_hier[rnks[ind]])[0]} not '
                               f'{taxon_hier[rnks[ind]]}')

        # CHECK THAT KINGDOM VS SUPERKINGDOM STORAGE IS WORKING OKAY
        # MOSTLY DONE, BUT STILL NEED TO KEEP IT IN MIND
        # SUPERSEEDED CASE GOING TO BE TRICKY WHEN THERE ARE MULTIPLE AKA IDS
