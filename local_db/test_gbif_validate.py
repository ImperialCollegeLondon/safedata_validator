from safe_dataset_checker import *

test_tax = [("Crematogaster borneensis", "Species", None),
            ("Dolichoderus", "Genus", None),
            ("Ponerinae", "Subfamily", None),
            ("Formicidae", "Family", None),
            ("Crematogaster ormei", "Species", None),
            ("Varanus salvator macromaculatus", "Subspecies", None),
            ("Alsomitra simplex", "Species", None),
            ("Morus", "Genus", 2480962L),
            ("Zenicomus photuroides", "Species", None),
            ("Cicada sanguinolenta", "Species", None),
            ("Melittia oedippus", "Species", None),
            ("Melaphorus potteri", "Species", None),
            ("Goniopholis tenuidens", "Species", None),
            ("Solenopsis abdita", "Species", None),
            ("Biarmosuchus tagax", "Species", None),
            ("Camponotites kraussei", "Species", None)]

web_tax_check = [web_gbif_validate(tax, rnk, id) for tax, rnk, id in test_tax]

con = sqlite3.connect('backbone-current-simple.sqlite')
local_tax_check = [local_gbif_validate(con, tax, rnk, gbif_id) for tax, rnk, gbif_id in test_tax]

print [w == t for w, t in zip(web_tax_check, local_tax_check)]