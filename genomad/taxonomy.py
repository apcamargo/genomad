from collections import defaultdict

import taxopy
from genomad import utils


def write_taxonomic_assignment(taxonomy_output, genes_output, database_obj):
    taxdb = database_obj.get_taxdb()
    contig_taxid_dict = defaultdict(lambda: ([], []))
    for line in utils.read_file(genes_output, skip_header=True):
        gene, *_, bitscore, _, _, _, taxid, _, _, _, _, _ = line.split("\t")
        contig = gene.rsplit("_", 1)[0]
        if taxid != "1":
            bitscore, taxid = int(bitscore), int(taxid)
            contig_taxid_dict[contig][0].append(taxid)
            contig_taxid_dict[contig][1].append(bitscore)
    with open(taxonomy_output, "w") as fout:
        fout.write("seq_name\tn_genes_with_taxonomy\tagreement\ttaxid\tlineage\n")
        for contig, (taxids, bitscores) in contig_taxid_dict.items():
            taxon_list = [taxopy.Taxon(i, taxdb) for i in taxids]
            if len(taxon_list) > 1:
                majority_taxon = taxopy.find_majority_vote(
                    taxon_list, taxdb, weights=bitscores, fraction=0.5
                )
                agreement = majority_taxon.agreement
                # If classification is at the family level, be more conservative
                if (majority_taxon.rank == "family") and (agreement < 0.7):
                    majority_taxon = taxopy.find_majority_vote(
                        taxon_list, taxdb, weights=bitscores, fraction=0.7
                    )
                    agreement = majority_taxon.agreement
                # If the contig was assigned to Nucleocytoviricota but contains at least
                # one Caudoviricetes marker, be more conservative
                elif (
                    (
                        majority_taxon.rank_name_dictionary.get("phylum")
                        == "Nucleocytoviricota"
                    )
                    and (agreement < 0.6)
                    and any(
                        t.rank_name_dictionary.get("class") == "Caudoviricetes"
                        for t in taxon_list
                    )
                ):
                    majority_taxon = taxopy.find_majority_vote(
                        taxon_list, taxdb, weights=bitscores, fraction=0.6
                    )
                    agreement = majority_taxon.agreement
            else:
                majority_taxon = taxon_list[0]
                agreement = 1
                if majority_taxon.rank == "family":
                    majority_taxon = taxopy.Taxon(
                        majority_taxon.taxid_lineage[1], taxdb
                    )
            lineage = ";".join(reversed(majority_taxon.name_lineage))
            if lineage.startswith("root"):
                lineage = lineage.replace("root", "Viruses", 1)
            fout.write(
                f"{contig}\t{len(taxon_list)}\t{agreement:.4f}\t{majority_taxon.taxid}\t{lineage}\n"
            )
