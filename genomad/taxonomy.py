from collections import defaultdict

import taxopy

from genomad import utils


def get_conservative_taxon(taxon, taxdb):
    for rank, taxid in taxon.ranked_taxid_lineage:
        if rank not in {"subfamily", "genus", "subgenus", "species"}:
            return taxopy.Taxon(taxid, taxdb)
    return taxon


def write_taxonomic_assignment(
    taxonomy_output,
    genes_output,
    database_obj,
    lenient_taxonomy=False,
    full_ictv_lineage=False,
):
    if full_ictv_lineage:
        output_ranks = [
            "realm",
            "subrealm",
            "kingdom",
            "subkingdom",
            "phylum",
            "subphylum",
            "class",
            "subclass",
            "order",
            "suborder",
            "family",
        ]
        if lenient_taxonomy:
            output_ranks += ["subfamily", "genus", "subgenus", "species"]
    else:
        output_ranks = ["realm", "kingdom", "phylum", "class", "order", "family"]
        if lenient_taxonomy:
            output_ranks += ["genus", "species"]
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
                # If the contig was assigned to Nucleocytoviricota but contains
                # at least one Caudoviricetes marker, be more conservative
                if (
                    majority_taxon.rank_name_dictionary.get("phylum")
                    == "Nucleocytoviricota"
                    and agreement < 0.6
                    and any(
                        t.rank_name_dictionary.get("class") == "Caudoviricetes"
                        for t in taxon_list
                    )
                ):
                    majority_taxon = taxopy.find_majority_vote(
                        taxon_list, taxdb, weights=bitscores, fraction=0.6
                    )
                    agreement = majority_taxon.agreement
                # If classification is below family level and `lenient_taxonomy`
                # is False,
                if not lenient_taxonomy and majority_taxon.rank in {
                    "subfamily",
                    "genus",
                    "subgenus",
                    "species",
                }:
                    majority_taxon = get_conservative_taxon(majority_taxon, taxdb)
                    agreement = 0
                    for t, w in zip(taxon_list, bitscores):
                        if (
                            t.rank_taxid_dictionary.get(majority_taxon.rank)
                            == majority_taxon.taxid
                        ):
                            agreement += w / sum(bitscores)
            else:
                majority_taxon = taxon_list[0]
                agreement = 1
                if not lenient_taxonomy and majority_taxon.rank in {
                    "subfamily",
                    "genus",
                    "subgenus",
                    "species",
                }:
                    majority_taxon = get_conservative_taxon(majority_taxon, taxdb)
            lineage = [
                majority_taxon.rank_name_dictionary.get(i, "") for i in output_ranks
            ]
            lineage = ";".join(["Viruses"] + lineage)
            fout.write(
                f"{contig}\t{len(taxon_list)}\t{agreement:.4f}\t"
                f"{majority_taxon.taxid}\t{lineage}\n"
            )
