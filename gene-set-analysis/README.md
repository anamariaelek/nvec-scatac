# gene-set-analysis

Functions to analyse gene sets (functional enrichment, etc)

Examples:

1. `topGO` to test for hierarchical **enrichment of gene ontologies**, using `fisher`'s test and the `elim` hierarchy-parsing algorithm (see `topGO` manual).

```R
source("gene-set-analysis.R")

# load GO annotations
go_annotations = gsa_topgo_load_emapper(emapper_fn = "test/Tadh_eggnog.Metazoa.emapper.annotations")

# get a list of genes of interest
your_list = read.table("test/your_list.txt")[,1]
head(your_list)

# enrichment test
go_enrichment_table = gsa_topgo_enrichment(
  annotation = go_annotations,
  genes_fg = your_list,
  output_prefix = "test/test_GO_enrichments",
  name_fg = "your_list",
  ontologyset = c("BP","MF","CC"),
  tg_test = "fisher",
  tg_algorithm = "elim",
  top_markers = 30,
  nodesize = 10,
  printfile = TRUE
)
```

2. Hypergeometric test to test for **enrichment of pfam domains** (or any other "flat" annotation, i.e. where each gene is mapped to a given annotation).

```R
source("gene-set-analysis.R")

# load pfam annotations
pfam_annotations = gsa_enrichment_load_pfam_list(pfam_architecture_file = "test/Tadh_long.pep.pfamscan_archs.csv")

# get a list of genes of interest
your_list = read.table("test/your_list.txt")[,1]

enrichment_table = gsa_enrichment_hypergeometric(
  annotation = pfam_annotations,
  genes_fg = your_list,
  output_prefix = "test/test_pfam_enrichments",
  name_fg = "your_list"
)
```

3. Venn and Euler diagrams, to visualise set overlaps (two or three categories). It also produces a **list of intersecting and unique terms** per set.

```R
source("gene-set-analysis.R")

# get a list of genes of interest, and get a (random) second one to compare it with
your_list = read.table("test/your_list.txt")[,1]
another_list = c(your_list [ sample(50) ], rep("randomrandom",50))

pdf("test/venn_output.pdf")
# euler diagram (size-scaled)
venn_output = venn.two(
  list1 = your_list,
  list2 = another_list,
  catname1 = "your list",
  catname2 = "a second list",
  main = "my random Euler diagram",
  eulerbool = TRUE
)

# venn diagram (all circles same size)
venn_output = venn.two(
  list1 = your_list,
  list2 = another_list,
  catname1 = "your list",
  catname2 = "a second list",
  main = "my random Venn diagram",
  eulerbool = FALSE
)
dev.off()
```
