# Change log since 2017/05
#### Updated by **Bofei Zhang** (zhang.5675@osu.edu)

## Build History
* 2017/11/06 The first correct version of RaMP was built. All source IDs from HMDB, KEGG, WikiPathways, and Reactome are successfully mapped with RaMP IDs with reduced duplicated items.


Table Names | Number of Items
-----------|--------------
analyte | 90493
analytehasontology | 180464
analytehaspathway | 201592
analytesynonym | 357662
catalyzed | 861526
ontology | 160
pathway | 3420
source | 158879


* 2017/11/30 KEGG genes have multiple synonyms in their source file. Previously, only the first one is retrieved. This update fully included all gene synonyms to the `analytesynonym` table. Also, Some IDs does not have correct ID mapping e.g. Glucose and ATP. This rebuild fix these errors from the algorithm, so the program pick up more items in mutiple tables.

Table Names | Number of Items
-----------|--------------
analyte | 91797
analytehasontology | 181689
analytehaspathway | 227940
analytesynonym | 393088
catalyzed | 861526
ontology | 160
pathway | 3559
source | 177804

* 17/12/06 We have parsed GPML file from WP to RaMP. It creates a problem that each analyte in GPML file only has one source ID instead of ID mapping. To reduces duplicated items in WP, CAS ID from HMDB is parsed to reduce 300 duplicated analytes.

Table Names | Number of Items
-----------|--------------
analyte | **91464**
analytehasontology | 181689
analytehaspathway | 238535
analytesynonym | 393583
catalyzed | 861526
ontology | 160
pathway | 3559
source | **181062**

* 17/12/16 HMDB is updated from HMDB 3.6 to HMDB 4.0. Therefore, HMDB has much much more SMPDB pathway IDs than before. Also, the number of metabolites jumped greatly.

Table Names | Number of Items
-----------|--------------
analyte | **130906**
analytehasontology | **245000**
analytehaspathway | **349916**
analytesynonym | **404312**
catalyzed | **865832**
ontology | 160
pathway | **36530**
source | **235702**

* 18/1/28 The algorithm to mapping ID across different sources has some problem due to update of their sources, which is fixed at this update. HMDB started to keep monthly update.

Table Names | Number of Items
-----------|--------------
analyte | **136802**
analytehasontology | 245000
analytehaspathway | **654109**
analytesynonym | **406512**
catalyzed | 865832
ontology | 160
pathway | **51526**
source | 217154

## Update History
* 18/3/2 Thanks to WP people's advice, we replace GPML file with their RDF files which have more comprehensive data inside. The way to update RaMP is implemented. Norm form for relational database is enforced by sqlalchemy e.g. no duplicate primary key, and well-defined relations between tables. All source ID is prepended with their source.

* 18/4/11 RaMP had all information from WP RDF file. This update reduced duplicated items, so some of table have less items. RaMP also have more IDs from more source e.g. WikiData and LIPIDMAPS.

Table Names | Number of Items
-----------|--------------
analyte | 130901
analytehasontology | 245000
analytehaspathway | 548950
analytesynonym | 406525
catalyzed | 865832
ontology | 160
pathway | 51994
source | 226897

* 18/4/30 Import lipidmaps id from HMDB to RaMP. Lipidmaps ID increases from 654 to **5218**.
