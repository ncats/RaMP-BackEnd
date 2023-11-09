-- MySQL dump 10.13  Distrib 8.0.21, for Win64 (x86_64)
--
-- Host: localhost    Database: ramp
-- ------------------------------------------------------
-- Server version	8.0.28

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `analyte`
--

DROP TABLE IF EXISTS `analyte`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `analyte` (
  `rampId` varchar(30) NOT NULL,
  `type` varchar(30) DEFAULT NULL,
  PRIMARY KEY (`rampId`),
  KEY `analyte_rampId_RampID_IDX` (`rampId`) USING BTREE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `analytehasontology`
--

DROP TABLE IF EXISTS `analytehasontology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `analytehasontology` (
  `rampCompoundId` varchar(30) DEFAULT NULL,
  `rampOntologyId` varchar(30) DEFAULT NULL,
  KEY `analyte_ont_ramp_id_idx` (`rampCompoundId`),
  KEY `analyte_ont_id_idx` (`rampOntologyId`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `analytehaspathway`
--

DROP TABLE IF EXISTS `analytehaspathway`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `analytehaspathway` (
  `rampId` varchar(30) DEFAULT NULL,
  `pathwayRampId` varchar(30) DEFAULT NULL,
  `pathwaySource` varchar(30) DEFAULT NULL,
  KEY `pathwayRampID_IDX` (`pathwayRampId`) USING BTREE,
  KEY `ahp_RampID_IDX` (`rampId`) USING BTREE,
  KEY `ahp_path_source_IDX` (`pathwaySource`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `analytesynonym`
--

DROP TABLE IF EXISTS `analytesynonym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `analytesynonym` (
  `Synonym` varchar(500) DEFAULT NULL,
  `rampId` varchar(30) DEFAULT NULL,
  `geneOrCompound` varchar(30) DEFAULT NULL,
  `source` varchar(30) DEFAULT NULL,
  KEY `analSyn_RampID_IDX` (`rampId`) USING BTREE,
  KEY `analSyn_syn_IDX` (`Synonym`) USING BTREE,
  KEY `idx_analytesynonym_geneOrCompound` (`geneOrCompound`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `catalyzed`
--

DROP TABLE IF EXISTS `catalyzed`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `catalyzed` (
  `rampCompoundId` varchar(30) DEFAULT NULL,
  `rampGeneId` varchar(30) DEFAULT NULL,
  `proteinType` varchar(32) DEFAULT NULL,
  KEY `catal_gene_idx` (`rampGeneId`),
  KEY `catal_comp_idx` (`rampCompoundId`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chem_props`
--

DROP TABLE IF EXISTS `chem_props`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `chem_props` (
  `ramp_id` varchar(30) NOT NULL,
  `chem_data_source` varchar(32) DEFAULT NULL,
  `chem_source_id` varchar(45) DEFAULT NULL,
  `iso_smiles` varchar(256) DEFAULT NULL,
  `inchi_key_prefix` varchar(32) DEFAULT NULL,
  `inchi_key` varchar(32) DEFAULT NULL,
  `inchi` varchar(4096) DEFAULT NULL,
  `mw` float DEFAULT NULL,
  `monoisotop_mass` float DEFAULT NULL,
  `common_name` varchar(1024) DEFAULT NULL,
  `mol_formula` varchar(64) DEFAULT NULL,
  KEY `prop_source_idx` (`chem_data_source`) USING BTREE,
  KEY `inchi_key_idx` (`inchi_key`) USING BTREE,
  KEY `inchi_key_prefix_idx` (`inchi_key_prefix`) USING BTREE,
  KEY `ramp_id_idx` (`ramp_id`) USING BTREE
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='Holds metabolite properties for all ramp metabolites';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `db_version`
--

DROP TABLE IF EXISTS `db_version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `db_version` (
  `ramp_version` varchar(20) NOT NULL,
  `load_timestamp` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `version_notes` varchar(256) DEFAULT NULL,
  `met_intersects_json` varchar(10000) DEFAULT NULL,
  `gene_intersects_json` varchar(10000) DEFAULT NULL,
  `met_intersects_json_pw_mapped` varchar(10000) DEFAULT NULL,
  `gene_intersects_json_pw_mapped` varchar(10000) DEFAULT NULL,
  `db_sql_url` varchar(256) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `entity_status_info`
--

DROP TABLE IF EXISTS `entity_status_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `entity_status_info` (
  `status_category` varchar(64) NOT NULL,
  `entity_source_id` varchar(32) NOT NULL,
  `entity_source_name` varchar(45) NOT NULL,
  `entity_count` int NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='holds entity counts';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `metabolite_class`
--

DROP TABLE IF EXISTS `metabolite_class`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `metabolite_class` (
  `ramp_id` varchar(32) NOT NULL,
  `class_source_id` varchar(32) NOT NULL,
  `class_level_name` varchar(128) NOT NULL,
  `class_name` varchar(128) NOT NULL,
  `source` varchar(32) NOT NULL,
  KEY `ramp_id_metclass_idx` (`ramp_id`),
  KEY `class_source_id_metclass_idx` (`class_source_id`),
  KEY `class_name_metclass_idx` (`class_name`),
  KEY `class_source_metclass_idx` (`source`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='holds rampid and class source id to metabolic class levels and names';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ontology`
--

DROP TABLE IF EXISTS `ontology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ontology` (
  `rampOntologyId` varchar(30) DEFAULT NULL,
  `commonName` varchar(64) DEFAULT NULL,
  `HMDBOntologyType` varchar(30) DEFAULT NULL,
  `metCount` int DEFAULT '0',
  KEY `ontol_parent_idx` (`commonName`),
  KEY `ontol_ramp_id_idx` (`rampOntologyId`),
  KEY `ontol_term_idx` (`HMDBOntologyType`),
  KEY `ontol_metCount_idx` (`metCount`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pathway`
--

DROP TABLE IF EXISTS `pathway`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `pathway` (
  `pathwayRampId` varchar(30) NOT NULL,
  `sourceId` varchar(30) DEFAULT NULL,
  `type` varchar(30) DEFAULT NULL,
  `pathwayCategory` varchar(30) DEFAULT NULL,
  `pathwayName` varchar(250) DEFAULT NULL,
  PRIMARY KEY (`pathwayRampId`),
  KEY `pathway_RampID_IDX` (`pathwayRampId`) USING BTREE,
  KEY `idx_pathway_sourceId` (`sourceId`),
  KEY `idx_pathway_type` (`type`),
  KEY `idx_pathway_pathwayCategory` (`pathwayCategory`),
  FULLTEXT KEY `idx_pathway_pathwayName` (`pathwayName`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ramp_data_object`
--

DROP TABLE IF EXISTS `ramp_data_object`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ramp_data_object` (
  `data_key` varchar(32) NOT NULL,
  `data_blob` longblob NOT NULL,
  PRIMARY KEY (`data_key`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `reaction`
--

DROP TABLE IF EXISTS `reaction`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `reaction` (
  `ramp_rxn_id` varchar(16) NOT NULL,
  `rxn_source_id` varchar(16) NOT NULL,
  `status` int NOT NULL,
  `is_transport` int NOT NULL,
  `direction` varchar(8) NOT NULL,
  `label` varchar(256) NOT NULL,
  `equation` varchar(256) NOT NULL,
  `html_equation` varchar(256) NOT NULL,
  `ec_num` varchar(256) DEFAULT NULL,
  `has_human_prot` int NOT NULL,
  `only_human_mets` int NOT NULL,
  PRIMARY KEY (`ramp_rxn_id`),
  KEY `reaction_src_id_idx` (`rxn_source_id`),
  KEY `reaction_ec_num_idx` (`ec_num`),
  KEY `reaction_has_human_prot_idx` (`has_human_prot`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='ramp primary reaction annotation table.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `reaction2met`
--

DROP TABLE IF EXISTS `reaction2met`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `reaction2met` (
  `ramp_rxn_id` varchar(16) NOT NULL,
  `rxn_source_id` varchar(16) NOT NULL,
  `ramp_cmpd_id` varchar(16) NOT NULL,
  `substrate_product` int NOT NULL,
  `met_source_id` varchar(32) NOT NULL,
  `met_name` varchar(256) DEFAULT NULL,
  `is_cofactor` int NOT NULL DEFAULT '0',
  KEY `rxn2met_rxn_ramp_id_idx` (`ramp_rxn_id`),
  KEY `rxn2met_rxn_source_id_idx` (`rxn_source_id`),
  KEY `rxn2met_met_ramp_id_idx` (`ramp_cmpd_id`),
  KEY `rxn2met_subs_prod_idx` (`substrate_product`),
  KEY `rxn2met_met_source_id_idx` (`met_source_id`),
  KEY `rxn2met_iscofactor_idx` (`is_cofactor`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='holds reaction to metabolite mapping.';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `reaction2protein`
--

DROP TABLE IF EXISTS `reaction2protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `reaction2protein` (
  `ramp_rxn_id` varchar(16) NOT NULL,
  `rxn_source_id` varchar(16) NOT NULL,
  `ramp_gene_id` varchar(16) NOT NULL,
  `uniprot` varchar(16) NOT NULL,
  `protein_name` varchar(16) NOT NULL,
  KEY `rxn_prot_ramp_id_idx` (`ramp_rxn_id`),
  KEY `rxn2prot_source_id_idx` (`rxn_source_id`),
  KEY `rxn2prot_ramp_gene_id_idx` (`ramp_gene_id`),
  KEY `rxn2prot_uniprot_idx` (`uniprot`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 COMMENT='maps reaction ids to associated uniprot ids and their ramp ids';
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `reaction_protein2met`
--

DROP TABLE IF EXISTS `reaction_protein2met`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `reaction_protein2met` (
  `ramp_rxn_id` varchar(16) NOT NULL,
  `rxn_source_id` varchar(16) NOT NULL,
  `ramp_gene_id` varchar(16) NOT NULL,
  `gene_source_id` varchar(16) NOT NULL,
  `substrate_product` int NOT NULL,
  `ramp_cmpd_id` varchar(16) NOT NULL,
  `cmpd_source_id` varchar(45) NOT NULL,
  `cmpd_name` varchar(256) DEFAULT NULL,
  `is_cofactor` int NOT NULL DEFAULT '0',
  KEY `rxn_p2m_rxn_ramp_id_idx` (`ramp_rxn_id`),
  KEY `rxn_p2m_rxn_source_id_idx` (`rxn_source_id`),
  KEY `rxn_p2m_rxn_gene_ramp_id_idx` (`ramp_gene_id`),
  KEY `rxn_p2m_gene_source_id_idx` (`gene_source_id`),
  KEY `rxn_subprod_idx` (`substrate_product`),
  KEY `rxn_p2m_ramp_cmpd_id_idx` (`ramp_cmpd_id`),
  KEY `rxn_p2m_cmpd_source_id_idx` (`cmpd_source_id`),
  KEY `rxn_p2m_iscofactor_idx` (`is_cofactor`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `source`
--

DROP TABLE IF EXISTS `source`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `source` (
  `sourceId` varchar(30) CHARACTER SET utf8 COLLATE utf8_general_ci NOT NULL,
  `rampId` varchar(30) CHARACTER SET utf8 COLLATE utf8_general_ci DEFAULT NULL,
  `IDtype` varchar(30) CHARACTER SET utf8 COLLATE utf8_general_ci DEFAULT NULL,
  `geneOrCompound` varchar(30) CHARACTER SET utf8 COLLATE utf8_general_ci DEFAULT NULL,
  `commonName` varchar(256) CHARACTER SET utf8 COLLATE utf8_general_ci DEFAULT NULL,
  `priorityHMDBStatus` varchar(32) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `dataSource` varchar(32) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL,
  `pathwayCount` int NOT NULL DEFAULT '0',
  KEY `source_RampID_IDX` (`rampId`) USING BTREE,
  KEY `source_sid_RampID_IDX` (`sourceId`) USING BTREE,
  KEY `source_comName_RampID_IDX` (`commonName`) USING BTREE,
  KEY `source_datasrc_IDX` (`dataSource`),
  KEY `source_pathCount_IDX` (`pathwayCount`),
  KEY `idx_source_geneOrCompound` (`geneOrCompound`),
  KEY `idx_source_IDtype` (`IDtype`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3 COLLATE=utf8_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `version_info`
--

DROP TABLE IF EXISTS `version_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `version_info` (
  `ramp_db_version` varchar(16) NOT NULL,
  `db_mod_date` date NOT NULL,
  `status` varchar(16) NOT NULL,
  `data_source_id` varchar(32) NOT NULL,
  `data_source_name` varchar(128) NOT NULL,
  `data_source_url` varchar(128) NOT NULL,
  `data_source_version` varchar(128) NOT NULL,
  KEY `status_index` (`status`),
  KEY `data_source_index` (`data_source_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2023-10-30 14:02:17
