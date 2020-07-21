using System;
using System.Text;
using MySql.Data;
using MySql.Data.MySqlClient;

namespace mysqlBatch
{
	public class Program
	{
        public String sqlFilesLocation = "/Users/pati13/Downloads/RaMP-BackEnd-manju/misc/sql/";

		static void Main(string[] args)
		{
            Program prog = new Program();
            prog.createDatabaseAndAllTables();
            //source
            prog.BulkImport("hmdbsource.sql", "source", "SELECT sourceId, rampId, IDtype, geneOrCompound FROM source");
            prog.BulkImport("keggsource.sql", "source", "SELECT sourceId, rampId, IDtype, geneOrCompound FROM source");
            prog.BulkImport("reactomesource.sql", "source", "SELECT sourceId, rampId, IDtype, geneOrCompound FROM source");
            prog.BulkImport("wikisource.sql", "source", "SELECT sourceId, rampId, IDtype, geneOrCompound FROM source");
            //analyte
            prog.BulkImport("wikianalyte.sql", "analyte", "SELECT rampId, type FROM analyte");
            prog.BulkImport("reactomeanalyte.sql", "analyte", "SELECT rampId, type FROM analyte");
            prog.BulkImport("kegganalyte.sql", "analyte", "SELECT rampId, type FROM analyte");
            prog.BulkImport("hmdbanalyte.sql", "analyte", "SELECT rampId, type FROM analyte");
            //pathway
            prog.BulkImport("wikipathway.sql", "pathway", "SELECT pathwayRampId, sourceId, type FROM pathway");
            prog.BulkImport("reactomepathway.sql", "pathway", "SELECT pathwayRampId, sourceId, type FROM pathway");
            prog.BulkImport("keggpathway.sql", "pathway", "SELECT pathwayRampId, sourceId, type FROM pathway");
            prog.BulkImport("hmdbpathway.sql", "pathway", "SELECT pathwayRampId, sourceId, type FROM pathway");
            //analyteHasPathway
            prog.BulkImport("kegganalyteHasPathway.sql", "analyteHasPathway", "SELECT rampId, pathwayRampId, pathwaySource FROM analyteHasPathway");
            prog.BulkImport("hmdbanalyteHasPathway.sql", "analyteHasPathway", "SELECT rampId, pathwayRampId, pathwaySource FROM analyteHasPathway");
            prog.BulkImport("wikianalyteHasPathway.sql", "analyteHasPathway", "SELECT rampId, pathwayRampId, pathwaySource FROM analyteHasPathway");
            prog.BulkImport("reactomeanalyteHasPathway.sql", "analyteHasPathway", "SELECT rampId, pathwayRampId, pathwaySource FROM analyteHasPathway");
            
            //analyteSynonym
            prog.BulkImport("hmdbanalyteSynonym.sql", "analyteSynonym", "SELECT Synonym, rampId, geneOrCompound FROM analyteSynonym");
            prog.BulkImport("kegganalyteSynonym.sql", "analyteSynonym", "SELECT Synonym, rampId, geneOrCompound FROM analyteSynonym");
            prog.BulkImport("reactomeanalyteSynonym.sql", "analyteSynonym", "SELECT Synonym, rampId, geneOrCompound FROM analyteSynonym");
            prog.BulkImport("wikianalyteSynonym.sql", "analyteSynonym", "SELECT Synonym, rampId, geneOrCompound FROM analyteSynonym");
            //catalyzed
            prog.BulkImport("reactomecatalyzed.sql", "catalyzed", "SELECT rampCompoundId, rampGeneId FROM catalyzed");
            prog.BulkImport("hmdbcatalyzed.sql", "catalyzed", "SELECT rampCompoundId, rampGeneId FROM catalyzed");
			prog.BulkImport("keggcatalyzed.sql", "catalyzed", "SELECT rampCompoundId, rampGeneId FROM catalyzed");
			prog.BulkImport("wikicatalyzed.sql", "catalyzed", "SELECT rampCompoundId, rampGeneId FROM catalyzed");
            //analyteHasOntology
            prog.BulkImport("hmdbanalyteHasOntologyLocation.sql", "analyteHasOntology", "SELECT rampId, rampOntologyIdLocation FROM analyteHasOntology");
            //ontology
            prog.BulkImport("hmdbOntologyLocation.sql", "ontology", "SELECT rampOntologyIdLocation, commonName, biofluidORcellular FROM ontology");
            // source Id + Origin
// prog.BulkImport("hmdbEndoExo.sql", "ontologyHasOrigin", "SELECT sourceId,origin FROM ontologyHasOrigins");
        }



        public void createDatabaseAndAllTables(){

			string connStr = "server=localhost;user=root;database=ramp;port=3306;password=Autumn6265";
			MySqlConnection conn = new MySqlConnection(connStr);

            //string stringCreateDatabase = "create database mathelabramp;";
            string stringCreateTable1 = "create table source(sourceId VARCHAR(30), rampId VARCHAR(30), IDtype VARCHAR(30), geneOrCompound VARCHAR(30),commonName VARCHAR(30), PRIMARY KEY (sourceId)) engine = InnoDB;";
            string stringCreateTable2 = "create table analyte (rampId VARCHAR(30), type VARCHAR(30), PRIMARY KEY (rampId)) engine = InnoDB;";
            string stringCreateTable3 = "create table pathway (pathwayRampId VARCHAR(30), sourceId VARCHAR(30),type VARCHAR(30), pathwayCategory VARCHAR(30), pathwayName VARCHAR(250), PRIMARY KEY (pathwayRampId)) engine = InnoDB;";
            string stringCreateTable4 = "create table analyteHasPathway (rampId VARCHAR(30), pathwayRampId VARCHAR(30), pathwaySource VARCHAR(30)) engine = InnoDB;";
            string stringCreateTable5 = "create table analyteSynonym (Synonym VARCHAR(500), rampId VARCHAR(30), geneOrCompound VARCHAR(30), source VARCHAR(30), sourceId VARCHAR(30)) engine = InnoDB;";
            string stringCreateTable6 = "create table catalyzed (rampId VARCHAR(30), rampGeneId VARCHAR(30), sourceId VARCHAR(30)) engine = InnoDB;";
            string stringCreateTable7 = "create table analyteHasOntology (rampId VARCHAR(30), OntologyRampId VARCHAR(30), sourceId VARCHAR(30)) engine = InnoDB;";
            string stringCreateTable8 = "create table ontology (rampOntologyIdLocation VARCHAR(30), commonName VARCHAR(30), biofluidORcellular VARCHAR(30)) engine = InnoDB;";
// string stringCreateTable9 = "create table ontologyHasOrigin (sourceId VARCHAR(30), origin VARCHAR(30)) engine = InnoDB;";

            //MySqlCommand CreateDatabase = new MySqlCommand(stringCreateDatabase, conn);
            MySqlCommand CreateTable1 = new MySqlCommand(stringCreateTable1, conn);
            MySqlCommand CreateTable2 = new MySqlCommand(stringCreateTable2, conn);
            MySqlCommand CreateTable3 = new MySqlCommand(stringCreateTable3, conn);
            MySqlCommand CreateTable4 = new MySqlCommand(stringCreateTable4, conn);
            MySqlCommand CreateTable5 = new MySqlCommand(stringCreateTable5, conn);
            MySqlCommand CreateTable6 = new MySqlCommand(stringCreateTable6, conn);
            MySqlCommand CreateTable7 = new MySqlCommand(stringCreateTable7, conn);
            MySqlCommand CreateTable8 = new MySqlCommand(stringCreateTable8, conn);
// MySqlCommand CreateTable9 = new MySqlCommand(stringCreateTable9, conn);
            try
			{
				Console.WriteLine("Connecting to MySQL...");
                conn.Open();
               // CreateDatabase.ExecuteNonQuery();
				CreateTable1.ExecuteNonQuery();
                CreateTable2.ExecuteNonQuery();
                CreateTable3.ExecuteNonQuery();
                CreateTable4.ExecuteNonQuery();
                CreateTable5.ExecuteNonQuery();
                CreateTable6.ExecuteNonQuery();
                CreateTable7.ExecuteNonQuery();
				CreateTable8.ExecuteNonQuery();
// CreateTable9.ExecuteNonQuery();
                Console.WriteLine("Tables Created Successfully...");
				conn.Close();
			}
			catch (Exception ex)
			{
				Console.WriteLine(ex.ToString());
			}


		}

        public void BulkImport(string fileName, string destinationTableName, string command)
		{

			string connStr = "server=localhost;user=root;database=ramp;port=3306;password=Autumn6265";
			MySqlConnection conn = new MySqlConnection(connStr);



			MySqlBulkLoader bulk1 = new MySqlBulkLoader(conn);
			bulk1.TableName = destinationTableName;
			bulk1.FieldTerminator = "\t";
			bulk1.LineTerminator = "\n";
            bulk1.FileName = sqlFilesLocation + fileName;
			bulk1.NumberOfLinesToSkip = 0;

			try
			{
                Console.WriteLine("Connecting to MySQL for " + destinationTableName + "...");
				conn.Open();

				//Upload data from file
				int count = bulk1.Load();
				

				string sql = command;
				MySqlCommand cmd = new MySqlCommand(sql, conn);
				MySqlDataReader rdr = cmd.ExecuteReader();
				rdr.Close();
				conn.Close();
                Console.WriteLine(count + " lines uploaded.");

			}
			catch (Exception ex)
			{
				Console.WriteLine(ex.ToString());
			}

        }


	}





}
