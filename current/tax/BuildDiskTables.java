package tax;

import java.io.IOException;
import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Builds disk-backed binary tables for the taxonomy server.
 * Reads NCBI accession-to-taxid files and GI tables, writes
 * memory-mappable binary files.
 *
 * Usage:
 *   BuildDiskTables accession=auto pattern=auto table=auto tree=auto \
 *     accout=accession_disk.bin giout=gi_disk.bin
 *
 * @author Brian Bushnell, Chloe
 */
public class BuildDiskTables {

	public static void main(String[] args){
		Timer t=new Timer();
		BuildDiskTables x=new BuildDiskTables(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public BuildDiskTables(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_UNPIGZ=true;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("accout") || a.equals("accessionout")){
				accOut=b;
			}else if(a.equals("giout")){
				giOut=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		if("auto".equalsIgnoreCase(patternFile)){patternFile=TaxTree.defaultPatternFile();}
	}

	void process(Timer t){
		if(taxTreeFile!=null){
			outstream.println("Loading tax tree.");
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, false, false);
		}

		if(patternFile!=null){
			outstream.println("Loading pattern table.");
			AnalyzeAccession.loadCodeMap(patternFile);
		}

		if(giTableFile!=null && giOut!=null){
			buildGiTable();
		}

		if(accessionFile!=null && accOut!=null){
			if(AnalyzeAccession.codeMap==null){
				throw new RuntimeException("Building a disk accession table requires pattern=...");
			}
			buildAccessionTable();
		}

		t.stop();
		outstream.println("Total time: "+t);
	}

	/*--------------------------------------------------------------*/

	private void buildGiTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building GI disk table ---");

		GiToTaxid.initialize(giTableFile);

		try{
			DiskGiTable.build(giOut, GiToTaxid.array, GiToTaxid.maxGiLoaded);
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		t.stop();
		outstream.println("GI table built in "+t);
		GiToTaxid.unload();
	}

	/*--------------------------------------------------------------*/

	private void buildAccessionTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building accession disk table ---");

		outstream.println("Pass 1: counting digitizable entries.");
		long digitizable=0, total=0, nonDigitizable=0;
		for(String fname : accessionFile.split(",")){
			ByteFile bf=ByteFile.makeByteFile(fname, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(skipAccessionLine(line)){continue;}
				total++;
				if(AnalyzeAccession.digitize(line)>0){digitizable++;}
				else{nonDigitizable++;}
			}
			bf.close();
		}
		outstream.println("Total lines: "+total);
		outstream.println("Digitizable: "+digitizable);
		outstream.println("Non-digitizable: "+nonDigitizable);

		outstream.println("Pass 2: building disk table.");
		long collected=0, added=0;
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(accOut, digitizable)){
			for(String fname : accessionFile.split(",")){
				ByteFile bf=ByteFile.makeByteFile(fname, true);
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(skipAccessionLine(line)){continue;}
					final long number=AnalyzeAccession.digitize(line);
					if(number>0){
						collected++;
						final int taxid=AccessionToTaxid.parseLineToTaxid(line, (byte)'\t');
						if(taxid>0 && builder.add(number, taxid)){added++;}
					}
				}
				bf.close();
			}
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		outstream.println("Collected "+collected+" entries.");
		outstream.println("Added "+added+" entries.");
		t.stop();
		outstream.println("Accession table built in "+t);
	}

	private static boolean skipAccessionLine(byte[] line){
		return line.length<1 || Tools.startsWith(line, "accession") || Tools.startsWith(line, "#");
	}

	/*--------------------------------------------------------------*/

	private String accessionFile=null;
	private String giTableFile=null;
	private String taxTreeFile="auto";
	private String patternFile=null;
	private String accOut=null;
	private String giOut=null;
	private TaxTree tree=null;
	private PrintStream outstream=System.err;

}
