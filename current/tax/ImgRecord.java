package tax;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import shared.Timer;
import shared.Tools;

/**
 * Represents an IMG (Integrated Microbial Genomes) database record.
 * Stores genome metadata including taxonomic ID, quality status, and file paths.
 * Used for processing IMG taxonomy dumps and managing genome collections.
 * @author Brian Bushnell
 */
public class ImgRecord implements Serializable {
	
	private static final long serialVersionUID = 6438551103300423985L;
	
	/**
	 * Program entry point for processing IMG dump files.
	 * Reads an IMG taxonomy dump and optionally writes output as text.
	 * @param args Command-line arguments: [input_file] [output_file]
	 */
	public static void main(String[] args){
		String in=args[0];
		String out=args.length>1 ? args[1] : null;

		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testOutputFiles(true, false, false, out)){
			throw new RuntimeException("\nCan't write to some output files.\n");
		}
		Timer t=new Timer();
		HashMap<Long, ImgRecord> map=toMap(in, TaxTree.IMG_HQ);
		t.stop();
		System.err.println(map.size()+"; "+t);
//		if(out!=null){ReadWrite.writeObjectInThread(map, out, false);}
		if(out!=null){writeAsText(map, out);}
	}
	
	/**
	 * Writes IMG records to text file format.
	 * Each record is written as a tab-separated line using toString().
	 * @param map Map of IMG ID to ImgRecord objects
	 * @param out Output file path
	 */
	private static void writeAsText(HashMap<Long, ImgRecord> map, String out){
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		for(Entry<Long, ImgRecord> e : map.entrySet()){
			tsw.println(e.toString());
		}
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(imgID);
		sb.append('\t').append(taxID);
		sb.append('\t').append(name);
		return sb.toString();
	}
	
	/**
	 * Creates a HashMap from IMG dump file with IMG ID as key.
	 * Optionally filters to include only high-quality genomes.
	 *
	 * @param fname Input file path containing IMG records
	 * @param highQuality If true, include only records marked as high quality
	 * @return HashMap mapping IMG IDs to ImgRecord objects
	 */
	public static HashMap<Long, ImgRecord> toMap(String fname, boolean highQuality){
		ImgRecord[] array=toArray(fname, highQuality);
		HashMap<Long, ImgRecord> map=new HashMap<Long, ImgRecord>((3+array.length*4)/3);
		for(ImgRecord ir : array){
			map.put(ir.imgID, ir);
		}
		return map;
	}
	
	/**
	 * Parses IMG dump file into an array of ImgRecord objects.
	 * Skips lines that don't start with digits and optionally filters by quality.
	 *
	 * @param fname Input file path containing IMG records
	 * @param highQuality If true, include only records marked as high quality
	 * @return Array of parsed ImgRecord objects
	 */
	public static ImgRecord[] toArray(String fname, boolean highQuality){
		TextFile tf=new TextFile(fname, false);
		ArrayList<ImgRecord> list=new ArrayList<ImgRecord>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()<1 || !Tools.isDigit(line.charAt(0))){
				//do nothing
			}else{
				ImgRecord record=new ImgRecord(line);
				if(!highQuality || record.highQuality){list.add(record);}
			}
		}
		tf.close();
		return list.toArray(new ImgRecord[0]);
	}
	
	/**
	 * Constructs ImgRecord from tab-separated line of IMG dump file.
	 * Parses IMG ID, name, tax ID, public status, obsolete status, genome type,
	 * and high quality flag from the input line.
	 *
	 * @param line Tab-separated line from IMG dump file
	 * @throws NumberFormatException If IMG ID or tax ID cannot be parsed
	 */
	public ImgRecord(String line){
		String[] split=line.split("\t");
		
		imgID=Long.parseLong(split[0]);
		name=(storeName ? split[1] : null);
		try {
			taxID=(split[2]==null || split[2].length()<1 ? -1 : Integer.parseInt(split[2]));
		} catch (NumberFormatException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			System.err.println(line);
			throw new RuntimeException();
		}
		isPublic=Parse.parseYesNo(split[3]);
		obsolete=Parse.parseYesNo(split[4]);
		genomeType=find(split[5], typeArray);
		boolean hq=false;
		if(split.length>7){
			try {
				hq=Parse.parseYesNo(split[7]);
			} catch (Exception e) {
				System.err.println(Arrays.toString(split));
				assert(false);
			}
		}
		highQuality=hq;
	}
	
	/** IMG database identifier for this genome */
	public final long imgID;
	/** NCBI taxonomy ID associated with this genome */
	public final int taxID;
	/** Genome type index: 0=isolate, 1=single_cell, 2=metagenome */
	public final int genomeType;
	/** Whether this genome is publicly available */
	public final boolean isPublic;
	/** Whether this genome record is marked as obsolete */
	public final boolean obsolete;
	/** Whether this genome meets high quality standards */
	public final boolean highQuality;
	/** Genome name/description (null if storeName is false) */
	public final String name;
	/**
	 * Returns the file system path to the genome FASTA file.
	 * Path format: /global/dna/projectdirs/microbial/img_web_data/taxon.fna/[imgID].fna
	 * @return Full path to the genome file
	 */
	public final String path(){return "/global/dna/projectdirs/microbial/img_web_data/taxon.fna/"+imgID+".fna";}
	
	final int ISOLATE=0, SINGLE_CELL=1, METAGENOME=2;
	final String[] typeArray={"isolate", "single_cell", "metagenome"};
	/**
	 * Finds the index of a string in an array using exact matching.
	 * Returns -1 if the string is not found.
	 *
	 * @param s String to search for
	 * @param array Array to search in
	 * @return Index of the string, or -1 if not found
	 */
	private static int find(String s, String[] array){
		for(int i=0; i<array.length; i++){
			if(array[i].equals(s)){return i;}
		}
		return -1;
	}
	
	/** Whether to store genome names in memory to save space */
	public static boolean storeName=true;
	/** Global map of IMG IDs to ImgRecord objects for lookup */
	public static HashMap<Long, ImgRecord> imgMap;
//	public static final String DefaultDumpFile="/global/cfs/cdirs/bbtools/tax/imgTaxDump.txt.gz";
	/** Default path to IMG taxonomy dump file */
	public static final String DefaultDumpFile="/global/u1/i/img/adhocDumps/taxonDumpForBrian.txt";
	
}
