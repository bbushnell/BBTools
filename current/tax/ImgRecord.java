package tax;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import shared.KillSwitch;
import shared.Timer;
import shared.Tools;

/**
 * Represents an IMG (Integrated Microbial Genomes) database record.
 * Stores genome metadata including taxonomic ID, quality status, and file paths.
 * Used for processing IMG taxonomy dumps and managing genome collections.
 * @author Brian Bushnell
 */
public class ImgRecord implements Serializable{

	/** Serialization version for persisted IMG record maps. */
	private static final long serialVersionUID=6438551103300423985L;

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point for processing IMG dump files.
	 * Reads {@code args[0]} as an IMG taxon dump, retains all genomes unless
	 * {@link TaxTree#IMG_HQ} is set, and optionally writes map-entry text rows.
	 * @param args input path, then optional output path
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
		if(out!=null){writeAsText(map, out);}
	}

	/**
	 * Writes IMG map entries to text file format.
	 * Each line uses {@code entry.toString()} and therefore includes the IMG ID
	 * key before the record: {@code IMGID=imgID	taxID	name}.
	 * @param map Map of IMG ID to ImgRecord objects
	 * @param out Output file path
	 */
	private static void writeAsText(HashMap<Long, ImgRecord> map, String out){
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		//FIXED [tax/ImgRecord#001]: TextStreamWriter must be started and closed.
		//Without this, tsw.println either crashes under -ea at addJob's started
		//assertion after the buffer fills or silently discards tiny outputs.
		//main() reaches this for every output argument; IMG dumps take the crash
		//path. start()+poisonAndWait() restores the writer lifecycle.
		tsw.start();
		for(Entry<Long, ImgRecord> e : map.entrySet()){
			tsw.println(e.toString());
		}
		final boolean error=tsw.poisonAndWait();
		if(error){
			final String message="Error writing IMG records to "+out;
			assert(false) : KillSwitch.assertDie(message);
			throw new RuntimeException(message);
		}
	}

	/**
	 * Returns a tab-delimited record body containing IMG ID, taxID, and name.
	 * @return IMG ID, taxonomy ID and name separated by tabs
	 */
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
	 * Parses IMG dump lines into records, closing the owned reader on success
	 * or failure.
	 * Skips nondata lines whose first character is not a digit and optionally
	 * filters records lacking a yes value in column 7.
	 * If reader cleanup also fails after a primary read or parse failure, the
	 * primary failure is preserved and the cleanup failure is suppressed.
	 *
	 * @param fname Input file path containing IMG records
	 * @param highQuality If true, include only records marked as high quality
	 * @return Array of parsed ImgRecord objects
	 * @throws RuntimeException if reading, parsing, or reported close status fails
	 */
	public static ImgRecord[] toArray(String fname, boolean highQuality){
		TextFile tf=new TextFile(fname, false);
		Throwable failure=null;
		try{
			ArrayList<ImgRecord> list=new ArrayList<ImgRecord>();
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(line.length()<1 || !Tools.isDigit(line.charAt(0))){
					//do nothing
				}else{
					ImgRecord record=new ImgRecord(line);
					if(!highQuality || record.highQuality){list.add(record);}
				}
			}
			return list.toArray(new ImgRecord[0]);
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			try{
				if(tf.close()){throw new RuntimeException("Error closing "+fname);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs ImgRecord from a tab-separated IMG dump line.
	 * Consumes columns 0=IMG ID, 1=name, 2=NCBI taxID, 3=public flag,
	 * 4=obsolete flag, 5=genome type, and optional column 7=high-quality flag.
	 * Empty taxID fields become -1; names are retained only when {@link #storeName}
	 * was true before parsing. A malformed IMG ID propagates
	 * {@link NumberFormatException}; a malformed nonempty taxID is printed with
	 * the input line, then discarded before a new {@link RuntimeException}.
	 * @param line Tab-separated line from IMG dump file
	 */
	public ImgRecord(String line){
		String[] split=line.split("\t");

		imgID=Long.parseLong(split[0]);
		name=(storeName ? split[1] : null);
		try{
			taxID=(split[2]==null || split[2].length()<1 ? -1 : Integer.parseInt(split[2]));
		}catch(NumberFormatException e1){
			e1.printStackTrace();
			System.err.println(line);
			throw new RuntimeException();
		}
		isPublic=Parse.parseYesNo(split[3]);
		obsolete=Parse.parseYesNo(split[4]);
		genomeType=find(split[5], typeArray);
		boolean hq=false;
		if(split.length>7){
			try{
				hq=Parse.parseYesNo(split[7]);
			}catch(Exception e){
				System.err.println(Arrays.toString(split));
				assert(false);
			}
		}
		highQuality=hq;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** IMG database identifier for this genome. */
	public final long imgID;
	/** NCBI taxonomy ID associated with this genome, or -1 when absent. */
	public final int taxID;
	/** Genome type index: 0=isolate, 1=single_cell, 2=metagenome, -1=unmatched. */
	public final int genomeType;
	/** Whether this genome is publicly available. */
	public final boolean isPublic;
	/** Whether this genome record is marked as obsolete. */
	public final boolean obsolete;
	/** Whether optional column 7 marks this genome as high quality. */
	public final boolean highQuality;
	/** Genome name, or null when {@link #storeName} is false during construction. */
	public final String name;

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the legacy IMG FASTA path under
	 * {@code /global/dna/projectdirs/microbial/img_web_data/taxon.fna/}.
	 * @return Full path formed from this record's IMG ID
	 */
	public final String path(){return "/global/dna/projectdirs/microbial/img_web_data/taxon.fna/"+imgID+".fna";}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Genome-type indices for isolate, single-cell and metagenome records. */
	final int ISOLATE=0, SINGLE_CELL=1, METAGENOME=2;
	/** IMG genome-type labels mapped to zero-based {@link #genomeType} values. */
	final String[] typeArray={"isolate", "single_cell", "metagenome"};

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Finds the index of a string in an array using exact matching.
	 * Returns -1 if the string is not found.
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

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether to retain genome names in newly parsed records. */
	public static boolean storeName=true;
	/** Global IMG ID lookup used by sketch metadata fallback. */
	public static HashMap<Long, ImgRecord> imgMap;
	/** Default path to IMG taxonomy dump file. */
	public static final String DefaultDumpFile="/global/u1/i/img/adhocDumps/taxonDumpForBrian.txt";

}
