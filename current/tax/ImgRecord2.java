package tax;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Timer;
import shared.Tools;

/**
 * Represents and parses IMG (Integrated Microbial Genomes) records.
 * Converts IMG record text files into HashMap or array containers for lookup.
 * Handles tab-delimited IMG record files with genome, taxonomy, and name information.
 *
 * @author Brian Bushnell
 */
public class ImgRecord2 implements Serializable{

	/** Serialization version identifier. */
	private static final long serialVersionUID=8596055205924485293L;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point for IMG record processing.
	 * Reads IMG records from input file and optionally writes serialized HashMap to output file.
	 * @param args Command-line arguments: input_file [output_file]
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
		HashMap<Long, ImgRecord2> map=toMap(in);
		t.stop();
		System.err.println(map.size()+"; "+t);
		if(out!=null){ReadWrite.writeObjectInThread(map, out, false);}
	}

	/**
	 * Converts an IMG record file to a HashMap for efficient lookup by IMG ID.
	 * @param fname Input file path containing IMG records
	 * @return HashMap mapping IMG IDs to ImgRecord2 objects
	 */
	public static HashMap<Long, ImgRecord2> toMap(String fname){
		ImgRecord2[] array=toArray(fname);
		HashMap<Long, ImgRecord2> map=new HashMap<Long, ImgRecord2>((3+array.length*4)/3);
		for(ImgRecord2 ir : array){
			map.put(ir.imgID, ir);
		}
		return map;
	}

	/**
	 * Parses IMG record file and returns array of ImgRecord2 objects.
	 * Skips lines that don't start with digits and processes valid IMG records.
	 * Closes its owned reader after successful and failed parses; parse errors remain
	 * primary if a distinct close error is reported during cleanup.
	 * @param fname Input file path containing tab-delimited IMG records
	 * @return Array of parsed ImgRecord2 objects
	 * @throws RuntimeException If reading fails, parsing fails, or closure reports an error
	 */
	public static ImgRecord2[] toArray(String fname){
		TextFile tf=new TextFile(fname, false);
		Throwable failure=null;
		try{
			ArrayList<ImgRecord2> list=new ArrayList<ImgRecord2>();
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(line.length()<1 || !Tools.isDigit(line.charAt(0))){
					//do nothing
				}else{
					ImgRecord2 record=new ImgRecord2(line);
					list.add(record);
				}
			}
			return list.toArray(new ImgRecord2[0]);
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

	/**
	 * Constructs an ImgRecord2 from a tab-delimited line.
	 * Parses genome ID, taxonomy ID, and, when {@link #storeName} is true, a third name field.
	 * Java's default split drops trailing empty fields, so a trailing empty name is missing.
	 * @param line Tab-delimited line beginning with imgID and taxID fields; a name field is
	 * required when {@link #storeName} is true
	 */
	public ImgRecord2(String line){
		String[] split=line.split("\t");

		imgID=Long.parseLong(split[0]);
		taxID=(split[1]==null || split[1].length()<1 ? 0 : Integer.parseInt(split[1]));
		name=(storeName ? split[2] : null);
	}

	/**
	 * Extracts IMG identifier from sequence header using delimiter detection.
	 * Searches for "img" followed by delimiter and parses the numeric ID.
	 *
	 * @param header Sequence header string containing an IMG identifier
	 * @param doAssertions Whether to assert missing delimiter, missing ID, and zero-ID cases
	 * @return Positive parsed IMG identifier, or -1 for missing delimiter, missing token or zero ID
	 * @throws AssertionError If guarded missing-token checks or the digit assertion fail
	 * @throws NumberFormatException If the numeric identifier exceeds Long.MAX_VALUE
	 */
	public static final long parseImgId(String header, boolean doAssertions){
		final char delimiter=TaxTree.ncbiHeaderDelimiter(header);
		assert(!doAssertions || delimiter!=' ') : header;
		if(delimiter==' '){return -1;}

		final int idx=header.indexOf(delimiter=='|' ? "img|" : "img~");
		assert(!doAssertions || idx>=0) : "Could not img id number from "+header;
		if(idx<0){return -1;}

		long img=0;
		for(int i=idx+4; i<header.length(); i++){
			final char c=header.charAt(i);
			if(c==delimiter || c==' '){break;}
			assert(Tools.isDigit(c)) : c+", '"+header+"'";
			final int digit=c-'0';
			if(img>Long.MAX_VALUE/10 || (img==Long.MAX_VALUE/10 && digit>Long.MAX_VALUE%10)){
				throw new NumberFormatException("IMG identifier exceeds Long.MAX_VALUE: "+header);
			}
			img=img*10+digit;
		}

		assert(!doAssertions || img>0) : "Could not img id number from "+header;
		return img>0 ? img : -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Controls whether to store name fields in ImgRecord2 objects. */
	public static boolean storeName=true;

	/** IMG genome identifier for the microbial genome. */
	public final long imgID;

	/** NCBI taxonomy identifier for the organism. */
	public final int taxID;

	/** Organism or genome name, or null if storeName is false. */
	public final String name;

}
