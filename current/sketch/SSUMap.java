package sketch;

import java.io.PrintStream;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;
import tax.GiToTaxid;
import tax.TaxTree;

/**
 * Manages Small Subunit ribosomal RNA (SSU) sequences for taxonomic identification.
 * Maintains maps of 16S and 18S rRNA sequences indexed by taxonomic ID for taxonomy assignment.
 * Loads sequences from FASTA files, keeping only the longest sequence per taxonomic ID.
 * @author Brian Bushnell
 */
public class SSUMap {
	
	/**
	 * Loads SSU sequence maps from configured files if not already loaded.
	 * Thread-safe initialization of both 16S and 18S rRNA sequence maps.
	 * @param outstream Print stream for status messages and debugging output
	 */
	public static synchronized void load(PrintStream outstream){
		if(r16SFile!=null && r16SMap==null){
			r16SMap=load(r16SFile, TaxTree.default16SFile(), outstream);
		}
		if(r18SFile!=null && r18SMap==null){
			r18SMap=load(r18SFile, TaxTree.default18SFile(), outstream);
		}
	}
	
	/**
	 * Loads SSU sequences from a specified file into a taxonomic ID-indexed map.
	 * Handles file format detection and maintains proper sequence input settings.
	 * Uses "auto" to select default file path if specified.
	 *
	 * @param ssuFile Path to SSU sequence file, or "auto" for default
	 * @param defaultFile Default file path to use when ssuFile is "auto"
	 * @param outstream Print stream for status messages
	 * @return Map of taxonomic IDs to SSU sequences, or null if no file specified
	 */
	private static synchronized HashMap<Integer, byte[]> load(String ssuFile, String defaultFile, PrintStream outstream){
		HashMap<Integer, byte[]> map=null;
		if(ssuFile!=null){
			final boolean oldAminoIn=Shared.AMINO_IN;
			final boolean oldInterleaved=FASTQ.FORCE_INTERLEAVED;
			Shared.AMINO_IN=false;//SSUs are nucleotide, which can cause a crash, esp. with IUPAC symbols
			FASTQ.FORCE_INTERLEAVED=false;

			String fname=ssuFile;
			if("auto".equalsIgnoreCase(fname)){fname=defaultFile;}
			final FileFormat ffssu=FileFormat.testInput(fname, FileFormat.FA, null, true, false);
			map=loadSSU(ffssu, outstream);

			Shared.AMINO_IN=oldAminoIn;
			FASTQ.FORCE_INTERLEAVED=oldInterleaved;
		}
		return map;
	}
	
	/**
	 * Reads SSU sequences from input stream and builds taxonomic ID-indexed map.
	 * Filters sequences by minimum length (>1000 bases) and keeps only the longest
	 * sequence per taxonomic ID to optimize memory usage and accuracy.
	 *
	 * @param ff File format specification for the input SSU file
	 * @param outstream Print stream for progress messages
	 * @return Map of taxonomic IDs to SSU sequence bytes
	 */
	private static HashMap<Integer, byte[]> loadSSU(FileFormat ff, PrintStream outstream){
		Streamer st=makeStreamer(ff, outstream);
		HashMap<Integer, byte[]> map=new HashMap<Integer, byte[]>(1000000);
		
		//As long as there is a nonempty read list...
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			if(verbose){outstream.println("Fetched "+ln.size()+" reads.");}
			for(Read r : ln){
				assert(r.mate==null);
				final int tid=GiToTaxid.getID(r.id);
				if(tid>=0 && r.length()>1000){
					byte[] old=map.get(tid);
					if(old==null || old.length<r.length()){map.put(tid, r.bases);}
				}
			}
		}
		errorState|=ReadWrite.closeStream(st);
		
		return map;
	}
	
	/**
	 * Creates and initializes a streamer for SSU file processing.
	 * Configures stream for unpaired sequence reading and starts the input pipeline.
	 *
	 * @param ff File format specification for the input file
	 * @param outstream Print stream for debug messages
	 * @return Initialized concurrent read input stream
	 */
	private static Streamer makeStreamer(FileFormat ff, PrintStream outstream){
		if(verbose){outstream.println("makeCris");}
		Streamer st=StreamerFactory.getReadInputStream(-1, false, ff, null, 1);
		st.start(); //Start the stream
		if(verbose){outstream.println("Loading "+ff.name());}
		boolean paired=st.paired();
		assert(!paired);
		return st;
	}

	/** Checks if any SSU sequence maps are loaded.
	 * @return true if either 16S or 18S maps are available, false otherwise */
	public static boolean hasMap(){return r16SMap!=null || r18SMap!=null;}
	/** Gets the number of 16S rRNA sequences currently loaded.
	 * @return Count of sequences in 16S map, or 0 if map is null */
	public static int r16SCount(){return r16SMap==null ? 0 : r16SMap.size();}
	/** Gets the number of 18S rRNA sequences currently loaded.
	 * @return Count of sequences in 18S map, or 0 if map is null */
	public static int r18SCount(){return r18SMap==null ? 0 : r18SMap.size();}
	
	/** Path to 16S rRNA sequence file for loading */
	public static String r16SFile=null;
	/** Path to 18S rRNA sequence file for loading */
	public static String r18SFile=null;
	/** Map of taxonomic IDs to 16S rRNA sequences */
	public static HashMap<Integer, byte[]> r16SMap=null;
	/** Map of taxonomic IDs to 18S rRNA sequences */
	public static HashMap<Integer, byte[]> r18SMap=null;
	/** Controls verbose output during SSU loading operations */
	static boolean verbose=false;
	/** Tracks error state during file processing operations */
	static boolean errorState=false;
	
}
