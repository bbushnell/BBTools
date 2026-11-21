package jgi;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;
import structures.IntList;

/**
 * Eagerly loads byte-oriented FASTA assembly statistics into mutable, unsynchronized state.
 * Empty reader lines are skipped; a first-byte {@code >} starts a header. Other
 * nonempty lines contribute every byte, including whitespace, to sequence lengths
 * and composition. Sequence before the first header is accepted; empty records
 * are omitted. This is not a strict FASTA validator.
 * An optional limit retains the first nonempty records in input order; the default
 * is unlimited. Successful loads retain positive contig lengths in descending order. Each fits
 * in an int; their total is a long. Sequence lookup requires ASCII bytes, while
 * header payloads are decoded with the platform default charset.
 * Reloads replace state without rollback after parsing begins. Queries require a
 * successful load; concurrent mutation, querying or lookup-table edits need coordination.
 * @author Brian Bushnell
 */
public class Assembly{
	
	/** Loads the named input before returning normally.
	 * @param fname_ Nonnull input name interpreted by FileFormat and ByteFile; retained for reloads */
	public Assembly(String fname_){
		this(fname_, -1);
	}

	/** Loads at most the requested number of nonempty sequence records in input order.
	 * Empty header records do not count; preheader sequence counts as a record.
	 * @param fname_ Nonnull input name retained for reloads
	 * @param maxRecords_ Per-load record limit; negative is unlimited, zero retains no records */
	public Assembly(String fname_, long maxRecords_){
		fname=fname_;
		maxRecords=maxRecords_;
		load();
	}
	
	/**
	 * Reloads the retained input name through FileFormat and the selected ByteFile backend.
	 * Acquisition precedes clearing this object's statistics, so acquisition failure
	 * can retain the previous state. Once acquired-reader processing starts, failure
	 * can leave partial state. Factory-internal partial acquisition is outside the
	 * acquired-reader cleanup contract.
	 */
	void load(){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		assert(ff.fasta());
		ByteFile bf=ByteFile.makeByteFile(ff);
		load(bf);
	}

	/** Replaces statistics from an acquired reader, then sorts and totals after checked closure.
	 * Explicitly attempts close in finally. A close error rejects normal completion;
	 * a distinct close error is suppressed on the original RuntimeException/Error.
	 * Borrowed underlying streams follow the backend's ownership rules. Replacement is not transactional.
	 * @param bf Nonnull reader whose closure is this helper's responsibility */
	private void load(final ByteFile bf){
		assert(bf!=null) : "Assembly parsing requires the reader acquired by ByteFile.makeByteFile";
		Throwable failure=null;
		try{
			clear();
			acgtnio=new long[8];
			int contigLen=0;
			for(byte[] line=(maxRecords==0 ? null : bf.nextLine()); line!=null; line=bf.nextLine()){
				//JGI-056: empty physical lines do not contribute sequence data.
				if(line.length==0){continue;}
				if(line[0]=='>'){
					if(contigLen>0){
						contigs.add(contigLen);
						contigLen=0;
						if(maxRecords>=0 && contigs.size>=maxRecords){break;}
					}
					headerLength+=(line.length-1);
					if(firstHeader==null && line.length>1){
						firstHeader=new String(line, 1, line.length-1);
					}
				}else{
					addToACGTNIO(line);
					contigLen=addContigLength(contigLen, line.length);
				}
			}
			//The final nonempty record has no following header to flush it (JGI-055).
			if(contigLen>0){contigs.add(contigLen);}
		}catch(RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{
			//JGI-057: close after rejection and preserve a primary parsing failure.
			try{
				if(bf.close()){throw new RuntimeException("Error reading assembly input "+fname);}
			}catch(RuntimeException|Error e){
				if(failure==null){throw e;}
				if(e!=failure){failure.addSuppressed(e);}
			}
		}
		contigs.sort();
		contigs.reverse();
		length=contigs.sumLong();
	}
	
	/** Adds nonnegative sequence-byte counts within IntList's per-contig int representation (JGI-061).
	 * @param current Nonnegative accumulated length from load's zero/reset/checked-return path
	 * @param added Nonnegative length of the current reader line
	 * @return Sum, at most Integer.MAX_VALUE
	 * @throws RuntimeException if the sum exceeds Integer.MAX_VALUE */
	private int addContigLength(final int current, final int added){
		assert(current>=0 && added>=0) : "Assembly.load starts at zero and adds nonnegative line lengths through this checked helper: "+current+" + "+added;
		if(added>Integer.MAX_VALUE-current){
			throw new RuntimeException("Contig exceeds 2147483647 bases in assembly input "+fname);
		}
		return current+added;
	}

	/**
	 * Counts supplied bytes through the shared lookup; the caller selects sequence lines.
	 * Requires the eight-slot composition array initialized by load. Bytes index the
	 * lookup directly, so non-ASCII sequence bytes are not supported by this path.
	 * @param line Nonnull ASCII sequence bytes; whitespace and other ASCII bytes are not removed
	 */
	void addToACGTNIO(byte[] line){
		for(byte b : line){
			byte x=baseToACGTNIO[b];
			acgtnio[x]++;
		}
	}
	
	/** Invalidates loaded statistics, retaining the contig list's backing storage for reuse.
	 * Resets lengths/header state and nulls composition; gc requires a subsequent load. */
	void clear(){
		contigs.clear();
		length=0;
		headerLength=0;
		firstHeader=null;
		acgtnio=null;
	}
	
	/** Computes an approximate float GC fraction from initialized composition counts.
	 * N, IUPAC and OTHER counts are excluded; U contributes to the denominator.
	 * @return (C+G)/max(1,A+T+U+C+G), including zero when no denominator bases exist */
	float gc(){
		float AT=acgtnio[A]+acgtnio[T]+acgtnio[U];
		float GC=acgtnio[G]+acgtnio[C];
		return GC/Tools.max(1, AT+GC);
	}
	
	/**
	 * Sums contigs meeting an inclusive length cutoff, using successful-load descending order.
	 * @param minimum Lower bound in counted sequence bytes; nonpositive includes all retained contigs
	 * @return Long sum of retained contig lengths at least minimum
	 */
	long lengthAtLeast(int minimum){
		long sum=0;
		for(int i=0; i<contigs.size; i++){
			int len=contigs.get(i);
			if(len<minimum){break;}
			sum+=len;
		}
		return sum;
	}
	
	/** Retained input name, assigned before the constructor invokes load; reload uses the same name. */
	final String fname;
	/** Per-load nonempty-record limit; negative is unlimited and zero skips reader iteration. */
	private final long maxRecords;
	/** Owned, reused list of positive sequence-byte lengths; sorted descending after successful load. */
	IntList contigs=new IntList();
	/** Long sum of retained contig lengths after success; zeroed before replacement parsing. */
	long length=0;
	/** Raw payload-byte total over included headers, excluding a next header that ends limited loading. */
	long headerLength=0;
	/** First header with at least one payload byte, decoded with the default charset; otherwise null.
	 * Empty headers do not set it; whitespace payload is retained. */
	String firstHeader=null;
	/**
	 * Counts in A,C,G,T,U,N,IUPAC,OTHER order. Allocated anew during load; null after clear.
	 */
	long[] acgtnio;
	
	/**
	 * Shared 128-entry ASCII-to-count-index lookup. The reference is final; array elements remain mutable.
	 */
	public static final byte[] baseToACGTNIO=makeBaseToACGTUNIO();
	/** Composition indices: A=0,C=1,G=2,T=3,U=4,N=5,IUPAC=6,OTHER=7. */
	private static final byte A=0, C=1, G=2, T=3, U=4, N=5, IUPAC=6, OTHER=7;
	
	/**
	 * Creates the ASCII lookup, initially OTHER, then assigns A,C,G,T,U,N in both cases.
	 * Remaining entries recognized by AminoAcid.baseToNumberExtended become IUPAC.
	 * @return Fresh 128-entry array of composition indices
	 */
	private static final byte[] makeBaseToACGTUNIO(){
		final byte[] array=new byte[128];
		Arrays.fill(array, OTHER);
		array['a']=array['A']=A;
		array['c']=array['C']=C;
		array['g']=array['G']=G;
		array['t']=array['T']=T;
		array['u']=array['U']=U;
		array['n']=array['N']=N;
		for(int i=0; i<array.length; i++){
			if(AminoAcid.baseToNumberExtended[i]>=0 && array[i]==OTHER){
				array[i]=IUPAC;
			}
		}
		return array;
	}
	
}
