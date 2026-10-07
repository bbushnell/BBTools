package stream.bam;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import map.IntObjectMap;
import shared.Timer;
import structures.BinaryByteWrapperLE;
import structures.LongList;

/**
 * Builds BAM index files (.bai) from coordinate-sorted BAM files.
 * 
 * <p>BAI format structure:
 * <ul>
 * <li>Magic bytes: "BAI\1"
 * <li>For each reference sequence:
 *   <ul>
 *   <li>Binning index: hierarchical bins (16kb to 512Mb) containing chunk lists
 *   <li>Linear index: 16kb-resolution array of file offsets
 *   <li>Optional pseudo-bin 37450 with summary statistics
 *   </ul>
 * <li>Count of unplaced unmapped reads
 * </ul>
 * 
 * <p>Implementation performs a single sequential pass over the BAM file,
 * building index structures on demand for references encountered in records.
 * Memory grows with reference slots, accumulated bin chunks, and linear-index windows.
 * 
 * <p>Historical performance estimate: about 10 seconds for a 10M-read BAM.
 * The original estimate does not identify its hardware or workload details.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 6, 2025
 */
public final class BamIndexWriter{

	/**
	 * Command-line entry point for indexing BAM files.
	 * Prints caught I/O exceptions and elapsed time; does not propagate those exceptions.
	 * 
	 * @param args At least one argument: input.bam, then optional output.bai
	 * (defaults to input.bam.bai); later arguments are ignored
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		try{
			if(args.length<2){writeIndex(args[0]);}else{
				writeIndex(args[0], args[1]);
			}
		//TODO: Possible bug [stream/bam/BamIndexWriter#003] - caught I/O errors are
		//printed but main returns normally, so callers cannot rely on failure status.
		}catch(IOException e){
			e.printStackTrace();
		}
		t.stopAndPrint();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this static utility. */
	private BamIndexWriter(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Write a .bai index next to the BAM file (adds ".bai" suffix).
	 * 
	 * @param bamPath Path to coordinate-sorted BAM file
	 * @throws IOException if BAM cannot be read or index cannot be written
	 */
	public static void writeIndex(String bamPath) throws IOException{
		writeIndex(bamPath, bamPath+".bai");
	}

	/**
	 * Write a .bai index to an explicit destination.
	 * 
	 * <p>Algorithm:
	 * <ol>
	 * <li>Validate BAM magic and check the declared sort order in the first at most 100 header bytes
	 * <li>Read reference dictionary to determine capacity hints
	 * <li>Stream through alignments, building bin and linear indices
	 * <li>Write completed index in BAI format
	 * </ol>
	 * 
	 * Input records must already be coordinate sorted; their order is not checked here.
	 * Placed-unmapped and zero-reference-span records occupy one base for binning and linear windows.
	 * For nonempty text, assertions require the inspected first-line prefix to start
	 * with @HD and contain SO:coordinate.
	 * @param bamPath Path to coordinate-sorted BAM file
	 * @param indexPath Output .bai path; opened for replacement before input validation
	 * @throws IOException If input cannot be read, a checked record field is invalid,
	 * or output cannot be written
	 * @throws AssertionError If the inspected header prefix fails the declared-sort
	 * check and assertions are enabled
	 */
	public static void writeIndex(String bamPath, String indexPath) throws IOException{
		try(FileInputStream fis=new FileInputStream(bamPath);
			BufferedInputStream bis=new BufferedInputStream(fis, 65536); //Buffer file I/O
			BgzfInputStream bgzf=new BgzfInputStream(bis); //BGZF decompression
			FileOutputStream fos0=new FileOutputStream(indexPath);
			BufferedOutputStream fos=new BufferedOutputStream(fos0); //Buffer index writes
			){

			BamWriterHelper writer=new BamWriterHelper(fos);
			BamReader reader=new BamReader(bgzf);

			//Validate BAM magic bytes
			byte[] magic=reader.readBytes(4);
			if(magic[0]!='B' || magic[1]!='A' || magic[2]!='M' || magic[3]!=1){throw new IOException("Input is not a BAM file: "+bamPath);}

			//Read the header prefix to check its declared sort order
			long lText=reader.readUint32();
			if(lText<0 || lText>Integer.MAX_VALUE){throw new IOException("Invalid BAM header length: "+lText);}
			if(lText>0){
				int checkLen=(int)Math.min(lText, 100); //Inspect at most 100 header bytes
				byte[] headerStart=reader.readBytes(checkLen);

				//Find end of first line (@HD line if present)
				int newline=0;
				while(newline<headerStart.length && headerStart[newline]!='\n'){newline++;}
				String firstLine=new String(headerStart, 0, newline, java.nio.charset.StandardCharsets.US_ASCII);

				//Require the inspected prefix to declare coordinate sorting
				assert(firstLine.startsWith("@HD") && firstLine.contains("SO:coordinate")) : 
					"BAM file must be coordinate-sorted (SO:coordinate) for indexing: "+bamPath
					+"\nAdd -da to override this warning.";

				//Skip rest of header text
				long remaining=lText-checkLen;
				if(remaining>0){reader.readBytes((int)remaining);}
			}

			//Read reference dictionary
			int nRef=reader.readInt32();
			if(nRef<0){throw new IOException("Negative reference count in BAM header");}

			//Calculate initial bin capacity based on reference count
			//More references = more fragmented data = smaller bins on average
			int binCapacity=256; //Default for small genomes
			if(nRef>8192){binCapacity=128;}
			if(nRef>32768){binCapacity=64;}
			if(nRef>131072){binCapacity=32;}

			//Read reference names/lengths but don't create indices yet (lazy allocation)
			ReferenceIndex[] references=new ReferenceIndex[nRef];
			for(int i=0; i<nRef; i++){
				long lName=reader.readUint32();
				if(lName<1 || lName>Integer.MAX_VALUE){throw new IOException("Invalid reference name length: "+lName);}
				reader.readBytes((int)lName); //Name (includes terminating NUL)
				reader.readUint32(); //Reference length (unused here)
				//ReferenceIndex created on-demand when first read for this reference is seen
			}

			long readsWithoutCoordinate=0L; //Count of unmapped reads with no RNAME

			//Reusable buffer for record parsing (grows as needed)
			BinaryByteWrapperLE bb=new BinaryByteWrapperLE(new byte[256]);

			//Stream through all alignment records
			while(true){
				long recordStart=bgzf.getVirtualOffset(); //BGZF virtual offset before record
				int blockSize;
				try{
					blockSize=reader.readInt32(); //Record size in bytes
				}catch(EOFException eof){
					break; //Reached EOF marker block
				}

				if(blockSize<0){throw new IOException("Negative BAM record block size");}

				byte[] recordData=reader.readBytes(blockSize);
				long recordEnd=bgzf.getVirtualOffset(); //BGZF virtual offset after record

				if(recordData.length<FIXED_RECORD_FIELDS){throw new IOException("Corrupted BAM record: truncated fixed fields");}

				//Wrap record for efficient little-endian parsing
				if(recordData.length>bb.array.length){
					bb.wrap(recordData); //Grows internal buffer
				}else{
					System.arraycopy(recordData, 0, bb.array, 0, recordData.length); //Reuse buffer
					bb.wrap(bb.array, 0, recordData.length);
				}

				//Parse fixed fields (32 bytes total)
				int refID=bb.getInt(); //Reference sequence ID (-1 when no reference is assigned)
				int pos=bb.getInt(); //0-based leftmost position (-1 if unavailable)
				int lReadName=bb.get()&0xFF; //Length of QNAME including NUL
				bb.get(); //MAPQ (unused for indexing)
				bb.getShort(); //stored bin field - IGNORED, recomputed via reg2bin below (#001 fix)
				int nCigar=bb.getShort()&0xFFFF; //Number of CIGAR operations
				int flag=bb.getShort()&0xFFFF; //SAM flags
				int lSeq=bb.getInt(); //Sequence length (unused here)
				bb.getInt(); //Next reference ID (unused)
				bb.getInt(); //Next position (unused)
				bb.getInt(); //Template length (unused)

				//Track unmapped reads without coordinates
				if(refID<0){
					readsWithoutCoordinate++;
					continue;
				}
				if(refID>=references.length){throw new IOException("Reference id out of bounds: "+refID+" >= "+references.length);}

				//Create ReferenceIndex on first read for this reference (lazy allocation)
				ReferenceIndex ref=references[refID];
				if(ref==null){
					ref=new ReferenceIndex(binCapacity);
					references[refID]=ref;
				}

				ref.incrementCounts(flag); //Track mapped/unmapped counts

				//Reads without valid positions don't contribute to spatial index
				if(pos<0){continue;}

				//Skip QNAME field
				if(lReadName>bb.remaining()){throw new IOException("Corrupted BAM record: read name exceeds record size");}
				bb.skip(lReadName);

				//Decode CIGAR to calculate reference span
				int cigarBytes=nCigar*4;
				if(cigarBytes>bb.remaining()){throw new IOException("Corrupted BAM record: CIGAR exceeds record size");}
				int refSpan=0;
				for(int c=0; c<nCigar; c++){
					int cigarOp=bb.getInt(); //Encoded as (length<<4)|op
					refSpan+=referenceSpanContribution(cigarOp);
				}

				//No need to parse SEQ/QUAL/AUX - all index info extracted

				//[stream/bam/BamIndexWriter#001] FIXED 2026-06-20 (greenlit by Brian): RECOMPUTE the bin via
				//reg2bin(pos, end) instead of trusting the record's stored `bin` field. A bin=0 / stale-bin
				//BAM (a legal placeholder - the SAM spec marks bin derivable) would otherwise get a silently
				//WRONG .bai (every chunk in bin 0 → region queries miss reads). samtools recomputes for the
				//same reason; refSpan is already in hand, so this is free robustness.
				//STR380: SAMv1 4.2.1 treats unmapped records as length one regardless of retained CIGAR.
				int alignmentEndExclusive=pos+((flag&BAM_FUNMAP)==0 ? Math.max(refSpan, 1) : 1);
				int computedBin=reg2bin(pos, alignmentEndExclusive);

				//Add alignment to bin index
				ref.addRecord(computedBin, recordStart, recordEnd);

				//Update linear index for this alignment's span
				ref.updateLinearIndex(pos, alignmentEndExclusive, recordStart);
			}

			//Write BAI file format
			writer.writeBytes(new byte[]{'B', 'A', 'I', 1}); //Magic
			writer.writeUint32(nRef); //Number of reference sequences

			//Write index for each reference
			for(int i=0; i<nRef; i++){
				ReferenceIndex ref=references[i];
				if(ref==null){
					ref=new ReferenceIndex(binCapacity); //Empty reference (no aligned reads)
				}

				//Write binning index
				int binCount=ref.binCount()+(ref.shouldEmitPseudoBin() ? 1 : 0);
				writer.writeUint32(binCount);

				//Write regular bins
				int[] binKeys=ref.bins.keys();
				for(int j=0; j<binKeys.length; j++){
					int binKey=binKeys[j];
					BinData data=ref.bins.get(binKey);
					if(data!=null){
						writer.writeUint32(binKey);
						writer.writeUint32(data.size()); //Number of chunks
						for(int k=0; k<data.size(); k++){
							writer.writeUint64(data.begList.get(k)); //Chunk start offset
							writer.writeUint64(data.endList.get(k)); //Chunk end offset
						}
					}
				}

				//Write pseudo-bin 37450 when an offset range was recorded
				if(ref.shouldEmitPseudoBin()){ref.writePseudoBin(writer);}

				//Write linear index
				LongList linear=ref.linear;
				int linearSize=linear.size();
				writer.writeUint32(linearSize);
				for(int j=0; j<linearSize; j++){
					long offset=linear.get(j);
					if(offset<0){offset=0;} //Unset entries become 0
					writer.writeUint64(offset);
				}
			}

			//Write count of unplaced unmapped reads
			writer.writeUint64(readsWithoutCoordinate);
		}
	}

	/**
	 * Calculate reference bases consumed by a CIGAR operation.
	 * 
	 * @param cigarEncoded CIGAR operation encoded as {@code (length<<4)|op}
	 * @return Encoded length for M, D, N, =, or X; zero for all other operation codes
	 */
	private static int referenceSpanContribution(int cigarEncoded){
		int op=cigarEncoded&0xF; //Bottom 4 bits = operation
		int len=cigarEncoded>>>4; //Top 28 bits = length
		switch(op){
			case 0: //M (match/mismatch)
			case 2: //D (deletion)
			case 3: //N (skipped region)
			case 7: //= (sequence match)
			case 8: //X (sequence mismatch)
				return len;
			default: //I, S, H, P (don't consume reference)
				return 0;
		}
	}

	/**
	 * Calculate the BAM bin for a 0-based region [beg, end).
	 * Original algorithm attribution: SAMv1.pdf section 5.3. Recomputes the bin per
	 * record (historical repair #001) instead of trusting its stored bin field.
	 * Uses the same arithmetic as SamToBamConverter.reg2bin.
	 * @param beg Inclusive start within the binning scheme
	 * @param end Exclusive end, decremented before comparing shifted endpoints
	 * @return First matching bin, or zero if no finer bin matches
	 */
	private static int reg2bin(int beg, int end){
		--end;
		if(beg>>14==end>>14){return ((1<<15)-1)/7+(beg>>14);}
		if(beg>>17==end>>17){return ((1<<12)-1)/7+(beg>>17);}
		if(beg>>20==end>>20){return ((1<<9)-1)/7+(beg>>20);}
		if(beg>>23==end>>23){return ((1<<6)-1)/7+(beg>>23);}
		if(beg>>26==end>>26){return ((1<<3)-1)/7+(beg>>26);}
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Index data for a single reference sequence.
	 * Contains binning index (hierarchical bins) and linear index (16kb windows).
	 */
	private static final class ReferenceIndex{
		
		/**
		 * Creates empty bin and linear indices with zero counts and unset offsets.
		 * @param binCapacity Initial capacity hint for the bin map
		 */
		ReferenceIndex(int binCapacity){
			this.linear=new LongList(16); //Grows as needed
			this.bins=new IntObjectMap<BinData>(binCapacity, BinData.class); //Sized based on genome fragmentation
		}

		/**
		 * Increment mapped/unmapped read counts based on SAM flags.
		 * @param flag SAM FLAG field
		 */
		void incrementCounts(int flag){
			if((flag&BAM_FUNMAP)==0){ //0x4 bit clear = mapped
				mappedReads++;
			}else{
				unmappedReads++;
			}
		}

		/**
		 * Add an alignment record to the bin index.
		 * Merges with the most recent chunk in that bin when adjacent or overlapping.
		 * Also expands the reference-wide minimum and maximum virtual offsets.
		 * 
		 * @param bin Recomputed bin number for the record
		 * @param start BGZF virtual offset at start of record
		 * @param end BGZF virtual offset after record
		 */
		void addRecord(int bin, long start, long end){
			BinData data=bins.get(bin);
			if(data==null){
				data=new BinData();
				bins.put(bin, data);
			}
			data.append(start, end); //May merge with previous chunk
			
			//Track global min/max offsets for pseudo-bin
			if(firstOffset<0 || start<firstOffset){firstOffset=start;}
			if(end>lastOffset){lastOffset=end;}
		}

		/**
		 * Update linear index for an alignment's reference span.
		 * Sets only previously unset 16kb windows overlapped by this alignment.
		 * Requires records to arrive in nondecreasing virtual-offset order.
		 * Returns without modification for a negative start position.
		 * 
		 * @param pos 0-based alignment start position
		 * @param endExclusive Alignment end position (exclusive)
		 * @param offset BGZF virtual offset at start of record
		 */
		void updateLinearIndex(int pos, int endExclusive, long offset){
			if(pos<0){return;}
			int linearBegin=pos>>LINEAR_INDEX_SHIFT; //Divide by 16384
			int linearEnd=Math.max(pos, endExclusive-1)>>LINEAR_INDEX_SHIFT;
			ensureLinearSize(linearEnd+1);
			//With nondecreasing supplied offsets, the first record touching a window
			//supplies its earliest offset. The header assertion checks the declared sort
			//order, not record order; coordinate-sorted input remains a caller precondition.
			for(int i=linearBegin; i<=linearEnd; i++){
				if(linear.get(i)==UNSET_OFFSET){ //Only set if unset (want earliest offset)
					linear.set(i, offset);
				}
			}
		}

		/** @return Number of bins with recorded chunks */
		int binCount(){return bins.size();}

		/** @return true if addRecord established a nonnegative, ordered offset range */
		//TODO: Possible bug [stream/bam/BamIndexWriter#002] - firstOffset is set only
		//by addRecord, which is skipped for pos<0. A reference with only refID>=0,
		//pos<0 records accumulates counts but emits no pseudo-bin. A positioned record
		//allows all accumulated counts to be emitted. Prior review proposed testing
		//counts instead and writing 0/0 offsets; format semantics and reader results
		//for that change remain unverified. No repair is made here.
		boolean shouldEmitPseudoBin(){
			return firstOffset>=0 && lastOffset>=firstOffset;
		}

		/**
		 * Write pseudo-bin 37450 containing summary statistics.
		 * Format: bin_id=37450, n_chunk=2, ref_beg, ref_end, n_mapped, n_unmapped.
		 * @param writer Destination for the six metadata fields; remains open
		 * @throws IOException If writing any field fails
		 */
		void writePseudoBin(BamWriterHelper writer) throws IOException{
			writer.writeUint32(PSEUDO_BIN); //Bin 37450
			writer.writeUint32(2); //Always 2 "chunks" (actually metadata fields)
			writer.writeUint64(firstOffset); //Earliest record offset
			writer.writeUint64(lastOffset); //Latest record offset
			writer.writeUint64(mappedReads); //Count of mapped reads
			writer.writeUint64(unmappedReads); //Count of unmapped reads
		}

		/**
		 * Appends unset sentinels until the linear index has at least size entries.
		 * @param size Minimum desired number of windows
		 */
		private void ensureLinearSize(int size){
			while(linear.size()<size){
				linear.add(UNSET_OFFSET);
			}
		}

		/** Bin number to accumulated chunk list. */
		private final IntObjectMap<BinData> bins;
		/** Earliest recorded offset per 16kb window, or UNSET_OFFSET. */
		private final LongList linear;
		/** Records counted with BAM_FUNMAP clear, including those without a position. */
		private long mappedReads=0L;
		/** Records counted with BAM_FUNMAP set, including those without a position. */
		private long unmappedReads=0L;
		/** Minimum start passed to addRecord; negative until a chunk is recorded. */
		private long firstOffset=-1L;
		/** Maximum end passed to addRecord; negative until a chunk is recorded. */
		private long lastOffset=-1L;
	}

	/**
	 * Chunk list for a single bin.
	 * Uses two parallel LongLists instead of {@code ArrayList<Chunk>} to reduce object overhead.
	 * Automatically merges adjacent/overlapping chunks on append.
	 */
	private static final class BinData{

		/**
		 * Appends a chunk, merging with the last chunk if they overlap or are adjacent.
		 * Chunks are added in file order, so we only check the last chunk for merging.
		 * 
		 * @param start BGZF virtual offset at start of chunk
		 * @param end BGZF virtual offset at end of chunk
		 */
		void append(long start, long end){
			if(begList.size()==0 || start>endList.get(endList.size()-1)){ //New non-overlapping chunk
				begList.add(start);
				endList.add(end);
			}else{ //Overlaps previous chunk - extend it
				int lastIdx=endList.size()-1;
				endList.set(lastIdx, Math.max(endList.get(lastIdx), end));
			}
		}

		/** @return Number of chunks in this bin */
		int size(){return begList.size();}

		/** Chunk start offsets, parallel to endList. */
		final LongList begList=new LongList(4);
		/** Chunk end offsets, parallel to begList. */
		final LongList endList=new LongList(4);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	/** Bytes in the fixed BAM record fields. */
	private static final int FIXED_RECORD_FIELDS=32;
	/** Bit shift selecting a 16384-base linear-index window. */
	private static final int LINEAR_INDEX_SHIFT=14;
	/** SAM flag bit identifying an unmapped record. */
	private static final int BAM_FUNMAP=0x4;
	/** Bin identifier used for per-reference summary metadata. */
	private static final int PSEUDO_BIN=37450;
	/** Sentinel for an unassigned linear-index entry. */
	private static final long UNSET_OFFSET=-1L;
}
