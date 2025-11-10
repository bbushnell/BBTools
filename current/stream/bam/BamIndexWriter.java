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
 * building index structures on-demand for references with aligned reads.
 * Memory usage scales with number of bins and reference sequences, not 
 * total read count.
 * 
 * <p>Performance: ~10 seconds for 10M read BAM on typical hardware.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 6, 2025
 */
public final class BamIndexWriter{

	/**
	 * Command-line entry point for indexing BAM files.
	 * 
	 * @param args [0]=input.bam, [1]=output.bai (optional, defaults to input.bam.bai)
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		try{
			if(args.length<2){
				writeIndex(args[0]);
			}else{
				writeIndex(args[0], args[1]);
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		t.stopAndPrint();
	}

	private BamIndexWriter(){} //Prevent instantiation

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
	 * <li>Validate BAM magic and sort order (SO:coordinate required)
	 * <li>Read reference dictionary to determine capacity hints
	 * <li>Stream through alignments, building bin and linear indices
	 * <li>Write completed index in BAI format
	 * </ol>
	 * 
	 * @param bamPath Path to coordinate-sorted BAM file
	 * @param indexPath Path for output .bai file
	 * @throws IOException if BAM cannot be read, is unsorted, or index cannot be written
	 * @throws AssertionError if BAM is not coordinate-sorted (disable with -da)
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
			if(magic[0]!='B' || magic[1]!='A' || magic[2]!='M' || magic[3]!=1){
				throw new IOException("Input is not a BAM file: "+bamPath);
			}

			//Read start of header to validate sort order
			long lText=reader.readUint32();
			if(lText<0 || lText>Integer.MAX_VALUE){
				throw new IOException("Invalid BAM header length: "+lText);
			}
			if(lText>0){
				int checkLen=(int)Math.min(lText, 100); //Only need first line
				byte[] headerStart=reader.readBytes(checkLen);

				//Find end of first line (@HD line if present)
				int newline=0;
				while(newline<headerStart.length && headerStart[newline]!='\n'){newline++;}
				String firstLine=new String(headerStart, 0, newline, java.nio.charset.StandardCharsets.US_ASCII);

				//Require coordinate sorting for indexing
				assert(firstLine.startsWith("@HD") && firstLine.contains("SO:coordinate")) : 
					"BAM file must be coordinate-sorted (SO:coordinate) for indexing: "+bamPath
					+"\nAdd -da to override this warning.";

				//Skip rest of header text
				long remaining=lText-checkLen;
				if(remaining>0){
					reader.readBytes((int)remaining);
				}
			}

			//Read reference dictionary
			int nRef=reader.readInt32();
			if(nRef<0){
				throw new IOException("Negative reference count in BAM header");
			}

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
				if(lName<1 || lName>Integer.MAX_VALUE){
					throw new IOException("Invalid reference name length: "+lName);
				}
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

				if(blockSize<0){
					throw new IOException("Negative BAM record block size");
				}

				byte[] recordData=reader.readBytes(blockSize);
				long recordEnd=bgzf.getVirtualOffset(); //BGZF virtual offset after record

				if(recordData.length<FIXED_RECORD_FIELDS){
					throw new IOException("Corrupted BAM record: truncated fixed fields");
				}

				//Wrap record for efficient little-endian parsing
				if(recordData.length>bb.array.length){
					bb.wrap(recordData); //Grows internal buffer
				}else{
					System.arraycopy(recordData, 0, bb.array, 0, recordData.length); //Reuse buffer
					bb.wrap(bb.array, 0, recordData.length);
				}

				//Parse fixed fields (32 bytes total)
				int refID=bb.getInt(); //Reference sequence ID (-1 for unmapped)
				int pos=bb.getInt(); //0-based leftmost position (-1 if unavailable)
				int lReadName=bb.get()&0xFF; //Length of QNAME including NUL
				bb.get(); //MAPQ (unused for indexing)
				int bin=bb.getShort()&0xFFFF; //Binning index bin number
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
				if(refID>=references.length){
					throw new IOException("Reference id out of bounds: "+refID+" >= "+references.length);
				}

				//Create ReferenceIndex on first read for this reference (lazy allocation)
				ReferenceIndex ref=references[refID];
				if(ref==null){
					ref=new ReferenceIndex(binCapacity);
					references[refID]=ref;
				}

				ref.incrementCounts(flag); //Track mapped/unmapped counts

				//Reads without valid positions don't contribute to spatial index
				if(pos<0){
					continue;
				}

				//Skip QNAME field
				if(lReadName>bb.remaining()){
					throw new IOException("Corrupted BAM record: read name exceeds record size");
				}
				bb.skip(lReadName);

				//Decode CIGAR to calculate reference span
				int cigarBytes=nCigar*4;
				if(cigarBytes>bb.remaining()){
					throw new IOException("Corrupted BAM record: CIGAR exceeds record size");
				}
				int refSpan=0;
				for(int c=0; c<nCigar; c++){
					int cigarOp=bb.getInt(); //Encoded as (length<<4)|op
					refSpan+=referenceSpanContribution(cigarOp);
				}

				//No need to parse SEQ/QUAL/AUX - all index info extracted

				//Add alignment to bin index
				ref.addRecord(bin, recordStart, recordEnd);

				//Update linear index for this alignment's span
				int alignmentEndExclusive=pos+Math.max(refSpan, 1);
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

				//Write pseudo-bin 37450 if reference has alignments
				if(ref.shouldEmitPseudoBin()){
					ref.writePseudoBin(writer);
				}

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
	 * @param cigarEncoded CIGAR operation encoded as (length<<4)|op
	 * @return Number of reference bases consumed (0 for insertions, soft clips, etc.)
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
	 * Index data for a single reference sequence.
	 * Contains binning index (hierarchical bins) and linear index (16kb windows).
	 */
	private static final class ReferenceIndex{
		
		ReferenceIndex(int binCapacity){
			this.linear=new LongList(16); //Grows as needed
			this.bins=new IntObjectMap<BinData>(binCapacity); //Sized based on genome fragmentation
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
		 * Merges adjacent or overlapping chunks within the same bin.
		 * 
		 * @param bin Bin number from BAM record
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
		 * Sets file offset for all 16kb windows overlapped by this alignment.
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
			for(int i=linearBegin; i<=linearEnd; i++){
				if(linear.get(i)==UNSET_OFFSET){ //Only set if unset (want earliest offset)
					linear.set(i, offset);
				}
			}
		}

		/** @return Number of bins with recorded chunks */
		int binCount(){
			return bins.size();
		}

		/** @return true if this reference has at least one aligned read */
		boolean shouldEmitPseudoBin(){
			return firstOffset>=0 && lastOffset>=firstOffset;
		}

		/**
		 * Write pseudo-bin 37450 containing summary statistics.
		 * Format: bin_id=37450, n_chunk=2, ref_beg, ref_end, n_mapped, n_unmapped
		 */
		void writePseudoBin(BamWriterHelper writer) throws IOException{
			writer.writeUint32(PSEUDO_BIN); //Bin 37450
			writer.writeUint32(2); //Always 2 "chunks" (actually metadata fields)
			writer.writeUint64(firstOffset); //Earliest record offset
			writer.writeUint64(lastOffset); //Latest record offset
			writer.writeUint64(mappedReads); //Count of mapped reads
			writer.writeUint64(unmappedReads); //Count of unmapped reads
		}

		/** Ensure linear index has at least 'size' entries */
		private void ensureLinearSize(int size){
			while(linear.size()<size){
				linear.add(UNSET_OFFSET);
			}
		}

		private final IntObjectMap<BinData> bins; //Map from bin number to chunk list
		private final LongList linear; //16kb-resolution file offset array
		private long mappedReads=0L; //Count for pseudo-bin
		private long unmappedReads=0L; //Count for pseudo-bin
		private long firstOffset=-1L; //Earliest record offset for pseudo-bin
		private long lastOffset=-1L; //Latest record offset for pseudo-bin
	}

	/**
	 * Chunk list for a single bin.
	 * Uses two parallel LongLists instead of ArrayList<Chunk> to reduce object overhead.
	 * Automatically merges adjacent/overlapping chunks on append.
	 */
	private static final class BinData{

		/**
		 * Append a new chunk, merging with the previous chunk if they overlap.
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
		int size(){
			return begList.size();
		}

		final LongList begList=new LongList(4); //Chunk start offsets
		final LongList endList=new LongList(4); //Chunk end offsets
	}

	//Constants
	private static final int FIXED_RECORD_FIELDS=32; //Size of fixed BAM record header
	private static final int LINEAR_INDEX_SHIFT=14; //log2(16384) for 16kb windows
	private static final int BAM_FUNMAP=0x4; //SAM FLAG bit for unmapped
	private static final int PSEUDO_BIN=37450; //Special bin for metadata
	private static final long UNSET_OFFSET=-1L; //Sentinel for unset linear index entries
}