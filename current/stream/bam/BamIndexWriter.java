package stream.bam;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import shared.Timer;
import structures.BinaryByteWrapperLE;
import structures.IntObjectMap;
import structures.LongList;

/**
 * Builds BAM index files (.bai) from an existing BAM.
 * This implementation performs a single sequential pass over the BAM file,
 * collecting bin and linear indexing information per reference.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 5, 2025
 */
public final class BamIndexWriter {

	public static void main(String[] args) {
		Timer t=new Timer();
		try{
			if(args.length<2) {
				writeIndex(args[0]);
			}else {
				writeIndex(args[0], args[1]);
			}
		}catch(IOException e){
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		t.stopAndPrint();
	}

	private BamIndexWriter(){}

	/**
	 * Write a .bai index next to the BAM file (adds ".bai" suffix).
	 */
	public static void writeIndex(String bamPath) throws IOException {
		writeIndex(bamPath, bamPath+".bai");
	}

	/**
	 * Write a .bai index to an explicit destination.
	 */
	public static void writeIndex(String bamPath, String indexPath) throws IOException {
		try(FileInputStream fis=new FileInputStream(bamPath);
			BufferedInputStream bis=new BufferedInputStream(fis, 65536);
			BgzfInputStream bgzf=new BgzfInputStream(bis);
			FileOutputStream fos0=new FileOutputStream(indexPath);
			BufferedOutputStream fos=new BufferedOutputStream(fos0);
			){

			BamWriterHelper writer=new BamWriterHelper(fos);
			BamReader reader=new BamReader(bgzf);

			// Validate BAM magic
			byte[] magic=reader.readBytes(4);
			if(magic[0]!='B' || magic[1]!='A' || magic[2]!='M' || magic[3]!=1){
				throw new IOException("Input is not a BAM file: "+bamPath);
			}

			// Skip header text
			long lText=reader.readUint32();
			if(lText<0 || lText>Integer.MAX_VALUE){
				throw new IOException("Invalid BAM header length: "+lText);
			}
			if(lText>0){
				reader.readBytes((int)lText);
			}

			// Read reference dictionary
			int nRef=reader.readInt32();
			if(nRef<0){
				throw new IOException("Negative reference count in BAM header");
			}

			ReferenceIndex[] references=new ReferenceIndex[nRef];
			for(int i=0; i<nRef; i++){
				long lName=reader.readUint32();
				if(lName<1 || lName>Integer.MAX_VALUE){
					throw new IOException("Invalid reference name length: "+lName);
				}
				reader.readBytes((int)lName); // includes terminating NUL
				reader.readUint32(); // reference length (unused, but read to advance)
				// Don't create ReferenceIndex here - wait until we see reads
			}

			long readsWithoutCoordinate=0L;

			// Reusable buffer for record parsing
			BinaryByteWrapperLE bb=new BinaryByteWrapperLE(new byte[256]);

			// Iterate over alignment records
			while(true){
				long recordStart=bgzf.getVirtualOffset();
				int blockSize;
				try{
					blockSize=reader.readInt32();
				}catch(EOFException eof){
					break; // reached BGZF EOF block
				}

				if(blockSize<0){
					throw new IOException("Negative BAM record block size");
				}

				byte[] recordData=reader.readBytes(blockSize);
				long recordEnd=bgzf.getVirtualOffset();

				if(recordData.length<FIXED_RECORD_FIELDS){
					throw new IOException("Corrupted BAM record: truncated fixed fields");
				}

				// Wrap record data for parsing
				if(recordData.length>bb.array.length){
					bb.wrap(recordData);
				}else{
					System.arraycopy(recordData, 0, bb.array, 0, recordData.length);
					bb.wrap(bb.array, 0, recordData.length);
				}

				int refID=bb.getInt();
				int pos=bb.getInt(); // 0-based, -1 if unavailable
				int lReadName=bb.get()&0xFF;
				bb.get(); // mapq (unused here)
				int bin=bb.getShort()&0xFFFF;
				int nCigar=bb.getShort()&0xFFFF;
				int flag=bb.getShort()&0xFFFF;
				int lSeq=bb.getInt();
				bb.getInt(); // nextRefID
				bb.getInt(); // pnext
				bb.getInt(); // tlen

				if(refID<0){
					readsWithoutCoordinate++;
					continue;
				}
				if(refID>=references.length){
					throw new IOException("Reference id out of bounds: "+refID+" >= "+references.length);
				}

				// Create ReferenceIndex on-demand
				ReferenceIndex ref=references[refID];
				if(ref==null){
					ref=new ReferenceIndex();
					references[refID]=ref;
				}

				ref.incrementCounts(flag);

				if(pos<0){
					// No coordinate available; nothing to index beyond count tracking
					continue;
				}

				// Skip QNAME
				if(lReadName>bb.remaining()){
					throw new IOException("Corrupted BAM record: read name exceeds record size");
				}
				bb.skip(lReadName);

				// Decode CIGAR to determine reference span
				int cigarBytes=nCigar*4;
				if(cigarBytes>bb.remaining()){
					throw new IOException("Corrupted BAM record: CIGAR exceeds record size");
				}
				int refSpan=0;
				for(int c=0; c<nCigar; c++){
					int cigarOp=bb.getInt();
					refSpan+=referenceSpanContribution(cigarOp);
				}

				// Don't need to parse SEQ/QUAL/AUX - we have all the info we need

				ref.addRecord(bin, recordStart, recordEnd);

				int alignmentEndExclusive=pos+Math.max(refSpan, 1);
				ref.updateLinearIndex(pos, alignmentEndExclusive, recordStart);
			}

			// Write BAI header
			writer.writeBytes(new byte[]{'B', 'A', 'I', 1});
			writer.writeUint32(nRef);

			for(int i=0; i<nRef; i++){
				ReferenceIndex ref=references[i];
				if(ref==null){
					ref=new ReferenceIndex(); // Empty reference
				}

				// Bins
				int binCount=ref.binCount()+(ref.shouldEmitPseudoBin() ? 1 : 0);
				writer.writeUint32(binCount);

				// Write regular bins
				int[] binKeys=ref.bins.keys();
				for(int j=0; j<binKeys.length; j++){
					int binKey=binKeys[j];
					BinData data=ref.bins.get(binKey);
					if(data!=null){
						List<Chunk> chunks=data.chunks;
						writer.writeUint32(binKey);
						writer.writeUint32(chunks.size());
						for(Chunk chunk : chunks){
							writer.writeUint64(chunk.beg);
							writer.writeUint64(chunk.end);
						}
					}
				}

				// Write pseudo-bin if needed
				if(ref.shouldEmitPseudoBin()){
					ref.writePseudoBin(writer);
				}

				// Linear index
				LongList linear=ref.linear;
				int linearSize=linear.size();
				writer.writeUint32(linearSize);
				for(int j=0; j<linearSize; j++){
					long offset=linear.get(j);
					if(offset<0){
						offset=0;
					}
					writer.writeUint64(offset);
				}
			}

			writer.writeUint64(readsWithoutCoordinate);
		}
	}

	private static int referenceSpanContribution(int cigarEncoded){
		int op=cigarEncoded&0xF;
		int len=cigarEncoded>>>4;
				switch(op){
					case 0: // M
					case 2: // D
					case 3: // N
					case 7: // =
					case 8: // X
						return len;
					default:
						return 0;
				}
	}

	private static final class ReferenceIndex {
		ReferenceIndex(){
			this.linear=new LongList(16);
			this.bins=new IntObjectMap<BinData>();
		}

		void incrementCounts(int flag){
			if((flag&BAM_FUNMAP)==0){
				mappedReads++;
			}else{
				unmappedReads++;
			}
		}

		void addRecord(int bin, long start, long end){
			BinData data=bins.get(bin);
			if(data==null){
				data=new BinData();
				bins.put(bin, data);
			}
			data.append(start, end);
			if(firstOffset<0 || start<firstOffset){
				firstOffset=start;
			}
			if(end>lastOffset){
				lastOffset=end;
			}
		}

		void updateLinearIndex(int pos, int endExclusive, long offset){
			if(pos<0){
				return;
			}
			int linearBegin=pos>>LINEAR_INDEX_SHIFT;
			int linearEnd=Math.max(pos, endExclusive-1)>>LINEAR_INDEX_SHIFT;
			ensureLinearSize(linearEnd+1);
			for(int i=linearBegin; i<=linearEnd; i++){
				if(linear.get(i)==UNSET_OFFSET){
					linear.set(i, offset);
				}
			}
		}

		int binCount(){
			return bins.size();
		}

		boolean shouldEmitPseudoBin(){
			return firstOffset>=0 && lastOffset>=firstOffset;
		}

		void writePseudoBin(BamWriterHelper writer) throws IOException {
			writer.writeUint32(PSEUDO_BIN);
			writer.writeUint32(2);
			writer.writeUint64(firstOffset);
			writer.writeUint64(lastOffset);
			writer.writeUint64(mappedReads);
			writer.writeUint64(unmappedReads);
		}

		private void ensureLinearSize(int size){
			while(linear.size()<size){
				linear.add(UNSET_OFFSET);
			}
		}

		private final IntObjectMap<BinData> bins;
		private final LongList linear;
		private long mappedReads=0L;
		private long unmappedReads=0L;
		private long firstOffset=-1L;
		private long lastOffset=-1L;
	}

	private static final class BinData {

		void append(long start, long end){
			if(current==null){
				current=new Chunk(start, end);
				chunks.add(current);
			}else if(start<=current.end){
				current.end=Math.max(current.end, end);
			}else{
				current=new Chunk(start, end);
				chunks.add(current);
			}
		}

		final ArrayList<Chunk> chunks=new ArrayList<>();
		private Chunk current;
	}

	private static final class Chunk {
		Chunk(long beg, long end){
			this.beg=beg;
			this.end=end;
		}

		final long beg;
		long end;
	}

	private static final int FIXED_RECORD_FIELDS=32;
	private static final int LINEAR_INDEX_SHIFT=14; // 16kb windows
	private static final int BAM_FUNMAP=0x4;
	private static final int PSEUDO_BIN=37450;
	private static final long UNSET_OFFSET=-1L;
}