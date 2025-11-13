package stream.bam;

import java.nio.charset.StandardCharsets;

import dna.AminoAcid;
import map.ObjectIntMap;
import shared.Parse;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Converts SAM text/SamLine to BAM binary format with zero allocation in hot path.
 * Handles encoding of CIGAR, SEQ (4-bit), QUAL, and auxiliary tags.
 *
 * @author Isla
 * @date November 1, 2025
 */
public class SamToBamConverter implements Cloneable {

//	private final Map<String, Integer> refMap;
	private final ObjectIntMap<String> refMap;

	//CIGAR operation lookup: direct array instead of HashMap
	private static final int[] CIGAR_OP_LOOKUP=new int[256];

	static {
		//Initialize with -1 for invalid ops
		for(int i=0; i<256; i++){
			CIGAR_OP_LOOKUP[i]=-1;
		}
		//CIGAR operations: MIDNSHP=X -> 0-8
		String ops="MIDNSHP=X";
		for(int i=0; i<ops.length(); i++){
			CIGAR_OP_LOOKUP[ops.charAt(i)]=i;
		}
	}
	
	public SamToBamConverter clone() {
		try{
			return (SamToBamConverter)super.clone();
		}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	public SamToBamConverter(String[] refNames){
		//Build reference name to ID map
		refMap=new ObjectIntMap<String>(512);
		for(int i=0; i<refNames.length; i++){
			refMap.put(refNames[i], i);
		}
	}
	
	/**
	 * Convert a SamLine to BAM binary format.
	 * @return Complete BAM record, including block_size prefix
	 */
	public byte[] convertAlignment(SamLine sl){
		ByteBuilder bb=new ByteBuilder(128);
		appendAlignment(sl, bb);
		bb.trimByAmount(4, 0);
		return bb.toBytes();
	}

	/**
	 * Convert a SamLine to BAM binary format, appending to ByteBuilder.
	 * Zero allocation hot path.
	 * @return ByteBuilder for chaining
	 */
	public ByteBuilder appendAlignment(final SamLine sl, final ByteBuilder bb){
		//Reserve 4 bytes for block_size (will patch at end)
		final int initialLength=bb.length();

		//Get reference IDs
		int refID=getRefID(sl.rnameS());
		int nextRefID=getNextRefID(sl.rnext(), refID);

		//Calculate bin and alignment length from CIGAR
		long binAndLen=calculateBinAndLength(sl);
		int bin=(int)(binAndLen>>32);
		int cigarOpCount=(int)binAndLen;

		//Get sequence length
		int seqLen=(sl.seq == null || sl.seq.length == 0) ? 0 : sl.seq.length;

		int estimatedSize=36+ 
			sl.qname.length()+1+
			cigarOpCount*4+
			(seqLen+1)/2+  // packed seq
			seqLen+              // qual
			(sl.optional==null ? 0 : 20*sl.optional.size());
		bb.expand(estimatedSize);

		bb.appendI32LE(0); //Placeholder
		//Write fixed-length fields (32 bytes)
		bb.appendI32LE(refID);
		bb.appendI32LE(sl.pos-1); //Convert 1-based to 0-based
		bb.appendU8(sl.qname.length()+1); //Include null terminator
		bb.appendU8(sl.mapq);
		bb.appendU16LE(bin);
		bb.appendU16LE(cigarOpCount);
		bb.appendU16LE(sl.flag);
		bb.appendU32LE(seqLen);
		bb.appendI32LE(nextRefID);
		bb.appendI32LE(sl.pnext-1); //Convert 1-based to 0-based
		bb.appendI32LE(sl.tlen);

		//QNAME with null terminator
		bb.append(sl.qname).append((byte)0);

		//CIGAR-encode directly to ByteBuilder
		appendCigar(bb, sl.cigar);
		
		//SEQ (4-bit encoded)-no temp arrays
		boolean mapped=(refID>=0);
		boolean reverseStrand=((sl.flag & 0x10)!=0);
		if(sl.seq==null || sl.seq.length==0) {
			//Do nothing
		}else if(mapped && reverseStrand && sl.seq!=null && sl.seq.length>0){
			appendSeqReverseComplement(bb, sl.seq);
		}else{appendSeq(bb, sl.seq);}
		
		//Qual is already 0-based
		assert(sl.qual==null || sl.qual.length==seqLen) : 
			"QUAL length mismatch: qual.length="+sl.qual.length+" != seqLen="+seqLen;
		if(sl.qual==null || sl.qual.length!=seqLen){
			appendSymbol(bb, (byte)0xFF, seqLen);
		}else if(mapped && reverseStrand){
			appendReversed(bb, sl.qual);
		}else{bb.append(sl.qual);}

		//Auxiliary tags-parse and append directly
		if(sl.optional!=null){
			for(String tag : sl.optional){appendTag(bb, tag);}
		}
		//Patch block_size at beginning (length excluding the 4-byte block_size itself)
		final int blockSize=bb.length()-initialLength-4;
		bb.setI32LE(blockSize, initialLength);
		return bb;
	}

	/**
	 * Get reference ID from reference name string.
	 */
	private int getRefID(String rname){
		if(rname == null || rname.equals("*")){
			return -1;
		}
		Integer id=refMap.get(rname);
		return (id != null) ? id : -1;
	}

	/**
	 * Get next reference ID from RNEXT.
	 */
	private int getNextRefID(byte[] rnext, int currentRefID){
		if(rnext == null || rnext.length == 0){
			return -1;
		}
		if(rnext.length == 1){
			if(rnext[0] == '*'){
				return -1;
			} else if(rnext[0] == '='){
				return currentRefID;
			}
		}
		String rnextStr=new String(rnext, StandardCharsets.US_ASCII);
		return getRefID(rnextStr);
	}

	/**
	 * Append CIGAR operations directly to ByteBuilder.
	 * Each operation: (length<<4) | op_code
	 */
	private void appendCigar(ByteBuilder bb, String cigar){
		if(cigar == null || cigar.equals("*")){return;}

		int len=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(c>='0' && c<='9'){
				len=len * 10+(c-'0');
			} else {
				int opCode=CIGAR_OP_LOOKUP[c];
				if(opCode<0){
					throw new RuntimeException("Unknown CIGAR operation: "+c);
				}
				bb.appendU32LE((len<<4) | opCode);
				len=0;
			}
		}
	}

	/**
	 * Calculate bin and CIGAR operation count in one pass.
	 * Returns packed long: (bin<<32) | cigarOpCount
	 */
	private long calculateBinAndLength(SamLine sl){
		if(sl.pos<=0 || sl.cigar == null || sl.cigar.equals("*")){
			return (4680L<<32); //Unmapped, 0 ops
		}

		String cigar=sl.cigar;
		int refLength=0;
		int cigarOpCount=0;
		int num=0;

		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(c>='0' && c<='9'){
				num=num * 10+(c-'0');
			} else {
				//Count operations
				cigarOpCount++;

				//Operations that consume reference: M, D, N, =, X
				if(c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X'){
					refLength += num;
				}
				num=0;
			}
		}

		int beg=sl.pos-1; //0-based
		int end=beg+refLength;
		int bin=reg2bin(beg, end);

		return ((long)bin<<32) | (cigarOpCount & 0xFFFFFFFFL);
	}

	/**
	 * Append 4-bit encoded sequence directly to ByteBuilder.
	 */
	private void appendSeq(ByteBuilder bb, byte[] seq){
		final byte[] array=bb.array;
		final int limit=(seq.length/2) * 2; //Even pairs
		int pos=bb.length;

		//Main loop-branchless
		for(int i=0; i<limit; i += 2){
			int hi=AminoAcid.baseToNumberExtended[seq[i]] & 0x0F;
			int lo=AminoAcid.baseToNumberExtended[seq[i+1]] & 0x0F;
			array[pos++]=(byte)((hi<<4) | lo);
		}

		//Handle odd length
		if((seq.length & 1) != 0){
			int hi=AminoAcid.baseToNumberExtended[seq[limit]] & 0x0F;
			array[pos++]=(byte)(hi<<4);
		}

		bb.length=pos;
	}

	/**
	 * Append 4-bit encoded reverse-complemented sequence directly to ByteBuilder.
	 */
	private void appendSeqReverseComplement(ByteBuilder bb, byte[] seq){
		final byte[] array=bb.array;
		int pos=bb.length;

		//Start from end, work backwards in pairs
		final int start=seq.length-1;
		final int limit=seq.length & 1; //Stop at 1 if odd, 0 if even

		//Main loop-branchless, iterate from end
		for(int i=start; i>=limit; i -= 2){
			int hi=AminoAcid.baseToComplementNumberExtended[seq[i]] & 0x0F;
			int lo=AminoAcid.baseToComplementNumberExtended[seq[i-1]] & 0x0F;
			array[pos++]=(byte)((hi<<4) | lo);
		}

		//Handle odd length (first base)
		if(limit != 0){
			int hi=AminoAcid.baseToComplementNumberExtended[seq[0]] & 0x0F;
			array[pos++]=(byte)(hi<<4);
		}

		bb.length=pos;
	}

	/**
	 * Append reversed quality scores directly to ByteBuilder.
	 */
	private int appendReversed(ByteBuilder bb, byte[] qual){
		final byte[] array=bb.array;
		int pos=bb.length;
		for(int i=qual.length-1; i>=0; i--){array[pos++]=qual[i];}
		return bb.length=pos;
	}
	
	private static int appendSymbol(final ByteBuilder bb, final byte symbol, final int amount){
		final byte[] array=bb.array;
		int pos=bb.length;
		for(int i=0; i<amount; i++){array[pos++]=symbol;}
		return bb.length=pos;
	}

	/**
	 * Append auxiliary tag directly to ByteBuilder.
	 * Format: TAG:TYPE:VALUE
	 * Zero allocation parsing and encoding.
	 */
	private void appendTag(ByteBuilder bb, String tagStr){
		if(tagStr.length()<5){
			throw new RuntimeException("Invalid tag format: "+tagStr);
		}

		//Write tag (2 bytes)
		bb.appendU8(tagStr.charAt(0));
		bb.appendU8(tagStr.charAt(1));

		char type=tagStr.charAt(3);

		//Write type and value based on type
		switch (type){
			case 'A': //Printable character
				bb.appendU8('A');
				bb.appendU8(tagStr.charAt(5));
				break;

			case 'i': { //Integer-choose smallest representation
				long intVal=Parse.parseLong(tagStr, 5, tagStr.length());
				if(intVal>=Byte.MIN_VALUE && intVal<=Byte.MAX_VALUE){
					bb.appendU8('c');
					bb.appendU8((int)intVal);
				} else if(intVal>=0 && intVal<=255){
					bb.appendU8('C');
					bb.appendU8((int)intVal);
				} else if(intVal>=Short.MIN_VALUE && intVal<=Short.MAX_VALUE){
					bb.appendU8('s');
					bb.appendU16LE((int)intVal);
				} else if(intVal>=0 && intVal<=65535){
					bb.appendU8('S');
					bb.appendU16LE((int)intVal);
				} else if(intVal>=Integer.MIN_VALUE && intVal<=Integer.MAX_VALUE){
					bb.appendU8('i');
					bb.appendI32LE((int)intVal);
				} else {
					bb.appendU8('I');
					bb.appendU32LE(intVal);
				}
				break;
			}

			case 'f': //Float
				bb.appendU8('f');
				float floatVal=Parse.parseFloat(tagStr, 5);
				bb.appendFloatLE(floatVal);
				break;

			case 'Z': //String
				bb.appendU8('Z');
				appendSubstring(bb, tagStr, 5, tagStr.length());
				bb.appendU8(0); //Null terminator
				break;

			case 'H': //Hex string
				bb.appendU8('H');
				appendSubstring(bb, tagStr, 5, tagStr.length());
				bb.appendU8(0); //Null terminator
				break;

			case 'B': //Array
				bb.appendU8('B');
				appendArrayTag(bb, tagStr, 5);
				break;

			default:
				throw new RuntimeException("Unknown tag type: "+type);
		}
	}

	/**
	 * Append array tag values directly to ByteBuilder.
	 * Format: type,val1,val2,...
	 */
	private void appendArrayTag(ByteBuilder bb, String value, int start){
		//Find array type (first char after TAG:TYPE:B:)
		char arrayType=value.charAt(start);
		bb.appendU8(arrayType);

		//Count elements and write count placeholder
		int countPos=bb.length();
		bb.appendI32LE(0); //Placeholder

		//Parse and write values
		int count=0;
		int i=start+2; //Skip type and comma

		while (i<value.length()){
			int commaPos=value.indexOf(',', i);
			if(commaPos<0) commaPos=value.length();

			switch (arrayType){
				case 'c':
				case 'C':
					bb.appendU8(Parse.parseInt(value, i, commaPos));
					break;
				case 's':
				case 'S':
					bb.appendU16LE(Parse.parseInt(value, i, commaPos));
					break;
				case 'i':
					bb.appendI32LE(Parse.parseInt(value, i, commaPos));
					break;
				case 'I':
					bb.appendU32LE(Parse.parseLong(value, i, commaPos));
					break;
				case 'f':
					bb.appendFloatLE(Parse.parseFloat(value, i, commaPos));
					break;
				default:
					throw new RuntimeException("Unknown array type: "+arrayType);
			}

			count++;
			i=commaPos+1;
		}

		//Patch count
		bb.setI32LE(count, countPos);
	}

	/**
	 * Append substring directly to ByteBuilder without allocation.
	 */
	private void appendSubstring(ByteBuilder bb, String s, int start, int end){
		for(int i=start; i<end; i++){
			bb.append((byte)s.charAt(i));
		}
	}

	/**
	 * Calculate BAM bin for a region [beg, end).
	 * Implementation from SAMv1.pdf page 20.
	 */
	private int reg2bin(int beg, int end){
		--end;
		if(beg>>14 == end>>14) return ((1<<15)-1)/7+(beg>>14);
		if(beg>>17 == end>>17) return ((1<<12)-1)/7+(beg>>17);
		if(beg>>20 == end>>20) return ((1<<9)-1)/7+(beg>>20);
		if(beg>>23 == end>>23) return ((1<<6)-1)/7+(beg>>23);
		if(beg>>26 == end>>26) return ((1<<3)-1)/7+(beg>>26);
		return 0;
	}
}