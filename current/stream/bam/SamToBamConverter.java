package stream.bam;

import java.nio.charset.StandardCharsets;

import dna.AminoAcid;
import map.ObjectIntMap;
import parse.Parse;
import shared.Tools;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * Encodes SamLine fields as BAM alignment records using a fixed reference dictionary.
 * appendAlignment retains the four-byte size prefix; convertAlignment returns only
 * the record body. Neither method emits a file header or performs complete validation.
 * Callers supply consistent fields representable by this encoder, including ASCII
 * names/tags and supported CIGAR/tag values. Existing CIGAR limitations are noted below.
 * Sequence and quality helpers write directly into reserved builder space. Builder
 * growth, named RNEXT lookup, and the byte-array convenience method can allocate.
 *
 * @author Isla
 * @date November 1, 2025
 */
public class SamToBamConverter implements Cloneable{

//	private final Map<String, Integer> refMap;
	/** Exact names and accepted aliases mapped to dictionary indices; populated at construction. */
	private final ObjectIntMap<String> refMap;

	/** CIGAR character to operation code, or -1; direct array instead of HashMap. */
	private static final int[] CIGAR_OP_LOOKUP=new int[256];

	static{
		//Initialize with -1 for invalid ops
		for(int i=0; i<256; i++){CIGAR_OP_LOOKUP[i]=-1;}
		//CIGAR operations: MIDNSHP=X -> 0-8
		String ops="MIDNSHP=X";
		for(int i=0; i<ops.length(); i++){CIGAR_OP_LOOKUP[ops.charAt(i)]=i;}
	}

	/**
	 * Makes a shallow copy sharing the reference map, which conversion only reads.
	 * @return New converter instance with the same dictionary and aliases
	 */
	public SamToBamConverter clone(){
		try{return (SamToBamConverter)super.clone();}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	/**
	 * Builds exact-name lookup from the caller's emitted dictionary order.
	 * @param refNames Non-null reference names in BAM dictionary order
	 * @throws IllegalArgumentException If one name identifies multiple indices
	 */
	public SamToBamConverter(String[] refNames){this(refNames, null);}

	/**
	 * Builds name lookup and adds aliases whose whitespace-trimmed names resolve.
	 * Unresolved aliases are ignored. Input arrays are not retained; indices must
	 * match the dictionary actually emitted by the caller.
	 * @param refNames Non-null reference names in BAM dictionary order
	 * @param refAliases Optional full names whose short names already resolve in the map
	 * @throws IllegalArgumentException If a name or alias identifies multiple indices
	 */
	public SamToBamConverter(String[] refNames, String[] refAliases){
		//Build reference name to ID map.  The emitted BAM dictionary is authoritative.
		refMap=new ObjectIntMap<String>(Math.max(512, refNames.length*2), String.class);
		for(int i=0; i<refNames.length; i++){putRefName(refNames[i], i, "BAM dictionary");}
		//Aliases are populated once so alignment conversion remains an exact hot-path lookup.
		if(refAliases!=null){
			for(String alias : refAliases){
				final String dictionaryName=Tools.trimToWhitespace(alias);
				final int id=refMap.get(dictionaryName);
				if(id>=0){putRefName(alias, id, "reference alias");}
			}
		}
	}

	/**
	 * Adds one name, accepting a repeated association with the same index.
	 * @param name Non-null lookup key
	 * @param id Nonnegative dictionary index
	 * @param source Description used in an ambiguity diagnostic
	 * @throws IllegalArgumentException If the name already identifies another index
	 */
	private void putRefName(String name, int id, String source){
		final int old=refMap.get(name);
		if(old>=0 && old!=id){
			throw new IllegalArgumentException("Ambiguous "+source+" name '"+name+
				"' identifies BAM reference IDs "+old+" and "+id+".");
		}
		if(old<0){refMap.put(name, id);}
	}

	/**
	 * Allocates a BAM record body, removing the size prefix written by appendAlignment.
	 * A caller writing a complete record must emit the returned array length first.
	 * @param sl Non-null alignment satisfying appendAlignment's input contract
	 * @return Newly allocated record body, excluding the four-byte block_size prefix
	 */
	public byte[] convertAlignment(SamLine sl){
		ByteBuilder bb=new ByteBuilder(128);
		appendAlignment(sl, bb);
		bb.trimByAmount(4, 0);
		return bb.toBytes();
	}

	/**
	 * Appends a complete alignment record and patches its four-byte size prefix.
	 * Existing builder contents are retained. The SamLine and its arrays are read,
	 * not modified; the caller exclusively owns the builder during this call.
	 * Reference IDs come from this converter's map. Sequence/quality reversal requires
	 * a resolved reference, reverse-strand flag, and SamLine.FLIP_ON_LOAD.
	 * Null qualities emit missing-quality bytes; enabled assertions check lengths.
	 * This is an encoder, not a complete consistency or range validator.
	 * @param sl Non-null alignment with a non-null QNAME, resolvable named references,
	 * consistent CIGAR/position/sequence fields, and representable lengths and tags
	 * @param bb Non-null destination, expanded as needed
	 * @return The supplied builder with the size-prefixed record appended
	 * @throws IllegalArgumentException If a named reference is absent from the map
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
		int seqLen=(sl.seq==null || sl.seq.length==0) ? 0 : sl.seq.length;

		int estimatedSize=36+
			sl.qname.length()+1+
			cigarOpCount*4+
			(seqLen+1)/2+  // packed seq
			seqLen+              // qual
			(sl.optional==null ? 0 : 20*sl.optional.size());
		//Direct-array-write safety invariant: estimatedSize EXACTLY covers the fixed 36B + qname + cigar +
		//packed-seq + qual, so after bb.expand the subsequent appendSeq/appendSeqReverseComplement/
		//appendReversed/appendSymbol can write straight into bb.array (caching it) with NO bounds risk - those
		//4 helpers do NOT ensureExtra themselves. Tags are the ONLY part whose true size can exceed the 20B/tag
		//estimate, and appendTag re-ensures per tag (ensureExtra L308); since tags are written AFTER seq+qual,
		//a tag-driven realloc can never invalidate the seq/qual direct writes. (Holds iff ByteBuilder.expand
		//ensures `estimatedSize` ADDITIONAL bytes - the standard BBTools contract.)
		bb.expand(estimatedSize);

		bb.appendI32LE(0); //Placeholder
		//Write fixed-length fields (32 bytes)
		bb.appendI32LE(refID);
		bb.appendI32LE(sl.pos-1); //Convert 1-based to 0-based
		bb.appendU8(sl.qname.length()+1); //Include null terminator
		//l_read_name is u8 -> qname.length()+1 must be <=255 i.e. qname<=254 chars. The SAM spec caps QNAME
		//at 254, so valid input never overflows this u8 (format-bounded, not a finding).
		bb.appendU8(sl.mapq);
		bb.appendU16LE(bin);
		//TODO: Possible bug [stream/bam/SamToBamConverter#003] - n_cigar_op is u16 (max 65535) with no
		//CG-tag fallback. A read with >65535 CIGAR ops (pathological ultra-long ONT alignments) overflows
		//this field silently (the full cigar bytes are still written by appendCigar -> op-count field
		//mismatches payload -> corrupt record). This is the known BAM-format limitation htslib works around
		//via the "CG" aux tag. LOW (extreme input), but silent; crash-loud would prefer an explicit
		//assert(cigarOpCount<=65535) over a wrapped count.
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
		//Historical report follows; its unconditional-RC description predates the FIXED guards below.
		//TODO: Possible bug [stream/bam/SamToBamConverter#001] - the reverse-strand RC here is UNCONDITIONAL,
		//but the symmetric read-side BamToSamConverter guards the SAME un-flip with `&& SamLine.FLIP_ON_LOAD`
		//(BamToSamConverter:409/610/817; flag = `flipsam`, Parser:1452, default true). This is CORRECT under
		//the default (flipsam=t): SamLine's text ctor flips reverse-strand seq to ORIGINAL orientation on load
		//(SamLine:474), so RC here restores forward-reference for the BAM. But under `flipsam=f` the load does
		//NOT flip (seq stays forward-reference), so RC here DOUBLE-flips -> silently corrupt BAM SEQ/QUAL for
		//every mapped reverse-strand read. SamLine.toBytes:2360 has the identical unconditional-RC (so SAM->SAM
		//with flipsam=f corrupts too). flag-for-Brian: is flipsam=f supported with mapped SAM/BAM output? If yes
		//-> add `&& SamLine.FLIP_ON_LOAD` here + at toBytes:2360 (match the read side); if it's FastqScan-only
		//(no mapped output) -> document the constraint + crash loud. See bug_reports/stream/bam/SamToBamConverter.md.
		//Note: `mapped` is derived from refID>=0 (rname resolution), NOT sl.mapped() (the 0x4 flag) that the
		//load-flip used - on inconsistent rname/flag input these diverge (#002 family).
		boolean mapped=(refID>=0);
		boolean reverseStrand=((sl.flag&0x10)!=0);
		if(sl.seq==null || sl.seq.length==0){
			//Do nothing
		}else if(mapped && reverseStrand && SamLine.FLIP_ON_LOAD && sl.seq!=null && sl.seq.length>0){
			//#001 FIXED 2026-06-20 (greenlit): +`&& SamLine.FLIP_ON_LOAD` mirrors the read-side
			//BamToSamConverter. With flipsam=f the SamLine seq is already forward-reference (load didn't flip),
			//so we take the else branch (appendSeq, as-is) instead of double-flipping. VALIDATED 2026-06-20 (Furina):
			//no-flip BAM decode confirms BOTH flipsam settings store reverse-strand SEQ in forward-reference orientation
			//(1002 phix rev reads, 0 mismatch, byte-identical BAMs); pre-fix flipsam=f stored RC(seq) — the BAM does not.
			appendSeqReverseComplement(bb, sl.seq);
		}else{appendSeq(bb, sl.seq);}

		//Qual is already 0-based
		assert(sl.qual==null || sl.qual.length==seqLen) :
			"QUAL length mismatch: qual.length="+sl.qual.length+" != seqLen="+seqLen;
		if(sl.qual==null || sl.qual.length!=seqLen){appendSymbol(bb, (byte)0xFF, seqLen);}else if(mapped && reverseStrand && SamLine.FLIP_ON_LOAD){ //#001 FIXED: +FLIP_ON_LOAD guard (mirror read side)
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
	 * Resolves an exact reference name without changing the dictionary.
	 * @param rname Reference name, null, or "*"
	 * @return Dictionary index, or -1 for null/"*"
	 * @throws IllegalArgumentException If a named reference is absent
	 */
	private int getRefID(String rname){
		if(rname==null || rname.equals("*")){return -1;}
		final int id=refMap.get(rname);
		if(id<0){
			throw new IllegalArgumentException("Reference name '"+rname+
				"' is absent from the BAM header dictionary.");
		}
		return id;
	}

	/**
	 * Resolves RNEXT, allocating an ASCII String for an ordinary reference name.
	 * @param rnext Null, empty, "*", "=", or an ASCII reference name
	 * @param currentRefID Index to reuse for "="
	 * @return Resolved index, or -1 for null/empty/"*"
	 * @throws IllegalArgumentException If a named reference is absent
	 */
	private int getNextRefID(byte[] rnext, int currentRefID){
		if(rnext==null || rnext.length==0){return -1;}
		if(rnext.length==1){
			if(rnext[0]=='*'){return -1;}else if(rnext[0]=='='){
				return currentRefID;
			}
		}
		String rnextStr=new String(rnext, StandardCharsets.US_ASCII);
		return getRefID(rnextStr);
	}

	/**
	 * Appends packed CIGAR words as {@code (length<<4)|op_code}.
	 * Null/"*" emits nothing. The caller reserves capacity and supplies valid,
	 * representable operations; this scan does not validate the complete grammar.
	 * @param bb Destination with space for four bytes per operation
	 * @param cigar CIGAR text, null, or "*"
	 */
	private void appendCigar(ByteBuilder bb, String cigar){
		if(cigar==null || cigar.equals("*")){return;}

		int len=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(c>='0' && c<='9'){len=len*10+(c-'0');}else{
				int opCode=CIGAR_OP_LOOKUP[c];
				if(opCode<0){throw new RuntimeException("Unknown CIGAR operation: "+c);}
				bb.appendU32LE((len<<4)|opCode);
				len=0;
			}
		}
	}

	/**
	 * Calculates the bin and operation count without changing the alignment.
	 * Nonpositive position or null/"*" CIGAR returns bin 4680 and zero operations.
	 * Otherwise M, D, N, =, and X contribute to the reference span.
	 * @param sl Alignment with representable position and CIGAR values
	 * @return Packed {@code (bin<<32)|cigarOpCount}, with the count in the low 32 bits
	 */
	private long calculateBinAndLength(SamLine sl){
		//The historical report's specification/severity judgments below are not input validation.
		//TODO: Possible bug [stream/bam/SamToBamConverter#002] - this forces cigarOpCount=0 when pos<=0,
		//but appendCigar() (called separately at L112) encodes ops from sl.cigar REGARDLESS of pos. So a
		//read with pos<=0 AND a non-"*" cigar writes n_cigar_op=0 into the record while appendCigar emits
		//real cigar bytes -> the BAM record's op count mismatches its cigar payload -> corrupt record (the
		//reader parses SEQ starting inside the cigar). Only diverges on INCONSISTENT input (a mapped cigar
		//at pos<=0; valid SAM has pos>0 <=> cigar!="*"), so LOW/latent - but crash-loud would prefer
		//catching it. Fix: count ops from the cigar string itself (decoupled from pos), or assert the
		//pos>0 <=> cigar!="*" invariant loud. The pos>0 + cigar="*" case is consistent (both 0 ops).
		if(sl.pos<=0 || sl.cigar==null || sl.cigar.equals("*")){
			return (4680L<<32); //Unmapped, 0 ops
		}

		String cigar=sl.cigar;
		int refLength=0;
		int cigarOpCount=0;
		int num=0;

		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(c>='0' && c<='9'){num=num*10+(c-'0');}else{
				//Count operations
				cigarOpCount++;

				//Operations that consume reference: M, D, N, =, X
				if(c=='M' || c=='D' || c=='N' || c=='=' || c=='X'){refLength+=num;}
				num=0;
			}
		}

		int beg=sl.pos-1; //0-based
		int end=beg+refLength;
		int bin=reg2bin(beg, end);

		return ((long)bin<<32)|(cigarOpCount&0xFFFFFFFFL);
	}

	/**
	 * Appends two four-bit base codes per byte without a temporary sequence array.
	 * An odd final base occupies the high nibble; the low nibble is zero.
	 * @param bb Destination with capacity for (seq.length+1)/2 additional bytes
	 * @param seq Non-null bases suitable for the AminoAcid lookup table
	 */
	private void appendSeq(ByteBuilder bb, byte[] seq){
		final byte[] array=bb.array;
		final int limit=(seq.length/2)*2; //Even pairs
		int pos=bb.length;

		//Main loop-branchless
		for(int i=0; i<limit; i+=2){
			int hi=AminoAcid.baseToNumberExtended[seq[i]]&0x0F;
			int lo=AminoAcid.baseToNumberExtended[seq[i+1]]&0x0F;
			array[pos++]=(byte)((hi<<4)|lo);
		}

		//Handle odd length
		if((seq.length&1)!=0){
			int hi=AminoAcid.baseToNumberExtended[seq[limit]]&0x0F;
			array[pos++]=(byte)(hi<<4);
		}

		bb.length=pos;
	}

	/**
	 * Appends reverse-complemented four-bit base codes without modifying the input.
	 * An odd final encoded base occupies the high nibble; the low nibble is zero.
	 * @param bb Destination with capacity for (seq.length+1)/2 additional bytes
	 * @param seq Non-null bases suitable for the AminoAcid complement lookup table
	 */
	private void appendSeqReverseComplement(ByteBuilder bb, byte[] seq){
		final byte[] array=bb.array;
		int pos=bb.length;

		//Start from end, work backwards in pairs
		final int start=seq.length-1;
		final int limit=seq.length&1; //Stop at 1 if odd, 0 if even

		//Main loop-branchless, iterate from end
		for(int i=start; i>=limit; i-=2){
			int hi=AminoAcid.baseToComplementNumberExtended[seq[i]]&0x0F;
			int lo=AminoAcid.baseToComplementNumberExtended[seq[i-1]]&0x0F;
			array[pos++]=(byte)((hi<<4)|lo);
		}

		//Handle odd length (first base)
		if(limit!=0){
			int hi=AminoAcid.baseToComplementNumberExtended[seq[0]]&0x0F;
			array[pos++]=(byte)(hi<<4);
		}

		bb.length=pos;
	}

	/**
	 * Copies quality bytes in reverse order without modifying the input.
	 * @param bb Destination with capacity for qual.length additional bytes
	 * @param qual Non-null numeric quality bytes
	 * @return Updated total builder length
	 */
	private int appendReversed(ByteBuilder bb, byte[] qual){
		final byte[] array=bb.array;
		int pos=bb.length;
		for(int i=qual.length-1; i>=0; i--){array[pos++]=qual[i];}
		return bb.length=pos;
	}

	/**
	 * Appends a repeated byte directly into reserved builder space.
	 * @param bb Destination with capacity for amount additional bytes
	 * @param symbol Byte to repeat
	 * @param amount Nonnegative repetition count
	 * @return Updated total builder length
	 */
	private static int appendSymbol(final ByteBuilder bb, final byte symbol, final int amount){
		final byte[] array=bb.array;
		int pos=bb.length;
		for(int i=0; i<amount; i++){array[pos++]=symbol;}
		return bb.length=pos;
	}

	/**
	 * Appends a SAM auxiliary field in TAG:TYPE:VALUE form, reserving its output space.
	 * Supports A, i, f, Z, H, and B; scalar integers choose a compact storage width.
	 * Numeric parsing delegates to Parse. String/hex payload characters are narrowed
	 * to bytes and terminated with zero; hex syntax is not checked here.
	 * @param bb Destination builder
	 * @param tagStr Non-null, well-formed ASCII tag with values representable by its type
	 */
	private void appendTag(ByteBuilder bb, String tagStr){
		if(tagStr.length()<5){throw new RuntimeException("Invalid tag format: "+tagStr);}
		bb.ensureExtra(20+4*tagStr.length());

		//Write tag (2 bytes)
		bb.appendU8(tagStr.charAt(0));
		bb.appendU8(tagStr.charAt(1));

		char type=tagStr.charAt(3);

		//Write type and value based on type
		switch(type){
			case 'A': //Printable character
				bb.appendU8('A');
				bb.appendU8(tagStr.charAt(5));
				break;

			case 'i':{ //Integer-choose smallest representation
				long intVal=Parse.parseLong(tagStr, 5, tagStr.length());
				if(intVal>=Byte.MIN_VALUE && intVal<=Byte.MAX_VALUE){
					bb.appendU8('c');
					bb.appendU8((int)intVal);
				}else if(intVal>=0 && intVal<=255){
					bb.appendU8('C');
					bb.appendU8((int)intVal);
				}else if(intVal>=Short.MIN_VALUE && intVal<=Short.MAX_VALUE){
					bb.appendU8('s');
					bb.appendU16LE((int)intVal);
				}else if(intVal>=0 && intVal<=65535){
					bb.appendU8('S');
					bb.appendU16LE((int)intVal);
				}else if(intVal>=Integer.MIN_VALUE && intVal<=Integer.MAX_VALUE){
					bb.appendU8('i');
					bb.appendI32LE((int)intVal);
				}else{
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
	 * Appends a B-tag subtype, element count, and comma-separated numeric values.
	 * Supports c, C, s, S, i, I, and f; patches the count after writing elements.
	 * The caller reserves capacity and supplies correctly delimited, representable values.
	 * @param bb Destination with enough space for the subtype, count, and elements
	 * @param value Full tag text containing a subtype followed by comma-separated values
	 * @param start Index of the subtype character
	 */
	private void appendArrayTag(ByteBuilder bb, String value, int start){
		//Find array subtype (first char after TAG:B:)
		char arrayType=value.charAt(start);
		bb.appendU8(arrayType);

		//Count elements and write count placeholder
		int countPos=bb.length();
		bb.appendI32LE(0); //Placeholder

		//Parse and write values
		int count=0;
		int i=start+2; //Skip type and comma

		while(i<value.length()){
			int commaPos=value.indexOf(',', i);
			if(commaPos<0){commaPos=value.length();}

			switch(arrayType){
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
	 * Appends characters narrowed to bytes, without creating a substring or terminator.
	 * @param bb Destination builder
	 * @param s Non-null text, expected to be ASCII
	 * @param start Inclusive first character index
	 * @param end Exclusive last character index
	 */
	private void appendSubstring(ByteBuilder bb, String s, int start, int end){
		for(int i=start; i<end; i++){bb.append((byte)s.charAt(i));}
	}

	/**
	 * Calculates a bin from a zero-based half-open reference interval.
	 * Original implementation attribution: SAMv1.pdf page 20.
	 * @param beg Inclusive start coordinate within the binning scheme's range
	 * @param end Exclusive end coordinate; decremented before comparing shifted endpoints
	 * @return First matching bin, or zero when no finer bin matches
	 */
	private int reg2bin(int beg, int end){
		--end;
		if(beg>>14==end>>14){return ((1<<15)-1)/7+(beg>>14);}
		if(beg>>17==end>>17){return ((1<<12)-1)/7+(beg>>17);}
		if(beg>>20==end>>20){return ((1<<9)-1)/7+(beg>>20);}
		if(beg>>23==end>>23){return ((1<<6)-1)/7+(beg>>23);}
		if(beg>>26==end>>26){return ((1<<3)-1)/7+(beg>>26);}
		return 0;
	}
}
