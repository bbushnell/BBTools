package stream.bam;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import shared.Shared;
import simd.Vector;
import stream.SamLine;
import structures.BinaryByteWrapperLE;
import structures.ByteBuilder;

/**
 * Decodes individual BAM alignment bodies into SAM text or SamLine objects.
 * Inputs exclude block_size and any file header; callers handle record framing
 * and decompression. These methods expect well-formed, supported field encodings
 * and lengths that fit the arrays and arithmetic used here, rather than providing
 * complete BAM validation. Reference indices use the supplied dictionary order.
 * Float tags use Float.toString through ByteBuilder.appendSlow to retain finite values
 * (STR377); their decimal spelling is not forced to six places.
 * The text path retains stored sequence orientation. SamLine paths can reverse-complement
 * sequence and reverse quality according to the alignment flags and FLIP_ON_LOAD, then
 * apply SamLine name trimming. Their selective-parsing behavior differs below.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamToSamConverter{

	/** Borrowed reference names in BAM dictionary order; callers must keep them stable. */
	private final String[] refNames;

	/** Four-bit sequence code to printable base symbol. */
	private static final byte[] SEQ_LOOKUP_B="=ACMGRSVTWYHKDBN".getBytes();
	/** Packed CIGAR operation code to printable operation symbol. */
	private static final byte[] CIGAR_OPS_B="MIDNSHP=X".getBytes();

	/**
	 * Retains the caller's reference dictionary without copying the array.
	 * @param refNames_ Non-null reference names indexed by the input BAM dictionary
	 */
	public BamToSamConverter(String[] refNames_){refNames=refNames_;}

	/**
	 * Allocates SAM text directly from one alignment body, without SamLine normalization.
	 * Does not consult the SamLine PARSE flags, flip sequence/quality, or trim names.
	 * Unsupported/out-of-range reference indices are printed as "*"; matching mate
	 * indices use "=". Scalar integer tags are rendered with SAM type i.
	 * @param bamRecord Non-null complete record body, excluding the four-byte block_size
	 * @return Newly allocated tab-delimited SAM text without a trailing newline
	 */
	//Historical call-site/optimization note below. Current source also has the ST reader's
	//toSamLine caller; absence of other in-tree text callers does not rule out external users.
	//TEST-ONLY [stream/bam/BamToSamConverter#001]: this direct BAM→SAM-text fast path has only ONE caller,
	//TestSeqReverse (grep-confirmed). The LIVE BAM→SAM-text route goes BamStreamer→toSamLine→SamLine→
	//SamLine.toBytes, NOT through here. Correctly does NOT reverse-complement (SAM text SEQ and BAM SEQ are
	//both forward-reference, so a direct passthrough is right — the RC is only for building a SamLine whose
	//internal seq is original-orientation). Kept as a potential future fast path; currently test-only.
	public byte[] convertAlignment(byte[] bamRecord){
		ByteBuffer bb=ByteBuffer.wrap(bamRecord).order(ByteOrder.LITTLE_ENDIAN);
		ByteBuilder sam=new ByteBuilder(16+bamRecord.length*2); //Estimate SAM is ~2x BAM size

		//Read fixed-length fields
		int refID=bb.getInt();
		int pos=bb.getInt();
		int l_read_name=bb.get()&0xFF;
		int mapq=bb.get()&0xFF;
		int bin=bb.getShort()&0xFFFF; //Ignore bin
		int n_cigar_op=bb.getShort()&0xFFFF;
		int flag=bb.getShort()&0xFFFF;
		long l_seq=bb.getInt()&0xFFFFFFFFL;
		int next_refID=bb.getInt();
		int next_pos=bb.getInt();
		int tlen=bb.getInt();

		//Read variable-length fields
		byte[] readNameBytes=new byte[l_read_name];
		bb.get(readNameBytes);

		//QNAME (exclude NUL terminator)
		sam.append(readNameBytes, 0, l_read_name-1).tab();

		//FLAG
		sam.append(flag).tab();

		//RNAME
		if(refID<0 || refID>=refNames.length){sam.append('*');}else{sam.append(refNames[refID]);}
		sam.tab();

		//POS (BAM is 0-based, SAM is 1-based)
		sam.append(pos+1).tab();

		//MAPQ
		sam.append(mapq).tab();

		//CIGAR
		if(n_cigar_op==0){sam.append('*');}else{
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				sam.append(opLen).append(CIGAR_OPS_B[op]);
			}
		}
		sam.tab();

		//RNEXT
		if(next_refID<0){sam.append('*');}else if(next_refID==refID){sam.append('=');}else if(next_refID<refNames.length){sam.append(refNames[next_refID]);}else{sam.append('*');}
		sam.tab();

		//PNEXT (BAM is 0-based, SAM is 1-based)
		sam.append(next_pos+1).tab();

		//TLEN
		sam.append(tlen).tab();

		//SEQ (4-bit encoded, 2 bases per byte)
		if(l_seq==0){sam.append('*');}else{
			int pairs=(int)(l_seq/2);  // Number of complete pairs
			for(int i=0; i<pairs; i++){
				int packed=bb.get()&0xFF;
				sam.append(SEQ_LOOKUP_B[packed>>>4]);
				sam.append(SEQ_LOOKUP_B[packed&0xF]);
			}
			// Handle odd length - last nibble
			if((l_seq&1)==1){
				int packed=bb.get()&0xFF;
				sam.append(SEQ_LOOKUP_B[packed>>>4]);
			}
		}
		sam.tab();

		//QUAL (raw phred scores, add 33 for SAM)
		if(l_seq==0){sam.append('*');}else{
			//Peek first byte to check if QUAL is missing (all 0xFF)
			byte firstByte=bb.get();
			if(firstByte==(byte)0xFF){ // Comparing -1 to -1
				//Skip remaining bytes and output '*'
				bb.position(bb.position()+(int)l_seq-1);
				sam.append('*');
			}else{
				//Read remaining bytes into array
				byte[] qualBytes=new byte[(int)l_seq];
				qualBytes[0]=(byte)firstByte;
				bb.get(qualBytes, 1, (int)l_seq-1);

				//Convert to ASCII (phred+33) and append directly
				for(int i=0; i<qualBytes.length; i++){sam.append((byte)(qualBytes[i]+33));}
			}
		}

		//Decode auxiliary tags
		boolean firstTag=true;
		while(bb.hasRemaining()){
			if(firstTag){
				sam.tab();
				firstTag=false;
			}else{sam.tab();}

			//Tag name (2 bytes)
			sam.append(bb.get()).appendColon(bb.get());

			char type=(char)(bb.get()&0xFF);
			switch(type){
				case 'A': //Printable character
					sam.appendColon('A').append((char)(bb.get()&0xFF));
					break;
				case 'c': //int8_t
					sam.appendColon('i').append((int)bb.get());
					break;
				case 'C': //uint8_t
					sam.appendColon('i').append(bb.get()&0xFF);
					break;
				case 's': //int16_t
					sam.appendColon('i').append((int)bb.getShort());
					break;
				case 'S': //uint16_t
					sam.appendColon('i').append(bb.getShort()&0xFFFF);
					break;
				case 'i': //int32_t
					sam.appendColon('i').append(bb.getInt());
					break;
				case 'I': //uint32_t
					long uintVal=bb.getInt()&0xFFFFFFFFL;
					sam.appendColon('i').append(uintVal);
					break;
				case 'f': //float
					sam.appendColon('f').appendSlow(bb.getFloat());
					break;
				case 'Z': //Null-terminated string
					sam.appendColon('Z');
					byte b;
					while((b=bb.get())!=0){sam.append(b);}
					break;
				case 'H': //Hex string
					sam.appendColon('H');
					while((b=bb.get())!=0){sam.append(b);}
					break;
				case 'B': //Array
					char arrayType=(char)(bb.get()&0xFF);
					int count=bb.getInt();
					sam.appendColon('B').append(arrayType);
					for(int i=0; i<count; i++){
						sam.comma();
						switch(arrayType){
							case 'c':
								sam.append((int)bb.get()); break;
							case 'C':
								sam.append(bb.get()&0xFF); break;
							case 's':
								sam.append((int)bb.getShort()); break;
							case 'S':
								sam.append(bb.getShort()&0xFFFF); break;
							case 'i':
								sam.append(bb.getInt()); break;
							case 'I':
								sam.append(bb.getInt()&0xFFFFFFFFL); break;
							case 'f':
								sam.appendSlow(bb.getFloat()); break;
							default:
								throw new RuntimeException("Unknown array type: "+arrayType);
						}
					}
					break;
				default:
					throw new RuntimeException("Unknown tag type: "+type);
			}
		}
		return sam.toBytes();
	}

	/**
	 * Decodes all fields into a new SamLine without gating them on the PARSE flags.
	 * RNAME_AS_BYTES still controls reference-name storage. Mapped reverse-strand
	 * sequence is reverse-complemented and quality reversed when FLIP_ON_LOAD is enabled, followed
	 * by trimNames. This method expects an alignment body, not a header.
	 * @param bamRecord Non-null complete record body without block_size
	 * @return Newly allocated SamLine; this method has no header/null-result branch
	 */
	//Historical cleanup proposal below. No in-tree calls were found in the current source;
	//the public method is retained. The earlier "functionally correct" claim is not a new audit result.
	//TODO: Possible bug [stream/bam/BamToSamConverter#001] - DEAD METHOD (zero callers tree-wide; grep
	//confirmed). Superseded by toSamLine(byte[],ByteBuilder) (the live BamStreamer path). Functionally
	//correct, but ~191 lines of dead duplicated decode. Cleanup candidate (delete, or consolidate the
	//4-way method duplication). LOW/latent. See bug_reports/stream/bam/BamToSamConverter.md.
	public SamLine toSamLine_slim(byte[] bamRecord){
//		ByteBuffer bb=ByteBuffer.wrap(bamRecord).order(ByteOrder.LITTLE_ENDIAN);
		BinaryByteWrapperLE bb=new BinaryByteWrapperLE(bamRecord);
		SamLine sl=new SamLine();

		//Read fixed-length fields
		int refID=bb.getInt();
		int pos=bb.getInt();
		int l_read_name=bb.get()&0xFF;
		int mapq=bb.get()&0xFF;
		int bin=bb.getShort()&0xFFFF; //Ignore
		int n_cigar_op=bb.getShort()&0xFFFF;
		int flag=bb.getShort()&0xFFFF;
		long l_seq=bb.getInt()&0xFFFFFFFFL;
		int next_refID=bb.getInt();
		int next_pos=bb.getInt();
		int tlen=bb.getInt();

		//QNAME - SamLine.PARSE_0
		byte[] readNameBytes=new byte[l_read_name];
		bb.get(readNameBytes);
		sl.qname=new String(readNameBytes, 0, l_read_name-1, StandardCharsets.US_ASCII); //Exclude NUL

		//FLAG
		sl.flag=flag;

		//RNAME - SamLine.PARSE_2
		if(refID<0 || refID>=refNames.length){
			//do nothing
		}else if(SamLine.RNAME_AS_BYTES){sl.setRname(refNames[refID].getBytes());}else{sl.setRnameS(refNames[refID]);}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){sl.setCigar(null);}else{
			StringBuilder cigar=new StringBuilder(n_cigar_op*4);
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
			}
			sl.setCigar(cigar.toString());
		}

		//RNEXT
		if(next_refID<0){sl.setRnext(null);}else if(next_refID==refID){sl.setRnext(byteequals);}else if(next_refID<refNames.length){sl.setRnext(refNames[next_refID].getBytes());}else{
			sl.setRnext(null);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		if(l_seq==0){sl.setSeq(null);}else{
			byte[] seq=new byte[(int)l_seq];
			int numBytes=(int)((l_seq+1)/2);
			int seqIdx=0;
			for(int i=0; i<numBytes; i++){
				int packed=bb.get()&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];}
			}
			sl.setSeq(seq);
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){sl.setQual(null);}else{
			byte firstByte=bb.get();
			if(firstByte==(byte)0xFF){
				bb.position(bb.position()+(int)l_seq-1);
				sl.setQual(null);
			}else{
				byte[] qual=new byte[(int)l_seq];
				qual[0]=firstByte;
				bb.get(qual, 1, (int)l_seq-1);
				//Don't add 33 - keep as phred scores
				sl.setQual(qual);
			}
		}

		//Auxiliary tags - SamLine.PARSE_OPTIONAL
		if(bb.hasRemaining()){
			sl.optional=new ArrayList<String>();
			ByteBuilder tag=new ByteBuilder(64);

			while(bb.hasRemaining()){
				tag.clear();

				//Tag name (2 bytes)
				tag.append(bb.get()).appendColon(bb.get());

				char type=(char)(bb.get()&0xFF);
				switch(type){
					case 'A':
						tag.appendColon('A').append((char)(bb.get()&0xFF));
						break;
					case 'c':
						tag.appendColon('i').append((int)bb.get());
						break;
					case 'C':
						tag.appendColon('i').append(bb.get()&0xFF);
						break;
					case 's':
						tag.appendColon('i').append((int)bb.getShort());
						break;
					case 'S':
						tag.appendColon('i').append(bb.getShort()&0xFFFF);
						break;
					case 'i':
						tag.appendColon('i').append(bb.getInt());
						break;
					case 'I':
						tag.appendColon('i').append(bb.getInt()&0xFFFFFFFFL);
						break;
					case 'f':
						tag.appendColon('f').appendSlow(bb.getFloat());
						break;
					case 'Z':
						tag.appendColon('Z');
						byte b;
						while((b=bb.get())!=0){tag.append(b);}
						break;
					case 'H':
						tag.appendColon('H');
						while((b=bb.get())!=0){tag.append(b);}
						break;
					case 'B':
						char arrayType=(char)(bb.get()&0xFF);
						int count=bb.getInt();
						tag.appendColon('B').append(arrayType);
						for(int i=0; i<count; i++){
							tag.comma();
							switch(arrayType){
								case 'c': tag.append((int)bb.get()); break;
								case 'C': tag.append(bb.get()&0xFF); break;
								case 's': tag.append((int)bb.getShort()); break;
								case 'S': tag.append(bb.getShort()&0xFFFF); break;
								case 'i': tag.append(bb.getInt()); break;
								case 'I': tag.append(bb.getInt()&0xFFFFFFFFL); break;
								case 'f': tag.appendSlow(bb.getFloat()); break;
								default:
									throw new RuntimeException("Unknown array type: "+arrayType);
							}
						}
						break;
					default:
						throw new RuntimeException("Unknown tag type: "+type);
				}

				String tagString=tag.toString();
				sl.optional.add(tagString);

				if(tagString.startsWith("MD:")){sl.mdTag=tagString.substring(5).getBytes();}
			}
		}

		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=null){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=null){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}

	/**
	 * Decodes an alignment using the older copied-buffer implementation.
	 * PARSE_0, PARSE_2, PARSE_5, PARSE_10, and PARSE_OPTIONAL gate their retained
	 * fields; PARSE_6 gates named RNEXT values but not the matching-reference "=".
	 * Position, mate position, template length, flags, mapping quality, and sequence
	 * are still populated. Zero CIGAR operations normalize to null through setCigar.
	 * Mapped reverse-strand sequence/quality follow FLIP_ON_LOAD; trimNames runs last.
	 * @param bamRecord Non-null complete record body without block_size
	 * @return Newly allocated SamLine; this method has no header/null-result branch
	 */
	//Historical cleanup note below. setCigar("*") canonicalizes to null, so the claimed
	//zero-CIGAR difference does not survive that setter. No in-tree calls were found, but
	//this public method is retained; external callers were not assessed.
	//DEAD METHOD [stream/bam/BamToSamConverter#001] - zero callers tree-wide (the "Old" version). Note the
	//latent inconsistency vs the live path: this sets cigar="*" for n_cigar_op==0 (L466), while toSamLine_slim
	//and the live toSamLine set null - harmless only because it's dead. Cleanup candidate (delete). LOW.
	public SamLine toSamLineOld(byte[] bamRecord){
//		ByteBuffer bb=ByteBuffer.wrap(bamRecord).order(ByteOrder.LITTLE_ENDIAN);
		BinaryByteWrapperLE bb=new BinaryByteWrapperLE(bamRecord);
		SamLine sl=new SamLine();

		//Read fixed-length fields
		int refID=bb.getInt();
		int pos=bb.getInt();
		int l_read_name=bb.get()&0xFF;
		int mapq=bb.get()&0xFF;
		int bin=bb.getShort()&0xFFFF; //Ignore
		int n_cigar_op=bb.getShort()&0xFFFF;
		int flag=bb.getShort()&0xFFFF;
		long l_seq=bb.getInt()&0xFFFFFFFFL;
		int next_refID=bb.getInt();
		int next_pos=bb.getInt();
		int tlen=bb.getInt();

		//QNAME - SamLine.PARSE_0
		byte[] readNameBytes=new byte[l_read_name];
		bb.get(readNameBytes);
		if(SamLine.PARSE_0){sl.qname=new String(readNameBytes, 0, l_read_name-1, StandardCharsets.US_ASCII);} //Exclude NUL

		//FLAG
		sl.flag=flag;

		//RNAME - SamLine.PARSE_2
		if(refID<0 || refID>=refNames.length || !SamLine.PARSE_2){
			//do nothing
		}else if(SamLine.RNAME_AS_BYTES){sl.setRname(refNames[refID].getBytes());}else{sl.setRnameS(refNames[refID]);}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){sl.setCigar("*");}else{
			StringBuilder cigar=new StringBuilder(n_cigar_op*4);
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
			}
			if(SamLine.PARSE_5){sl.setCigar(cigar.toString());}
		}

		//TODO: Use bb.skip() for skipped fields, ByteBuilder instead of StringBuilder
		//TODO: Use direct array access, not array copies to byte[]

		//RNEXT
		if(next_refID<0){sl.setRnext(null);}else if(next_refID==refID){sl.setRnext(byteequals);}else if(next_refID<refNames.length && SamLine.PARSE_6){sl.setRnext(refNames[next_refID].getBytes());}else{
			sl.setRnext(null);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		if(l_seq==0){sl.setSeq(null);}else{
			byte[] seq=new byte[(int)l_seq];
			int numBytes=(int)((l_seq+1)/2);
			int seqIdx=0;
			for(int i=0; i<numBytes; i++){
				int packed=bb.get()&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];}
			}
			sl.setSeq(seq);
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){sl.setQual(null);}else{
			byte firstByte=bb.get();
			if(firstByte==(byte)0xFF){
				bb.position(bb.position()+(int)l_seq-1);
				sl.setQual(null);
			}else{
				byte[] qual=new byte[(int)l_seq];
				qual[0]=firstByte;
				bb.get(qual, 1, (int)l_seq-1);
				//Don't add 33 - keep as phred scores
				if(SamLine.PARSE_10){sl.setQual(qual);}
			}
		}

		//Auxiliary tags - SamLine.PARSE_OPTIONAL
		if(SamLine.PARSE_OPTIONAL && bb.hasRemaining()){
			sl.optional=new ArrayList<String>();
			ByteBuilder tag=new ByteBuilder(64);

			while(bb.hasRemaining()){
				tag.clear();

				//Tag name (2 bytes)
				tag.append(bb.get()).appendColon(bb.get());

				char type=(char)(bb.get()&0xFF);
				switch(type){
					case 'A':
						tag.appendColon('A').append((char)(bb.get()&0xFF));
						break;
					case 'c':
						tag.appendColon('i').append((int)bb.get());
						break;
					case 'C':
						tag.appendColon('i').append(bb.get()&0xFF);
						break;
					case 's':
						tag.appendColon('i').append((int)bb.getShort());
						break;
					case 'S':
						tag.appendColon('i').append(bb.getShort()&0xFFFF);
						break;
					case 'i':
						tag.appendColon('i').append(bb.getInt());
						break;
					case 'I':
						tag.appendColon('i').append(bb.getInt()&0xFFFFFFFFL);
						break;
					case 'f':
						tag.appendColon('f').appendSlow(bb.getFloat());
						break;
					case 'Z':
						tag.appendColon('Z');
						byte b;
						while((b=bb.get())!=0){tag.append(b);}
						break;
					case 'H':
						tag.appendColon('H');
						while((b=bb.get())!=0){tag.append(b);}
						break;
					case 'B':
						char arrayType=(char)(bb.get()&0xFF);
						int count=bb.getInt();
						tag.appendColon('B').append(arrayType);
						for(int i=0; i<count; i++){
							tag.comma();
							switch(arrayType){
								case 'c': tag.append((int)bb.get()); break;
								case 'C': tag.append(bb.get()&0xFF); break;
								case 's': tag.append((int)bb.getShort()); break;
								case 'S': tag.append(bb.getShort()&0xFFFF); break;
								case 'i': tag.append(bb.getInt()); break;
								case 'I': tag.append(bb.getInt()&0xFFFFFFFFL); break;
								case 'f': tag.appendSlow(bb.getFloat()); break;
								default:
									throw new RuntimeException("Unknown array type: "+arrayType);
							}
						}
						break;
					default:
						throw new RuntimeException("Unknown tag type: "+type);
				}

				String tagString=tag.toString();
				sl.optional.add(tagString);

				if(tagString.startsWith("MD:")){sl.mdTag=tagString.substring(5).getBytes();}
			}
		}

		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=null){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=null){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}

	/**
	 * Decodes one alignment while reusing optional CIGAR scratch space.
	 * PARSE_0 and PARSE_5 skip unwanted QNAME/CIGAR payloads. PARSE_2 gates RNAME;
	 * PARSE_6 gates named RNEXT values but not "=" for matching reference indices.
	 * PARSE_10 and PARSE_OPTIONAL gate retained qualities and auxiliary tags.
	 * Flags, coordinates, template length, mapping quality, and sequence are populated
	 * independently of the other PARSE flags. Mapped reverse-strand sequence/quality
	 * follow FLIP_ON_LOAD, and trimNames applies the current description-trimming setting.
	 * The record array is read without modification or retention by the returned object.
	 * The caller must keep the dictionary and relevant global options stable during use.
	 * @param bamRecord Non-null complete record body, excluding block_size
	 * @param cigar Optional caller-owned builder, empty when nonzero CIGAR operations
	 * are parsed; capacity is estimated from the operation count, with the limitation
	 * noted below, and the builder is cleared after successful CIGAR conversion
	 * @return Newly allocated SamLine; no file-header recognition or null result
	 */
	//Historical optimization note below. BamStreamer and BamReadInputStreamST both use this
	//entry. The write-side FLIP_ON_LOAD guard is now present in SamToBamConverter; the old
	//"MISSING" and "Traced clean" wording is retained as history, not a current correctness claim.
	//Scratch reuse limits allocation; returned strings/arrays and buffer growth can still allocate.
	//THE LIVE METHOD: BamStreamer:362 (the default BAM read path) calls this per record; the other three
	//parse methods are dead/test-only (see their headers). Optimized: skip-based (BinaryByteWrapperLE.skip)
	//+ DIRECT bamRecord[] reads for qname/seq/qual (no intermediate ByteBuffer copies), and a caller-supplied
	//reusable `cigar` ByteBuilder (zero-alloc; assert(cigar.isEmpty()) guards a dirty handoff). Traced clean:
	//4-bit SEQ map "=ACMGRSVTWYHKDBN", CIGAR "MIDNSHP=X", per-type tag decode, MD:Z→sl.mdTag, and the
	//FLIP_ON_LOAD-guarded reverse-strand RC below = the CORRECT read-side mirror of the write side (this is
	//exactly the `&& SamLine.FLIP_ON_LOAD` guard SamToBamConverter#001 is MISSING).
	public SamLine toSamLine(byte[] bamRecord, ByteBuilder cigar){
		BinaryByteWrapperLE bbw=new BinaryByteWrapperLE(bamRecord);
		SamLine sl=new SamLine();

		//Read fixed-length fields
		int refID=bbw.getInt();
		int pos=bbw.getInt();
		int l_read_name=bbw.get()&0xFF;
		int mapq=bbw.get()&0xFF;
		int bin=bbw.getShort()&0xFFFF; //Ignore
		int n_cigar_op=bbw.getShort()&0xFFFF;
		int flag=bbw.getShort()&0xFFFF;
		long l_seq=bbw.getInt()&0xFFFFFFFFL;
		int next_refID=bbw.getInt();
		int next_pos=bbw.getInt();
		int tlen=bbw.getInt();

		//QNAME - SamLine.PARSE_0
		if(SamLine.PARSE_0){
			int qnameStart=bbw.position();
			bbw.skip(l_read_name);
			sl.qname=new String(bamRecord, qnameStart, l_read_name-1, StandardCharsets.US_ASCII); //Exclude NUL
		}else{bbw.skip(l_read_name);}

		//FLAG
		sl.flag=flag;

		//RNAME - SamLine.PARSE_2
		if(refID<0 || refID>=refNames.length || !SamLine.PARSE_2){
			//do nothing
		}else if(SamLine.RNAME_AS_BYTES){sl.setRname(refNames[refID].getBytes());}else{sl.setRnameS(refNames[refID]);}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){
			sl.setCigar(null);
			//No skip needed - 0 cigar ops
		}else if(SamLine.PARSE_5){
			//TODO: Probable bug [stream/bam/BamToSamConverter#002] - the estimate reserves
			//8 characters per operation plus 8 total, but a 28-bit length can require 9 digits
			//plus its operation letter. appendOpUnsafe does not grow the backing array, so
			//this estimate is not sufficient for every representable CIGAR. Existing capacity
			//can mask the shortfall. Review sizing separately; source finding STR-187, not a runtime test.
			if(cigar==null){cigar=new ByteBuilder(n_cigar_op*8+8);}else{cigar.expand(n_cigar_op*8+8);}
			assert(cigar.isEmpty());
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bbw.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
//				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
				cigar.appendOpUnsafe(opLen, CIGAR_OPS_B[op]);
			}
			sl.setCigar(cigar.toString());
			cigar.clear();
		}else{bbw.skip(n_cigar_op*4);}

		//RNEXT
		if(next_refID<0){sl.setRnext(null);}else if(next_refID==refID){sl.setRnext(byteequals);}else if(next_refID<refNames.length && SamLine.PARSE_6){sl.setRnext(refNames[next_refID].getBytes());}else{
			sl.setRnext(null);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		int numSeqBytes=(int)((l_seq+1)/2);
		if(l_seq==0){sl.setSeq(null);}else{
			int seqStart=bbw.position();
			byte[] seq=new byte[(int)l_seq];
			int seqIdx=0;
			for(int i=0; i<numSeqBytes; i++){
				int packed=bamRecord[seqStart+i]&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];}
			}
			sl.setSeq(seq);
			bbw.skip(numSeqBytes);
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){sl.setQual(null);}else{
			int qualStart=bbw.position();
			byte firstByte=bamRecord[qualStart];
			//Missing-QUAL detection by the FIRST byte only is valid: BAM's "qual unavailable" convention is
			//ALL bytes==0xFF, and a real raw phred score is 0-93 (never 255), so firstByte==0xFF ⟺ qual
			//absent. (A malformed BAM with qual[0]=0xFF but later real bytes would misread → '*', but that's
			//format-impossible for valid input.)
			if(firstByte==(byte)0xFF){
				bbw.skip((int)l_seq);
				sl.setQual(null);
			}else if(SamLine.PARSE_10){
				byte[] qual=new byte[(int)l_seq];
				System.arraycopy(bamRecord, qualStart, qual, 0, (int)l_seq);
				bbw.skip((int)l_seq);
				sl.setQual(qual);
			}else{bbw.skip((int)l_seq);}
		}

		//Auxiliary tags - SamLine.PARSE_OPTIONAL
		if(SamLine.PARSE_OPTIONAL && bbw.hasRemaining()){
			sl.optional=new ArrayList<String>();
			ByteBuilder tag=new ByteBuilder(64);

			while(bbw.hasRemaining()){
				tag.clear();

				//Tag name (2 bytes)
				tag.append(bbw.get()).appendColon(bbw.get());

				char type=(char)(bbw.get()&0xFF);
				switch(type){
					case 'A':
						tag.appendColon('A').append((char)(bbw.get()&0xFF));
						break;
					case 'c':
						tag.appendColon('i').append((int)bbw.get());
						break;
					case 'C':
						tag.appendColon('i').append(bbw.get()&0xFF);
						break;
					case 's':
						tag.appendColon('i').append((int)bbw.getShort());
						break;
					case 'S':
						tag.appendColon('i').append(bbw.getShort()&0xFFFF);
						break;
					case 'i':
						tag.appendColon('i').append(bbw.getInt());
						break;
					case 'I':
						tag.appendColon('i').append(bbw.getInt()&0xFFFFFFFFL);
						break;
					case 'f':
						tag.appendColon('f').appendSlow(bbw.getFloat());
						break;
					case 'Z':
						tag.appendColon('Z');
						byte b;
						while((b=bbw.get())!=0){tag.append(b);}
						break;
					case 'H':
						tag.appendColon('H');
						while((b=bbw.get())!=0){tag.append(b);}
						break;
					case 'B':
						char arrayType=(char)(bbw.get()&0xFF);
						int count=bbw.getInt();
						tag.appendColon('B').append(arrayType);
						for(int i=0; i<count; i++){
							tag.comma();
							switch(arrayType){
								case 'c': tag.append((int)bbw.get()); break;
								case 'C': tag.append(bbw.get()&0xFF); break;
								case 's': tag.append((int)bbw.getShort()); break;
								case 'S': tag.append(bbw.getShort()&0xFFFF); break;
								case 'i': tag.append(bbw.getInt()); break;
								case 'I': tag.append(bbw.getInt()&0xFFFFFFFFL); break;
								case 'f': tag.appendSlow(bbw.getFloat()); break;
								default:
									throw new RuntimeException("Unknown array type: "+arrayType);
							}
						}
						break;
					default:
						throw new RuntimeException("Unknown tag type: "+type);
				}

				String tagString=tag.toString();
				sl.optional.add(tagString);

				if(tagString.startsWith("MD:")){sl.mdTag=tagString.substring(5).getBytes();}
			}
		}

		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=null){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=null){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}

//	private static final String stringstar=SamLine.stringstar;
//	private static final String stringequals=SamLine.stringequals;
//	private static final byte[] bytestar=SamLine.bytestar;
	/** Shared canonical RNEXT marker for a mate on the same reference. */
	private static final byte[] byteequals=SamLine.byteequals;

}
