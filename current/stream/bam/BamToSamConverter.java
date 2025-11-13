package stream.bam;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import shared.Shared;
import shared.Vector;
import stream.SamLine;
import structures.BinaryByteWrapperLE;
import structures.ByteBuilder;

/**
 * Converts BAM binary alignment records to SAM text format.
 * Handles all BAM field encodings including 4-bit SEQ, packed CIGAR, and auxiliary tags.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamToSamConverter {

	private final String[] refNames;

	//Lookup tables for decoding
	private static final byte[] SEQ_LOOKUP_B="=ACMGRSVTWYHKDBN".getBytes();
	private static final byte[] CIGAR_OPS_B="MIDNSHP=X".getBytes();

	public BamToSamConverter(String[] refNames_){
		refNames=refNames_;
	}

	/**
	 * Convert a BAM alignment record to SAM text format.
	 * @param bamRecord The raw BAM record bytes (not including block_size field)
	 * @return Tab-delimited SAM line as byte array
	 */
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
		if(refID<0 || refID>=refNames.length){sam.append('*');}
		else{sam.append(refNames[refID]);}
		sam.tab();

		//POS (BAM is 0-based, SAM is 1-based)
		sam.append(pos+1).tab();

		//MAPQ
		sam.append(mapq).tab();

		//CIGAR
		if(n_cigar_op==0){sam.append('*');}
		else{
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				sam.append(opLen).append(CIGAR_OPS_B[op]);
			}
		}
		sam.tab();

		//RNEXT
		if(next_refID<0){sam.append('*');}
		else if(next_refID==refID){sam.append('=');}
		else if(next_refID<refNames.length){sam.append(refNames[next_refID]);}
		else{sam.append('*');}
		sam.tab();

		//PNEXT (BAM is 0-based, SAM is 1-based)
		sam.append(next_pos+1).tab();

		//TLEN
		sam.append(tlen).tab();

		//SEQ (4-bit encoded, 2 bases per byte)
		if(l_seq==0){sam.append('*');}
		else{
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
		if(l_seq==0){sam.append('*');}
		else{
			//Peek first byte to check if QUAL is missing (all 0xFF)
			byte firstByte=bb.get();
			if(firstByte == (byte)0xFF){ // Comparing -1 to -1
				//Skip remaining bytes and output '*'
				bb.position(bb.position()+(int)l_seq-1);
				sam.append('*');
			}else{
				//Read remaining bytes into array
				byte[] qualBytes=new byte[(int)l_seq];
				qualBytes[0]=(byte)firstByte;
				bb.get(qualBytes, 1, (int)l_seq-1);

				//Convert to ASCII (phred+33) and append directly
				for(int i=0; i<qualBytes.length; i++){
					sam.append((byte)(qualBytes[i]+33));
				}
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
			switch (type){
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
					sam.appendColon('f').append(bb.getFloat(), 6);
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
						switch (arrayType){
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
								sam.append(bb.getFloat(), 6); break;
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
	 * Convert a BAM alignment record directly to a SamLine object.
	 * @param bamRecord The raw BAM record bytes (not including block_size field)
	 * @return Populated SamLine object, or null if header
	 */
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
		}else if(SamLine.RNAME_AS_BYTES){
			sl.setRname(refNames[refID].getBytes());
		}else {
			sl.setRname(refNames[refID]);
		}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){
			sl.cigar="*";
		}else{
			StringBuilder cigar=new StringBuilder(n_cigar_op*4);
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
			}
			sl.cigar=cigar.toString();
		}

		//RNEXT
		if(next_refID<0){
			sl.setRnext(bytestar);
		}else if(next_refID==refID){
			sl.setRnext(byteequals);
		}else if(next_refID<refNames.length){
			sl.setRnext(refNames[next_refID].getBytes());
		}else{
			sl.setRnext(bytestar);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		if(l_seq==0){
			sl.seq=bytestar;
		}else{
			byte[] seq=new byte[(int)l_seq];
			int numBytes=(int)((l_seq+1)/2);
			int seqIdx=0;
			for(int i=0; i<numBytes; i++){
				int packed=bb.get()&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){
					seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];
				}
			}
			sl.seq=seq;
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){
			sl.qual=bytestar;
		}else{
			byte firstByte=bb.get();
			if(firstByte==(byte)0xFF){
				bb.position(bb.position()+(int)l_seq-1);
				sl.qual=bytestar;
			}else{
				byte[] qual=new byte[(int)l_seq];
				qual[0]=firstByte;
				bb.get(qual, 1, (int)l_seq-1);
				//Don't add 33 - keep as phred scores
				sl.qual=qual;
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
						tag.appendColon('f').append(bb.getFloat(), 6);
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
								case 'f': tag.append(bb.getFloat(), 6); break;
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

				if(tagString.startsWith("MD:")){
					sl.mdTag=tagString.substring(5).getBytes();
				}
			}
		}
		
		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=bytestar){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=bytestar){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}
	
	/**
	 * Convert a BAM alignment record directly to a SamLine object.
	 * @param bamRecord The raw BAM record bytes (not including block_size field)
	 * @return Populated SamLine object, or null if header
	 */
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
		if(SamLine.PARSE_0) {sl.qname=new String(readNameBytes, 0, l_read_name-1, StandardCharsets.US_ASCII);} //Exclude NUL

		//FLAG
		sl.flag=flag;

		//RNAME - SamLine.PARSE_2
		if(refID<0 || refID>=refNames.length || !SamLine.PARSE_2){
			//do nothing
		}else if(SamLine.RNAME_AS_BYTES){
			sl.setRname(refNames[refID].getBytes());
		}else {
			sl.setRname(refNames[refID]);
		}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){
			sl.cigar="*";
		}else{
			StringBuilder cigar=new StringBuilder(n_cigar_op*4);
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bb.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
			}
			if(SamLine.PARSE_5) {sl.cigar=cigar.toString();}
		}
		
		//TODO: Use bb.skip() for skipped fields, ByteBuilder instead of StringBuilder
		//TODO: Use direct array access, not array copies to byte[]

		//RNEXT
		if(next_refID<0){
			sl.setRnext(bytestar);
		}else if(next_refID==refID){
			sl.setRnext(byteequals);
		}else if(next_refID<refNames.length && SamLine.PARSE_6){
			sl.setRnext(refNames[next_refID].getBytes());
		}else{
			sl.setRnext(bytestar);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		if(l_seq==0){
			sl.seq=bytestar;
		}else{
			byte[] seq=new byte[(int)l_seq];
			int numBytes=(int)((l_seq+1)/2);
			int seqIdx=0;
			for(int i=0; i<numBytes; i++){
				int packed=bb.get()&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){
					seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];
				}
			}
			sl.seq=seq;
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){
			sl.qual=bytestar;
		}else{
			byte firstByte=bb.get();
			if(firstByte==(byte)0xFF){
				bb.position(bb.position()+(int)l_seq-1);
				sl.qual=bytestar;
			}else{
				byte[] qual=new byte[(int)l_seq];
				qual[0]=firstByte;
				bb.get(qual, 1, (int)l_seq-1);
				//Don't add 33 - keep as phred scores
				if(SamLine.PARSE_10) {sl.qual=qual;}
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
						tag.appendColon('f').append(bb.getFloat(), 6);
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
								case 'f': tag.append(bb.getFloat(), 6); break;
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

				if(tagString.startsWith("MD:")){
					sl.mdTag=tagString.substring(5).getBytes();
				}
			}
		}
		
		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=bytestar){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=bytestar){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}
	
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
		}else{
			bbw.skip(l_read_name);
		}

		//FLAG
		sl.flag=flag;

		//RNAME - SamLine.PARSE_2
		if(refID<0 || refID>=refNames.length || !SamLine.PARSE_2){
			//do nothing
		}else if(SamLine.RNAME_AS_BYTES){
			sl.setRname(refNames[refID].getBytes());
		}else{
			sl.setRname(refNames[refID]);
		}

		//POS (BAM is 0-based, SAM is 1-based)
		sl.pos=pos+1;

		//MAPQ
		sl.mapq=mapq;

		//CIGAR - SamLine.PARSE_5
		if(n_cigar_op==0){
			sl.cigar="*";
			//No skip needed - 0 cigar ops
		}else if(SamLine.PARSE_5){
			if(cigar==null) {cigar=new ByteBuilder(n_cigar_op*8+8);}
			else {cigar.expand(n_cigar_op*8+8);}
			assert(cigar.isEmpty());
			for(int i=0; i<n_cigar_op; i++){
				int cigOp=bbw.getInt();
				int opLen=cigOp>>>4;
				int op=cigOp&0xF;
//				cigar.append(opLen).append((char)CIGAR_OPS_B[op]);
				cigar.appendOpUnsafe(opLen, CIGAR_OPS_B[op]);
			}
			sl.cigar=cigar.toString();
			cigar.clear();
		}else{
			bbw.skip(n_cigar_op*4);
		}

		//RNEXT
		if(next_refID<0){
			sl.setRnext(bytestar);
		}else if(next_refID==refID){
			sl.setRnext(byteequals);
		}else if(next_refID<refNames.length && SamLine.PARSE_6){
			sl.setRnext(refNames[next_refID].getBytes());
		}else{
			sl.setRnext(bytestar);
		}

		//PNEXT (BAM is 0-based, SAM is 1-based) - SamLine.PARSE_7
		sl.pnext=next_pos+1;

		//TLEN - SamLine.PARSE_8
		sl.tlen=tlen;

		//SEQ
		int numSeqBytes=(int)((l_seq+1)/2);
		if(l_seq==0){
			sl.seq=bytestar;
		}else{
			int seqStart=bbw.position();
			byte[] seq=new byte[(int)l_seq];
			int seqIdx=0;
			for(int i=0; i<numSeqBytes; i++){
				int packed=bamRecord[seqStart+i]&0xFF;
				seq[seqIdx++]=SEQ_LOOKUP_B[packed>>>4];
				if(seqIdx<l_seq){
					seq[seqIdx++]=SEQ_LOOKUP_B[packed&0xF];
				}
			}
			sl.seq=seq;
			bbw.skip(numSeqBytes);
		}

		//QUAL - SamLine.PARSE_10
		if(l_seq==0){
			sl.qual=bytestar;
		}else{
			int qualStart=bbw.position();
			byte firstByte=bamRecord[qualStart];
			if(firstByte==(byte)0xFF){
				bbw.skip((int)l_seq);
				sl.qual=bytestar;
			}else if(SamLine.PARSE_10){
				byte[] qual=new byte[(int)l_seq];
				System.arraycopy(bamRecord, qualStart, qual, 0, (int)l_seq);
				bbw.skip((int)l_seq);
				sl.qual=qual;
			}else{
				bbw.skip((int)l_seq);
			}
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
						tag.appendColon('f').append(bbw.getFloat(), 6);
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
								case 'f': tag.append(bbw.getFloat(), 6); break;
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

				if(tagString.startsWith("MD:")){
					sl.mdTag=tagString.substring(5).getBytes();
				}
			}
		}

		if(sl.mapped() && sl.strand()==Shared.MINUS && SamLine.FLIP_ON_LOAD){
			if(sl.seq!=bytestar){Vector.reverseComplementInPlaceFast(sl.seq);}
			if(sl.qual!=bytestar){Vector.reverseInPlace(sl.qual);}
		}

		sl.trimNames();
		return sl;
	}
	
	private static final byte[] bytestar=new byte[] {(byte)'*'};
	private static final byte[] byteequals=new byte[] {(byte)'='};
	
}
