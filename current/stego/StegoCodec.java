package stego;

import java.util.Arrays;
import java.util.zip.CRC32;

/** Sequence-only framing for StegoFast version 1. No header or quality dependency.
 * @author Yaoyao */
public final class StegoCodec{

	/** Twenty-base signature; version, integers and payload follow it. */
	public static final String MAGIC="ACGTTGCACTGATCGATGCA";
	public static final int VERSION=1;
	private static final byte[] BASES={'A','C','G','T'};
	private StegoCodec(){}

	/** Integers use a three-base digit-count-minus-one prefix, then big-endian base 4. */
	public static int integerBases(long value){
		if(value<0){throw new IllegalArgumentException("Negative StegoFast integer");}
		return 3+Math.max(1, (64-Long.numberOfLeadingZeros(value)+1)/2);
	}

	public static int capacity(int length, long group, long part){
		final int fixed=MAGIC.length()+2+integerBases(group)+integerBases(part)+16;
		int n=Math.max(0, (length-fixed-4)/4);
		while(n>0 && fixed+integerBases(n)+4*n>length){n--;}
		return n;
	}

	public static byte[] encode(int length, long group, long part, byte[] data, int count){
		if(count<=0 || count>data.length || count>capacity(length, group, part)){
			throw new IllegalArgumentException("Payload does not fit metadata read");
		}
		byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'A');
		for(int i=0; i<MAGIC.length(); i++){bases[i]=(byte)MAGIC.charAt(i);}
		Cursor c=new Cursor(bases, MAGIC.length(), length-16, false);
		c.put(VERSION, 2);
		c.putInteger(group);
		c.putInteger(part);
		c.putInteger(count);
		for(int i=0; i<count; i++){c.put(data[i]&255, 4);}
		CRC32 crc=new CRC32();
		crc.update(bases, 0, length-16);
		new Cursor(bases, length-16, length, false).put(crc.getValue(), 16);
		return bases;
	}

	/** Returns null for ordinary sequence; a recognized but invalid envelope throws. */
	public static Chunk decode(byte[] bases){
		if(bases==null || bases.length<MAGIC.length()){return null;}
		final boolean forward=magic(bases, false), reverse=magic(bases, true);
		if(!forward && !reverse){return null;}
		if(forward){
			try{return decode(bases, false);}catch(IllegalArgumentException e){if(!reverse){throw e;}}
		}
		return decode(bases, true);
	}

	private static Chunk decode(byte[] bases, boolean reverse){
		if(bases.length<48){throw bad("Truncated metadata envelope");}
		final int end=bases.length-16;
		CRC32 crc=new CRC32();
		for(int i=0; i<end; i++){crc.update(BASES[digit(bases, i, reverse)]);}
		long expected=new Cursor(bases, end, bases.length, reverse).get(16);
		if(crc.getValue()!=expected){throw bad("Metadata read CRC mismatch");}
		Cursor c=new Cursor(bases, MAGIC.length(), end, reverse);
		if(c.get(2)!=VERSION){throw bad("Unsupported metadata version");}
		long group=c.getInteger(), part=c.getInteger(), count=c.getInteger();
		if(count<=0 || count>(end-c.pos)/4){throw bad("Invalid metadata payload length");}
		byte[] data=new byte[(int)count];
		for(int i=0; i<data.length; i++){data[i]=(byte)c.get(4);}
		while(c.pos<end){if(c.get(1)!=0){throw bad("Nonzero metadata padding");}}
		return new Chunk(group, part, data, bases.length);
	}

	private static boolean magic(byte[] bases, boolean reverse){
		for(int i=0; i<MAGIC.length(); i++){
			int x=code(bases[reverse ? bases.length-1-i : i]);
			if(x<0 || BASES[reverse ? 3-x : x]!=MAGIC.charAt(i)){return false;}
		}
		return true;
	}

	private static int code(byte b){
		switch(b){
			case 'A': case 'a': return 0;
			case 'C': case 'c': return 1;
			case 'G': case 'g': return 2;
			case 'T': case 't': return 3;
			default: return -1;
		}
	}

	private static int digit(byte[] bases, int pos, boolean reverse){
		int x=code(bases[reverse ? bases.length-1-pos : pos]);
		if(x<0){throw bad("Non-ACGT base in metadata read");}
		return reverse ? 3-x : x;
	}

	private static IllegalArgumentException bad(String message){return new IllegalArgumentException(message);}

	private static final class Cursor{
		Cursor(byte[] b, int p, int e, boolean r){bases=b; pos=p; end=e; reverse=r;}
		void put(long value, int digits){
			if(pos+digits>end){throw bad("Metadata envelope overflow");}
			for(int i=digits-1; i>=0; i--){bases[pos++]=BASES[(int)(value>>>(2*i))&3];}
		}
		long get(int digits){
			if(pos+digits>end){throw bad("Truncated metadata integer");}
			long value=0;
			for(int i=0; i<digits; i++){value=(value<<2)|digit(bases, pos++, reverse);}
			return value;
		}
		void putInteger(long value){
			int digits=integerBases(value)-3;
			put(digits-1, 3);
			put(value, digits);
		}
		long getInteger(){
			int digits=(int)get(3)+1;
			if(digits>32){throw bad("Metadata integer too wide");}
			long value=get(digits);
			if(value<0 || integerBases(value)!=digits+3){throw bad("Noncanonical metadata integer");}
			return value;
		}
		final byte[] bases;
		final int end;
		final boolean reverse;
		int pos;
	}

	public static final class Chunk{
		Chunk(long g, long p, byte[] d, int l){group=g; part=p; data=d; length=l;}
		public final long group, part;
		public final byte[] data;
		public final int length;
	}
}
