package test.clump2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;

/**
 * Canonicalizes bounded four-line FASTQ fixtures as unsigned-byte-sorted,
 * NUL-terminated complete records. Preserves duplicates, names, bases, separators
 * and qualities. All records are retained in heap; this is a verification helper,
 * not a streaming production sorter. Its output matches LC_ALL=C sort -z after
 * joining each four-line record with newlines and replacing its final newline by NUL.
 * @author Shinobu
 */
public final class CanonicalFastq{

	/** Takes in=, out= and expected=; rejects malformed or unexpectedly sized data. */
	public static void main(String[] args){
		String in=null, out=null;
		int expected=-1;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("in")){in=value;}
			else if(key.equalsIgnoreCase("out")){out=value;}
			else if(key.equalsIgnoreCase("expected")){expected=Parse.parseIntKMG(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(in==null || out==null || expected<0){throw new IllegalArgumentException("Required: in= out= expected=");}
		assert(COMPARE.compare(new byte[] {(byte)128}, new byte[] {127})>0) : "Canonical order must use unsigned bytes";
		assert(COMPARE.compare(new byte[] {65, 0}, new byte[] {65, 10, 0})<0) : "NUL terminates a shorter canonical record";
		final ArrayList<byte[]> records=new ArrayList<byte[]>(expected);
		final ByteFile reader=ByteFile.makeByteFile1(in, false);
		try{
			for(byte[] name=reader.nextLine(); name!=null; name=reader.nextLine()){
				final byte[] bases=reader.nextLine(), plus=reader.nextLine(), qualities=reader.nextLine();
				if(bases==null || plus==null || qualities==null || name.length<1 || name[0]!='@'
						|| plus.length<1 || plus[0]!='+' || bases.length!=qualities.length){
					throw new IllegalArgumentException("Malformed four-line FASTQ record "+records.size());
				}
				if(records.size()>=expected){throw new IllegalArgumentException("More than "+expected+" records");}
				final long length=(long)name.length+bases.length+plus.length+qualities.length+4;
				if(length>Integer.MAX_VALUE){throw new IllegalArgumentException("Canonical record too large: "+length);}
				final byte[] record=new byte[(int)length];
				int offset=append(record, 0, name);
				offset=append(record, offset, bases);
				offset=append(record, offset, plus);
				offset=append(record, offset, qualities);
				record[offset-1]=0;
				records.add(record);
			}
		}finally{
			if(reader.close()){throw new IllegalStateException("Failed to close input "+in);}
		}
		if(records.size()!=expected){throw new IllegalArgumentException("Expected "+expected+" records, got "+records.size());}
		Collections.sort(records, COMPARE);
		final ByteStreamWriter writer=new ByteStreamWriter(out, false, false, false, FileFormat.TEXT);
		writer.start();
		for(byte[] record : records){writer.print(record);}
		if(writer.poisonAndWait()){throw new IllegalStateException("Failed to write "+out);}
		System.out.println("records\t"+records.size());
	}

	/** Copies one complete field and its newline, rejecting embedded record delimiters. */
	private static int append(final byte[] dest, final int offset, final byte[] field){
		for(byte b : field){if(b==0){throw new IllegalArgumentException("Embedded NUL in FASTQ field");}}
		System.arraycopy(field, 0, dest, offset, field.length);
		dest[offset+field.length]='\n';
		return offset+field.length+1;
	}

	/** Lexicographic unsigned-byte comparison, including the terminating NUL. */
	private static final Comparator<byte[]> COMPARE=new Comparator<byte[]>(){
		@Override
		public int compare(byte[] a, byte[] b){
			final int lim=Math.min(a.length, b.length);
			for(int i=0; i<lim; i++){
				final int difference=(a[i]&255)-(b[i]&255);
				if(difference!=0){return difference;}
			}
			return Integer.compare(a.length, b.length);
		}
	};
}
