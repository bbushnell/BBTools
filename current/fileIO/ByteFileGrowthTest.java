package fileIO;

import java.io.FileOutputStream;
import java.util.Arrays;
import shared.Shared;
import structures.IntList;

/** Regression inputs and exact sequence checks for the reader growth repair. @author Raiden */
public final class ByteFileGrowthTest {

	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Expected mode=generate|check|controls and path=FILE";
		final String mode=args[0].substring("mode=".length()), path=args[1].substring("path=".length());
		assert(args[0].startsWith("mode=") && args[1].startsWith("path=")) : "Use mode= and path=";
		Shared.setThreads(1);
		if(mode.equals("generate")){
			final byte[] block=new byte[65536];
			Arrays.fill(block, (byte)'A');
			try(FileOutputStream out=new FileOutputStream(path)){
				out.write(">repro_gt1g\n".getBytes("US-ASCII"));
				for(long written=0; written<PREFIX;){
					final int n=(int)Math.min(block.length, PREFIX-written);
					out.write(block, 0, n);
					written+=n;
				}
				out.write((TAIL+"\n").getBytes("US-ASCII"));
			}
		}else if(mode.equals("check")){
			long bases=0;
			int records=0;
			final ByteFile1 in=new ByteFile1(path, true);
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				if(line.length>0 && line[0]=='>'){records++; continue;}
				for(byte b : line){
					assert(bases<PREFIX+TAIL.length()) : "Extra sequence at "+bases;
					final byte expected=bases<PREFIX ? (byte)'A' : (byte)TAIL.charAt((int)(bases-PREFIX));
					assert(b==expected) : "Wrong base at "+bases+": "+(char)b+" != "+(char)expected;
					bases++;
				}
			}
			assert(!in.close()) : "Input read failed";
			assert(records==1 && bases==PREFIX+TAIL.length()) : "Expected one intact record; records="+records+", bases="+bases;
		}else if(mode.equals("controls")){
			final ByteFile1Fc in=new ByteFile1Fc(path, true);
			final IntList ends=new IntList();
			final StringBuilder observed=new StringBuilder();
			for(byte[] block=in.nextLine(ends); block!=null; block=in.nextLine(ends)){
				assert(ends.size>0) : "Record block lacks boundaries";
				observed.append(new String(block, 0, ends.get(ends.size-1)+1, "US-ASCII"));
			}
			assert(!in.close()) : "Control input read failed";
			assert(observed.toString().equals(">a\nACGT\n>b\nACGTACGT\n")) : "EOF control changed: "+observed;
		}else{throw new IllegalArgumentException("Unknown mode "+mode);}
		System.out.println("BYTEFILE_GROWTH_TEST_PASS mode="+mode);
	}

	private static final long PREFIX=1100000000L;
	private static final String TAIL="CGTACGTACGTACGTTAG";
}
