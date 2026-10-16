package test.structures;

import fileIO.ReadWrite;
import structures.CoverageArray;
import structures.CoverageArray2;

/** Old-writer to fixed-reader serialization fixture for CoverageArray.devSum.
 * @author Sucrose
 */
public final class CoverageArraySerializationCase{

	/** Writes an old CoverageArray or reads it after the keyed devSum repair. */
	public static void main(final String[] args){
		String mode=null, path=null;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			check(split.length==2, "Expected flag=value, got "+arg);
			if("mode".equals(split[0])){mode=split[1];}
			else if("file".equals(split[0])){path=split[1];}
			else{throw new AssertionError("Unknown argument: "+arg);}
		}
		check(mode!=null, "Missing mode=write or mode=read");
		check(path!=null, "Missing file=path.ca");
		if("write".equals(mode)){
			writeLegacyObject(path);
		}else if("read".equals(mode)){
			readFixedObject(path);
		}else{
			throw new AssertionError("Expected write or read, got "+mode);
		}
		System.out.println("PASS\t"+mode+'\t'+path);
	}

	/** Writes a CoverageArray2 carrying a stale legacy unkeyed devSum value. */
	private static void writeLegacyObject(final String path){
		final CoverageArray2 ca=new CoverageArray2(3, 4);
		for(int i=0; i<4; i++){
			ca.set(i, 2);
		}
		checkEquals("legacy mean-4 devSum", 16, ca.devSum(4));
		ReadWrite.write(ca, path, false);
	}

	/** Reads the old object and verifies fixed devSum ignores its stale cache. */
	private static void readFixedObject(final String path){
		final CoverageArray ca=CoverageArray.read(path);
		check(ca instanceof CoverageArray2, "Expected CoverageArray2, got "+ca.getClass());
		check(ca.chromosome==3, "chromosome expected 3, got "+ca.chromosome);
		check(ca.length()==4, "length expected 4, got "+ca.length());
		check(ca.minIndex==0, "minIndex expected 0, got "+ca.minIndex);
		check(ca.maxIndex==3, "maxIndex expected 3, got "+ca.maxIndex);
		for(int i=0; i<4; i++){
			check(ca.get(i)==2, "depth "+i+" expected 2, got "+ca.get(i));
		}
		checkEquals("fixed mean-2 devSum", 0, ca.devSum(2));
		checkEquals("fixed mean-2 cached devSum", 0, ca.devSum(2));
	}

	/** Checks an exact double oracle. */
	private static void checkEquals(final String name, final double expected,
			final double observed){
		check(Double.compare(expected, observed)==0,
				name+" expected "+expected+", got "+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
