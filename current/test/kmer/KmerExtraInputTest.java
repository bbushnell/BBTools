package test.kmer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;

/**
 * Finite extra-file expansion checks through construction and real input loading.
 * Run one case per fresh JVM, with a separate dir for each assertion mode/case.
 * Existing input files are never overwritten. No production internals are replaced.
 * @author Sucrose
 */
public final class KmerExtraInputTest {

	/** Accepts case=expanded|explicit|literal-hash and dir= for retained test inputs. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(dir==null || !("expanded".equals(mode) || "explicit".equals(mode) || "literal-hash".equals(mode))){
			throw new IllegalArgumentException("Expected case=expanded|explicit|literal-hash and dir=<new input directory>");
		}
		Files.createDirectories(dir);
		final Path main=dir.resolve("main.fq"), first=dir.resolve("extra1.fq"), second=dir.resolve("extra2.fq");
		write(main, "main", "AAA");
		write(first, "extra1", "ACG");
		write(second, "extra2", "GAA");
		final boolean literal=mode.equals("literal-hash");
		final String extra;
		if(literal){
			final Path path=dir.resolve("literal#.fq");
			write(path, "literal", "ACG");
			extra=path.toString();
		}else if(mode.equals("expanded")){
			extra=dir.resolve("extra#.fq").toString();
		}else{extra=first+","+second;}

		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		final KmerTableSet tables=new KmerTableSet(new String[] {
			"in="+main, "extra="+extra, "k=3", "ways=1", "threads=1", "prealloc=f",
			"prefilter=0", "onepass=f", "rcomp=t", "qin=33", "interleaved=f",
			"minprob=0", "minavgquality=0", "qtrim=f", "ecco=f", "merge=f"
		}, 12);
		final long added=tables.processInput();
		final int expected=literal ? 2 : 3;
		check(!tables.errorState, "Input loading reported an error");
		check(tables.readsIn==expected && tables.basesIn==3L*expected, "Wrong input read/base totals");
		check(added==expected && tables.kmersLoaded==expected && tables.kmersIn==expected,
			"Every three-base read must contribute one distinct k-mer");
		// Literal canonical two-bit encodings: AAA -> TTT=63, ACG -> CGT=27, GAA -> TTC=61.
		for(int key=0; key<64; key++){
			final int count=(key==63 || key==27 || (!literal && key==61)) ? 1 : 0;
			check(Math.max(0, tables.getCount(key))==count, "Wrong count for canonical key "+key);
		}
		System.out.println("PASS case="+mode+" reads="+expected+" distinctKmers="+expected);
	}

	/** Writes one high-quality unpaired read, retaining exact input bytes for reproduction. */
	private static void write(final Path path, final String id, final String bases) throws Exception{
		Files.write(path, ("@"+id+"\n"+bases+"\n+\nIII\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Keeps the oracle active with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
